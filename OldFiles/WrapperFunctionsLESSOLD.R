GenerateData = function(){
	covMat <<- DoCovMat(p=p, corVal=corVal, stdVect=stdVect)
	doDatOut <<- DoDat(n=n, meanVect=meanVect, covMat=covMat, trueW=trueW, trueW0=trueW0, noiseLevel=noiseLevel) # sample data
	doMissOut <<- DoMiss(dat=doDatOut, missingY=missingY, missingVarProb=missingVarProb) # inject missing values
	doDataSplitOutOuter = DoOriginalMissingDataSplit(doDatOut=doDatOut, doMissOut = doMissOut, numFolds=numFolds) #partition data (with missing values) in multiple training/testing sets using STRATIFIED x-fold validation -- CONSIDER NON-STRATIFIED???(same in inner below)
	
	for(i in 1:length(doDataSplitOutOuter)){
		currTrain = doDataSplitOutOuter[[i]]$inDat$missing #current training portion
		currDoDataSplitOutInner = DoOriginalMissingDataSplit(doDatOut=currTrain, doMissOut=currTrain, numFolds=numFolds) # partition training into multiple training/validation sets
		doDataSplitOutOuter[[i]]$innerSplit = currDoDataSplitOutInner
	}
	doDataSplitOutOuter<<-doDataSplitOutOuter
}

MultiplyImpute = function(){
	doDataSplitOutOuter =  DoMultipleImputationFolds(doMissOut = doDataSplitOutOuter,	method=method, maxIter=maxIter, 
		numImput=numImput)
	for(i in 1:length(doDataSplitOutOuter)){
		doMultipleImputationFoldsOutInner = DoMultipleImputationFolds(doMissOut = doDataSplitOutOuter[[i]]$innerSplit,
			method=method, maxIter=maxIter, numImput=numImput)
		doDataSplitOutOuter[[i]]$innerSplit = doMultipleImputationFoldsOutInner
	}
	doDataSplitOutOuter<<-doDataSplitOutOuter
	#doDataSplitOutOuter[[i]]$inDat : i-th training data ($imputed/$missing/$original)
	#doDataSplitOutOuter[[i]]$outDat : i-th testing data ($imputed/$missing/$original)
	#doDataSplitOutOuter[[i]]$innerSplit : set of training/validation data + imputations
	#doDataSplitOutOuter[[i]]$innerSplit[[j]]$inDat: j-th training data within i-th training data ($imputed/$missing/$original)
	#doDataSplitOutOuter[[i]]$innerSplit[[j]]$outDat: j-th validation data within i-th training data ($imputed/$missing/$original)
	
}

CalculateValidationErrors = function(){
	if(approach=="doMedian"){
		parValuesListNoUncert = parValuesList
		parValuesListNoUncert$Cuncertain = 0 # these are 
		parValuesListNoUncert$extraEpsilonUncertain = 0 # actually irrelevant
		doParListGridOut = DoParListGrid(parValuesList = parValuesListNoUncert)
		pcPropValuesVect = "irrelevantForDoMedian"
	} else	if(approach=="doPCbb"){
		doParListGridOut = DoParListGrid(parValuesList = parValuesList)
		pcPropValuesVect = pcPropValues
	} else{
		stop("Don't know what to do\n")		
	}
	doErrorFoldOutInnerList  = list()
	for(i in 1:length(doDataSplitOutOuter)){
		#inner cross-validation
		doErrorFoldOutInnerListTmp = list()
		for(k in 1:length(pcPropValuesVect)){
			doPolyListFoldsOutInner = DoPolyListFolds(doMultipleImputationFoldsOut=doDataSplitOutOuter[[i]]$innerSplit, pcProp=pcPropValues[k], scaleData=scaleData, maxUncertainDims=maxUncertainDims, doMedian=approach=="doMedian") # for each training/testing pair, adding the polytope representation of the training set
		
			currPcPropErrFold = DoErrorFold(doPolyListFoldsOut=doPolyListFoldsOutInner, doParListGridOut=doParListGridOut, replaceImputedWithTrueY = replaceImputedWithTrueY)
			for(kk in 1:length(currPcPropErrFold)){
				currPcPropErrFold[[kk]]$parList$pcProp = pcPropValues[k]
			}
			doErrorFoldOutInnerListTmp = c(doErrorFoldOutInnerListTmp, currPcPropErrFold) # doErrorFoldOutInnerListTmp[[j]] contains results for the j-th model parameter combination over each of the folds
		}
		doErrorFoldOutInnerList[[i]] = DoExtractErrMat(doErrorFoldOut = doErrorFoldOutInnerListTmp) # these are cross validation results for the i-th training data set (divided into numFolds training/validation); doErrorFoldOutInnerList[[i]][[j]] are the results for the j-th model parameter combination
	}
	doErrorFoldOutInnerList<<-doErrorFoldOutInnerList
}

SetUpTest = function(){
	getTrainResReadyForTest = list()

	for(i in 1:length(doErrorFoldOutInnerList)){
		getTrainResReadyForTest[[i]] = list()
		tmp = DoBestParList(doExtractErrMatOut=doErrorFoldOutInnerList[[i]], errMeasure=errMeasure) # best parameter combination and corresponding error
		getTrainResReadyForTest[[i]]$bestParList = tmp$bestParList
		getTrainResReadyForTest[[i]]$bestAvgErr = tmp$bestParListAvgError
	
		 # best to keep this outside of the former loop, so we can't change errMeasure without having to recompute everything else, which should be fixed
	
		currTrain = doDataSplitOutOuter[[i]]$inDat$missing
		currTest = doDataSplitOutOuter[[i]]$outDat$missing
		missingTestDataLogical = apply(currTest, 1, function(x)any(is.na(x)))
	
		currTrainImput = doDataSplitOutOuter[[i]]$inDat$imputed
		currTestImput = doDataSplitOutOuter[[i]]$outDat$imputed
		
		#currTrainMedianImput = doDataSplitOutOuter[[i]]$inDat$imputed$medianImputDat
		#currTestMedianImput = doDataSplitOutOuter[[i]]$outDat$imputed$medianImputDat
	
		getTrainResReadyForTest[[i]]$currTrain = currTrain
		getTrainResReadyForTest[[i]]$currTest = currTest
		getTrainResReadyForTest[[i]]$missingTestDataLogical = missingTestDataLogical
		getTrainResReadyForTest[[i]]$currTrainImput = currTrainImput
		getTrainResReadyForTest[[i]]$currTestImput = currTestImput
		#getTrainResReadyForTest[[i]]$currTrainMedianImput = currTrainMedianImput
		#getTrainResReadyForTest[[i]]$currTestMedianImput = currTestMedianImput

	}
	
	getTrainResReadyForTest<<-getTrainResReadyForTest
}

CalculateTestErrors = function(){
	errVect = NULL

	for(i in 1:length(getTrainResReadyForTest)){

		 # best to keep this outside of the former loop, so we can't change 'what' without having to recompute everything else, which should be fixed
	 
		bestParList = getTrainResReadyForTest[[i]]$bestParList
		currTrain = getTrainResReadyForTest[[i]]$currTrain
		currTest = getTrainResReadyForTest[[i]]$currTest
		missingTestDataLogical = getTrainResReadyForTest[[i]]$missingTestDataLogical
	
		currTrainImput = getTrainResReadyForTest[[i]]$currTrainImput
		currTestImput = getTrainResReadyForTest[[i]]$currTestImput
		

		currTrainPolyList = DoPolyList(missDat=currTrain, imputDatList = currTrainImput$imputDatList,
			medianImputDat=currTrainImput$medianImputDat, pcProp=bestParList$pcProp,  
			scaleData=scaleData, maxUncertainDims=maxUncertainDims, doMedian=approach=="doMedian")	
		
	
		if(approach == "doSquarebb"){
			stop("don't know what to do\n")
		}
	
		#if(scaleData)
		#	currTestMedianImput = ScaleCenter(currTestMedianImput, currTrainPolyList$scaleInfo$mean, currTrainPolyList$scaleInfo$std)
		# NO! this is taken care of by doErrorList, via attributes of currTrainPolyList
		currModel = DoTrainModel(polyList=currTrainPolyList, parList = bestParList)
		currError = DoErrorList(doTrainModelOut=currModel, medianImputOut=currTestImput$medianImputDat, 
			doPolyListOut=currTrainPolyList, missingDatOutLogical = missingTestDataLogical)[[errMeasure]]
		errVect = c(errVect, currError)
	}
	errVect <<-errVect
}

CollectTestPerformance=function(errMes, appr){
	testPerf = NULL
	for(i in 1:length(testRes)){
		testPerf[i] = testRes[[i]][[errMes]][[appr]]
	}
	res = AggregateTestError(testPerf)
	return(res)
}