GenerateData = function(){
	covMat <<- DoCovMat(p=p, corVal=corVal, stdVect=stdVect)
	doDatOut <<- DoDat(n=n, meanVect=meanVect, covMat=covMat, trueW=trueW, trueW0=trueW0, theoRsq=theoRsq) # sample data
	doMissOut <<- DoMiss(dat=doDatOut, missingY=missingY, missingObsProp = missingObsProp, missingVarProp=missingVarProp) # inject missing values
	doDataSplitOutOuter = DoOriginalMissingDataSplit(doDatOut=doDatOut, doMissOut = doMissOut, numFolds=numFolds) #partition data (with missing values) in multiple training/testing sets using STRATIFIED x-fold validation -- CONSIDER NON-STRATIFIED???(same in inner below)
	
	for(i in 1:length(doDataSplitOutOuter)){
		currTrain = doDataSplitOutOuter[[i]]$inDat$missing #current training portion
		currDoDataSplitOutInner = DoOriginalMissingDataSplit(doDatOut=currTrain, doMissOut=currTrain, numFolds=numFolds) # partition training into multiple training/validation sets
		doDataSplitOutOuter[[i]]$innerSplit = currDoDataSplitOutInner
	}
	doDataSplitOutOuter<<-doDataSplitOutOuter
}

MultiplyImpute = function(){
	doDataSplitOutOuter <<- DoMultipleImputationFolds(doMissOut = doDataSplitOutOuter,	method=method, maxIter=maxIter, 
		numImput=numImput)
		cat("Outer multiple imputations done\n")
	for(i in 1:length(doDataSplitOutOuter)){
		doMultipleImputationFoldsOutInner = DoMultipleImputationFolds(doMissOut = doDataSplitOutOuter[[i]]$innerSplit,
			method=method, maxIter=maxIter, numImput=numImput)
		doDataSplitOutOuter[[i]]$innerSplit = doMultipleImputationFoldsOutInner
		cat(i, "-th inner multiple imputations done\n")
	}
	doDataSplitOutOuter<<-doDataSplitOutOuter
	#doDataSplitOutOuter[[i]]$inDat : i-th training data ($imputed/$missing/$original)
	#doDataSplitOutOuter[[i]]$outDat : i-th testing data ($imputed/$missing/$original)
	#doDataSplitOutOuter[[i]]$innerSplit : set of training/validation data + imputations
	#doDataSplitOutOuter[[i]]$innerSplit[[j]]$inDat: j-th training data within i-th training data ($imputed/$missing/$original)
	#doDataSplitOutOuter[[i]]$innerSplit[[j]]$outDat: j-th validation data within i-th training data ($imputed/$missing/$original)
	
}

CalculateValidationErrors = function(){
	if(approach=="doMedian" | approach=="doNoMiss"){
		parValuesListNoUncert = parValuesList
		parValuesListNoUncert$Cuncertain = ifelse(approach=="doMedian", "irrelevantForDoMedian", "irrelevantForDoMiss")
		parValuesListNoUncert$extraEpsilonUncertain = ifelse(approach=="doMedian", "irrelevantForDoMedian", "irrelevantForDoMiss")
		parValuesListSquaredBB$uncertaintySpecialTreatment = F
		doParListGridOut = DoParListGrid(parValuesList = parValuesListNoUncert)
		quantOrSdPropValuesVect = ifelse(approach=="doMedian", "irrelevantForDoMedian", "irrelevantForDoMiss")
	} else if(approach=="doSquarebbSd" | approach=="doSquarebbQuant"){
		parValuesListSquaredBB = parValuesList
		parValuesListSquaredBB$Cuncertain = ifelse(approach=="doSquarebbSd", "irrelevantForDoSquarebbSd", "irrelevantForDoSquarebbQuant")
		parValuesListSquaredBB$extraEpsilonUncertain = ifelse(approach=="doSquarebbSd", "irrelevantForDoSquarebbSd", "irrelevantForDoSquarebbQuant")
		parValuesListSquaredBB$uncertaintySpecialTreatment = F
		doParListGridOut = DoParListGrid(parValuesList = parValuesListSquaredBB)
		quantOrSdPropValuesVect = quantOrSdPropValues
	} else	if(approach=="doPCbb"){
		doParListGridOut = DoParListGrid(parValuesList = parValuesList)
		quantOrSdPropValuesVect = quantOrSdPropValues
	} else{
		stop("Don't know what to do\n")		
	}
	doErrorFoldOutInnerList  = list()
	for(i in 1:length(doDataSplitOutOuter)){
		#inner cross-validation
		doErrorFoldOutInnerListTmp = list()
		for(k in 1:length(quantOrSdPropValuesVect)){
			doPolyListFoldsOutInner = DoPolyListFolds(doMultipleImputationFoldsOut=doDataSplitOutOuter[[i]]$innerSplit, quantOrSdProp=quantOrSdPropValuesVect[k], scaleData=scaleData, maxUncertainDims=maxUncertainDims, doMedian=approach=="doMedian", doNoMiss=approach=="doNoMiss", doSquarebbSd=approach=="doSquarebbSd", doSquarebbQuant=approach=="doSquarebbQuant") # for each training/testing pair, adding the polytope representation of the training set
			currPcPropErrFold = DoErrorFold(doPolyListFoldsOut=doPolyListFoldsOutInner, doParListGridOut=doParListGridOut, replaceImputedWithTrueY = replaceImputedWithTrueY)
			for(kk in 1:length(currPcPropErrFold)){
				currPcPropErrFold[[kk]]$parList$quantOrSdProp = quantOrSdPropValues[k]
			}
			doErrorFoldOutInnerListTmp = c(doErrorFoldOutInnerListTmp, currPcPropErrFold) # doErrorFoldOutInnerListTmp[[j]] contains results for the j-th model parameter combination over each of the folds
		}
		doErrorFoldOutInnerList[[i]] = DoExtractErrMat(doErrorFoldOut = doErrorFoldOutInnerListTmp) # these are cross validation results for the i-th training data set (divided into numFolds training/validation); doErrorFoldOutInnerList[[i]][[j]][[k]] are the results for the j-th model parameter combination, in the k-th training/validation inner partition of the i-th outer training/testing fold
		cat("Inner cross-validation", i, "out of ", length(doDataSplitOutOuter), "done\n")
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
			medianImputDat=currTrainImput$medianImputDat, quantOrSdProp=bestParList$quantOrSdProp,  
			scaleData=scaleData, maxUncertainDims=maxUncertainDims, doMedian=approach=="doMedian", doNoMiss=approach=="doNoMiss",
			doSquarebbSd=approach=="doSquarebbSd", doSquarebbQuant=approach=="doSquarebbQuant")	
		
	
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