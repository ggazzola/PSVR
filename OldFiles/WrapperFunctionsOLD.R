GenerateData = function(){
	covMat <<- DoCovMat(p=p, corVal=corVal, stdVect=stdVect)
	doDatOut <<- DoDat(n=n, meanVect=meanVect, covMat=covMat, trueW=trueW, trueW0=trueW0, noiseLevel=noiseLevel) # sample data
	doMissOut <<- DoMiss(dat=doDatOut, missingY=missingY, missingVarProb=missingVarProb) # inject missing values
	doDataSplitOutOuter <<- DoOriginalMissingDataSplit(doDatOut=doDatOut, doMissOut = doMissOut, numFolds=numFolds) #partition data (with missing values) in multiple training/testing sets using STRATIFIED x-fold validation -- CONSIDER NON-STRATIFIED???(same in inner below)
}


CalculateValidationErrors = function(){
	doParListGridOut = DoParListGrid(parValuesList = parValuesList)
	doErrorFoldOutInnerList  = list()
	for(i in 1:length(doDataSplitOutOuter)){
		#inner cross-validation
		currTrain = doDataSplitOutOuter[[i]]$inDat$missing #current training portion
		currTrainOriginal = doDataSplitOutOuter[[i]]$inDat$original # use for debugging, then eliminate
		currDoDataSplitOutInner = DoOriginalMissingDataSplit(doDatOut=currTrain, doMissOut=currTrain, numFolds=numFolds) # partition training into multiple training/validation sets

		doMultipleImputationFoldsOutInner = DoMultipleImputationFolds(doMissOut = currDoDataSplitOutInner, method=method, maxIter=maxIter, numImput=numImput)
		# for each training/validation set, run multiple imputation; training set by itself and testing set using both training and testing, and then disregarding training
		doErrorFoldOutInnerListTmp = list()
		for(k in 1:length(pcPropValues)){
			doPolyListFoldsOutInner = DoPolyListFolds(doMultipleImputationFoldsOut=doMultipleImputationFoldsOutInner, pcProp=pcPropValues[k], scaleData=scaleData, maxUncertainDims=maxUncertainDims) # for each training/testing pair, adding the polytope representation of the training set
		
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


TrainOnBestModel = function(){
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
	
		currTrainImput = DoMultipleImputation(missDat=currTrain, method=method, numImput=numImput, maxIter=maxIter)
		currTestImput = DoMultipleImputation(missDat=currTest, method=method, numImput=numImput, maxIter=maxIter, 
			extraPossiblyMissDat=currTrain)
		
		currTrainMedianImput = currTrainImput$medianImputDat
		currTestMedianImput = currTestImput$medianImputDat
	
		getTrainResReadyForTest[[i]]$currTrain = currTrain
		getTrainResReadyForTest[[i]]$currTest = currTest
		getTrainResReadyForTest[[i]]$missingTestDataLogical = missingTestDataLogical
		getTrainResReadyForTest[[i]]$currTrainImput = currTrainImput
		getTrainResReadyForTest[[i]]$currTestImput = currTestImput
		getTrainResReadyForTest[[i]]$currTrainMedianImput = currTrainMedianImput
		getTrainResReadyForTest[[i]]$currTestMedianImput = currTestMedianImput

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
		
		currTrainMedianImput = getTrainResReadyForTest[[i]]$currTestMedianImput
		currTestMedianImput = getTrainResReadyForTest[[i]]$currTestMedianImput
	
		if(approach =="doMedian"){
	   	 # compare with median imput-training
			
			currTrainMedianFakeImput = DoMultipleImputation(missDat=currTrainMedianImput, method="irrelevant", numImput=2, maxIter=0) 
			#the 'multiple imputation' in $imputDatList won't be used, we get it just so that the stopifnot functions in DoPolyList
			# don't complain; the $
		
			currTrainPolyList = DoPolyList(missDat=currTrainMedianImput, imputDatList=currTrainMedianFakeImput$imputDatList, 
				medianImputDat=currTrainMedianImput, pcProp=bestParList$pcProp,  
				scaleData=scaleData, maxUncertainDims=maxUncertainDims) 
		}	
	
		if(approach == "doPCbb"){
	   	 # compare with traditional-box-training
		
			currTrainPolyList = DoPolyList(missDat=currTrain, imputDatList=currTrainImput$imputDatList, 
				medianImputDat=currTrainImput$medianImputDat, pcProp=bestParList$pcProp,  
				scaleData=scaleData, maxUncertainDims=maxUncertainDims)
		}	
	
		if(approach == "doSquarebb"){
			stop("don't know what to do\n")
		}
	
		#if(scaleData)
		#	currTestMedianImput = ScaleCenter(currTestMedianImput, currTrainPolyList$scaleInfo$mean, currTrainPolyList$scaleInfo$std)
		# NO! this is taken care of by doErrorList, via attributes of currTrainPolyList
		currModel = DoTrainModel(polyList=currTrainPolyList, parList = bestParList)
		currError = DoErrorList(doTrainModelOut=currModel, medianImputOut=currTestMedianImput, 
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