GenerateData = function(){
	
	success = F
	givenUp = F
	generateDataCnt = 1
	
	if(realData){
		dd = load(paste0(dataFolder, realDataFileName))
		stopifnot("dat"%in%dd)
		stopifnot(is.matrix(dat))
		if(realDataFileName=="Boston.RData"){
			if(nSubSelect<Inf){
				stopifnot(nSubSelect<=nrow(dat))
				dat=dat[sample(nrow(dat), nSubSelect),]
			}
		}
		doDatOut <<- dat
	}
	
	while(!success){
		if(!realData){
			covMat <<- DoCovMat(p=p, corVal=corVal, stdVect=stdVect)
			doDatOut <<- DoDat(n=n, meanVect=meanVect, covMat=covMat, trueW=trueW, trueW0=trueW0, theoRsq=theoRsq) # sample data
		}
		if(injectMissingness & any(is.na(doDatOut)))
			stop("Should not inject missingness in data set that already contains missing values")
		doMissOut <<- DoMiss(dat=doDatOut, missingY=missingY, missingObsProp = missingObsProp, missingVarProp=missingVarProp, doMCAR) # inject missing values
		doDataSplitOutOuter = DoOriginalMissingDataSplit(doDatOut=doDatOut, doMissOut = doMissOut, numFolds=numFolds) #partition data (with missing values) in multiple training/testing sets using STRATIFIED x-fold validation -- CONSIDER NON-STRATIFIED???(same in inner below)
		stopifnot(colnames(doDatOut)[ncol(doDatOut)]=="Y")
	
		for(i in 1:length(doDataSplitOutOuter)){
			currTrain = doDataSplitOutOuter[[i]]$inDat$missing #current training portion
			currDoDataSplitOutInner = DoOriginalMissingDataSplit(doDatOut=currTrain, doMissOut=currTrain, numFolds=numFolds) # partition training into multiple training/validation sets
			doDataSplitOutOuter[[i]]$innerSplit = currDoDataSplitOutInner
		}
		
		doDataSplitOutOuter<<-doDataSplitOutOuter
		HasNA = function(x, opposite=F) {if(!opposite) any(is.na(x)) else all(!is.na(x))}
		for(i in 1:length(doDataSplitOutOuter)){
						
			n1=sum(apply(doDataSplitOutOuter[[i]]$inDat$missing, 1, HasNA))
			n2=sum(apply(doDataSplitOutOuter[[i]]$outDat$missing, 1, HasNA))
			n3=sum(apply(doDataSplitOutOuter[[i]]$inDat$missing, 1, HasNA, opposite=T))
			n4=sum(apply(doDataSplitOutOuter[[i]]$outDat$missing, 1, HasNA, opposite=T))
			if(n1<2 | n2==0 | n3<2 | n4==0){
				success = F				
				if(!rejectSilently){
					cat("Rejecting outer generated data. Starting over\n")
				}	
				break
			}
			
			breakOuter=F
			for(j in 1:length(doDataSplitOutOuter[[i]]$innerSplit)){
				n5 = sum(apply(doDataSplitOutOuter[[i]]$innerSplit[[j]]$inDat$missing, 1, HasNA))
				n6 = sum(apply(doDataSplitOutOuter[[i]]$innerSplit[[j]]$outDat$missing, 1, HasNA))
				n7 = sum(apply(doDataSplitOutOuter[[i]]$innerSplit[[j]]$inDat$missing, 1, HasNA, opposite=T))
				n8 = sum(apply(doDataSplitOutOuter[[i]]$innerSplit[[j]]$outDat$missing, 1, HasNA, opposite=T))
				if(n5<2 | n6==0 | n7<2 | n8==0){
					success = F	
					if(!rejectSilently){
						cat("Rejecting inner generated data. Starting over\n")
					}	
					breakOuter = T
					break
				}
			}
			if(breakOuter)
				break
			if (i == length(doDataSplitOutOuter))
				success = T
		}
		
		if(!success & (generateDataCnt>maxGenerateDataAttempts)){
			givenUp =T
			break
		}
		generateDataCnt = generateDataCnt+1
		if(!success)
			cat("Starting attempt", generateDataCnt, "\n")
	}
	givenUp <<- givenUp
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
		parValuesListNoUncert$Cuncertain = ifelse(approach=="doMedian", "irrelevantForDoMedian", "irrelevantForDoNoMiss")
		parValuesListNoUncert$extraEpsilonUncertain = ifelse(approach=="doMedian", "irrelevantForDoMedian", "irrelevantForDoNoMiss")
		parValuesListNoUncert$uncertaintySpecialTreatment = F
		parValuesListNoUncert$twoSlacks=F
		doParListGridOut = DoParListGrid(parValuesList = parValuesListNoUncert)
		quantOrSdPropValuesVect = ifelse(approach=="doMedian", "irrelevantForDoMedian", "irrelevantForDoNoMiss")
	} else if(approach=="doSquarebbSd" | approach=="doSquarebbQuant"){
		parValuesListSquaredBB = parValuesList
		parValuesListSquaredBB$Cuncertain = ifelse(approach=="doSquarebbSd", "irrelevantForDoSquarebbSd", "irrelevantForDoSquarebbQuant")
		parValuesListSquaredBB$extraEpsilonUncertain = ifelse(approach=="doSquarebbSd", "irrelevantForDoSquarebbSd", "irrelevantForDoSquarebbQuant")
		parValuesListSquaredBB$uncertaintySpecialTreatment = F
		parValuesListSquaredBB$twoSlacks=T
		doParListGridOut = DoParListGrid(parValuesList = parValuesListSquaredBB)
		quantOrSdPropValuesVect = quantOrSdPropValues
	} else	if(approach=="doPCbb"){
		parValuesList$twoSlacks=F
		doParListGridOut = DoParListGrid(parValuesList = parValuesList)
		quantOrSdPropValuesVect = quantOrSdPropValues
	} else{
		stop("Don't know what to do\n")		
	}
	doErrorFoldOutInnerList  = list()
	for(i in 1:length(doDataSplitOutOuter)){
		innerProgressOut = paste("Inner cross-validation", i, "out of ", length(doDataSplitOutOuter), "START at", date(),"\n")
		cat(innerProgressOut)
		write.table(innerProgressOut, quote=F, row.names=F, col.names=F, append=T, file=progressFile)
		#inner cross-validation
		doErrorFoldOutInnerListTmp = list()
		for(k in 1:length(quantOrSdPropValuesVect)){
			innerProgressOut = paste("quantOrSdProp value", k, "out of ", length(quantOrSdPropValuesVect), "START at", date(),"\n")
			cat(innerProgressOut)
			write.table(innerProgressOut, quote=F, row.names=F, col.names=F, append=T, file=progressFile)
	
			doPolyListFoldsOutInner = DoPolyListFolds(doMultipleImputationFoldsOut=doDataSplitOutOuter[[i]]$innerSplit, quantOrSdProp=quantOrSdPropValuesVect[k], scaleData=scaleData, maxUncertainDims=maxUncertainDims, doMedian=approach=="doMedian", doNoMiss=approach=="doNoMiss", doSquarebbSd=approach=="doSquarebbSd", doSquarebbQuant=approach=="doSquarebbQuant") # for each training/testing pair, adding the polytope representation of the training set
			currPcPropErrFold = DoErrorFold(doPolyListFoldsOut=doPolyListFoldsOutInner, doParListGridOut=doParListGridOut, 
				replaceImputedWithTrueY = replaceImputedWithTrueY, approach=approach) # this calculates errors according to all error measures (see DoError)
			for(kk in 1:length(currPcPropErrFold)){
				currPcPropErrFold[[kk]]$parList$quantOrSdProp = quantOrSdPropValuesVect[k]
			}
			doErrorFoldOutInnerListTmp = c(doErrorFoldOutInnerListTmp, currPcPropErrFold) # doErrorFoldOutInnerListTmp[[j]] contains results for the j-th model parameter combination over each of the folds
			innerProgressOut = paste("quantOrSdProp value", k, "out of ", length(quantOrSdPropValuesVect), "END at", date(),"\n")
			cat(innerProgressOut)
			write.table(innerProgressOut, quote=F, row.names=F, col.names=F, append=T, file=progressFile)
		
		}
		doErrorFoldOutInnerList[[i]] = DoExtractErrMat(doErrorFoldOut = doErrorFoldOutInnerListTmp) # these are cross validation results for the i-th training data set (divided into numFolds training/validation); doErrorFoldOutInnerList[[i]][[j]]$ElementName[[k]] are the results for the j-th model parameter combination, in the k-th training/validation inner partition of the i-th outer training/testing fold, with all error measures; ElementName is one of "avgError"  "errorList" "errorMat"  "model"     "parList"
		doErrorFoldOutInnerList<<-doErrorFoldOutInnerList # adding this here, so partial results can still be viewed
		innerProgressOut = paste("Inner cross-validation", i, "out of ", length(doDataSplitOutOuter), "DONE at", date(), "\n")
		cat(innerProgressOut)
		write.table(innerProgressOut, quote=F, row.names=F, col.names=F, append=T, file=progressFile)
		
	}
	doErrorFoldOutInnerList<<-doErrorFoldOutInnerList
}

SetUpTest = function(){
	getTrainResReadyForTest = list()

	for(i in 1:length(doErrorFoldOutInnerList)){
		getTrainResReadyForTest[[i]] = list()
		tmp = DoBestParList(doExtractErrMatOut=doErrorFoldOutInnerList[[i]], errMeasure=errMeasure) # best parameter combination and corresponding error for a given error measure
		getTrainResReadyForTest[[i]]$bestParList = tmp$bestParList
		getTrainResReadyForTest[[i]]$bestAvgErr = tmp$bestParListAvgError
	
		 # best to keep this outside of the former loop, so we can't change errMeasure without having to recompute everything else, which should be fixed
	
		currTrain = doDataSplitOutOuter[[i]]$inDat$missing
		currTest = doDataSplitOutOuter[[i]]$outDat$missing
		missingTestDataLogical = apply(currTest, 1, function(x)any(is.na(x)))
	
		currTrainImput = doDataSplitOutOuter[[i]]$inDat$imputed
		currTestImput = doDataSplitOutOuter[[i]]$outDat$imputed
	
	
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
	modelList = list()

	for(i in 1:length(getTrainResReadyForTest)){
		outerProgressOut = paste("Outer testing", i, "out of ", length(getTrainResReadyForTest), "START at", date(), "\n")
		cat(outerProgressOut)
		write.table(outerProgressOut, quote=F, row.names=F, col.names=F, append=T, file=progressFile)
		 # best to keep this outside of the former loop, so we can't change 'what' without having to recompute everything else, which should be fixed
	 
		bestParList = getTrainResReadyForTest[[i]]$bestParList
		currTrain = getTrainResReadyForTest[[i]]$currTrain
		currTest = getTrainResReadyForTest[[i]]$currTest
		missingTestDataLogical = getTrainResReadyForTest[[i]]$missingTestDataLogical
	
		currTrainImput = getTrainResReadyForTest[[i]]$currTrainImput
		currTestImput = getTrainResReadyForTest[[i]]$currTestImput
		

		currTrainPolyList = DoPolyList(missDat=currTrain, imputDatList = currTrainImput$imputDatList,
			medianOrientedOrNonOrientedImputDat=currTrainImput$medianOrientedBoxImputDat, quantOrSdProp=bestParList$quantOrSdProp,  
			scaleData=scaleData, maxUncertainDims=maxUncertainDims, doMedian=approach=="doMedian", doNoMiss=approach=="doNoMiss",
			doSquarebbSd=approach=="doSquarebbSd", doSquarebbQuant=approach=="doSquarebbQuant")	
		
		if(approach %in% c("doSquarebbQuant", "doSquarebbSd")){
			stopifnot(bestParList$twoSlacks==T)
			testMeanOrMedian = currTestImput$meanBoxImputDat
		} else{
			stopifnot(bestParList$twoSlacks==F)
			
			if(approach%in%c("doMedian", "doPCbb", "doNoMiss")){
				# relevant for doNoMiss too because this is for the prediction of testing points, not training the model
				testMeanOrMedian = currTestImput$medianOrientedBoxImputDat
			} else{
				stop("Don't know what to do")
			}
		}
	
		#if(scaleData)
		#	currTestMedianImput = ScaleCenter(currTestMedianImput, currTrainPolyList$scaleInfo$mean, currTrainPolyList$scaleInfo$std)
		# NO! this is taken care of by doErrorList, via attributes of currTrainPolyList
		currModel = DoTrainModel(polyList=currTrainPolyList, parList = bestParList) # for the current error measure, select  parameters
		#	that give the best error measure results in the validation set, train model on train + valid data, and then use same measure
		#	to calculate testing performance
		modelList[[i]] = currModel
		
		currError = DoErrorList(doTrainModelOut=currModel, medianOrMeanImputOut=testMeanOrMedian, 
			doPolyListOut=currTrainPolyList, missingDatOutLogical = missingTestDataLogical)[[errMeasure]]
		errVect = c(errVect, currError)
		outerProgressOut = paste("Outer testing", i, "out of ", length(getTrainResReadyForTest), "DONE at", date(), "\n")
		cat(outerProgressOut)
		write.table(outerProgressOut, quote=F, row.names=F, col.names=F, append=T, file=progressFile)
	}
	errVect <<-errVect
	modelList <<- modelList
	
}

CollectTestPerformance=function(errMes, appr){
	testPerf = NULL
	for(i in 1:length(testRes)){
		testPerf[i] = testRes[[i]][[errMes]][[appr]]
	}
	res = AggregateTestError(testPerf)
	return(res)
}