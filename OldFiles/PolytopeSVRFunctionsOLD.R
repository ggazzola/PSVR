DoCovMat = function(p, corVal, stdVect){
	# Returns a covariance matrix, with non-diagonal elements
	# derived from correlation value corVal and all elements
	# rescaled via stdVect
	stopifnot(p>1)
	stopifnot(corVal>=0 & corVal<=1)
	stopifnot(length(stdVect)==p)
	corMat  = diag(1, p, p)
	for(i in 1:(p-1))
		for(j in (i+1):p)
			corMat[i,j]=corMat[j,i]= corVal

	covMat = CorToCovMat(corMat, stdVect)
	return(covMat)
}	

DoDat = function(n, meanVect, covMat, trueW, trueW0, noiseLevel){
	#Generates multivariate normal x and corresponding y noisy linear combination
	
	#require(MASS)
	
	stopifnot(length(meanVect) == nrow(covMat))
	stopifnot(nrow(covMat)==ncol(covMat))
	stopifnot(length(meanVect)==length(trueW))
	stopifnot(length(trueW0)==1)
	stopifnot(length(noiseLevel)==1)
	
	stopifnot(n>1)
	
	datX=mvrnorm(n, meanVect, covMat)
	datY= trueW0+datX%*%trueW+noiseLevel*rnorm(n)
	dat = cbind(datX, datY)
	colnames(dat)=c(paste("X", 1:length(meanVect), sep=""), "Y")
	return(dat)
}

#doDatOut = DoDat(n=100, meanVect=c(0,1), covMat=matrix(c(.5, 1, .5, 1), nrow=2, byrow=T), trueW=1:2, trueW0=3, noiseLevel=1)

DoScale = function(dat){
	dat = scale(dat)  
	datMean =  as.numeric(attr(dat, "scaled:center"))
	datSd =  as.numeric(attr(dat, "scaled:scale"))
	#outDat = ScaleCenter(outDat, inDatMean, inDatSd)  	

	return(list(mean=datMean, std=datSd))
}

DoOriginalMissingDataSplit = function(doDatOut, doMissOut, numFolds){
	# Splits 'original' data set and same after application of DoMiss,
	# into train/test(validation) partitions
	# The inner CreateFolds function does stratified (on Y) partitioning
	# so dat should not have missing data along Y
	#source("MiscFunctions.R")
	stopifnot(is.matrix(doDatOut) | is.data.frame(doDatOut))
	stopifnot(is.matrix(doMissOut) | is.data.frame(doMissOut))
	
	stopifnot(length(numFolds)==1)
	stopifnot("Y"==colnames(doDatOut)[ncol(doDatOut)])
	stopifnot(all(!is.na(doDatOut[, ncol(doDatOut)])))
	stopifnot("Y"==colnames(doMissOut)[ncol(doMissOut)])
	stopifnot(all(dim(doDatOut)==dim(doMissOut)))
	
	inOutFolds = CreateFolds(doDatOut[,ncol(doDatOut)], numFolds)

	inOutFoldsList = list()
	
	for(i in 1:numFolds){
		inIdx =  inOutFolds$train[[i]]
		outIdx = inOutFolds$test[[i]]
		inDatOriginal = doDatOut[inIdx,]
		outDatOriginal = doDatOut[outIdx,]
		inDatMissing = doMissOut[inIdx, ]
		outDatMissing = doMissOut[outIdx, ]
		
		#if(scaleData){
		#	inDat = scale(inDat)  
		#	inDatMean =  as.numeric(attr(inDat, "scaled:center"))
		#	inDatSd =  as.numeric(attr(inDat, "scaled:scale"))
		#	outDat = ScaleCenter(outDat, inDatMean, inDatSd)  	
		#	attr(inDat, "scaled:center") = NULL
		#	attr(inDat, "scaled:scale") = NULL
		#}
		inOutFoldsList[[i]] = list()
		inOutFoldsList[[i]]$inDat = list(original = inDatOriginal, missing = inDatMissing)
		inOutFoldsList[[i]]$outDat = list(original = outDatOriginal, missing = outDatMissing)
	}
	
	return(inOutFoldsList)
}

DoDataSplit = function(doDatOut, numFolds, scaleData){
	# Splits data set into train/test(validation) partitions
	# The inner CreateFolds function does stratified (on Y) partitioning
	# so dat should not have missing data along Y
	#source("MiscFunctions.R")
	stopifnot(is.matrix(doDatOut) | is.data.frame(doDatOut))
	stopifnot(length(numFolds)==1)
	stopifnot(is.logical(scaleData))
	stopifnot("Y"==colnames(doDatOut)[ncol(doDatOut)])
	stopifnot(all(!is.na(doDatOut[, ncol(doDatOut)])))
	inOutFolds = CreateFolds(doDatOut[,ncol(doDatOut)], numFolds)

	inOutFoldsList = list()
	
	for(i in 1:numFolds){
		inIdx =  inOutFolds$train[[i]]
		outIdx = inOutFolds$test[[i]]
		inDat = doDatOut[inIdx,]
		outDat = doDatOut[outIdx,]
		if(scaleData){
			inDat = scale(inDat)  
			inDatMean =  as.numeric(attr(inDat, "scaled:center"))
			inDatSd =  as.numeric(attr(inDat, "scaled:scale"))
			outDat = ScaleCenter(outDat, inDatMean, inDatSd)  	
			attr(inDat, "scaled:center") = NULL
			attr(inDat, "scaled:scale") = NULL
		}
		inOutFoldsList[[i]] = list()
		inOutFoldsList[[i]]$inDat = list(original = inDat)
		inOutFoldsList[[i]]$outDat = list(original = outDat)
	}
	
	return(inOutFoldsList)
}

#doDataSplitOut=DoDataSplit(doDatOut,numFolds=ncol(doDatOut),scaleData=T)

DoMiss = function(dat, missingY, missingVarProb){
	#Injects missing data totally at random
	#MissingProb is the proportion of variables with missing values for any given observation
	#If missingY is true, then the Y variable can have missing values as well
	
	stopifnot(is.matrix(dat) | is.data.frame(dat))
	stopifnot(is.logical(missingY))
	stopifnot(is.numeric(missingVarProb) & missingVarProb>=0 & missingVarProb<1)
	stopifnot("Y"%in%colnames(dat))

	p = ncol(dat)
	n = nrow(dat)
	totPotentialMissVar = ifelse(missingY, p, p-1)

	for(i in 1:n){
		dat[i, sample(totPotentialMissVar, round(totPotentialMissVar*missingVarProb))] = NA # each row independently of the other has a certain proportion of missing columns
	}
	return(dat)	
}

DoMissFold = function(doDataSplitOut, missingY, missingVarProb){
	doDataSplitDoMissOut = doDataSplitOut
	for(i in 1:length(doDataSplitOut)){
		stopifnot(is.list(doDataSplitDoMissOut[[i]]$inDat))
		stopifnot(is.list(doDataSplitDoMissOut[[i]]$outDat))
		stopifnot(is.matrix(doDataSplitDoMissOut[[i]]$inDat$original) | 
			is.data.frame(doDataSplitDoMissOut[[i]]$inDat$original))
		stopifnot(is.matrix(doDataSplitDoMissOut[[i]]$outDat$original) | 
			is.data.frame(doDataSplitDoMissOut[[i]]$outDat$original))
		
	#  note this applies the missing process on training/testing independently 
		doDataSplitDoMissOut[[i]]$inDat$missing = DoMiss(doDataSplitOut[[i]]$inDat$original, missingY, missingVarProb)
		doDataSplitDoMissOut[[i]]$outDat$missing = DoMiss(doDataSplitOut[[i]]$outDat$original, missingY, missingVarProb) 
	}
	return(doDataSplitDoMissOut)
}

#doMissFoldOut = DoMissFold(doDataSplitOut, missingY=T, missingVarProb = .3)

DoMultipleImputation = function(missDat, method, numImput, maxIter, extraPossiblyMissDat=NULL){
	# Returns list of multiple-imputation matrices associated to a matrix with missing data
	# as well as the median of such multiple imputations; if the matrix has no missing data
	# then the multiple imputations will be simply copies of that matrix
	
	# missDat is the matrix with missing data to be imputed;
	# extraPossiblyMissDat is an optional extra data set to use as a base for the imputations
	# e.g., training data if missDat is testing data. extraPossiblyMissDat should NOT already 
	# include imputations --> maybe consider doing imputation of missDat one data point at a time + extraPossiblyMissDat
	# if missDat is a test set
	
	stopifnot(is.matrix(missDat) | is.data.frame(missDat))
	if(!is.null(extraPossiblyMissDat)){
		stopifnot(is.matrix(extraPossiblyMissDat) | is.data.frame(extraPossiblyMissDat))
		stopifnot(ncol(extraPossiblyMissDat)==ncol(missDat))
	}	

	anyMissingInMissDat = any(is.na(missDat))
	if(anyMissingInMissDat){
		stopifnot(numImput>=ncol(missDat)) # else can't run pca on it 
		stopifnot(maxIter>1)
		stopifnot(is.character(method))
		n = nrow(missDat)
		# there is no point in doing this if there are no missing values
		missDat = rbind(missDat, extraPossiblyMissDat)
		multiImputDat = mice(missDat, method = method, m=numImput, maxit=maxIter, print=F) # multiple imputations
	} else{
		if(!is.null(extraPossiblyMissDat)){
			stop("extraPossiblyMissDat should be NULL if missDat has no missing values.\n")
		}
	}
	imputDatList = list()
	for(j in 1:numImput){
		
		if(anyMissingInMissDat){
			tmpImpMat = as.matrix(complete(multiImputDat, j)[1:n,])
			if(any(is.na(tmpImpMat))){
				stop("Imputation didn't work along all variables (collinearity...)\n")
			}
		} else{
			tmpImpMat = as.matrix(missDat)
		}	
		imputDatList[[j]] = tmpImpMat
		row.names(imputDatList[[j]]) = NULL
	}
	medianImputDat = apply(simplify2array(imputDatList), 1:2, median) #1:2 --> rows and columns
	rownames(medianImputDat) = NULL
	res = list(imputDatList=imputDatList, medianImputDat = medianImputDat)
	return(res)
}


DoMultipleImputationFolds = function(doMissOut, method, maxIter, numImput){
	#returns list of multiple imputations and corresponding median imputation
	# of each training/testing partition of data, assuming it includes missing data
	
	#require(mice)
	stopifnot(is.character(method))
	stopifnot(maxIter>1)
	stopifnot(numImput>=ncol(doMissOut[[1]]$inDat$missing))
	stopifnot(is.list(doMissOut))

	resImputeFolds = doMissOut
	for(i in 1:length(doMissOut)){
		stopifnot("inDat"%in%names(doMissOut[[i]]))
		stopifnot("outDat"%in%names(doMissOut[[i]]))
		stopifnot("Y"%in%colnames(doMissOut[[i]]$inDat$missing))
		stopifnot("Y"%in%colnames(doMissOut[[i]]$outDat$missing))
		
		inDat = doMissOut[[i]]$inDat$missing
		outDat = doMissOut[[i]]$outDat$missing
		resImputeFolds[[i]]$inDat$imputed = DoMultipleImputation(inDat, method, numImput, maxIter)
		resImputeFolds[[i]]$outDat$imputed = DoMultipleImputation(outDat, method, numImput, maxIter, inDat) # assuming test data given in batch, not one point at a time... 
	}
	attr(resImputeFolds, "method")=method
	attr(resImputeFolds, "maxIter")=maxIter
	return(resImputeFolds)
}


#	doMultipleImputationFoldsOut = DoMultipleImputationFolds(doMissFoldOut, method="pmm", maxIter=2, numImput=4)

DoPolyList = function(missDat, imputDatList, medianImputDat, pcProp, scaleData, maxUncertainDims){
	# calculates Pc bounding box, and the corresponding uncertainty value, around data with 
	# missing values defined by missDat
	# and its corresponding multiple imputation matrix list
	# pcProp is the proportion of data to consider around the median in the 
	# pc space
	# scaleData, if TRUE, all data is scaled to mean and std of meanImputDat
	#source("ConstraintFunctions.R")
	stopifnot(is.matrix(medianImputDat) | is.data.frame(medianImputDat))
	stopifnot(is.matrix(missDat) | is.data.frame(missDat))
	
	stopifnot(all(dim(medianImputDat)==dim(missDat)))
	stopifnot(all(dim(medianImputDat)==dim(imputDatList[[1]])))
	stopifnot("Y"%in%colnames(missDat))
	stopifnot("Y"%in%colnames(medianImputDat))
	stopifnot(is.list(imputDatList))
	stopifnot(length(imputDatList)>1)
	stopifnot(is.numeric(pcProp))
	stopifnot(pcProp>0)
	stopifnot(pcProp<=1)
	stopifnot(is.matrix(missDat)| is.data.frame(missDat))
	stopifnot(is.logical(scaleData))
	#stopifnot(any(is.na(missDat)))
	
	polyList=list()
	overallProjPc = list()
	projDimHash = NULL
	n = nrow(missDat)
	pPlusOne = ncol(missDat)
	
	if(scaleData){
		scaleInfo  = DoScale(medianImputDat)
		medianImputDat = ScaleCenter(medianImputDat, scaleInfo$mean, scaleInfo$std)
		for(i in 1:length(imputDatList)){
			imputDatList[[i]] = ScaleCenter(imputDatList[[i]], scaleInfo$mean, scaleInfo$std)
		}
	} 
	for(i in 1:n){
		currPoint = missDat[i,]
		if(any(is.na(currPoint))){
			singleMissDatMultImput = NULL
			for(j in 1:length(imputDatList)){
				singleMissDatMultImput = rbind(singleMissDatMultImput, imputDatList[[j]][i,]) # stacking multiple imputations for i-th point
			}
			singleMissDatMultImput = as.matrix(singleMissDatMultImput)
			###polyBuiltFrom = singleMissDatMultImput
			singleDatRes = PCConstraintsFull(singleMissDatMultImput, pcProp) # extracting bounding box for i-th point's multiple imputations
			singleDatResDimHash = PasteFast(singleDatRes$volumeDims) # this string contains the indices of the uncertain dimensions
			if(!singleDatResDimHash%in%projDimHash){
				projDimHash = c(projDimHash, singleDatResDimHash) # adding the string if not added before to "database" projDimHash
				overallProjPc[[singleDatResDimHash]] = PCConstraintsFull(medianImputDat, propIn=1, projectDimsLogical = (1:pPlusOne)%in%singleDatRes$volumeDims)
				# calculating pc bounding box on subspace defined by redDimHash (doing this only once per subset -- <= once per point)
				# considering the median imputation of all the data; this bounding box corresponds to the 'overall uncertainty', to compare the the 
				# individual point uncertainty
				# note we set propIn = 1, even if imputDatList may have been calculated with propIn<1 (so that acting on the latter,
				# actually has an effect in determining currUncertainty)
			} 
			currUncertainty = PcBBUncertainty(overallProjPc[[singleDatResDimHash]], singleDatRes, missDat, maxUncertainDims) 
			# a [0,1] uncertainty level -- note 'ALL'
			
		} else{
			if(scaleData)
				currPoint = ScaleCenter(currPoint, scaleInfo$mean, scaleInfo$std)
			
			###polyBuiltFrom=currPoint
			singleDatRes = FixedPointConstraints(currPoint) # fixed-point bounding box for fully known points
			currUncertainty = 0
		}
		tmp=PolytopeDescription(singleDatRes)
		tmp$uncertainty = currUncertainty
		###tmp$polyBuiltFrom = polyBuiltFrom
		polyList[[i]] = tmp
	}
	attr(polyList, "pcProp") = pcProp
	
	if(scaleData){
		attr(polyList, "scaled") = T
		attr(polyList, "scaleInfo") = scaleInfo
	} else{
		attr(polyList, "scaled") = F
	}
	return(polyList)
}

DoPolyListFolds = function(doMultipleImputationFoldsOut, pcProp, scaleData,maxUncertainDims){
	# returns polyList object (polytope representation of data) for each data point
	# in each inDat component of doMultipleImputationFoldsOut,
	polyListFolds=list()
	for(i in 1:length(doMultipleImputationFoldsOut)){
		polyListFolds[[i]] = list()
		stopifnot(!is.null(doMultipleImputationFoldsOut[[i]]$inDat$missing))
		stopifnot(!is.null(doMultipleImputationFoldsOut[[i]]$inDat$imputed$imputDatList))
		
		polyListFolds[[i]]$polyList = DoPolyList(doMultipleImputationFoldsOut[[i]]$inDat$missing, 
			doMultipleImputationFoldsOut[[i]]$inDat$imputed$imputDatList, 
			doMultipleImputationFoldsOut[[i]]$inDat$imputed$medianImputDat, 
			pcProp, scaleData, maxUncertainDims)
		polyListFolds[[i]]$inDat= doMultipleImputationFoldsOut[[i]]$inDat
		polyListFolds[[i]]$outDat= doMultipleImputationFoldsOut[[i]]$outDat		
	}
	return(polyListFolds)	
}

#doPolyListFoldsOut = DoPolyListFolds(doMultipleImputationFoldsOut, pcProp=0.8)


DoParListGrid = function(parValuesList){
	# calculates a grid of parameter-value combinations and stores it into a list,
	# each element of which corresponds to a different combination
	# parValuesList should be a list with parameter names and corresponding parameter
	# values to try; Parameter names should be all specified, even if a parameter
	# to be kept fixed
	stopifnot(is.list(parValuesList))
	requiredParams = names(formals(PolytopeSVR))
	requiredParams = requiredParams[!requiredParams%in%c("polyList", "...")]
	parNames = names(parValuesList)
	
	stopifnot(all(sort(parNames)==sort(requiredParams)))
	
	parGrid = expand.grid(parValuesList)
	parListGrid = list()
	for(i in 1:nrow(parGrid)){
		parListGrid[[i]] = list()
		for(j in 1:ncol(parGrid)){
			parListGrid[[i]][[parNames[j]]] = parGrid[i,j]
		}
	}
	
	return(parListGrid)
}

#parValuesList = list(Ccertain=c(1,5), Cuncertain=c(1,5), epsilonCertain=c(0,1), 
#	extraEpsilonUncertain = c(0,1), uncertaintySpecialTreatment = T)
# doParListGridOut = DoParListGrid(parValuesList)

DoTrainModel = function(polyList, parList){
	# polyList: polytope representation of a data set
	# parList: list of parameters to feed to PolytopeSVR
	# returns list with components w and w0  resulting from model training
	parList$polyList = polyList
	model = do.call(PolytopeSVR, parList)
	res = gurobi(model)
	w = GetSolution(model, res, "w")
	w0 = GetSolution(model, res, "w0")
	res = list(w=w, w0=w0)
	class(res) = "PolytopeSVRWW0"
	return(res)
}

DoError = function(predY, trueY, missingDatOutLogical){
	# returns a list with different prediction error measures, 
	# given vector of predictions and vector of true values
	errList = list()
	errList$mae = MAE(predY, trueY)
	errList$rmse = RMSE(predY, trueY)
	errList$Maxae = MaxAE(predY, trueY)
	errList$maeUncert = MAE(predY, trueY, missingDatOutLogical)
	errList$rmseUncert = RMSE(predY, trueY, missingDatOutLogical)
	errList$MaxaeUncert = MaxAE(predY, trueY, missingDatOutLogical)
	errList$maeCert = MAE(predY, trueY, !missingDatOutLogical)
	errList$rmseCert = RMSE(predY, trueY, !missingDatOutLogical)
	errList$MaxaeCert = MaxAE(predY, trueY, !missingDatOutLogical)
	return(errList)
}

DoErrorList = function(doTrainModelOut, medianImputOut, doPolyListOut, missingDatOutLogical){
	# returns list with different out-of-sample prediction error measures 
	# using model doTrainModelOut and ***out-of sample*** data medianImputOut
	# which is supposed to come from data with uncertainties, if any, only 
	# on the x variables!!! (otherwise can't compute errors correctly)
	#doPolyListOut is the polyList corresponding to doTrainModelOut; used 
	# for scaling medianImputOut
	# missingDatOutLogical: logical, True for rows of medianImputOut corresponding to
	# observations with missing values
	cat("Are you sure medianImput does not come from uncertain-Y out-of-sample data?\n")
	stopifnot(class(doTrainModelOut)=="PolytopeSVRWW0")
	stopifnot(is.list(doPolyListOut))
	stopifnot(is.logical(missingDatOutLogical))
	stopifnot(length(missingDatOutLogical)==nrow(medianImputOut))
	
	w= doTrainModelOut$w
	w0 = doTrainModelOut$w0
	# medianImputOut: median imputation of testing or validation data
	stopifnot("Y"%in%colnames(medianImputOut))
	p = ncol(medianImputOut)-1
	
	if(attr(doPolyListOut, "scaled")){
		scaleInfo = attr(doPolyListOut, "scaleInfo")
		medianImputOut = ScaleCenter(medianImputOut, scaleInfo$mean, scaleInfo$std)
	}	
	
	predY = w0+medianImputOut[, 1:p]%*%w
	trueY =  medianImputOut[,p+1]
	errList  = DoError(predY, trueY, missingDatOutLogical)
	return(errList)
}

DoErrorFold = function(doPolyListFoldsOut, doParListGridOut, replaceImputedWithTrueY = F){
	stopifnot(is.list(doParListGridOut))
	stopifnot(is.list(doPolyListFoldsOut))
	stopifnot(is.logical(replaceImputedWithTrueY))
	
	# for each training/testing fold, and for each parameter combination, train model on in-sample,
	# calculate its error measure list, using the corresponding imputed out-of-sample. 
	# polyListWasScaled
	# Use replaceImputedWithTrueY to replace uncertain Y values in the imputed out-of-sample with true
	# values
	parListRes  = list()
	for(j in 1:length(doParListGridOut)){
		parListRes[[j]] = list()
		currParList = parListRes[[j]]$parList = doParListGridOut[[j]]
		parListRes[[j]]$errorList = list()
		parListRes[[j]]$model = list() 
		
		for(i in 1:length(doPolyListFoldsOut)){
			currFoldPolyList = doPolyListFoldsOut[[i]]$polyList
			
			
			parListRes[[j]]$model[[i]] = DoTrainModel(currFoldPolyList, currParList)
			currFoldMissingDatOut = doPolyListFoldsOut[[i]]$outDat$missing
			currFoldImputedOutDat = doPolyListFoldsOut[[i]]$outDat$imputed$medianImputDat
			currFoldMissingOutDatLogical = apply(currFoldMissingDatOut, 1, function(x){any(is.na(x))})  # T for out-of-sample points with missing values
			stopifnot("Y"==colnames(currFoldImputedOutDat)[ncol(currFoldImputedOutDat)])
			outDatMissingY = any(is.na(currFoldMissingDatOut[,ncol(currFoldMissingDatOut)]))
			if(outDatMissingY){
				stop("No missing Y allowed")
				# NOTE: if replacing, we must rescale!!!!!!!!!!
				#if(replaceImputedWithTrueY){
				#	cat("Warning: replacing imputed Y with true Y.\n")
				#	currFoldImputedOutDat[, ncol(currFoldImputedOutDat)] = doPolyListFoldsOut[[i]]$outDat$original[, ncol(currFoldImputedOutDat)]
				#} else{
				#	stop("Can't compute performance with uncertain Y in out-of-sample data\n")
				#}
			}
			parListRes[[j]]$errorList[[i]] = DoErrorList(parListRes[[j]]$model[[i]], currFoldImputedOutDat, 
				currFoldPolyList, currFoldMissingOutDatLogical)
		}		
	}
	#parListRes[[j]]errorList[[i]]: errorList of parameter combination j on fold i
	return(parListRes)	
}


DoExtractErrMat = function(doErrorFoldOut, average=mean){
	stopifnot(is.list(doErrorFoldOut))
	stopifnot(is.function(average))
	errMatRes = doErrorFoldOut
	errMeasureVect = names(errMatRes[[1]]$errorList[[1]])
	for(i in  1:length(errMatRes)){
		# for each parameter combination
		errMat = NULL
		for (errMeasure in errMeasureVect){
			#for each error measure
			errVect = NULL
			for(j in 1:length(errMatRes[[i]])){
				# for each fold
				currErr = errMatRes[[i]]$errorList[[j]][[errMeasure]]
				stopifnot(!is.null(currErr))
				errVect = c(errVect, currErr)
			}
			errMat = cbind(errMat, errVect) # each row is a different parameter combination, each column is a fold
		}
		colnames(errMat)=errMeasureVect
		errMat=as.data.frame(errMat)
		errMatRes[[i]]$errorMat=errMat
		avgErr = apply(errMat,2, average)
		names(avgErr)= colnames(errMat)
		errMatRes[[i]]$avgError = data.frame(t(avgErr))
		attr(errMatRes[[i]]$avgError, "average")=average
	}
	return(errMatRes)
	
}

DoBestParList = function(doExtractErrMatOut, errMeasure){
	stopifnot(is.list(doExtractErrMatOut))
	stopifnot(is.character(errMeasure))
	
	averageErrVect = NULL
	for(i in 1:length(doExtractErrMatOut)){
		currAvgError = doExtractErrMatOut[[i]]$avgError
		if(is.null(currAvgError[[errMeasure]])){
			stop("Error measure '", errMeasure, "' does not exist in doExtractErrMatOut.\n")
		}
		averageErrVect[i] = currAvgError[[errMeasure]]
	}

	minAverageErrIdx =Which.min.tie(averageErrVect)
	bestParList = doExtractErrMatOut[[minAverageErrIdx]]$parList
	bestAvgError = averageErrVect[minAverageErrIdx]
	if(bestAvgError==Inf)
		stop("Error measure '", errMeasure, "' can't be calculated because there is no data for it.\n")
	res = list(errorMeasure = errMeasure, bestParList = bestParList, bestParListAvgError = bestAvgError,
		 bestParListIdx = minAverageErrIdx)
	return(res)
}		

