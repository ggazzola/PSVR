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

DoDat = function(n, meanVect, covMat, trueW, trueW0, theoRsq){
	#Generates multivariate normal x and corresponding y noisy linear combination
	
	#require(MASS)
	
	stopifnot(length(meanVect) == nrow(covMat))
	stopifnot(nrow(covMat)==ncol(covMat))
	stopifnot(length(meanVect)==length(trueW))
	stopifnot(length(trueW0)==1)
	stopifnot(length(theoRsq)==1)
	
	stopifnot(n>1)
	
	datX=mvrnorm(n, meanVect, covMat)
	
	noiseLevel = NoiseLevel(trueW, covMat, theoRsq)
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
		inDatOriginal = VectToMat(doDatOut[inIdx,])
		outDatOriginal = VectToMat(doDatOut[outIdx,])
		inDatMissing = VectToMat(doMissOut[inIdx, ])
		outDatMissing = VectToMat(doMissOut[outIdx, ])
		
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

DoMiss = function(dat, missingY, missingObsProp, missingVarProp){
	#Injects missing data totally at random
	#MissingProb is the proportion of variables with missing values for any given observation
	#If missingY is true, then the Y variable can have missing values as well
	
	stopifnot(is.matrix(dat) | is.data.frame(dat))
	stopifnot(is.logical(missingY))
	stopifnot(is.numeric(missingVarProp) & missingVarProp>=0 & ifelse(missingY, missingVarProp<1, missingVarProp<=1))
	stopifnot(is.numeric(missingObsProp) & missingObsProp>=0 & missingObsProp<1)
	
	stopifnot("Y"%in%colnames(dat))
	p = ncol(dat)
	n = nrow(dat)
	totPotentialMissVar = ifelse(missingY, p, p-1)
	numMissVar = round(totPotentialMissVar*missingVarProp)
	numMissObs = round(n*missingObsProp)
	stopifnot(numMissVar<totPotentialMissVar)
	stopifnot(numMissObs<n)
	
	whichMissObs = sample(n, numMissObs)
	for(i in whichMissObs){
		dat[i, sample(totPotentialMissVar, numMissVar)] = NA # each row independently of the other has a certain proportion of missing columns
	}
	return(dat)	
}

DoMissFold = function(doDataSplitOut, missingY, missingObsProp, missingVarProp){
	doDataSplitDoMissOut = doDataSplitOut
	for(i in 1:length(doDataSplitOut)){
		stopifnot(is.list(doDataSplitDoMissOut[[i]]$inDat))
		stopifnot(is.list(doDataSplitDoMissOut[[i]]$outDat))
		stopifnot(is.matrix(doDataSplitDoMissOut[[i]]$inDat$original) | 
			is.data.frame(doDataSplitDoMissOut[[i]]$inDat$original))
		stopifnot(is.matrix(doDataSplitDoMissOut[[i]]$outDat$original) | 
			is.data.frame(doDataSplitDoMissOut[[i]]$outDat$original))
		
	#  note this applies the missing process on training/testing independently 
		doDataSplitDoMissOut[[i]]$inDat$missing = DoMiss(doDataSplitOut[[i]]$inDat$original, missingY, missingObsProp, missingVarProp)
		doDataSplitDoMissOut[[i]]$outDat$missing = DoMiss(doDataSplitOut[[i]]$outDat$original, missingY, missingObsProp, missingVarProp) 
	}
	return(doDataSplitDoMissOut)
}

#doMissFoldOut = DoMissFold(doDataSplitOut, missingY=T, missingVarProp = .3)

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

	stopifnot(numImput>=ncol(missDat)) # else can't run pca on it 
	stopifnot(maxIter>1)
	stopifnot(is.character(method))
	n = nrow(missDat)
	# there is no point in doing this if there are no missing values
	missDat = rbind(missDat, extraPossiblyMissDat)

	imputDatList = list()
	
	if(!any(is.na(missDat))){
		cat("Warning: no missing data. Using multiple copies of the data as 'multiple imputations'\n")
		for(j in 1:numImput){
			imputDatList[[j]] = as.matrix(missDat[1:n,]) # the [1:n,] is because of the rbind to missDat above
		}
	} else{
		
		multiImputDat = mice(missDat, method = method, m=numImput, maxit=maxIter, print=F) # multiple imputations
		
		for(j in 1:numImput){
			tmpImpMat = as.matrix(complete(multiImputDat, j)[1:n,])
			if(any(is.na(tmpImpMat))){
				stop("Imputation didn't work along all variables (collinearity...)\n")
			}
			imputDatList[[j]] = tmpImpMat
			row.names(imputDatList[[j]]) = NULL
		}
	}
	meanImput = apply(missDat, 2, mean, na.rm=T)
	meanBoxImputDat = missDat[1:n,]
	for(i in 1:n){
		missCol = which(is.na(meanBoxImputDat[i,]))
		meanBoxImputDat[i,missCol] = meanImput[missCol]
	}
	
	medianImputDat = apply(simplify2array(imputDatList), 1:2, median) #1:2 --> rows and columns
	rownames(medianImputDat) = NULL
	
	res = list(imputDatList=imputDatList, medianImputDat = medianImputDat, meanBoxImputDat=meanBoxImputDat)
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

DoPolyList = function(missDat, imputDatList, medianImputDat, quantOrSdProp, scaleData, maxUncertainDims, 
	doMedian, doNoMiss, doSquarebbSd, doSquarebbQuant){
	# calculates Pc bounding box, and the corresponding uncertainty value, around data with 
	# missing values defined by missDat
	# and its corresponding multiple imputation matrix list
	# quantOrSdProp is the proportion of data to consider around the median in the 
	# pc space
	# scaleData, if TRUE, all data is scaled to mean and std of meanImputDat
	# if missDat has no missing points, then imputDatList is irrelevant
	#source("ConstraintFunctions.R")
	#doMedian: if TRUE, then the polytope representation will be based on the median imputation of the missing
	# data, with no uncertainty around it
	#doNoMiss: if TRUE, then the (degenerate) polytope representation will be based on the only fully-available data points;
	#	this approach is therefore equivalent to complete-case analysis SVR
	stopifnot(is.matrix(medianImputDat) | is.data.frame(medianImputDat))
	stopifnot(is.matrix(missDat) | is.data.frame(missDat))
	stopifnot(is.logical(doMedian)& is.logical(doNoMiss) &is.logical(doSquarebbSd) & is.logical(doSquarebbQuant))
	
	sumAltMethods = doMedian+ doNoMiss+ doSquarebbSd +doSquarebbQuant
	doOrientedbb = sumAltMethods==0
	
	if(sumAltMethods>1)
		stop("Only one between doMedian,  doMiss, doSquarebbSd, and doSquarebbQuant can be TRUE\n")
	
	stopifnot(all(dim(medianImputDat)==dim(missDat)))
	stopifnot("Y"%in%colnames(missDat))
	stopifnot("Y"%in%colnames(medianImputDat))
	
	if(doMedian){
		missDat = medianImputDat
		imputDatList = quantOrSdProp = "irrelevantForDoMedian"
	} else if(doNoMiss){
		noMissIdx = which(apply(missDat, 1, function(x){all(!is.na(x))}))
		if(length(noMissIdx)==0)
			stop("No complete observations available in the data\n")
		cat(round(length(noMissIdx)/nrow(missDat)*100), "% of the data is complete for doNoMiss\n")
		missDat = VectToMat(missDat[noMissIdx,])
		medianImputDat = missDat
		imputDatList = quantOrSdProp = "irrelevantForDoMiss"
	} else if(doOrientedbb){
		missDatHasNAs = any(is.na(missDat))
		if(missDatHasNAs){
			stopifnot(is.list(imputDatList))
			stopifnot(length(imputDatList)>1)
			stopifnot(all(dim(medianImputDat)==dim(imputDatList[[1]])))
		}
		overallProjPc = list()
		projDimHash = NULL
	}
	
	stopifnot(is.matrix(missDat)| is.data.frame(missDat))
	stopifnot(is.logical(scaleData))
	if(!doMedian & ! doNoMiss){
		stopifnot(is.numeric(quantOrSdProp))
		stopifnot(quantOrSdProp>=0)
		stopifnot(quantOrSdProp<=1)
	}
	polyList=list()
	n = nrow(missDat)
	pPlusOne = ncol(missDat)
	
	if(scaleData){ 
		if(doSquarebbSd|doSquarebbQuant){
			scaleInfo = DoScale(missDat) # we don't want to do scale on meanImputDat b/c artificially small sd along imputed columns (mean is the same, though)
			stopifnot(all(scaleInfo$std>0))

		} else{
			scaleInfo  = DoScale(medianImputDat)
			# note that this is scaling is different from that of doSquareBB above (here we consider also the imputed values, above we don't)
			if(doOrientedbb){
				medianImputDat = ScaleCenter(medianImputDat, scaleInfo$mean, scaleInfo$std)
				if(missDatHasNAs){
					for(i in 1:length(imputDatList)){
						imputDatList[[i]] = ScaleCenter(imputDatList[[i]], scaleInfo$mean, scaleInfo$std)
					}
				}
			}
		}
	} 
	
	if(doSquarebbQuant){
		lowerQuantile = (1-quantOrSdProp)/2
		upperQuantile = 1-lowerQuantile
		lowerQuantMissDat = apply(missDat, 2, quantile, probs=lowerQuantile, type=1, na.rm=T) 
		upperQuantMissDat = apply(missDat, 2, quantile, probs=upperQuantile, type=1, na.rm=T) 
	}
	for(i in 1:n){
		currPoint = missDat[i,]
		if(any(is.na(currPoint))){
			if(doOrientedbb){
				stopifnot(missDatHasNAs) # just for debugging
				singleMissDatMultImput = NULL
				for(j in 1:length(imputDatList)){
					singleMissDatMultImput = rbind(singleMissDatMultImput, imputDatList[[j]][i,]) # stacking multiple imputations for i-th point
				}
				singleMissDatMultImput = as.matrix(singleMissDatMultImput)
				if(all(apply(singleMissDatMultImput, 2, sd)==0)){
					cat("WARNING: All imputations are identical\n")
					currPointConstantImput = singleMissDatMultImput[1,] # any row is ok, since they are all equal
					singleDatRes = FixedPointConstraints(currPointConstantImput) # fixed-point bounding box around constant imputation
					currUncertainty = 0
				} else{
					singleDatRes = PCConstraintsFull(dat=singleMissDatMultImput, propIn=quantOrSdProp) # extracting bounding box for i-th point's multiple imputations
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
					currUncertainty = PcBBUncertainty(overallProjPc[[singleDatResDimHash]], singleDatRes, missDat, maxUncertainDims, propIn=quantOrSdProp) 
				}	# a [0,1] uncertainty level -- note 'ALL'
			} else if(doSquarebbSd){
				currPointMeanImput = currPoint
				currPointMeanImput[is.na(currPoint)] = scaleInfo$mean[is.na(currPoint)]
				sdVect = rep(0, pPlusOne)
				if(scaleData){
					currPointMeanImput = ScaleCenter(currPointMeanImput, scaleInfo$mean, scaleInfo$std)
					sdVect[is.na(currPoint)] = 1
				} else{
					sdVect[is.na(currPoint)] = scaleInfo$std[is.na(currPoint)]
				}
				singleDatRes = BoxConstraints(currPointMeanImput, sdVect, quantOrSdProp)
				currUncertainty = "irrelevantFordoSquarebbSd"
				
			} else if(doSquarebbQuant){
				currPointLowerQuantileImput= currPointUpperQuantileImput = currPoint
				currPointLowerQuantileImput[is.na(currPoint)] = lowerQuantMissDat[is.na(currPoint)]
				currPointUpperQuantileImput[is.na(currPoint)] = upperQuantMissDat[is.na(currPoint)]
				if(scaleData){
					currPointLowerQuantileImput = ScaleCenter(currPointLowerQuantileImput, scaleInfo$mean, scaleInfo$std)
					currPointUpperQuantileImput = ScaleCenter(currPointUpperQuantileImput, scaleInfo$mean, scaleInfo$std)
					#scaling currPoint is irrelevant for doSquarebbQuant
				}
				singleDatRes = BoxConstraintsFromUpperAndLowerBound(currPointLowerQuantileImput, currPointUpperQuantileImput)
				currUncertainty = "irrelevantFordoSquarebbQuant"
				
			} else{
				stop("Should never get here")
			}
				
		} else{
			if(scaleData)
				currPoint = ScaleCenter(currPoint, scaleInfo$mean, scaleInfo$std)
			
			#for(i in 1:nrow(medianImputDat)){
			#	currRow = medianImputDat[i, ]
			#	missingVarIdx = is.na(currRow)
			#	currRow[missingVarIdx] = missDatMean[missingVarIdx]
			#	medianImputDat[i, ] = currRow
			#} #for squareBB
			
			###polyBuiltFrom=currPoint
			singleDatRes = FixedPointConstraints(currPoint) # fixed-point bounding box for fully known points
			currUncertainty = 0
		}
		tmp=PolytopeDescription(singleDatRes)
		tmp$uncertainty = currUncertainty
		###tmp$polyBuiltFrom = polyBuiltFrom
		polyList[[i]] = tmp
	}
	attr(polyList, "quantOrSdProp") = quantOrSdProp
	
	if(scaleData){
		attr(polyList, "scaled") = T
		attr(polyList, "scaleInfo") = scaleInfo
	} else{
		attr(polyList, "scaled") = F
	}
	return(polyList)
}

DoPolyListFolds = function(doMultipleImputationFoldsOut, quantOrSdProp, scaleData, maxUncertainDims, doMedian, doNoMiss, doSquarebbSd, doSquarebbQuant){
	# returns polyList object (polytope representation of data) for each data point
	# in each inDat component of doMultipleImputationFoldsOut,

	polyListFolds=list()
	for(i in 1:length(doMultipleImputationFoldsOut)){
		polyListFolds[[i]] = list()
		stopifnot(!is.null(doMultipleImputationFoldsOut[[i]]$inDat$missing))
		stopifnot(!is.null(doMultipleImputationFoldsOut[[i]]$inDat$imputed$imputDatList))
		polyListFolds[[i]]$polyList = DoPolyList(missDat=doMultipleImputationFoldsOut[[i]]$inDat$missing, 
			imputDatList=doMultipleImputationFoldsOut[[i]]$inDat$imputed$imputDatList, 
			medianImputDat=doMultipleImputationFoldsOut[[i]]$inDat$imputed$medianImputDat, 
			quantOrSdProp=quantOrSdProp, scaleData=scaleData, maxUncertainDims=maxUncertainDims, doMedian=doMedian, doNoMiss=doNoMiss, 
				doSquarebbSd=doSquarebbSd, doSquarebbQuant=doSquarebbQuant)
	
		polyListFolds[[i]]$inDat= doMultipleImputationFoldsOut[[i]]$inDat
		polyListFolds[[i]]$outDat= doMultipleImputationFoldsOut[[i]]$outDat		
	}
	return(polyListFolds)	
}

#doPolyListFoldsOut = DoPolyListFolds(doMultipleImputationFoldsOut, quantOrSdProp=0.8)


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
	model <<- do.call(PolytopeSVR, parList)
	gurobiParams = list(LogFile="", OutputFlag=0)
	resGurobi <<- gurobi(model, params=gurobiParams)
	w = GetSolution(model, resGurobi, "w")
	w0 = GetSolution(model, resGurobi, "w0")
	u = GetSolution(model, resGurobi, "u")
	v = GetSolution(model, resGurobi, "v")
	if(is.logical(parList$twoSlacks)){
		if(parList$twoSlacks){
			csiPlus = GetSolution(model, resGurobi, "allCsiPlus")
			csiMinus = GetSolution(model, resGurobi, "allCsiMinus")
			res = list(w=w, w0=w0, u=u, v=v, csiPlus=csiPlus, csiMinus=csiMinus)
		} else{
			csi = GetSolution(model, resGurobi, "allCsi")
			res = list(w=w, w0=w0, u=u, v=v, csi=csi)
		}
	}
	class(res) = "PolytopeSVRWW0"
	return(res)
}

DoError = function(predY, trueY, missingDatOutLogical){
	# returns a list with different prediction error measures, 
	# given vector of predictions and vector of true values
	errList = list()
	errList$mae = MAE(trueY, predY)
	errList$rmse = RMSE(trueY, predY)
	errList$Maxae = MaxAE(trueY, predY)
	errList$quantNineAe = QuantNineAE(trueY, predY)
	errList$quantEightAe = QuantEightAE(trueY, predY)
	errList$quantSevenAe = QuantSevenAE(trueY, predY)

	errList$maeUncert = MAE(trueY, predY, missingDatOutLogical)	
	errList$rmseUncert = RMSE(trueY, predY, missingDatOutLogical)
	errList$MaxaeUncert = MaxAE(trueY, predY, missingDatOutLogical)
	errList$quantNineAeUncert = QuantNineAE(trueY, predY, missingDatOutLogical)
	errList$quantEightAeUncert = QuantEightAE(trueY, predY, missingDatOutLogical)
	errList$quantSevenAeUncert = QuantSevenAE(trueY, predY, missingDatOutLogical)
	
	errList$maeCert = MAE(trueY, predY, !missingDatOutLogical)
	errList$rmseCert = RMSE(trueY, predY, !missingDatOutLogical)
	errList$MaxaeCert = MaxAE(trueY, predY, !missingDatOutLogical)
	errList$quantNineAeCert = QuantNineAE(trueY, predY, !missingDatOutLogical)
	errList$quantEightAeCert = QuantEightAE(trueY, predY, !missingDatOutLogical)
	errList$quantSevenAeCert = QuantSevenAE(trueY, predY, !missingDatOutLogical)
	
	errList$predRes = list(predY=predY, trueY=trueY, missingDatOutLogical=missingDatOutLogical)
	return(errList)
}

DoErrorList = function(doTrainModelOut, medianOrMeanImputOut, doPolyListOut, missingDatOutLogical){
	# returns list with different out-of-sample prediction error measures 
	# using model doTrainModelOut and ***out-of sample*** data medianOrMeanImputOut
	# which is supposed to come from data with uncertainties, if any, only 
	# on the x variables!!! (otherwise can't compute errors correctly)
	#doPolyListOut is the polyList corresponding to doTrainModelOut; used 
	# for scaling medianOrMeanImputOut
	# missingDatOutLogical: logical, True for rows of medianOrMeanImputOut corresponding to
	# observations with missing values
	#cat("Are you sure medianImput does not come from uncertain-Y out-of-sample data?\n")
	stopifnot(class(doTrainModelOut)=="PolytopeSVRWW0")
	stopifnot(is.list(doPolyListOut))
	stopifnot(is.logical(missingDatOutLogical))

	stopifnot(length(missingDatOutLogical)==nrow(medianOrMeanImputOut))
	
	w= doTrainModelOut$w
	w0 = doTrainModelOut$w0
	# medianOrMeanImputOut: median/mean imputation of testing or validation data
	stopifnot("Y"%in%colnames(medianOrMeanImputOut))	
	p = ncol(medianOrMeanImputOut)-1
	stopifnot(which(colnames(medianOrMeanImputOut)=="Y")==(p+1))
	
	if(attr(doPolyListOut, "scaled")){
		scaleInfo = attr(doPolyListOut, "scaleInfo")
		scaleInfoMeanX = scaleInfo$mean[1:p]
		scaleInfoSdX = scaleInfo$std[1:p]
		scaleInfoMeanY = scaleInfo$mean[p+1]
		scaleInfoSdY = scaleInfo$std[p+1]
		#medianOrMeanImputOut = ScaleCenter(medianOrMeanImputOut, scaleInfo$mean, scaleInfo$std)
		scaledMedianImputOutX = ScaleCenter(medianOrMeanImputOut[, 1:p], scaleInfoMeanX, scaleInfoSdX)		
	} else{
		scaledMedianImputOutX = medianOrMeanImputOut[, 1:p]
		scaleInfoMeanY = 0
		scaleInfoSdY = 1
		#nothing=100000 # remove
	}	
	
	#predY = w0+medianOrMeanImputOut[, 1:p]%*%w
	predY = as.numeric((w0+scaledMedianImputOutX%*%w)*scaleInfoSdY+scaleInfoMeanY)
	trueY = medianOrMeanImputOut[, p+1]
	errList  = DoError(predY, trueY, missingDatOutLogical)
	return(errList)
}

DoErrorFold = function(doPolyListFoldsOut, doParListGridOut, replaceImputedWithTrueY = F, approach){
	stopifnot(is.list(doParListGridOut))
	stopifnot(is.list(doPolyListFoldsOut))
	stopifnot(is.logical(replaceImputedWithTrueY))
	stopifnot(is.character(approach))
	
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
			if(approach %in%c("doSquarebbQuant", "doSquarebbSd")){
				currFoldImputedOutDat = doPolyListFoldsOut[[i]]$outDat$imputed$meanBoxImputDat
			} else{
				cat("Assuming non-box polyhedra: extracting median of multiple imputations for out of sample predictions\n")
				currFoldImputedOutDat = doPolyListFoldsOut[[i]]$outDat$imputed$medianImputDat
			}
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


DoExtractErrMat = function(doErrorFoldOut, average=median){ ### median helps avoid issues when have a few Infs...
	stopifnot(is.list(doErrorFoldOut))
	stopifnot(is.function(average))
	errMatRes = doErrorFoldOut
	errMeasureVect = names(errMatRes[[1]]$errorList[[1]])
	errMeasureVect = errMeasureVect[errMeasureVect!="predRes"]
	
	for(i in  1:length(errMatRes)){
		# for each parameter combination
		errMat = NULL
		predYMat = NULL
		trueYMat = NULL
		missingDatOutLogicalMat = NULL
		
		for (errMeasure in errMeasureVect){
			#for each error measure
			errVect = NULL
			for(j in 1:length(errMatRes[[i]]$errorList)){
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
	bestParList = doExtractErrMatOut[[minAverageErrIdx]]$parList # for a given error measure, this is the combination of parameters that give
																# best (minimum) error performance
	bestAvgError = averageErrVect[minAverageErrIdx]
	if(bestAvgError==Inf){
		cat("Warning: Error measure '", errMeasure, "' can't be calculated because there is no data for it.\n")
		cat("Arbitrarily picking a parameter combination as 'best'\n")
	}
	res = list(errorMeasure = errMeasure, bestParList = bestParList, bestParListAvgError = bestAvgError,
		 bestParListIdx = minAverageErrIdx)
	return(res)
}		

DoMinMaxPrediction = function(polyListOutIndiv,doTrainModelOut, w=NULL, w0=NULL){
	if(any(is.null(c(w, w0)))){
		w=doTrainModelOut$w
		w0=doTrainModelOut$w0
	}
	p=length(w)
	stopifnot(p==ncol(polyListOutIndiv$A)-1)
	pPlusOne = p+1
	
	pPlusOneCol = polyListOutIndiv$A[,pPlusOne]
	nonZeroPlusOneCoeff = which(pPlusOneCol!=0)
	whichIsOne = which(pPlusOneCol==1)
	stopifnot(sum(pPlusOneCol[nonZeroPlusOneCoeff])==0)
	stopifnot(length(nonZeroPlusOneCoeff)==2)
	stopifnot(all(abs(pPlusOneCol[nonZeroPlusOneCoeff])==1))
	stopifnot(sum(polyListOutIndiv$a[nonZeroPlusOneCoeff])==0)
	
	trueY = polyListOutIndiv$a[whichIsOne]
	modelMinMax = list()
	modelMinMax$A = polyListOutIndiv$A
	modelMinMax$rhs = polyListOutIndiv$a
	modelMinMax$sense = polyListOutIndiv$dir
	modelMinMax$lb = rep(-Inf, p+1)
	modelMinMax$obj = c(w,0)
	gurobiParams = list(LogFile="", OutputFlag=0)
	
	resU = gurobi(modelMinMax, params=gurobiParams)
	xU = resU$x[1:p]
	yPredU = resU$objval + w0
	errU = trueY - yPredU 
	
	modelMinMax$obj = -c(w,0)
	resV = gurobi(modelMinMax, params=gurobiParams)
	xV = resV$x[1:p]
	yPredV = -resV$objval + w0
	errV = yPredV - trueY
	

	res = list(u=list(x=xU, yPred=yPredU, err=errU), v=list(x=xV, yPred=yPredV, err=errV))
	return(res)

}

