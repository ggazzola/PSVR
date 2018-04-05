BoxConstraintsFromUpperAndLowerBound = function(lowerBoundVect, upperBoundVect){
	
	stopifnot(is.vector(lowerBoundVect) & is.numeric(lowerBoundVect))
	stopifnot(is.vector(upperBoundVect) & is.numeric(upperBoundVect))
	stopifnot(length(upperBoundVect)==length(upperBoundVect))
	stopifnot(all(lowerBoundVect<=upperBoundVect))
	
	numBoxed = sum(lowerBoundVect< upperBoundVect)
	boxedLogical = lowerBoundVect< upperBoundVect
	numVars = length(lowerBoundVect)
		
	upperRHS = upperBoundVect
	lowerRHS = lowerBoundVect

	dimLength = upperRHS - lowerRHS
	if(numBoxed>0){
		volumeBoxed =  prod(dimLength[boxedLogical])
	} else{
		volumeBoxed = 0
	}
	
	lhs = rbind(diag(1, numVars, numVars), diag(-1, numVars, numVars))
	rhs = c(upperRHS, -lowerRHS)
	direction = c(rep("<=", length(lowerRHS)), rep("<=", length(upperRHS)))
	res = list(lhs=lhs, rhs=rhs, dir=direction, volume = volumeBoxed, length = dimLength, numUncertain = numBoxed, volumeDims = boxedLogical)
	return(res)
}

BoxConstraints = function(centerVect, halfSideVect, halfSideMultiplier) {
	# Returns lhs matrix, rhs vector and corresponding inequality operators to define 
	#	a box around a center
	# centerVect: middle point coordinates of the box
	# halfSideVect: half length of the sides of the box
	# halfSideMultiplier: non-negative scaler that stretches (>1) or shrinks (<1) each component of halfSideVect
	
	stopifnot(is.vector(centerVect) & is.numeric(centerVect))
	stopifnot(is.vector(halfSideVect) & is.numeric(halfSideVect) & all(halfSideVect>=0))
	stopifnot(length(centerVect)==length(halfSideVect))
	stopifnot(is.numeric(halfSideMultiplier))
	stopifnot(length(halfSideMultiplier)==1)
	stopifnot(halfSideMultiplier>=0)
	
	numBoxed = sum(halfSideVect>0 & halfSideMultiplier>0)
	boxedLogical = halfSideVect>0 & halfSideMultiplier>0
	numVars = length(halfSideVect)
		
	upperRHS = lowerRHS = centerVect
	upperRHS = upperRHS+halfSideVect*halfSideMultiplier
	lowerRHS = lowerRHS-halfSideVect*halfSideMultiplier
	dimLength = upperRHS - lowerRHS
	if(numBoxed>0){
		volumeBoxed =  prod(dimLength[boxedLogical])
	} else{
		volumeBoxed = 0
	}

	lhs = rbind(diag(1, numVars, numVars), diag(-1, numVars, numVars))
	
	rhs = c(upperRHS, -lowerRHS)
	direction = c(rep("<=", length(lowerRHS)), rep("<=", length(upperRHS)))
	res = list(lhs=lhs, rhs=rhs, dir=direction, volume = volumeBoxed, numUncertain = numBoxed, length = dimLength, volumeDims = which(boxedLogical))
	return(res)
	
}

DiagonalConstraints = function(centerVect, halfSideVect, halfSideMultiplier, slopeMat, interceptMultiplier){
	# Given a box, calculating enclosing 2D hyperplanes with given slope and intercept no larger/smaller than that of
	# of such hyperplanes intersecting extremal points of the box
	#centerVect, halfSideVect,  halfSideMultiplier: same as BoxConstraints
	# slopeMat: a symmetric matrix whose i,j entry is the slope of the hyperplane for pair i,j
	# interceptMultiplier: non-negative scaler that stretches (>1) or shrinks (<1) the intercepts of the 2D hyperplanes
	
	stopifnot(is.vector(centerVect) & is.numeric(centerVect))
	stopifnot(is.vector(halfSideVect) & is.numeric(halfSideVect))
	stopifnot(is.matrix(slopeMat) & is.numeric(slopeMat))	
	stopifnot(is.numeric(interceptMultiplier))
	stopifnot(is.numeric(halfSideMultiplier))
	stopifnot(length(centerVect)==length(halfSideVect))
	stopifnot(nrow(slopeMat)==ncol(slopeMat))
	stopifnot(length(centerVect)==nrow(slopeMat))
	stopifnot(length(interceptMultiplier)==1)
	stopifnot(length(halfSideMultiplier)==1)
	stopifnot(halfSideMultiplier>=0)
	
	if(all(halfSideVect*halfSideMultiplier==0)){
		cat("No diagonal constraints to be added, since the box is a point\n")
		return(NULL)
	}
	if(all(slopeMat==0)){
		cat("No diagonal constraints to be added, since all slopes are zero\n")
		return(NULL)
	}
	
	upperVect = lowerVect = centerVect
	upperVect = upperVect+halfSideVect*halfSideMultiplier #  box upper bound
	lowerVect = lowerVect-halfSideVect*halfSideMultiplier #  box lower bound
	
	numVars = length(lowerVect)
	whichBoxed = which(upperVect>lowerVect)
	numBoxed = length(whichBoxed)
	largerDirectionInterceptMat = smallerDirectionInterceptMat = matrix( , numVars, numVars)
	
	lhsTemp = rep(0, numVars)
	lhs = rhs = direction = NULL
	for(ii in 1:(numBoxed-1)){
		for(jj in (ii+1):numBoxed){
			# for each pair of variables, consider the corresponding 2d box
			i = whichBoxed[ii]
			j = whichBoxed[jj]
			upperUpper = upperVect[c(i,j)] #upper right corner
			lowerLower = lowerVect[c(i,j)] # lower left corner
				
			upperLower = c(upperVect[i], lowerVect[j]) #lower right corner
			lowerUpper = c(lowerVect[i], upperVect[j]) #upper left corner
	
			stopifnot(slopeMat[i,j]==slopeMat[j,i])

			if(slopeMat[i,j]!=0){
				if(slopeMat[i,j]<0){
					pointForSmallerDirection = upperUpper # if slope is negative, these are the extremal points of the box to consider
					pointForLargerDirection = lowerLower
				} else{
					pointForSmallerDirection = lowerUpper # if slope is positive, these are the extremal points of the box to consider
					pointForLargerDirection = upperLower
				}
				largerDirectionIntercept = pointForLargerDirection[2]-slopeMat[i,j]*pointForLargerDirection[1] #intercept (with the cartesian axes) of the intersecting "higher" line 
				smallerDirectionIntercept = pointForSmallerDirection[2]-slopeMat[i,j]*pointForSmallerDirection[1] #intercept (with the cartesian axes) of the intersecting "lower" line 
				largerDirectionInterceptMat[i,j] = largerDirectionInterceptMat[j,i] = largerDirectionIntercept*interceptMultiplier # shrinking/enlarging distance two of two lines from the center (acting on intercepts)
				smallerDirectionInterceptMat[i,j] = smallerDirectionInterceptMat[j,i] = smallerDirectionIntercept*interceptMultiplier
				lhsIJ = lhsTemp
				lhsIJ[c(j,i)] = c(1, -1) #x_j = mx_i +q
				rhsIJLarger =  largerDirectionInterceptMat[i,j]
				rhsIJSmaller =  smallerDirectionInterceptMat[i,j]
				lhs = rbind(lhs, lhsIJ, lhsIJ)
				rhs = rbind(rhs, rhsIJSmaller, rhsIJLarger)
				direction = c(direction, c("<=", ">="))
			}
		}
	}
	res = list(lhs=lhs, rhs=rhs, dir=direction)
	return(res)	
}

PCConstraintsFull = function(dat, propIn = .8,  projectDimsLogical = NULL){
	#, preCalcPc = NULL
	
	#projectDimsLogical: T for variables along which the pc should be calculated
	# essentially the same as PCConstraintsObsolete, except the former removes constant dimensions prior to calculating PCs 
	# to speed up computation. Perhaps use this to make sure not making mistakes, unless cpu time really too high
	# Note: propIn is a quantile calculated relative to the PC, not the original axes! so propIn =0 will correspond
	# to the median of the data along the PCs, not along the original axes
	# RECALL: no standardization is carried out on the data prior to calculating the PCs here!!!
	# dat could, e.g., be repeated "measurements" of the same observation
	stopifnot(is.matrix(dat)) # even if single data point
	p = ncol(dat)
	n = nrow(dat)
	stopifnot(n>=p)
	
	stopifnot(propIn>=0 & propIn<=1)
	if(!is.null(projectDimsLogical)){
		stopifnot(is.logical(projectDimsLogical)) # false for dimensions to disregard
	#	stopifnot(!is.null(preCalcPc))
	#	stopifnot(is.matrix(preCalcPc))
	#	stopifnot(ncol(preCalcPc)==nrow(preCalcPc))
		stopifnot(length(projectDimsLogical)==ncol(dat))
		stopifnot(sum(projectDimsLogical)>0)
		stopifnot(propIn==1)
		
		stopifnot(all(apply(dat, 2, function(x)diff(range(x)))>0))
		dat[, !projectDimsLogical] = 10^6  # useful to calculate numCertainDims correctly (can be any finite constant value)
		project = T
	} else{
	#	stopifnot(is.null(preCalcPc))
		project = F
		projectDimsLogical = rep(F, p)
	}
		
	#if(project){
	#	dat[, !projectDimsLogical] = 1000  # useful to calculate numCertainDims correctly (can be any constant value)
	#}

	
	minDat = apply(dat, 2, min)
	maxDat = apply(dat, 2, max)

	uncertainDimsLogical = ((maxDat-minDat)>0 | projectDimsLogical) 
	numUncertainDims = sum(uncertainDimsLogical)
	numCertainDims = sum(!uncertainDimsLogical)
		
#	if(!project){
		datTmp = dat[, uncertainDimsLogical]
		#ee<<-datTmp
		#dd<<-dat
		#uu<<-uncertainDimsLogical
		pc = prcomp(datTmp, center=T)
		fullDimPc = matrix(0 , nrow=p, ncol=p)

		fullDimPc[uncertainDimsLogical, uncertainDimsLogical] = pc$rotation
		fullDimPc[!uncertainDimsLogical, !uncertainDimsLogical] = diag(1, nrow=numCertainDims, ncol=numCertainDims)
		fullDimPc = cbind(fullDimPc[, uncertainDimsLogical], fullDimPc[, !uncertainDimsLogical])
	#}
#	} else{
#		fullDimPc = preCalcPc
#		fullDimPc[!projectDimsLogical,] = 0 
#	}	
	
	projDat = dat%*%fullDimPc
	lowerQuantile = (1-propIn)/2
	upperQuantile = 1-lowerQuantile
	
	projDatSd = as.numeric(apply(projDat, 2, sd))
	nonConstantProjDatDimsLogical = projDatSd>0
	constantProjDatDimsLogical = projDatSd==0
	
	minProjDat = as.numeric(apply(projDat, 2, quantile, probs=lowerQuantile)) #  smallest-coordinate point along each pc box 
	maxProjDat = as.numeric(apply(projDat, 2, quantile, probs=upperQuantile)) #  largest-coordinate point along each pc box
	if(any(constantProjDatDimsLogical))
		minProjDat[constantProjDatDimsLogical] =  maxProjDat[constantProjDatDimsLogical] #  else min/max may come out of quantile numerically different 
	maxNoLargerThanMinNonConstantProjDatDimsLogical = (maxProjDat - minProjDat<=0) & nonConstantProjDatDimsLogical
		 # this could still happen along non-constant dims b/c of quantile (although with propIn>0 it should never happen)
		 minProjDatPrima = minProjDat
	if(any(maxNoLargerThanMinNonConstantProjDatDimsLogical)){
		#print(ifelse(propIn>0, 1e-16, 0))
		#print(minProjDat[maxNoLargerThanMinNonConstantProjDatDimsLogical])
		#aa<<-minProjDat[maxNoLargerThanMinNonConstantProjDatDimsLogical]
		minProjDat[maxNoLargerThanMinNonConstantProjDatDimsLogical] =  
			maxProjDat[maxNoLargerThanMinNonConstantProjDatDimsLogical]-.Machine$double.eps*2# min will forced to be at least an epsilon smaller than max along non constant dims
		if(any(minProjDat[maxNoLargerThanMinNonConstantProjDatDimsLogical]-maxProjDat[maxNoLargerThanMinNonConstantProjDatDimsLogical]==0)){
			stop("The epsilon trick didn't work\n")
		}
			#ifelse(propIn>0, .Machine$double.neg.eps, 0) 
		#bb<<-minProjDat[maxNoLargerThanMinNonConstantProjDatDimsLogical]	
			
	}
	#print(minProjDat)
	#print(maxProjDat)
	#print(projDat)
	lengthProjDat = maxProjDat - minProjDat
	#if(any(lengthProjDat<0)){
	#	projDat<<-projDat
	#	whichNonConstantProjDatDims<<-whichNonConstantProjDatDims
	#	lowerQuantile<<-lowerQuantile
	#	upperQuantile<<-upperQuantile
	#	stop()
	#}
	
	numLenghtNonZeroProjDat = sum(lengthProjDat!=0) # this could be < numUncertainDims if the top and bottom quantile happen to be the same along uncertain dimensions (e.g., because propIn is small or because all values along a certain projected dimension are  equal except for one, etc. ), which is why we make sure the difference is at least 1e-16
			
	stopifnot(numLenghtNonZeroProjDat==numUncertainDims)
		
	if(numUncertainDims>0 & propIn>0){
		volumeProjDat = prod(lengthProjDat[lengthProjDat!=0])
	} else{
		volumeProjDat = 0
	}
	
	lhs = rhs = direction = NULL
	for(i in 1:p){
		#if(sum(fullDimPc[,i]==1)!=1 | sum(fullDimPc[,i]==0)!=(p-1)){
			lhs = rbind(lhs, -fullDimPc[,i], fullDimPc[,i]) # minus signs for having <= inequality direction 
			rhs = c(rhs, -minProjDat[i], maxProjDat[i])# minus signs for having <= inequality direction 
			direction = c(direction, c("<=", "<="))
			#}
	}
	
	#volumeDims are the dimensions of dat that have uncertain data/we project to
	res = list(lhs=lhs, rhs=rhs, dir=direction, pc = fullDimPc, propIn = propIn, length = lengthProjDat, volume = volumeProjDat, volumeDims = which(uncertainDimsLogical)) #pc = fullDimPc
	return(res)
}

FixedPointConstraints = function(datPoint){
	# writing an individual p-dimensional point as the intersection of 2*p inequalities
	stopifnot(is.vector(datPoint))
	p = length(datPoint)
	lhs = rbind(diag(1, p, p), diag(-1, p, p))
	rhs = c(datPoint, -datPoint)
	direction = rep("<=", 2*p)
	res = list(lhs=lhs, rhs=rhs, dir = direction, length = rep(0, p), volume = 0, volumeDims = integer(0), numUncertainDims = 0)
	return(res)
}


MonteCarloPolytopeFilling = function(missingPoint, imputationMat, lhsMat, rhsVect, pertN, pertExtraSize, dimToPlot, mainTitle=NULL){
	
	#plot points inside pc-bb of a given point
	# missingPoint is the point with NAs
	# imputationMat is the matrix of imputations on missingPoint
	#lhsMat/rhsVect are the left/right hand side of the linear inequalities describing the pc-bb
	#pertN is the number of 'perturbations' to sample within the ranges of imputationMat +/ pertExtraSize (a constant across all variables)
	#dimToPlot are the two dimensions along which the results should be plotted
	p = length(missingPoint) # normally it actually corresponds to pPlusOne...
	
	sdImputMat = apply(imputationMat,2, sd)
	stopifnot(which(sdImputMat>0)==which(is.na(missingPoint)))
	#stopifnot(any(is.na(missingPoint)))
	
	stopifnot(is.matrix(imputationMat))
	stopifnot(all(!is.na(imputationMat)))
	stopifnot(p==ncol(imputationMat))
	stopifnot(p==ncol(lhsMat))
	stopifnot(nrow(lhsMat)==length(rhsVect))
	stopifnot(pertN>0)
	stopifnot(pertExtraSize>=0)
	stopifnot(dimToPlot%in%(1:p))
	colNames = colnames(imputationMat)
	stopifnot("Y"%in%colNames)
	
	
	perturbMat = matrix( , nrow=pertN, ncol=p)
	minImput = apply(imputationMat, 2, min)
	maxImput = apply(imputationMat, 2, max)
	missingDims = which(is.na(missingPoint))
	perturbMat = NULL
	for(j in 1:p){
		smallestPerturbedVal = minImput[j]
		largestPerturbedVal = maxImput[j]
		if(j %in%missingDims){
			smallestPerturbedVal = smallestPerturbedVal - pertExtraSize
			largestPerturbedVal = largestPerturbedVal + pertExtraSize
			
		}
	
		perturbMat = cbind(perturbMat, runif(pertN, smallestPerturbedVal, largestPerturbedVal))
	}
	colorVect = numeric(n)
	
	lhsPerturb = perturbMat%*%t(lhsMat)
	
	for(i in 1:pertN){
		if(all(lhsPerturb[i,]<=rhsVect)){
			colorVect[i]=4 # inside polytope
		} else{
			colorVect[i]=1 # outside polytope
		}
	}
	perturbMatFirstDim = perturbMat[, dimToPlot[1]]
	perturbMatSecondDim = perturbMat[, dimToPlot[2]]
	
	firstDim = imputationMat[, dimToPlot[1]]
	firstDimLab = colNames[dimToPlot[1]]
	secondDim = imputationMat[, dimToPlot[2]]
	secondDimLab = colNames[dimToPlot[2]]
	firstDimRange = range(firstDim+pertExtraSize, firstDim-pertExtraSize)
	secondDimRange = range(secondDim+pertExtraSize, secondDim-pertExtraSize)
	
	
	plot(firstDim, secondDim, xlab=firstDimLab, ylab=secondDimLab, xlim=firstDimRange, ylim=secondDimRange, col=2, lwd=2, main = mainTitle) #imputations
	points(perturbMatFirstDim, perturbMatSecondDim, col=colorVect)
	
	
	if(F){ # just to plot something...
		pertExtraSize=.5
		pertN=10000
		qualeInnerFold = 2
		qualePoint=3
		originalPoint = currDoDataSplitOutInnerOriginal[[qualeInnerFold]]$inDat$original[qualePoint,]
		missingPoint = currDoDataSplitOutInnerOriginal[[qualeInnerFold]]$inDat$missing[qualePoint,]

		dimToPlot = which(is.na(missingPoint))[1:2] # need to change this if the point has less than two missing dimensions

		imputationMat = NULL
		for(k in 1:length(doMultipleImputationFoldsOutInner[[qualeInnerFold]]$inDat$imputed$imputDatList)){
			imputationMat=rbind(imputationMat, doMultipleImputationFoldsOutInner[[qualeInnerFold]]$inDat$imputed$imputDatList[[k]][qualePoint,])
		}

		polyPoint = doPolyListFoldsOutInner[[qualeInnerFold]]$polyList[[qualePoint]]
		lhsMat = polyPoint$A
		rhsVect = polyPoint$a

		if(attr(doPolyListFoldsOutInner[[qualeInnerFold]]$polyList, "scaled")){
			scaleTitle = "scale = T"
			mm = attr(doPolyListFoldsOutInner[[qualeInnerFold]]$polyList, "scaleInfo")$mean
			ss = attr(doPolyListFoldsOutInner[[qualeInnerFold]]$polyList, "scaleInfo")$std
			imputationMat = ScaleCenter(imputationMat,mm,ss)
			originalPoint = ScaleCenter(originalPoint, mm, ss)
		} else{
			scaleTitle = "scale = F"
		}
		quantOrSdPropTitle = paste("quantOrSdProp = ", attr(doPolyListFoldsOutInner[[qualeInnerFold]]$polyList, "quantOrSdProp"))
		
		
		MonteCarloPolytopeFilling (missingPoint, imputationMat, lhsMat, rhsVect, pertN, pertExtraSize, dimToPlot, paste(quantOrSdPropTitle, scaleTitle))

	}		
	
}

PCConstraintsObsolete = function(dat, propIn = .8){
	

	stopifnot(is.matrix(dat))
	stopifnot(propIn>0 & propIn<=1)
	p = ncol(dat)
	#scaleCenterDat = scale(dat)
	
	#datMean = attr(scaleCenterDat, "scaled:center")
	#datSd = attr(scaleCenterDat, "scaled:scale")
	#uncertainDimsLogical = datSd!=0
	minDat = apply(dat, 2, min)
	maxDat = apply(dat, 2, max)
	
	uncertainDimsLogical = (maxDat-minDat)>0
	numUncertainDims = sum(uncertainDimsLogical)
	
	if(numUncertainDims==0){
		cat("No uncertainty: returning NULL\n")
		return(NULL)
	}
	#datMean = datMean[uncertainDimsLogical]
	#datSd = datSd[uncertainDimsLogical]
	#sdMat = matrix(0, nrow=numUncertainDims, ncol=numUncertainDims)
	#diag(sdMat) = datSd	
	#scaleCenterDat = scaleCenterDat[, uncertainDimsLogical]
	#pcScaleCenter = prcomp(scaleCenterDat, scale=T, center=T)$rotation
	pc = prcomp(dat[, uncertainDimsLogical], center=T)$rotation 	# if working on standardized data, end up with non-orthogonal components when un-standardizing

	#pc = sdMat%*%pcScaleCenter 
	fullDimPc = matrix(0 , nrow=p, ncol=p)
	for(i in 1:p){
		if (uncertainDimsLogical[i]){
			fullDimPc[uncertainDimsLogical,i] = pc[,i]
		} else{
			fullDimPc[i, ] = 1
		}
	}
	#projDat = dat%*%fullDimPc
	projDat = dat[,uncertainDimsLogical]%*%pc
	
	lowerQuantile = (1-propIn)/2
	upperQuantile = 1-lowerQuantile
	minProjDat = apply(projDat, 2, quantile, probs=lowerQuantile, type=1) # index of smallest-coordinate point along each pc
	maxProjDat = apply(projDat, 2, quantile, probs=upperQuantile, type=1) # index of largest-coordinate point along each pc
	lhs = rhs = direction = NULL
	for(i in 1:numUncertainDims){
		lhs = rbind(lhs, -fullDimPc[,i], fullDimPc[,i]) # minus signs for having <= inequality direction 
		rhs = c(rhs, -minProjDat[i], maxProjDat[i])# minus signs for having <= inequality direction 
		direction = c(direction, c("<=", "<="))
	}
	
	res = list(lhs=lhs, rhs=rhs, dir=direction, pc = fullDimPc, volume) #pc = fullDimPc
	return(res)
}
