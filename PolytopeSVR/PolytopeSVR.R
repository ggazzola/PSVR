PolytopeSVR = function(polyList, Ccertain, Cuncertain, epsilonCertain, extraEpsilonUncertain, uncertaintySpecialTreatment=T, twoSlacks = F,
	 linear = F, returnAllInfo = F, ...){
	# Creates SVR model object for Gurobi
	# Equivalent to PolytopeSVRNoUncertainSpecialTreatment if uncertaintySpecialTreatment=F, 
	#	provided Ccertain =C and epsilonCertain=epsilon
	#polyList: a list of n lists, one for each training point
	# polyList[[i]]$A and polyList[[i]]$a are the lhs and rhs of the inequalities defining the uncertainty on the i-th point
	# 	ACROSS ALL p+1 VARIABLES
	# The direction of all inequalities (polyList[[i]]$dir) should always be '<='
	# If there is no uncertainty on any of the p+1 variables, 
	#	e.g., x_1=3, then A will contain two rows corresponding to x_1<=3 and -x_1<=-3)
	# Cuncertain, extraEpsilonUncertain same as c^m and \epsilon^m in formulation
	stopifnot(is.list(polyList))
	stopifnot(is.numeric(epsilonCertain) & length(epsilonCertain)==1)
		
	stopifnot(is.numeric(Ccertain) & length(Ccertain)==1)
	stopifnot(is.logical(uncertaintySpecialTreatment) & is.logical(twoSlacks))
	if(uncertaintySpecialTreatment)
		stopifnot(!twoSlacks) 
	
	if(uncertaintySpecialTreatment & !twoSlacks){
		stopifnot(is.numeric(Cuncertain) & length(Cuncertain)==1)
		stopifnot(is.numeric(extraEpsilonUncertain) & length(extraEpsilonUncertain)==1)
		uncertaintyVect = 	sapply(polyList, function(x) {x$uncertainty})
		if(is.character(uncertaintyVect))
			stop("Can't proceed with uncertaintySpecialTreatment because uncertaintyVect[1] is ", uncertaintyVect[1])
		stopifnot(is.numeric(uncertaintyVect) & min(uncertaintyVect)>=0 & max(uncertaintyVect)<=1)
		uncertainPoints = uncertaintyVect>0
	}
	n = length(polyList)
	stopifnot(n>1)
	
	twoTimesN = 2*n
	oneThroughN = 1:n
	if(twoSlacks){
		costVect = rep(0, twoTimesN)
	} else{
		costVect = rep(0, n)
	}
	
	if(uncertaintySpecialTreatment & !twoSlacks){
		if(sum(uncertainPoints)>0){
			costVect[uncertainPoints] = Cuncertain*(1-uncertaintyVect[uncertainPoints])
			costVect[!uncertainPoints] = Ccertain
		} else{
			costVect[oneThroughN]=Ccertain
		}
	} else{
		if(twoSlacks){
			costVect[1:(2*n)]=Ccertain
		} else{
			costVect[oneThroughN]=Ccertain
		}
	}
	totNumRowA = 0
	
	pPlusOne = ncol(polyList[[1]]$A)
	twoTimesPplusOne = pPlusOne*2
	p = pPlusOne-1	
	twoTimesP = 2*p
	oneThroughP = 1:p
	
	wBlock = rbind(matrix(0, nrow=2, ncol=p), diag(1, nrow=p, ncol=p), matrix(0, nrow=1, ncol=p), diag(-1, nrow=p, ncol=p), matrix(0, nrow=1, ncol=p)) 
	wZeroBlock = matrix(c(-1, 1, rep(0, twoTimesPplusOne)), ncol=1)
	if(twoSlacks){
		block1 = c(-1, rep(0, 1+twoTimesPplusOne))
		block2 = c(0, -1, rep(0, twoTimesPplusOne))
		csiColBlock = matrix(c(block1, block2), ncol=2)
	} else {
		csiColBlock = matrix(c(-1, -1, rep(0, twoTimesPplusOne)), ncol=1)
	}
	
	for(i in 1:n){
		A = polyList[[i]]$A
		a = polyList[[i]]$a
		dir = polyList[[i]]$dir
		uvSize = length(a)
		stopifnot(is.matrix(A))
		stopifnot(is.numeric(a))
		stopifnot(nrow(A)==uvSize)
		stopifnot(ncol(A)==pPlusOne)
		stopifnot(all(dir=="<="))
		aTranspose = t(a)
		ATranspose = t(A)
		totNumRowA = totNumRowA + uvSize
		uBlock = rbind(aTranspose, matrix(0, ncol = uvSize), ATranspose, matrix(0, nrow=pPlusOne, ncol=uvSize))
		vBlock = rbind(matrix(0, ncol = uvSize), aTranspose, matrix(0, nrow=pPlusOne, ncol=uvSize),ATranspose)
		polyList[[i]]$uvBlock = cbind(uBlock, vBlock)
	}

	twoTimesTotNumRowA =  2*totNumRowA
	oneThroughUvSize =  1:uvSize

	if(twoSlacks){
		lbVect = c(rep(-Inf, pPlusOne), rep(0, twoTimesTotNumRowA+twoTimesN))  
	} else{
		lbVect = c(rep(-Inf, pPlusOne), rep(0, twoTimesTotNumRowA+n))  
	}
	
	blockRowSize = 2+twoTimesPplusOne
	nTimesBlockRowSize = n*blockRowSize
	
	if(twoSlacks){
		ncolLhsMatrix = pPlusOne+twoTimesTotNumRowA+twoTimesN
	} else{
		ncolLhsMatrix = pPlusOne+twoTimesTotNumRowA+n		
	}
	lhsMatrix = matrix(0 , nrow = nTimesBlockRowSize, ncol = ncolLhsMatrix)
	
	
	if(returnAllInfo){
		uVNames = NULL
		for(i in oneThroughN){
			uVNames = c(uVNames, c(paste0(paste0(rep("u", uvSize), i), oneThroughUvSize), paste0(paste0(rep("v", uvSize), i), oneThroughUvSize)))
		}	
	
		if(twoSlacks){
			varNames = c(paste0("w", oneThroughP), "w0", uVNames, paste0("csiPlus", oneThroughN), paste0("csiMinus", oneThroughN))
		} else{
			varNames = c(paste0("w", oneThroughP), "w0", uVNames, paste0("csi", oneThroughN))
		}
	}
	
	senseVect = rep(c("<=", "<=", rep("=", twoTimesPplusOne)), n)
	
	coefVect = c(rep(0, pPlusOne+twoTimesTotNumRowA), costVect)

	if(twoSlacks){
		QMatrix = diag(0, pPlusOne+twoTimesTotNumRowA+twoTimesN, pPlusOne+twoTimesTotNumRowA+twoTimesN)
	} else{
		QMatrix = diag(0, pPlusOne+twoTimesTotNumRowA+n, pPlusOne+twoTimesTotNumRowA+n)
	}	
	
	QMatrix[oneThroughP, oneThroughP]= diag(1, p, p)
	startUvCol = pPlusOne+1

	uIdx = vIdx = list() 
	if(twoSlacks){
		csiPlusIdx = csiMinusIdx = list()
	} else{
		csiIdx = list()
	}
	
	rhsVect = NULL
	for(i in oneThroughN){
		startRow = (i-1)*blockRowSize+1
		endRow = i*blockRowSize
		blockRows = startRow:endRow
		endUvCol = startUvCol + ncol(polyList[[i]]$uvBlock)-1
		uvBlockCol = startUvCol:endUvCol
		numUvCol = length(uvBlockCol)
		lhsMatrix[blockRows,oneThroughP]=wBlock
		lhsMatrix[blockRows,pPlusOne]=wZeroBlock
		lhsMatrix[blockRows, uvBlockCol] = polyList[[i]]$uvBlock
		
		currEpsilon = epsilonCertain
		if(uncertaintySpecialTreatment & !twoSlacks)
			if(uncertainPoints[i])
				currEpsilon = currEpsilon + extraEpsilonUncertain*uncertaintyVect[i]
		
		rhsVect = c(rhsVect, c(currEpsilon, currEpsilon, rep(0, p), 1, rep(0,p), -1))		
		
		uIdx[[i]] = uvBlockCol[1:(numUvCol/2)] 
		vIdx[[i]] = uvBlockCol[(numUvCol/2+1):numUvCol]
		if(twoSlacks){
			csiPlusIdx[[i]] = pPlusOne+twoTimesTotNumRowA + 2*i-1
			csiMinusIdx[[i]] = pPlusOne+twoTimesTotNumRowA + 2*i
		} else{
			csiIdx[[i]] = pPlusOne+twoTimesTotNumRowA + i
		}
			
		startUvCol = startUvCol + ncol(polyList[[i]]$uvBlock)
		if(twoSlacks){
			lhsMatrix[blockRows, c(csiPlusIdx[[i]], csiMinusIdx[[i]])] = csiColBlock
		} else{
			lhsMatrix[blockRows, csiIdx[[i]]] = csiColBlock
		}
	}
	
	if(linear){
		extendedLhsMatrixNcol = ncolLhsMatrix+p
		extendedLhsMatrixNrow = nTimesBlockRowSize+twoTimesP
		extendedLhsMatrix = matrix(0 , nrow = extendedLhsMatrixNrow, ncol = extendedLhsMatrixNcol)
		
		extraTopBlockLhs = matrix(0 , nrow=twoTimesP, ncol=extendedLhsMatrixNcol)
		extraLeftBlockLhs = matrix(0, nrow = nTimesBlockRowSize, ncol=p)
		cntRow = 1
		for(i in oneThroughP){
			extraTopBlockLhs[cntRow, c(i, i+p)] = -1
			extraTopBlockLhs[cntRow+1, c(i, i+p)] = c(-1,1)
			cntRow = cntRow+2
		}
		
		extraTopBlockRhs = rep(0, twoTimesP)
		
		extraTopSenseVect = rep("<=", twoTimesP)
		
		extraLbVect = rep(-Inf, p)
		extraCoef = rep(1, p)
		extraColumns = p
		
		#lhsMatrix = rbind(extraTopBlockLhs, cbind(extraLeftBlockLhs,lhsMatrix))
		
		extendedLhsMatrix[1:twoTimesP, ] = extraTopBlockLhs
		extendedLhsMatrix[(twoTimesP+1):(extendedLhsMatrixNrow), oneThroughP] = extraLeftBlockLhs
		extendedLhsMatrix[(twoTimesP+1):(extendedLhsMatrixNrow), pPlusOne:(extendedLhsMatrixNcol)] = lhsMatrix
		lhsMatrix = extendedLhsMatrix
		
		rhsVect = c(extraTopBlockRhs, rhsVect)
		senseVect = c(extraTopSenseVect, senseVect)
		lbVect = c(extraLbVect, lbVect)
		coefVect = c(extraCoef, coefVect)
		
		if(returnAllInfo){
			extraNames = paste0("b", oneThroughP)
			varNames = c(extraNames, varNames)
		}
		
	} else{
		extraColumns = 0
	}

	
	model = list()
	model$A = lhsMatrix
	model$rhs = rhsVect
	model$sense = senseVect
	model$lb = lbVect
	model$obj = coefVect
	if(!linear)
		model$Q=QMatrix
	varIdx = list()

	varIdx$w = extraColumns+ oneThroughP
	varIdx$w0 = extraColumns+pPlusOne
	if(returnAllInfo){
		if(linear)
			varIdx$b=oneThroughP
		varIdx$u = lapply(uIdx, function(x) {x+ extraColumns})
		varIdx$v = lapply(vIdx, function(x) {x+ extraColumns})
		varIdx$allU = unlist(varIdx$u)
		varIdx$allV = unlist(varIdx$v)
		if(twoSlacks){
			varIdx$csiPlus = lapply(csiPlusIdx, function(x) {x+ extraColumns})
			varIdx$csiMinus = lapply(csiMinusIdx, function(x) {x+ extraColumns})
			varIdx$allCsiPlus = unlist(varIdx$csiPlus)
			varIdx$allCsiMinus = unlist(varIdx$csiMinus)
		} else{
			varIdx$csi = lapply(csiIdx, function(x) {x+ extraColumns})
			varIdx$allCsi = unlist(varIdx$csi)
		}
		model$varNames = varNames
		
	}
	model$varIdx = varIdx
	class(model)="polytopeSVR"
	return(model)
}	

PolytopeSVRNoUncertainSpecialTreatment = function(polyList, C = 10, epsilon = 0){
	#polyList: a list of n lists, one for each training point
	# polyList[[i]]$A and polyList[[i]]$a are the lhs and rhs of the inequalities defining the uncertainty on the i-th point
	# 	ACROSS ALL p+1 VARIABLES
	# The direction of all inequalities (polyList[[i]]$dir) should always be '<='
	# If there is no uncertainty on any of the p+1 variables, 
	#	e.g., x_1=3, then A will contain two rows corresponding to x_1<=3 and -x_1<=-3)
	
	stopifnot(is.list(polyList))
	stopifnot(is.numeric(epsilon) & length(epsilon)==1)
	stopifnot(is.numeric(C) & length(C)==1)

	totNumRowA = 0
	pPlusOne = ncol(polyList[[1]]$A)
	p = pPlusOne-1	
	wBlock = rbind(matrix(0, nrow=2, ncol=p), diag(-1, nrow=p, ncol=p), matrix(0, nrow=1, ncol=p), diag(1, nrow=p, ncol=p), matrix(0, nrow=1, ncol=p)) 
	wZeroBlock = matrix(c(1, -1, rep(0, pPlusOne*2)), ncol=1)
	csiColBlock = matrix(c(-1, -1, rep(0, pPlusOne*2)), ncol=1)

	n = length(polyList)
	for(i in 1:n){
		A = polyList[[i]]$A
		a = polyList[[i]]$a
		dir = polyList[[i]]$dir
		uvSize = length(a)
		stopifnot(is.matrix(A))
		stopifnot(is.numeric(a))
		stopifnot(nrow(A)==uvSize)
		stopifnot(ncol(A)==pPlusOne)
		stopifnot(all(dir=="<="))
		totNumRowA = totNumRowA + uvSize
		uBlock = rbind(t(a), matrix(0, ncol = uvSize), t(A), matrix(0, nrow=pPlusOne, ncol=uvSize))
		vBlock = rbind(matrix(0, ncol = uvSize), t(a), matrix(0, nrow=pPlusOne, ncol=uvSize),t(A))
		polyList[[i]]$uvBlock = cbind(uBlock, vBlock)
	}

	lbVect = c(rep(-Inf, pPlusOne), rep(0, 2*totNumRowA+n))  

	blockRowSize = 2+2*pPlusOne
	lhsMatrix = matrix(0 , nrow = n*blockRowSize, ncol = pPlusOne+2*totNumRowA+n)
	rhsVect = rep(c(epsilon, epsilon, rep(0, p), -1, rep(0,p), 1), n)
	senseVect = rep(c("<=", "<=", rep("=", pPlusOne*2)), n)
	coefVect = c(rep(0, pPlusOne+2*totNumRowA), rep(C, n))
	QMatrix = diag(0, pPlusOne+2*totNumRowA+n, pPlusOne+2*totNumRowA+n)
	QMatrix[1:p, 1:p]= diag(1, p, p)
	startUvCol = pPlusOne+1

	uIdx = vIdx = csiIdx = list()
	for(i in 1:n){
		startRow = (i-1)*blockRowSize+1
		endRow = i*blockRowSize
		blockRows = startRow:endRow
		endUvCol = startUvCol + ncol(polyList[[i]]$uvBlock)-1
		uvBlockCol = startUvCol:endUvCol
		numUvCol = length(uvBlockCol)
		lhsMatrix[blockRows,1:p]=wBlock
		lhsMatrix[blockRows,pPlusOne]=wZeroBlock
		lhsMatrix[blockRows, uvBlockCol] = polyList[[i]]$uvBlock
		uIdx[[i]] = uvBlockCol[1:(numUvCol/2)]
		vIdx[[i]] = uvBlockCol[(numUvCol/2+1):numUvCol]
		csiIdx[[i]] = pPlusOne+2*totNumRowA + i
		startUvCol = startUvCol + ncol(polyList[[i]]$uvBlock)
		lhsMatrix[blockRows, pPlusOne+2*totNumRowA+i] = csiColBlock
	}

	model = list()
	model$A = lhsMatrix
	model$rhs = rhsVect
	model$sense = senseVect
	model$lb = lbVect
	model$obj = coefVect
	model$Q=QMatrix
	varIdx = list()
	varIdx$w = 1:p
	varIdx$w0 = pPlusOne
	varIdx$u = uIdx
	varIdx$v = vIdx
	varIdx$csi = csiIdx
	model$varIdx = varIdx
	class(model)="polytopeSVR"
	return(model)
}	

