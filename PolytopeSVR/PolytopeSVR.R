PolytopeSVR = function(polyList, Ccertain, Cuncertain, epsilonCertain, extraEpsilonUncertain, uncertaintySpecialTreatment=T, twoSlacks = F, ...){
	# Creates SVR model object for Gurobi
	# Equivalent to PolytopeSVRNoUncertainSpecialTreatment if uncertaintySpecialTreatment=F, 
	#	provided Ccertain =C and epsilonCertain=epsilon
	#polyList: a list of n lists, one for each training point
	# polyList[[i]]$A and polyList[[i]]$a are the lhs and rhs of the inequalities defining the uncertainty on the i-th point
	# 	ACROSS ALL p+1 VARIABLES
	# The direction of all inequalities (polyList[[i]]$dir) should always be '<='
	# If there is no uncertainty on any of the p+1 variables, 
	#	e.g., x_1=3, then A will contain two rows corresponding to x_1<=3 and -x_1<=-3)
	# Cuncertain, extraEpsilonUncertain same as alpha and beta in formulation
	stopifnot(is.list(polyList))
	stopifnot(is.numeric(epsilonCertain) & length(epsilonCertain)==1)
		
	stopifnot(is.numeric(Ccertain) & length(Ccertain)==1)
	stopifnot(is.logical(uncertaintySpecialTreatment) & is.logical(twoSlacks))
	if(uncertaintySpecialTreatment)
		stopifnot(!twoSlacks) 
	
	if(uncertaintySpecialTreatment & !twoSlacks){
		stopifnot(is.numeric(Cuncertain) & length(Cuncertain)==1)
		stopifnot(is.numeric(extraEpsilonUncertain) & length(extraEpsilonUncertain)==1)
		uncertaintyVect = 	sapply(polyList, function(x) x$uncertainty)
		if(is.character(uncertaintyVect))
			stop("Can't proceed with uncertaintySpecialTreatment because uncertaintyVect[1] is ", uncertaintyVect[1])
		stopifnot(is.numeric(uncertaintyVect) & min(uncertaintyVect)>=0 & max(uncertaintyVect)<=1)
		uncertainPoints = uncertaintyVect>0
	}
	n = length(polyList)
	stopifnot(n>1)
	
	if(twoSlacks){
		costVect = rep(0, 2*n)
	} else{
		costVect = rep(0, n)
	}
	
	if(uncertaintySpecialTreatment & !twoSlacks){
		if(sum(uncertainPoints)>0){
			costVect[uncertainPoints] = Cuncertain*(1-uncertaintyVect[uncertainPoints])
			costVect[!uncertainPoints] = Ccertain
		} else{
			costVect[1:n]=Ccertain
		}
	} else{
		if(twoSlacks){
			costVect[1:(2*n)]=Ccertain
		} else{
			costVect[1:n]=Ccertain
		}
	}
	totNumRowA = 0
	pPlusOne = ncol(polyList[[1]]$A)
	p = pPlusOne-1	
	wBlock = rbind(matrix(0, nrow=2, ncol=p), diag(1, nrow=p, ncol=p), matrix(0, nrow=1, ncol=p), diag(-1, nrow=p, ncol=p), matrix(0, nrow=1, ncol=p)) 
	wZeroBlock = matrix(c(-1, 1, rep(0, pPlusOne*2)), ncol=1)
	if(twoSlacks){
		block1 = c(-1, rep(0, 1+pPlusOne*2))
		block2 = c(0, -1, rep(0, pPlusOne*2))
		csiColBlock = matrix(c(block1, block2), ncol=2)
	} else {
		csiColBlock = matrix(c(-1, -1, rep(0, pPlusOne*2)), ncol=1)
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
		totNumRowA = totNumRowA + uvSize
		uBlock = rbind(t(a), matrix(0, ncol = uvSize), t(A), matrix(0, nrow=pPlusOne, ncol=uvSize))
		vBlock = rbind(matrix(0, ncol = uvSize), t(a), matrix(0, nrow=pPlusOne, ncol=uvSize),t(A))
		polyList[[i]]$uvBlock = cbind(uBlock, vBlock)
	}

	if(twoSlacks){
		lbVect = c(rep(-Inf, pPlusOne), rep(0, 2*totNumRowA+2*n))  
	} else{
		lbVect = c(rep(-Inf, pPlusOne), rep(0, 2*totNumRowA+n))  
	}
	
	blockRowSize = 2+2*pPlusOne
	
	if(twoSlacks){
		lhsMatrix = matrix(0 , nrow = n*blockRowSize, ncol = pPlusOne+2*totNumRowA+2*n)
	} else{
		lhsMatrix = matrix(0 , nrow = n*blockRowSize, ncol = pPlusOne+2*totNumRowA+n)
	}
	
	uVNames = NULL
	for(i in 1:n){
		uVNames = c(uVNames, c(paste0(paste0(rep("u", uvSize), i), 1:uvSize), paste0(paste0(rep("v", uvSize), i), 1:uvSize)))
	}	
	
	if(twoSlacks){
		varNames = c(paste0("w", 1:p), "w0", uVNames, paste0("csiPlus", 1:n), paste0("csiMinus", 1:n))
	} else{
		varNames = c(paste0("w", 1:p), "w0", uVNames, paste0("csi", 1:n))
	}
	
	
	senseVect = rep(c("<=", "<=", rep("=", pPlusOne*2)), n)
	
	coefVect = c(rep(0, pPlusOne+2*totNumRowA), costVect)

	if(twoSlacks){
		QMatrix = diag(0, pPlusOne+2*totNumRowA+2*n, pPlusOne+2*totNumRowA+2*n)
	} else{
		QMatrix = diag(0, pPlusOne+2*totNumRowA+n, pPlusOne+2*totNumRowA+n)
	}	
	QMatrix[1:p, 1:p]= diag(1, p, p)
	startUvCol = pPlusOne+1

	uIdx = vIdx = list() 
	if(twoSlacks){
		csiPlusIdx = csiMinusIdx = list()
	} else{
		csiIdx = list()
	}
	
	rhsVect = NULL
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
		
		currEpsilon = epsilonCertain
		if(uncertaintySpecialTreatment & !twoSlacks)
			if(uncertainPoints[i])
				currEpsilon = currEpsilon + extraEpsilonUncertain*uncertaintyVect[i]
		
		rhsVect = c(rhsVect, c(currEpsilon, currEpsilon, rep(0, p), 1, rep(0,p), -1))
		
		uIdx[[i]] = uvBlockCol[1:(numUvCol/2)] 
		vIdx[[i]] = uvBlockCol[(numUvCol/2+1):numUvCol]
		if(twoSlacks){
			csiPlusIdx[[i]] = pPlusOne+2*totNumRowA + 2*i-1
			csiMinusIdx[[i]] = pPlusOne+2*totNumRowA + 2*i
		} else{
			csiIdx[[i]] = pPlusOne+2*totNumRowA + i
		}
		startUvCol = startUvCol + ncol(polyList[[i]]$uvBlock)
		if(twoSlacks){
			lhsMatrix[blockRows, c(csiPlusIdx[[i]], csiMinusIdx[[i]])] = csiColBlock
		} else{
			lhsMatrix[blockRows, csiIdx[[i]]] = csiColBlock
		}
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
	varIdx$allU = unlist(varIdx$u)
	varIdx$allV = unlist(varIdx$v)
	if(twoSlacks){
		varIdx$csiPlus = csiPlusIdx
		varIdx$csiMinus = csiMinusIdx
		varIdx$allCsiPlus = unlist(varIdx$csiPlus)
		varIdx$allCsiMinus = unlist(varIdx$csiMinus)
	} else{
		varIdx$csi = csiIdx
		varIdx$allCsi = unlist(varIdx$csi)
	}
	model$varIdx = varIdx
	model$varNames = varNames
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

