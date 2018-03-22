OneSlackSVR = function(dat, C = 10, epsilon = 0){
	stopifnot(is.matrix(dat))
	#dat should be a ***p+1*** matrix!!!
	n=nrow(dat)
	p=ncol(dat)-1
	datX = dat[, 1:(ncol(dat)-1)]
	datY = dat[, ncol(dat)]
	rowOne = 1:n
	rowTwo = (n+1):(2*n)

	matOne = matrix(0, nrow=2*n, ncol=p)
	matOne[rowOne, ] = datX
	matOne[rowTwo,] = -datX

	matTwo = matrix(0, nrow=2*n, ncol=1)
	matTwo[rowOne,]=1
	matTwo[rowTwo,]=-1

	matThree  = matrix(0, nrow=2*n, ncol = n)
	matThree[rowOne, ] = matThree[rowTwo, ] =  diag(-1, n, n)

	lhs = cbind(matOne, matTwo, matThree)

	rhs = rep(0, 2*n)
	rhs[rowOne] = datY+epsilon
	rhs[rowTwo] = -datY+epsilon

	sense =rep("<=", 2*n)

	Q = diag(0, p+1+n, p+1+n)
	Q[1:(p), 1:(p)] = diag(1, p, p)
	coeff = rep(0, p+1+n)
	coeff[(p+2):length(coeff)] = C
	lowerBound = c(rep(-Inf, p+1), rep(0, n))


	model = list()
	model$A = lhs
	model$rhs = rhs
	model$sense = sense
	model$Q = Q
	model$obj = coeff
	model$lb = lowerBound
	varIdx = list()
	varIdx$w = 1:p
	varIdx$w0 = p+1
	varIdx$csi = list()
	for(i in 1:n)
		varIdx$csi[[i]] = p+1+i
	model$varIdx = varIdx
	
	return(model)

}

TwoSlackSVR = function(dat, C=10, epsilon = 0) {
	#dat should be a ***p+1*** matrix!!!
	stopifnot(is.matrix(dat))
	datX = dat[, 1:(ncol(dat)-1)]
	datY = dat[, ncol(dat)]
	n=nrow(dat)
	p=ncol(dat)-1
	rowOne = 1:n
	rowTwo = (n+1):(2*n)

	matOne = matrix(0, nrow=2*n, ncol=p)
	matOne[rowOne, ] = datX
	matOne[rowTwo,] = -datX

	matTwo = matrix(0, nrow=2*n, ncol=1)
	matTwo[rowOne,]=1
	matTwo[rowTwo,]=-1
	matThree = matFour = matrix(0, nrow=2*n, ncol = n)
	matThree[rowOne, ]  = diag(-1, n, n)
	matFour[rowTwo,] = diag(-1, n, n)

	lhs = cbind(matOne, matTwo, matThree, matFour)

	rhs = rep(0, 2*n)
	rhs[rowOne] = datY+epsilon
	rhs[rowTwo] = -datY+epsilon

	sense =rep("<=", 2*n)

	Q = diag(0, p+1+2*n, p+1+2*n)
	Q[1:(p), 1:(p)] = diag(1, p, p)
	coeff = rep(0, p+1+2*n)
	coeff[(p+2):length(coeff)] = C
	lowerBound = c(rep(-Inf, p+1), rep(0, 2*n))

	model = list()
	model$A = lhs
	model$rhs = rhs
	model$sense = sense
	model$Q = Q
	model$obj = coeff
	model$lb = lowerBound
	
	varIdx = list()
	varIdx$w = 1:p
	varIdx$w0 = p+1
	varIdx$csiPlus = varIdx$csiMinus = list()
	for(i in 1:n){
		varIdx$csiPlus[[i]] = p+1+i
		varIdx$csiMinus[[i]] = p+1+n+i
	}	
	model$varIdx = varIdx
	return(model)	
}

OneSlackSVROLD = function(dat, C = 2){
	# has non-negativity constraints for eps but no lb component
	# so also w and w_0 are (incorrectly) assumed forced to be non-negative
	
	n=nrow(dat)
	p=ncol(dat)-1
	datX = dat[, 1:(ncol(dat)-1)]
	datY = dat[, ncol(dat)]
	rowOne = 1:n
	rowTwo = (n+1):(2*n)
	rowThree = ((2*n)+1):(3*n)

	matOne = matrix(0, nrow=3*n, ncol=p)
	matOne[rowOne, ] = datX
	matOne[rowTwo,] = -datX

	matTwo = matrix(0, nrow=3*n, ncol=1)
	matTwo[rowOne,]=1
	matTwo[rowTwo,]=-1

	matThree  = matrix(0, nrow=3*n, ncol = n)
	matThree[rowOne, ] = matThree[rowTwo, ] = matThree[rowThree, ] =  diag(-1, n, n)

	lhs = cbind(matOne, matTwo, matThree)

	rhs = rep(0, 3*n)
	rhs[rowOne] = datY
	rhs[rowTwo] = -datY

	sense =rep("<=", 3*n)

	Q = diag(0, p+1+n, p+1+n)
	Q[1:(p), 1:(p)] = diag(1, p, p)
	coeff = rep(0, p+1+n)
	coeff[(p+2):length(coeff)] = C


	model = list()
	model$A = lhs
	model$rhs = rhs
	model$sense = sense
	model$Q = Q
	model$obj = coeff
	
	return(model)

}

TwoSlackSVROLD = function(dat, C=2) {
	# has non-negativity constraints for eps_+ and eps_-, but no lb component
	# so also w and w_0 are (incorrectly) assumed forced to be non-negative
	datX = dat[, 1:(ncol(dat)-1)]
	datY = dat[, ncol(dat)]
	n=nrow(dat)
	p=ncol(dat)-1
	rowOne = 1:n
	rowTwo = (n+1):(2*n)
	rowThree = ((2*n)+1):(3*n)
	rowFour = ((3*n)+1):(4*n)

	matOne = matrix(0, nrow=4*n, ncol=p)
	matOne[rowOne, ] = datX
	matOne[rowTwo,] = -datX

	matTwo = matrix(0, nrow=4*n, ncol=1)
	matTwo[rowOne,]=1
	matTwo[rowTwo,]=-1

	matThree = matFour = matrix(0, nrow=4*n, ncol = n)
	matThree[rowOne, ] = matThree[rowThree, ] = diag(-1, n, n)

	matFour[rowTwo,] = matFour[rowFour,] = diag(-1, n, n)

	lhs = cbind(matOne, matTwo, matThree, matFour)

	rhs = rep(0, 4*n)
	rhs[rowOne] = datY
	rhs[rowTwo] = -datY

	sense =rep("<=", 4*n)

	Q = diag(0, p+1+2*n, p+1+2*n)
	Q[1:(p), 1:(p)] = diag(1, p, p)
	coeff = rep(0, p+1+2*n)
	coeff[(p+2):length(coeff)] = C


	model = list()
	model$A = lhs
	model$rhs = rhs
	model$sense = sense
	model$Q = Q
	model$obj = coeff
	return(model)	
}
