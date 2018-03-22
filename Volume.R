Volume = function(A, a, path="./PolyVest/", nPts = 1600, removeFiles = T){
	startingWd = getwd()
	stopifnot(is.matrix(A))
	stopifnot(is.vector(a))
	stopifnot(is.numeric(A))
	stopifnot(is.numeric(a))
	stopifnot(nrow(A)==length(a))
	
	M = nrow(A)
	N = ncol(A)
	aA = cbind(a, A)
	input = c(M, N, as.vector(t(aA)))
	inputFileName = "volumeIn.tmp"
	outputFileName = "volumeOut.tmp"
	setwd(path)
	write.table(input, inputFileName, quote=F, row.names=F, col.names=F)
	polyVestCall = paste("PolyVest", inputFileName, nPts, outputFileName)
	print(polyVestCall)
	system(polyVestCall)
	volume = as.numeric(read.table(outputFileName))
	if(removeFiles){
		system(paste("rm -rf",inputFileName))
		system(paste("rm -rf",outputFileName))
	}
	setwd(startingWd)	
	return(volume)
}
