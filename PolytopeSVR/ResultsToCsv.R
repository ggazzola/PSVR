#!/usr/bin/env Rscript

nameSim = commandArgs(trailingOnly=TRUE)
if (length(nameSim)!=1) {
  stop("Problem with provided argument", call.=FALSE)
}

load(nameSim)

approachVectWanted = c("doPCbb", "doSquarebbSd", "doSquarebbQuant", "doNoMiss","doMedian") 
stopifnot(all(approachVectWanted%in%approachVect))
approachVectNew = c("PSVR", "CarrizosaSD", "CarrizosaQ", "CSVR", "MSVR")

errMeasureVectWanted=c("mae", "rmse", "Maxae",  "quantEightAe", "maeCert", "rmseCert", "MaxaeCert", "quantEightAeCert",
	"maeUncert", "rmseUncert", "MaxaeUncert", "quantEightAeUncert") #maeCert #maeUncert, ...
	stopifnot(all(errMeasureVectWanted%in%errMeasureVect))
errMeasureVectNew = c(
	 "\\bs{e_{ma}^a}",            "\\bs{e_rms^a}"   ,        "\\bs{e_{maxa}^a}"  ,     
	    "\\bs{e_{q.8a}^a}",
	 "\\bs{e_{ma}^c}",            "\\bs{e_rms^c}"   ,        "\\bs{e_{maxa}^c}"  ,     
   "\\bs{e_{q.8a}^c}", "\\bs{e_{ma}^u}",            "\\bs{e_rms^u}"   ,        "\\bs{e_{maxa}^u}",
 "\\bs{e_{q.8a}^u}" )      

numRep = length(testRes)

resMatList =  list()
for(nRep in 1:numRep){ 
	resMat =  matrix(, nrow=length(approachVectWanted), ncol=length(errMeasureVectWanted))
 
	for(j in 1:length(errMeasureVectWanted)) {
		for(i in 1:length(approachVectWanted)) {
			errMeasure = errMeasureVectWanted[j]
			approach = approachVectWanted[i]
			resMat[i,j]=testRes[[nRep]][[errMeasure]][[approach]]$testErrorsAggregate
		}
	}	
	resMatList[[nRep]] = resMat

}

meanResMat = seResMat = matrix(, nrow=length(approachVectWanted), ncol=length(errMeasureVectWanted))


for(j in 1:length(errMeasureVectWanted)) {
	for(i in 1:length(approachVectWanted)) {
		vect = NULL
		for(nRep in 1:numRep){
			vect = c(vect, resMatList[[nRep]][i,j])
		}	
		meanResMat[i,j] = round(mean(vect),3)
		seResMat[i,j] = round(sd(vect)/sqrt(numRep),3)
	}
}	

meanResMat = as.data.frame(meanResMat)
meanResMat = cbind(approachVectNew, meanResMat)

seResMat = as.data.frame(seResMat)
seResMat = cbind(approachVectNew, seResMat)

colnames(meanResMat)=c("\\text{\\bfseries{Method}}", errMeasureVectNew)
write.csv(meanResMat, file=paste(nameSim, "MeanRes.csv", sep=""), quote=F, row.names=F)

colnames(seResMat)=c("\\text{\\bfseries{Method}}", errMeasureVectNew)
write.csv(seResMat, file=paste(nameSim, "SeRes.csv", sep=""), quote=F, row.names=F)

newMeanMat = matrix(, nrow=nrow(meanResMat), ncol=ncol(meanResMat))

for(i in 1:nrow(meanResMat)){
	for(j in 2:ncol(meanResMat)){
		newMeanMat[i, j] = paste(meanResMat[i,j], seResMat[i,j], sep="\\pm")
	}
	currMeth = meanResMat[i,1] 
	currMeth = paste0("\\text{\\textsf{", currMeth, "}}")
	newMeanMat[i,1]  = currMeth
}	

colnames(newMeanMat)=c("\\text{\\bfseries{Method}}", errMeasureVectNew)

write.csv(newMeanMat, file=paste(nameSim, "MeanSeRes.csv", sep=""), quote=F, row.names=F)
cat(nameSim, "done\n")
