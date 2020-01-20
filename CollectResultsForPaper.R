#!/usr/bin/env Rscript

fileList = list.files()
fileName = fileList[grep("SqbbQuant", fileList)] # will load this file, which is supposed to contain results from all methods if SqbbQuant is the last method in the loop
stopifnot(length(fileName)==1)
cat("Loading", fileName, "...\n")

load(fileName)

numRep = length(testRes)

simulatedRepNum = 10

stopifnot(simulatedRepNum>=numRep)
if(numRep==1 & simulatedRepNum>1)
	stop("The actual number of repeats is 1. Aborting.\n")

if(simulatedRepNum>numRep)
	cat("WARNING: simulating", simulatedRepNum, "repeats from", numRep, "actual repeats.\n")

roundingFactor = 3
errMeasureVectToShow = c("rmse", "mae", "Maxae",  
						"rmseCert", "maeCert",  "MaxaeCert",
						"rmseUncert", "maeUncert",  "MaxaeUncert")
						
errMeasureVectToShowRenamed = errMeasureVectToShow

errMeasureVectToShowRenamed[errMeasureVectToShowRenamed=="rmse"] = "\\bs{\\bar{e}_{rms}}"
errMeasureVectToShowRenamed[errMeasureVectToShowRenamed=="mae"] = "\\bs{\\bar{e}_{ma}}"
errMeasureVectToShowRenamed[errMeasureVectToShowRenamed=="quantEightAe"] = "\\bs{\\bar{e}_{q.8a}}"
errMeasureVectToShowRenamed[errMeasureVectToShowRenamed=="Maxae"] = "\\bs{\\bar{e}_{maxa}}"
errMeasureVectToShowRenamed[errMeasureVectToShowRenamed=="rmseCert"] = "\\bs{\\bar{e}_{rms}^c}"
errMeasureVectToShowRenamed[errMeasureVectToShowRenamed=="maeCert"] = "\\bs{\\bar{e}_{ma}^c}"
errMeasureVectToShowRenamed[errMeasureVectToShowRenamed=="quantEightAeCert"] = "\\bs{\\bar{e}_{q.8a}^c}"
errMeasureVectToShowRenamed[errMeasureVectToShowRenamed=="MaxaeCert"] = "\\bs{\\bar{e}_{maxa}^c}"
errMeasureVectToShowRenamed[errMeasureVectToShowRenamed=="rmseUncert"] = "\\bs{\\bar{e}_{rms}^u}"
errMeasureVectToShowRenamed[errMeasureVectToShowRenamed=="maeUncert"] = "\\bs{\\bar{e}_{ma}^u}"
errMeasureVectToShowRenamed[errMeasureVectToShowRenamed=="quantEightAeUncert"] = "\\bs{\\bar{e}_{q.8a}^u}"
errMeasureVectToShowRenamed[errMeasureVectToShowRenamed=="MaxaeUncert"] = "\\bs{\\bar{e}_{maxa}^u}"


approachVectToShow = c("doPCbb", "doNoMiss", "doSquarebbSd", "doSquarebbQuant")

approachVectToShowRenamed = approachVectToShow

approachVectToShowRenamed[approachVectToShowRenamed=="doPCbb"] = "PSVR"
approachVectToShowRenamed[approachVectToShowRenamed=="doMedian"] = "MSVR"
approachVectToShowRenamed[approachVectToShowRenamed=="doNoMiss"] = "CSVR"
approachVectToShowRenamed[approachVectToShowRenamed=="doSquarebbSd"] = "BoxSD"
approachVectToShowRenamed[approachVectToShowRenamed=="doSquarebbQuant"] = "BoxQ"


numErrMeas = length(errMeasureVectToShow)
numApproach = length(approachVectToShow)
meanMat = seMat = matrix(, nrow=numApproach , ncol=numErrMeas)

colnames(meanMat) = colnames(seMat) = errMeasureVectToShowRenamed
rownames(meanMat) = rownames(seMat) = approachVectToShowRenamed

for(j in 1:numErrMeas){
       errMeasure = errMeasureVectToShow[j]
       for(i in 1:numApproach){
               approach = approachVectToShow[i]
               res = NULL
               for(k in 1:numRep)
                       res = c(res, testRes[[k]][[errMeasure]][[approach]]$testErrorsAggregate)
               stopifnot(length(res)==numRep)
               meanMat[i,j] = round(mean(res), roundingFactor)
               seMat[i,j] = sd(res)/sqrt(simulatedRepNum) # CAREFUL WITH THIS!!!!!
       }
}

meanMat <- t(meanMat)
approachVectToShowRenamed2=approachVectToShowRenamed

for(i in 1:nrow(meanMat))
	approachVectToShowRenamed2[i] = paste0("\\text{\\textsf{", approachVectToShowRenamed[i], "}}")

newMeanMat = meanMat

rownames(newMeanMat) = approachVectToShowRenamed2

for(i in 1:nrow(meanMat)){
	for(j in 1:ncol(meanMat)){
		newMeanMat[i, j] = paste(round(meanMat[i,j], roundingFactor), round(seMat[i,j], roundingFactor), sep="\\pm")
	}
}	

newMeanMat = t(newMeanMat)

newMeanMat = cbind(rownames(newMeanMat), newMeanMat)
row.names(newMeanMat) = NULL
colnames(newMeanMat)[1]="\\text{\\bfseries{Error Measure}}"


#finalFileName = paste0("NormalMissObs",missingObsProp, "MissVar", missingVarProp, "Cor", corVal)
 finalFileName = resultsFolderName ##### TYPICALLY, THIS IS THE LINE TO USE

write.csv(meanMat, file=paste0(finalFileName, "Mean.csv"), quote=F)
write.csv(seMat, file=paste0(finalFileName, "SE.csv"), quote=F)
write.csv(newMeanMat, file=paste0(finalFileName, "MeanSE.csv"), quote=F, row.names=F)



