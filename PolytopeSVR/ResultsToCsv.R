#!/usr/bin/env Rscript

nameSim = commandArgs(trailingOnly=TRUE)
if (length(nameSim)!=1) {
  stop("Problem with provided argument", call.=FALSE)
}

load(nameSim)

res = testRes[[1]] 
approachVect = c("doPCbb" ,  "doMedian", "doNoMiss")

approachVectNew = c("PolySVR", "MedSVR", "ComplSVR")


#errMeasureVect = c("mae", "rmse", "Maxae", "quantNineAe", "quantEightAe", "quantSevenAe",
#				"maeCert", "rmseCert", "MaxaeCert", "quantNineAeCert", "quantEightAeCert", "quantSevenAeCert",
#				"maeUncert", "rmseUncert", "MaxaeUncert", "quantNineAeUncert", "quantEightAeUncert", "quantSevenAeUncert") #maeCert #maeUncert, ...

errMeasureVect = c(
				"maeCert", "rmseCert", "MaxaeCert", "quantNineAeCert", "quantEightAeCert", "quantSevenAeCert",
				"maeUncert", "rmseUncert", "MaxaeUncert", "quantNineAeUncert", "quantEightAeUncert", "quantSevenAeUncert") #maeCert #maeUncert, ...

errMeasureVectNew=errMeasureVect
#errMeasureVect = c(  "mae", "rmse", "Maxae",
#"maeCert",            "rmseCert"   ,        "MaxaeCert"  ,         
# "maeUncert"    ,      "rmseUncert"      ,   "MaxaeUncert" )   
 
 errMeasureVectNew = c(  "\\bs{e_{ma}^c}",            "\\bs{e_rms^c}"   ,        "\\bs{e_{maxa}^c}"  ,     
 "\\bs{e_{q.9a}^c}",     "\\bs{e_{q.8a}^c}",  "\\bs{e_{q.7a}^c}",
"\\bs{e_{ma}^u}",            "\\bs{e_rms^u}"   ,        "\\bs{e_{maxa}^u}",
"\\bs{e_{q.9a}^u}",     "\\bs{e_{q.8a}^u}",  "\\bs{e_{q.7a}^u}"
 )      
 
resMat = resMatSE =  matrix(, nrow=length(approachVect), ncol=length(errMeasureVect))
 
for(j in 1:length(errMeasureVect)) {
	for(i in 1:length(approachVect)) {
		errMeasure = errMeasureVect[j]
		approach = approachVect[i]
		resMat[i,j]=round(res[[errMeasure]][[approach]]$testErrorsAggregate, 2)
		resMatSE[i,j]=round(sd(res[[errMeasure]][[approach]]$testErrors)/sqrt(length(res[[j]][[i]]$testErrors)), 2) # redo, SE should be calculated across repeats, for comparability with other results
	}
}		

resMat = as.data.frame(resMat)

resMat = cbind(approachVectNew, resMat)

colnames(resMat)=c("\\text{\\bfseries{Method}}", errMeasureVectNew)

write.csv(resMat, file=paste(nameSim, "Res.csv", sep=""), quote=F, row.names=F)

cat(nameSim, "done\n")

