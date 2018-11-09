#load results .RData file first

source("MiscFunctions.R")


errMeasureVectToShow = c("rmse", "mae", "quantEightAe", "Maxae",  
						"rmseCert", "maeCert", "quantEightAeCert", "MaxaeCert",
						"rmseUncert", "maeUncert", "quantEightAeUncert",  "MaxaeUncert")
						
performanceByParameterValueOut = PerformanceByParameterValue(doErrorFoldOutInnerList)

for(foldIdx in 1:5){
	system(paste("mkdir", paste0("Fold", foldIdx)))
	PlotBestPerformanceByParameterValueOverlap(performanceByParameterValueOut, foldIdx=foldIdx)
	system(paste("mv *pdf", paste0("Fold", foldIdx)))	
	
}	





for(foldIdx in 1:5){
	system(paste("mkdir", paste0("Fold", foldIdx)))
	for(errMeasureName in errMeasureVectToShow){
		PlotBestPerformanceByParameterValue(performanceByParameterValueOut, foldIdx=foldIdx, errMeasureName=errMeasureName, parName= "all")
		system(paste("mv *pdf", paste0("Fold", foldIdx)))	
	}
}	
