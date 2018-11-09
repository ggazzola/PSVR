# to compare coefficients of different models

simulatedRepNum = 40

numRep = length(testRes)
stopifnot(simulatedRepNum>=numRep)
if(numRep==1 & simulatedRepNum>1)
	stop("The actual number of repeats is 1. Aborting.\n")

if(simulatedRepNum>numRep)
	cat("WARNING: simulating", simulatedRepNum, "repeats from", numRep, "actual repeats.\n")

roundingFactor = 3

errMeasureVectToShow = c("rmse", "mae", "quantEightAe", "Maxae",  
						"rmseCert", "maeCert", "quantEightAeCert", "MaxaeCert",
						"rmseUncert", "maeUncert", "quantEightAeUncert",  "MaxaeUncert")
						
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
	
approachVectToShow = c("doPCbb", "doMedian" , "doNoMiss", "doSquarebbSd", "doSquarebbQuant")

approachVectToShowRenamed = approachVectToShow

approachVectToShowRenamed[approachVectToShowRenamed=="doPCbb"] = "PSVR"
approachVectToShowRenamed[approachVectToShowRenamed=="doMedian"] = "MSVR"
approachVectToShowRenamed[approachVectToShowRenamed=="doNoMiss"] = "CSVR"
approachVectToShowRenamed[approachVectToShowRenamed=="doSquarebbSd"] = "CarrizosaSD"
approachVectToShowRenamed[approachVectToShowRenamed=="doSquarebbQuant"] = "CarrizosaQ"



for(errMeasure in errMeasureVectToShow){
	matMeanTmp = matSETmp = matrix(, ncol=length(approachVectToShow), nrow = 11) # note: will be overwritten at each cycle of errMeasure!
	cntApproach = 1
	for(approach in approachVectToShow){
		innerMatTmp = matrix(, ncol=length(repVect), nrow = 11)
		
		for(repIdx in repVect){
			currDat = testRes[[repIdx]][[errMeasure]][[approach]]$testModels
			wVect = rep(0, 10)
			w0Val = 0
			for(innerRep in 1:5){
				wVect = wVect + currDat[[innerRep]]$w/10
				w0Val = w0Val + currDat[[innerRep]]$w0/10
			}
			innerMatTmp[, repIdx] = c(w0Val, wVect)
		}
		
		matMeanTmp[, cntApproach] = apply(innerMatTmp, 1, mean)
		matSETmp[, cntApproach] = apply(innerMatTmp, 1, sd)/sqrt(simulatedRepNum)
		#write.csv ...
		cntApproach = cntApproach + 1
	}
}

#PSVR and CSVR usually quite similar, if different, no discernible patterns (sometimes PSVR uses bigger intercept in absolute value)

#True coeff after scaling (but if there are correlations, SVM may use some predictors as substitutes for others, etc.)
#  (Intercept)            X1            X2            X3            X4            X5            X6            X7            X8            X9           X10 
#-2.534677e-17  2.879704e-02  5.796170e-02  8.759637e-02  1.148619e-01  1.466080e-01  1.740065e-01  2.019095e-01  2.296192e-01  2.601165e-01  2.908092e-01 