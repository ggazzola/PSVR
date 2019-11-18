#!/usr/bin/env Rscript

require(MASS)
require(mice)
require(gurobi)
source("../MiscFunctions.R")
source("../ConstraintFunctions.R")
source("../PolytopeSVRFunctions.R")
source("../OneVsTwoSlackSVR.R")
source("PolytopeSVR.R")
source("../WrapperFunctions.R")
dataFolder="../Data/"

realData = F
injectMissingness = T
doMCAR = F
missingVarProp = 0.2 #####PAPERCHANGE
missingObsProp = 0.2 #####PAPERCHANGE
corVal = 0.5
# do .9, .9, X; .3 .3 X, with X = 0.3, .9
# MAR only?


#######noMissDat --> the best testing parameters for uncertain data should be chosen based on the best for certain predictions: WHY? you can still calculate validation performance according to uncertain data only...

#### are the extraEpsilonUncertain and Cuncertain parameters playing a role?

# we don't need to have a model that works well on both certain and uncertain; we need one that does better at least one of those; for the other thing
# we can use a standard model, etc.


maxGenerateDataAttempts = 20
repVect=1 # ######PAPERCHANGE?       MUST DO MORE REPEATS, the results don't seem stable
numFolds = 2#### #####PAPERCHANGE?
scaleData = T
method = "pmm" #norm, cart, rf
maxIter = 10 #####PAPERCHANGE? ################ try with small numbers of these two, to verify if imputation is possible first
numImput = 20 #####PAPERCHANGE?

AggregateTestError = mean
replaceImputedWithTrueY = F

maxUncertainDims = "all" # NULL #("all" considers the p+1 dims; NULL considers the actual max number of missing dims in all the data )

if(T){
	missingY = F # Do not modify this
}

if(!realData){
	n = 50 ######PAPERCHANGE do 120 for 10 folds is the least to have at least one miss/non miss point if missingObsPropVect = 0.1 or 0.9
	p = 5 #####PAPERCHANGE
	meanVect = rep(0,p) 
	stdVect = rep(1, p)
	trueW = 1:p
	trueW0 = p/2
	theoRsq = 0.9
	stopifnot(injectMissingness)
} else{
	#"Automobile.RData" # kept numerical variables, removed 4 obs with NA Y; has natural NAs
	#Boston.RData --boston corrected: kept numerical variables (removed boolean); has no natural NAs
	#Communities.RData -- kept all; has natural NAs
	#Ozone.RData -- removed V2, V3 (day 1-31, day of the week Mon-Sun); removed few observations with missing Y
	#Mammalsleep -- removed 'species' (useless factor with n levels), used 'bwt' as output (mammalsleep in mice)
	#Dutch -- removed 'reg' (factor); 'gen', 'phb' (ordered factors) converted to numeric;used age (only one with no NAs) as output (boys in mice) 
	nSubSelect = 100 #Inf
	realDataFileName = "Boston.RData" # 
	corVal = "irrelevant"
	theoRsq = "irrelevant"
	if(realDataFileName%in%c("Communities.RData", "Ozone.RData", "Automobile.RData", "Mammalsleep.RData", "Dutch.RData"))
		stopifnot(injectMissingness==F)
	if(realDataFileName=="Boston.RData")
		stopifnot(injectMissingness==T)
}

if(!injectMissingness){
	cat("Ignoring missingVarProp and missingObsProp (setting to 0), and doMCAR, since injectMissingness=F\n")
	missingVarProp=0	
	missingObsProp=0
	doMCAR = "irrelevant"
}

parValuesList = list( #####PAPERCHANGE?
	Ccertain=c(0, 2,  5) # add one
	Cuncertain=c(0, 2, 5) # add one
	epsilonCertain=c(0,  .5, 1) # # add one
	extraEpsilonUncertain = c(0, .5, 1) # add one
	#Ccertain=c(0, .05, .1, .5, 1, 2, 5)
	#Cuncertain=c(0, .05, .1, .5, 1, 2, 5)
	#epsilonCertain=c(0, 0.25, .5, 1) # no sense having these large if standardizing output (so to magnitude within 1 or so..)
	#extraEpsilonUncertain = c(0, 0.25, .5, 1), for the two UNCERTAIN METAPARAMETERS, GO BACK TO THE DEFINITIONS TO CHECK IF THIS SCALE IS OK
	uncertaintySpecialTreatment = T,
	linear =T
	)	
	#####PAPERCHANGE?
quantOrSdPropValues = c(0,  0.5, 1)#c(0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1) # 
	
#quantOrSdPropValues = c(0, 0.0001, 0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 1)#c(0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1) # 

errMeasureVect=c("mae", "rmse", "Maxae", "maeCert", "rmseCert", "MaxaeCert", 
	"maeUncert", "rmseUncert", "MaxaeUncert") #maeCert #maeUncert, ...
approachVect = c("doPCbb")#, "doMedian", "doNoMiss", "doSquarebbSd", "doSquarebbQuant")  ####################



for(repIdx in repVect){
	rejectSilently=F	
	set.seed(repIdx)
	GenerateData() # inefficient, because redundant with the below, but useful to do prescreening of generated data
	if(realData){
		n = nrow(doDatOut)
		p = ncol(doDatOut)-1
	}
	if(givenUp)
		stop("Couldn't generate data partitions containing at least one missing point and one non-missing point")
	
}

#innerTrainingN = floor(n*(1-1/numFolds)^2) # this "trimming" may actually be counter productive (e.g., even if n points, can obtain >n different quantilic box sizes out of quantile())
#uselessQuanOrSdPropValuesIdx = duplicated(round(innerTrainingN* quantOrSdPropValues))
#quantOrSdPropValues = quantOrSdPropValues[!uselessQuanOrSdPropValuesIdx]

cat("Data partitions containing at least one missing point and one non-missing point can be generated\n")
cat("Proceeding to actual experiments...\n\n\n")


currDate = system('date +%Y%m%d-%H%M%S', intern=T)

if(realData){
	fileNameRoot = paste0(realDataFileName)
} else{
	fileNameRoot = paste0("NormalN", n, "P", p)
}		

if(is.logical(doMCAR)){
	missMechString = ifelse(doMCAR, "MCAR", "MAR")
} else{
	stopifnot(doMCAR=="irrelevant")
	missMechString = "Irrel" 
}

fileNameRoot = paste0(fileNameRoot, "MissObs", missingObsProp, "MissVar", missingVarProp, "MissMech", missMechString, "Meth", method)

resultsFolderName = paste0(fileNameRoot, "Date", currDate)
system(paste("mkdir", resultsFolderName))
system(paste("cp *.R *sh", resultsFolderName))
system(paste("cp ../*.R", resultsFolderName))
progressFile = paste0(resultsFolderName, "/Progress.txt")
				
progressOut=paste("SAVING results in", resultsFolderName, "\n")
cat(progressOut)
write.table(progressOut, quote=F, row.names=F, col.names=F, append=T, file=progressFile)
	
totComb = length(repVect)*length(approachVect)*length(errMeasureVect)
	
progressOut=paste("STARTING a total of", totComb, "combinations at", date(), "\n")
cat(progressOut)
write.table(progressOut, quote=F, row.names=F, col.names=F, append=T, file=progressFile)
	
cnt = 1
				

if(!realData){
	fileNameRoot = paste0(fileNameRoot, "Cor", corVal, "Rsq", theoRsq, "N", n, "P", p)
} else{
	fileNameRoot = paste0(fileNameRoot, "N", nSubSelect)
}

testRes = list()
numFailedImputReps = 0
maxNumReps = length(repVect)*2
times <- NULL
times$start = Sys.time()
for(repIdx in repVect){
	testRes[[repIdx]]=list()
	rejectSilently=T	
	set.seed(repIdx)
	GenerateData()
	times$generatedata = Sys.time()
	MultiplyImpute()
	times$multiplyimpute = Sys.time()
	
	if(imputationsFailed){
		#imputationsFailed is declared by MultiplyImpute()
		numFailedImputReps = numFailedImputReps+1
		if(numFailedImputReps>maxNumReps)
			stop("Too many failed imputations\n")
			
		remainingRepNum = max(repVect)-repIdx 
		repVect = (repIdx+1):(repIdx+1+remainingRepNum) # BUGBUGBUG WON'T WORK --> repVect above isn't redeclared within the loop
		
	}
	gc()

	for(approach in approachVect){
		appShort = sub("do", "", approach)
		appShort = sub("Square", "Sq", appShort)
		fileName = paste0(fileNameRoot, "Appr", appShort, "Date", currDate, ".RData", sep="")
			
		progressOut = paste("STARTING cross-validation of", fileName, "Rep", repIdx, "at", date(), "\n")
		cat(progressOut)
		write.table(progressOut, quote=F, row.names=F, col.names=F, append=T, file=progressFile)
		
		CalculateValidationErrors() # errors calculated with all possible error measures
		times$calculatevalidationerrors = Sys.time()
				
		gc()
		
		
		for(errMeasure in errMeasureVect){

			if(is.null(testRes[[repIdx]][[errMeasure]]))
				testRes[[repIdx]][[errMeasure]]=list()

			testRes[[repIdx]][[errMeasure]][[approach]] = list()
			#testRes[[repIdx]][[errMeasure]][[approach]]$valid = doErrorFoldOutInnerList #to save memory
			SetUpTest()
			testRes[[repIdx]][[errMeasure]][[approach]]$testSetup = getTrainResReadyForTest 
			
			testRes[[repIdx]][[errMeasure]][[approach]]$testSetup$currTest = 
				testRes[[repIdx]][[errMeasure]][[approach]]$testSetup$currTest = 
				testRes[[repIdx]][[errMeasure]][[approach]]$testSetup$currTestImput = 
				testRes[[repIdx]][[errMeasure]][[approach]]$testSetup$currTrain = 
				testRes[[repIdx]][[errMeasure]][[approach]]$testSetup$currTrainImput = NULL #to save memory
				 
			progressOut = paste("STARTING testing on", fileName, "Rep", repIdx, "Error Measure", errMeasure, "at", date(), "\n")
			cat(progressOut)
			write.table(progressOut, quote=F, row.names=F, col.names=F, append=T, file=progressFile)
		
			CalculateTestErrors() # choose best model based on errMeasure, train new model on outer training data, and test it
			gc()
			testRes[[repIdx]][[errMeasure]][[approach]]$testErrors = errVect 
			testRes[[repIdx]][[errMeasure]][[approach]]$testModels = modelList 

			stopifnot(length(errVect)==numFolds)
			testRes[[repIdx]][[errMeasure]][[approach]]$testErrorsAggregate=AggregateTestError(errVect) 
			times$aggregatetesterrors[[errMeasure]] = Sys.time()				
			
			testRes[[repIdx]][[errMeasure]][[approach]]$testErrorsSd=sd(errVect)

			progressOut=paste("Combination", cnt, "out of", totComb, "DONE at", date(), "\n")
			cat(progressOut)
			write.table(progressOut, quote=F, row.names=F, col.names=F, append=T, file=progressFile)
			save.image(file=paste0(resultsFolderName, "/", fileName)) # --> note: we save ONLY the inner cross-validation results of the last repeat
																		# since they get overwritten	

			cnt = cnt+1
			
		}
		times$calculatetesterrors = Sys.time()				
		
	}
}


numHyperPar = prod(sapply(parValuesList, length))*length(quantOrSdPropValues)
numErrMeasure = length(errMeasureVect)
save(times,missingVarProp,missingObsProp,repVect,numFolds,maxIter,numImput, n, p, approachVect, numHyperPar, numErrMeasure, 
	file=paste0(resultsFolderName, "/", paste0("Essentials", fileName)))


progressOut=paste("All results collected and saved in", resultsFolderName, "\n")
cat(progressOut)
write.table(progressOut, quote=F, row.names=F, col.names=F, append=T, file=progressFile)
	
if(F){
       numErrMeas = length(errMeasureVect)
       numApproach = length(approachVect)
       numRep = length(testRes)
       meanMat = seMat = matrix(, nrow=numApproach , ncol=numErrMeas)
       colnames(meanMat) = colnames(seMat) = errMeasureVect
       rownames(meanMat) = rownames(seMat) = approachVect
       
       for(j in 1:numErrMeas){
               errMeasure = errMeasureVect[j]
               for(i in 1:numApproach){
                       approach = approachVect[i]
                       res = NULL
                       for(k in 1:numRep)
                               res = c(res, testRes[[k]][[errMeasure]][[approach]]$testErrorsAggregate)
                       stopifnot(length(res)==numRep)
                       meanMat[i,j] = mean(res)
                       seMat[i,j] = sd(res)/sqrt(numRep)
               }
       }
       
}


if(F){
	
	#load("") must be the correct PC bb file! # note -- this will plot the validation error, calculated as a mean across the 5 inner folds of each outer fold i (as given by doErrorFoldOutInnerList[[i]]),  and a subsequent grand mean over all outer folds i --> note that this grand mean is such that the min value of the curve across all hyper-parameters may differ (it wouldn't if we plotted
	# one curve per outer fold separately)
	res = PerformanceByParameterValue(doErrorFoldOutInnerList)	
	PlotBestPerformanceByParameterValue (res, foldIdx=1:length(res), errMeasureName=errMeasureVect[1], parName= "all")
	
}

if(F){
	#load("") must be the last file (SquareQuant?)
	for (errMeasure in errMeasureVect){
		for(approach in approachVect){
			cat(approach, errMeasure, testRes[[1]][[errMeasure]][[approach]]$testErrorsAggregate, "\n")
		}
	}
}


if(F){
	testResMedian = testRes
	for (errMeasure in errMeasureVect){
		for(approach in approachVect){
			testResMedian[[1]][[errMeasure]][[approach]]$testErrorsAggregate = median(testResMedian[[1]][[errMeasure]][[approach]]$testErrors)
		}
	}
}