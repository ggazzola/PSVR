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
doMCAR = T 
missingVarProp = 0.9 
missingObsProp = 0.9
corVal = 0.9
# do .9, .9, X; .3 .3 X, with X = 0.3, .9



#######noMissDat --> the best testing parameters for uncertain data should be chosen based on the best for certain predictions: WHY? you can still calculate validation performance according to uncertain data only...

#### are the extraEpsilonUncertain and Cuncertain parameters playing a role?

# we don't need to have a model that works well on both certain and uncertain; we need one that does better at least one of those; for the other thing
# we can use a standard model, etc.


maxGenerateDataAttempts = 20
numFolds = 5#####################
scaleData = T
method = "pmm" #norm, cart, rf
maxIter = 20 ################ try with small numbers of these two, to verify if imputation is possible first
numImput = 40 ###############

AggregateTestError = mean
replaceImputedWithTrueY = F

maxUncertainDims = "all" # NULL #("all" considers the p+1 dims; NULL considers the actual max number of missing dims in all the data )

if(T){
	missingY = F # Do not modify this
}

if(!realData){
	n = 100 # do 120 for 10 folds is the least to have at least one miss/non miss point if missingObsPropVect = 0.1 or 0.9
	p = 10
	meanVect = rep(0,p) 
	stdVect = rep(1, p)
	trueW = 1:p
	trueW0 = p/2
	theoRsq = 0.9
	stopifnot(injectMissingness)
} else{
	#"Automobile.RData" # kept numerical variables, removed 4 obs with NA Y; has natural NAs
	#Boston.RData --boston corrected: kept numerical variables (removed boolean); has no natural NAs
	realDataFileName = "Automobile.RData" # kept
	corVal = "irrelevant"
	theoRsq = "irrelevant" 
}

if(!injectMissingness){
	cat("Ignoring missingVarProp and missingObsProp (setting to 0), and doMCAR, since injectMissingness=F\n")
	missingVarProp=0	
	missingObsProp=0
	doMCAR = "irrelevant"
}

parValuesList = list(
	Ccertain=c(0, .05, .1, .5, 1, 2, 5),#c(0,10^(-2:1)),   ##################
	Cuncertain=c(0, .05, .1, .5, 1, 2, 5),#c(0,10^(-2:1)), ##################
	epsilonCertain=c(0, 0.25, .5, 1),#c(0,10^(-2:1)),  ################## no sense having these large if standardizing output (so to magnitude within 1 or so..)
	extraEpsilonUncertain = c(0, 0.25, .5, 1),# c(0,10^(-2:1)),  ################# for the two UNCERTAIN METAPARAMETERS, GO BACK TO THE DEFINITIONS TO CHECK IF THIS SCALE IS OK
	uncertaintySpecialTreatment = T,
	linear =T
	)	


quantOrSdPropValues = c(0, 0.0001, 0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 1)#c(0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1) # 
#errMeasureVect=c("mae", "rmse", "Maxae", "cor", "quantNineAe", "quantEightAe", "quantSevenAe",
#"maeCert", "rmseCert", "MaxaeCert", "quantNineAeCert", "quantEightAeCert", "quantSevenAeCert", "corCert",
#"maeUncert", "rmseUncert", "MaxaeUncert", "quantNineAeUncert", "quantEightAeUncert", "quantSevenAeUncert", "corUncert") #maeCert #maeUncert, ...
errMeasureVect=c("mae", "rmse", "Maxae", "cor",  "quantEightAe", "maeCert", "rmseCert", "MaxaeCert",  "corCert", "quantEightAeCert",
	"maeUncert", "rmseUncert", "MaxaeUncert", "corUncert", "quantEightAeUncert") #maeCert #maeUncert, ...
#approachVect = c("doPCbb", "doSquarebbSd", "doSquarebbQuant", "doMedian", "doNoMiss")  ####################
approachVect = c("doPCbb", "doMedian", "doNoMiss", "doSquarebbSd", "doSquarebbQuant")  ####################


repVect=1:4 # # MUST DO MORE REPEATS, the results don't seem stable

for(repIdx in repVect){
	rejectSilently=F	
	set.seed(repIdx)
	GenerateData() # inefficient, because redundant with the below, but useful to do prescreening of generated data
	if(realData){
		n = nrow(dat)
		p = ncol(dat)-1
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
	fileNameRoot = paste0(fileNameRoot, "Cor", corVal, "Rsq", theoRsq)
}		

testRes = list()

for(repIdx in repVect){
	testRes[[repIdx]]=list()
	rejectSilently=T	
	set.seed(repIdx)
	GenerateData()
	MultiplyImpute()
	gc()

	for(approach in approachVect){
		appShort = sub("do", "", approach)
		appShort = sub("Square", "Sq", appShort)
		fileName = paste0(fileNameRoot, "Appr", appShort, "Date", currDate, ".RData", sep="")
			
		progressOut = paste("STARTING cross-validation of", fileName, "Rep", repIdx, "at", date(), "\n")
		cat(progressOut)
		write.table(progressOut, quote=F, row.names=F, col.names=F, append=T, file=progressFile)
		
		CalculateValidationErrors() # errors calculated with all possible error measures
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

			stopifnot(length(errVect)==numFolds)
			testRes[[repIdx]][[errMeasure]][[approach]]$testErrorsAggregate=AggregateTestError(errVect)
			testRes[[repIdx]][[errMeasure]][[approach]]$testErrorsSd=sd(errVect)

			progressOut=paste("Combination", cnt, "out of", totComb, "DONE at", date(), "\n")
			cat(progressOut)
			write.table(progressOut, quote=F, row.names=F, col.names=F, append=T, file=progressFile)
			save.image(file=paste0(resultsFolderName, "/", fileName)) # save only testRes?

			cnt = cnt+1
			
		}
	}
}
				
				
				



progressOut=paste("All results collected and saved in", resultsFolderName, "\n")
cat(progressOut)
write.table(progressOut, quote=F, row.names=F, col.names=F, append=T, file=progressFile)
	

if(F){
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