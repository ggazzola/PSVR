require(MASS)
require(mice)
require(gurobi)
source("../MiscFunctions.R")
source("../ConstraintFunctions.R")
source("../PolytopeSVRFunctions.R")
source("../OneVsTwoSlackSVR.R")
source("PolytopeSVR.R")
source("../WrapperFunctions.R")
#source("../Volume.R")

#######MissObs0.9Var0.25Cor0.9.RData: why are results not as good?

#######noMissDat --> the best testing parameters for uncertain data should be chosen based on the best for certain predictions: WHY? you can still calculate validation performance according to uncertain data only...

#### new un-standardization of prediction working ok?
#### are the extraEpsilonUncertain and Cuncertain parameters playing a role?

# Miss prop should be applied separately to train/valid/test, so that they all have that proportion; else, e.g., if n is small and missingObsProp is large can easily not have any certain observations;change back average=median to average=mean in DoExtractErrMat after that  --> too complicated, instead, just checking that there is at least one fully available obs, and one obs with some missing value in each train/valid/test

# we don't need to have a model that works well on both certain and uncertain; we need one that does better at least one of those; for the other thing
# we can use a standard model, etc.

maxGenerateDataAttempts = 10
numFolds = 5#####################
scaleData = T
method = "pmm"
maxIter = 20 ################ try with small numbers of these two, to verify if imputation is possible first
numImput = 50 ###############

n = 100
p = 10
meanVect = rep(0,p) 
stdVect = rep(1, p)
	
trueW = 1:p
trueW0 = p/2

parValuesList = list(
	Ccertain=10^(-1),   ##################
	Cuncertain=10^(-1), ##################
	epsilonCertain=10^(-1),  ################## no sense having these large if standardizing output (so to magnitude within 1 or so..)
	extraEpsilonUncertain = 10^(-1:0),  ################# for the two UNCERTAIN METAPARAMETERS, GO BACK TO THE DEFINITIONS TO CHECK IF THIS SCALE IS OK
	uncertaintySpecialTreatment = T
)	

missingVarPropVect = c(0.2, 0.9)########
missingObsPropVect = c(0.2, 0.9) ############
corValVect = c(0.9) ####################
theoRsqVect = c(0.8) ####################
quantOrSdPropValues = c(0.1) ####################
errMeasureVect=c("mae", "rmse", "Maxae", "cor", "quantNineAe", "quantEightAe", "quantSevenAe",
				"maeCert", "rmseCert", "MaxaeCert", "quantNineAeCert", "quantEightAeCert", "quantSevenAeCert", "corCert",
				"maeUncert", "rmseUncert", "MaxaeUncert", "quantNineAeUncert", "quantEightAeUncert", "quantSevenAeUncert", "corUncert") #maeCert #maeUncert, ...
approachVect = c("doPCbb", "doMedian", "doNoMiss", "doSquarebbSd", "doSquarebbQuant") 
AggregateTestError = mean
replaceImputedWithTrueY = F

maxUncertainDims = "all" # NULL #("all" considers the p+1 dims; NULL considers the actual max number of missing dims in all the data )
nRep=1 # # MUST DO MORE REPEATS, the results don't seem stable

if(T){
	missingY = F # Do not modify this
}



for(missingVarProp in missingVarPropVect){#############
for(missingObsProp in missingObsPropVect){############  
for(corVal in corValVect){ #################
for(theoRsq in theoRsqVect){
for(repIdx in 1:nRep){
rejectSilently=F	
set.seed(repIdx)
GenerateData() # inefficient, because redundant with the below, but useful to do prescreening of generated data
}
if(givenUp)
stop("Couldn't generate data partitions containing at least one missing point and one non-missing point")
}
}
}
}

cat("Data partitions containing at least one missing point and one non-missing point can be generated\n")
cat("Proceeding to actual experiments...\n\n\n")
	
for(missingVarProp in missingVarPropVect){#############
for(missingObsProp in missingObsPropVect){############  
for(corVal in corValVect){ #################
for(theoRsq in theoRsqVect){

cat("MissObs", missingObsProp, "Var", missingVarProp, "Cor", corVal, "Rsq", theoRsq, "starting\n")

totComb = nRep*length(errMeasureVect)*length(approachVect)
testRes = list()
cnt = 0

		
for(repIdx in 1:nRep){
	testRes[[repIdx]]=list()
	cat("Starting data generation, outer/inner splitting, and imputation, rep =", repIdx, "...\n")
	rejectSilently=T	
	set.seed(repIdx)
	GenerateData()
	MultiplyImpute()
	gc()
#	for(errMeasure in errMeasureVect){
#		testRes[[repIdx]][[errMeasure]]=list()
#		for(approach in approachVect){
	
	for(approach in approachVect){
		cat("Starting inner cross-validation, rep =", repIdx, ", approach =", approach,  "...\n")
		
		CalculateValidationErrors() # errors calculated with all possible error measures
		gc()
		for(errMeasure in errMeasureVect){
		
			if(is.null(testRes[[repIdx]][[errMeasure]]))
				testRes[[repIdx]][[errMeasure]]=list()
			
			testRes[[repIdx]][[errMeasure]][[approach]] = list()
			testRes[[repIdx]][[errMeasure]][[approach]]$valid = doErrorFoldOutInnerList 
			SetUpTest()
			testRes[[repIdx]][[errMeasure]][[approach]]$testSetup = getTrainResReadyForTest 
			
			cat("Starting inner testing, rep =", repIdx, ", approach =", approach, ", errMeasure =", errMeasure, "...\n")
			CalculateTestErrors() # choose best model based on errMeasure, train new model on outer training data, and test it
			gc()
			testRes[[repIdx]][[errMeasure]][[approach]]$testErrors = errVect 
			
			stopifnot(length(errVect)==numFolds)
			testRes[[repIdx]][[errMeasure]][[approach]]$testErrorsAggregate=AggregateTestError(errVect)
			testRes[[repIdx]][[errMeasure]][[approach]]$testErrorsSd=sd(errVect)
			
			cnt = cnt+1
			if(round(cnt/totComb*100)%%10==0){
				cat(cnt/totComb*100, "% done\n")
			}
		}
	}
}

if(F){
	for (errMeasure in errMeasureVect){
		for(approach in approachVect){
			cat(approach, errMeasure, testRes[[1]][[errMeasure]][[approach]]$testErrorsAggregate, "\n")
			cat("\n")
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

# for doMedian, the Cuncertain is irrelevant!
#set.seed(seed)

#totComb = nRep*length(errMeasureVect)*length(approachVect)
#testRes = list()
#cnt = 0
#for(rep in 1:nRep){
#	testRes[[rep]]=list()
#	cat("Starting data generation, outer/inner splitting, and imputation, rep =", rep, "...\n")
#	GenerateData()
#	MultiplyImpute()
#	for(errMeasure in errMeasureVect){
#		testRes[[rep]][[errMeasure]]=list()
#		for(approach in approachVect){
#			cat("Starting inner cross-validation, rep =", rep, ", errMeasure =", errMeasure,", approach =", approach,  "...\n")
#			CalculateValidationErrors()
#			SetUpTest()
#			cat("Starting inner testing, rep =", rep, ", errMeasure =", errMeasure, ", approach =", approach, "...\n")
#			CalculateTestErrors() # before this, can change 'approach'
#			stopifnot(length(errVect)==numFolds)
#			testRes[[rep]][[errMeasure]][[approach]]=AggregateTestError(errVect)
#			cnt = cnt+1
#			if(round(cnt/totComb*100)%%10==0){
#				cat(cnt/totComb*100, "% done\n")
#			}
#		}
#	}
#}


save.image(file=paste("MissObs", missingObsProp, "Var", missingVarProp, "Cor", corVal, "Rsq", theoRsq, ".RData", sep=""))
cat(paste("MissObs", missingObsProp, "Var", missingVarProp, "Cor", corVal, "Rsq", theoRsq, "done\n"))
}
}
}
}

