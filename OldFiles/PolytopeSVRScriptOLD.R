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

seed=1

numFolds = 3###########5
scaleData = T

missingVarProp = 0.25###
missingObsProp = .9 
method = "pmm"
maxIter = 10##############20
numImput = 50

n = 50#############100
p = 5
meanVect = rep(0,p) 
stdVect = rep(1, p)
corVal = .9 

trueW = 1:p
trueW0 = p/2
noiseLevel = 1

parValuesList = list(
	Ccertain=10,#10^(-1:1), ########## 
	Cuncertain=10,#10^(-1:1),##########
	epsilonCertain=10,#10^(-1:1),##########
	extraEpsilonUncertain =c(1,10),#10^(-1:1), ############
	uncertaintySpecialTreatment = T
)	

# Miss prop should be applied separately to train/valid/test, so that they all have that proportion

pcPropValues = c(0.1, .25, 0.5, 1)

errMeasureVect=c("mae", "rmse", "Maxae", "quantNineAe", "quantEightAE", "quantSevenAE",
				"maeCert", "rmseCert", "MaxaeCert", "quantNineAeCert", "quantEightAECert", "quantSevenAECert",
				"maeUncert", "rmseUncert", "MaxaeUncert", "quantNineAeUncert", "quantEightAEUncert", "quantSevenAEUncert") #maeCert #maeUncert, ...
approachVect = c("doPCbb", "doMedian", "doNoMiss") #doSquarebb # doNoMissData
AggregateTestError = mean
replaceImputedWithTrueY = F

maxUncertainDims = "all" # NULL #("all" considers the p+1 dims; NULL considers the actual max number of missing dims in all the data )


nRep=1

if(T){
	missingY = F # Do not modify this
}



totComb = nRep*length(errMeasureVect)*length(approachVect)
testRes = list()
cnt = 0
for(repIdx in 1:nRep){
	testRes[[repIdx]]=list()
	cat("Starting data generation, outer/inner splitting, and imputation, rep =", repIdx, "...\n")
	GenerateData()
	MultiplyImpute()
	for(errMeasure in errMeasureVect){
		testRes[[repIdx]][[errMeasure]]=list()
		for(approach in approachVect){
			cat("Starting inner cross-validation, rep =", repIdx, ", errMeasure =", errMeasure,", approach =", approach,  "...\n")
			testRes[[repIdx]][[errMeasure]][[approach]] = list()
			CalculateValidationErrors() 
			testRes[[repIdx]][[errMeasure]][[approach]]$valid = doErrorFoldOutInnerList 
			SetUpTest()
			testRes[[repIdx]][[errMeasure]][[approach]]$testSetup = getTrainResReadyForTest 
			
			cat("Starting inner testing, rep =", repIdx, ", errMeasure =", errMeasure, ", approach =", approach, "...\n")
			CalculateTestErrors() # before this, can change 'approach'
			testRes[[repIdx]][[errMeasure]][[approach]]$testErrors = errVect 
			
			stopifnot(length(errVect)==numFolds)
			testRes[[repIdx]][[errMeasure]][[approach]]$testErrorsAggregate=AggregateTestError(errVect)
			cnt = cnt+1
			if(round(cnt/totComb*100)%%10==0){
				cat(cnt/totComb*100, "% done\n")
			}
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


save.image(file=paste("Miss", missingVarProb, "Cor", corVal, ".RData", sep=""))


