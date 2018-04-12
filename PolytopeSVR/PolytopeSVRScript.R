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
#source("../Volume.R")

#######MissObs0.9Var0.25Cor0.9.RData: why are results not as good?

#######noMissDat --> the best testing parameters for uncertain data should be chosen based on the best for certain predictions: WHY? you can still calculate validation performance according to uncertain data only...

#### are the extraEpsilonUncertain and Cuncertain parameters playing a role?

# we don't need to have a model that works well on both certain and uncertain; we need one that does better at least one of those; for the other thing
# we can use a standard model, etc.


maxGenerateDataAttempts = 20
numFolds = 5#####################
scaleData = T
method = "pmm" #norm, cart, rf
maxIter = 20 ################ try with small numbers of these two, to verify if imputation is possible first
numImput = 20 ###############

realData = F
injectMissingness = T
doMCAR = T #######

if(!realData){
	n = 100 # do 120 for 10 folds is the least to have at least one miss/non miss point if missingObsPropVect = 0.1 or 0.9
	p = 10
	meanVect = rep(0,p) 
	stdVect = rep(1, p)

	trueW = 1:p
	trueW0 = p/2
	corValVect = 0.9#c(0.9, 0.5, 0) ####################
	theoRsqVect = 0.95#c(0.95, 0.9, 0.8) ####################

} else{
	#"Automobile.RData" # kept numerical variables, removed 4 obs with NA Y; has natural NAs
	#Boston.RData --boston corrected: kept numerical variables (removed boolean); has no natural NAs
	realDataFileName = "Automobile.RData" # kept
	corValVect = "irrelevant"
	theoRsqVect = "irrelevant" 
}

parValuesList = list(
	Ccertain=c(0, .05, .1, 1, 2, 10),#c(0,10^(-2:1)),   ##################
	Cuncertain=c(0, .05, .1, 1, 2, 10),#c(0,10^(-2:1)), ##################
	epsilonCertain=c(0, .5, 1, 5, 10),#c(0,10^(-2:1)),  ################## no sense having these large if standardizing output (so to magnitude within 1 or so..)
	extraEpsilonUncertain = c(0, .5, 1, 5, 10),# c(0,10^(-2:1)),  ################# for the two UNCERTAIN METAPARAMETERS, GO BACK TO THE DEFINITIONS TO CHECK IF THIS SCALE IS OK
	uncertaintySpecialTreatment = T,
	linear =T
	)	

missingVarPropVect = 0.9#c(0.9, 0.2)########
missingObsPropVect = 0.9# c(0.9, 0.2) ############
quantOrSdPropValues = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1) ####################
#errMeasureVect=c("mae", "rmse", "Maxae", "cor", "quantNineAe", "quantEightAe", "quantSevenAe",
#"maeCert", "rmseCert", "MaxaeCert", "quantNineAeCert", "quantEightAeCert", "quantSevenAeCert", "corCert",
#"maeUncert", "rmseUncert", "MaxaeUncert", "quantNineAeUncert", "quantEightAeUncert", "quantSevenAeUncert", "corUncert") #maeCert #maeUncert, ...
errMeasureVect=c("mae", "rmse", "Maxae", "cor",  "quantEightAe", "maeCert", "rmseCert", "MaxaeCert",  "corCert", "quantEightAeCert",
	"maeUncert", "rmseUncert", "MaxaeUncert", "corUncert", "quantEightAeUncert") #maeCert #maeUncert, ...
approachVect = c("doPCbb", "doSquarebbSd", "doSquarebbQuant", "doMedian", "doNoMiss")  ####################
AggregateTestError = mean
replaceImputedWithTrueY = F

maxUncertainDims = "all" # NULL #("all" considers the p+1 dims; NULL considers the actual max number of missing dims in all the data )
nRep=1 # # MUST DO MORE REPEATS, the results don't seem stable

if(T){
	missingY = F # Do not modify this
}


if(!injectMissingness){
	cat("Ignoring missingVarPropVect and missingObsPropVect (setting to 0), and doMCAR, since injectMissingness=F\n")
	missingVarPropVect=0	
	missingObsPropVect=0
	doMCAR = "irrelevant"
}

for(missingVarProp in missingVarPropVect){#############
	for(missingObsProp in missingObsPropVect){############  
		for(corVal in corValVect){ #################
			for(theoRsq in theoRsqVect){
				for(repIdx in 1:nRep){
					rejectSilently=F	
					set.seed(repIdx)
					GenerateData() # inefficient, because redundant with the below, but useful to do prescreening of generated data
					if(givenUp)
						stop("Couldn't generate data partitions containing at least one missing point and one non-missing point")
				}
			}
		}
	}
}

cat("Data partitions containing at least one missing point and one non-missing point can be generated\n")
cat("Proceeding to actual experiments...\n\n\n")


currDate = system('date +%Y%m%d-%H%M%S', intern=T)

if(realData){
	fileNameRoot = paste0(realDataFileName)
} else{
	fileNameRoot = paste0("NormalN", n, "P", p)
}		

resultsFolderName = paste0(fileNameRoot, "Date", currDate)
system(paste("mkdir", resultsFolderName))
system(paste("cp *.R *sh", resultsFolderName))
system(paste("cp ../*.R", resultsFolderName))
progressFile = paste0(resultsFolderName, "/Progress.txt")
				
totComb = length(missingVarPropVect)*length(missingObsPropVect)*length(corValVect)*length(theoRsqVect)*
	nRep*length(approachVect)*length(errMeasureVect)
	
progressOut=paste("STARTING a total of", totComb, "combinations at", date(), "\n")
cat(progressOut)
write.table(progressOut, quote=F, row.names=F, col.names=F, append=T, file=progressFile)
	
cnt = 1
				
for(missingVarProp in missingVarPropVect){#############
	for(missingObsProp in missingObsPropVect){############  
		for(corVal in corValVect){ #################
			for(theoRsq in theoRsqVect){
				if(realData){
					fileNameRoot = paste0(fileNameRoot, "Cor", corVal, "Rsq", theoRsq)
				}		
				
				testRes = list()

				for(repIdx in 1:nRep){
					testRes[[repIdx]]=list()
					rejectSilently=T	
					set.seed(repIdx)
					GenerateData()
					MultiplyImpute()
					gc()

					for(approach in approachVect){
						appShort = sub("do", "", approach)
						appShort = sub("Square", "Sq", appShort)
						fileName = paste0(fileNameRoot, "MissObs", missingObsProp, "MissVar", missingVarProp, ifelse(doMCAR, "MCAR", "MAR"), "Meth", method, 
							"Appr", appShort, "Date", currDate, ".RData", sep="")
							
						progressOut = paste("STARTING cross-validation of", fileName, "Rep", repIdx, "at", date(), "\n")
						cat(progressOut)
						write.table(progressOut, quote=F, row.names=F, col.names=F, append=T, file=progressFile)
						
						CalculateValidationErrors() # errors calculated with all possible error measures
						gc()
						for(errMeasure in errMeasureVect){

							if(is.null(testRes[[repIdx]][[errMeasure]]))
								testRes[[repIdx]][[errMeasure]]=list()

							testRes[[repIdx]][[errMeasure]][[approach]] = list()
							testRes[[repIdx]][[errMeasure]][[approach]]$valid = doErrorFoldOutInnerList 
							SetUpTest()
							testRes[[repIdx]][[errMeasure]][[approach]]$testSetup = getTrainResReadyForTest 

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

			}
		}
	}
}

