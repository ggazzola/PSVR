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
missingVarProp = 0.5
missingObsProp = 0.75
corVal = 0.2


#######noMissDat --> the best testing parameters for uncertain data should be chosen based on the best for certain predictions: WHY? you can still calculate validation performance according to uncertain data only...

#### are the extraEpsilonUncertain and Cuncertain parameters playing a role?

# we don't need to have a model that works well on both certain and uncertain; we need one that does better at least one of those; for the other thing
# we can use a standard model, etc.


maxGenerateDataAttempts = 20
numFolds = 3#####################
scaleData = T
method = "pmm" #norm, cart, rf
maxIter = 20 ################ try with small numbers of these two, to verify if imputation is possible first
numImput = 50 ###############

AggregateTestError = mean
replaceImputedWithTrueY = F

maxUncertainDims = "all" # NULL #("all" considers the p+1 dims; NULL considers the actual max number of missing dims in all the data )

if(T){
	missingY = F # Do not modify this
}

if(!realData){
	n = 30 # do 120 for 10 folds is the least to have at least one miss/non miss point if missingObsPropVect = 0.1 or 0.9
	p = 4
	meanVect = rep(0,p) 
	stdVect = rep(1, p)
	trueW = 1:p
	trueW0 = p/2
	theoRsq = 0.95
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
	Ccertain=c(0),#c(0,10^(-2:1)),   ##################
	Cuncertain=c(0),#c(0,10^(-2:1)), ##################
	epsilonCertain=c(0),#c(0,10^(-2:1)),  ################## no sense having these large if standardizing output (so to magnitude within 1 or so..)
	extraEpsilonUncertain = c(0, 10),# c(0,10^(-2:1)),  ################# for the two UNCERTAIN METAPARAMETERS, GO BACK TO THE DEFINITIONS TO CHECK IF THIS SCALE IS OK
	uncertaintySpecialTreatment = T,
	linear =T
	)	


quantOrSdPropValues = c(0) ####################
#errMeasureVect=c("mae", "rmse", "Maxae", "cor", "quantNineAe", "quantEightAe", "quantSevenAe",
#"maeCert", "rmseCert", "MaxaeCert", "quantNineAeCert", "quantEightAeCert", "quantSevenAeCert", "corCert",
#"maeUncert", "rmseUncert", "MaxaeUncert", "quantNineAeUncert", "quantEightAeUncert", "quantSevenAeUncert", "corUncert") #maeCert #maeUncert, ...
errMeasureVect=c("mae","MaxaeCert",  "quantEightAeUncert") #maeCert #maeUncert, ...
#approachVect = c("doPCbb", "doSquarebbSd", "doSquarebbQuant", "doMedian", "doNoMiss")  ####################
approachVect = c("doPCbb", "doMedian", "doNoMiss")  ####################


repVect=1:1 # # MUST DO MORE REPEATS, the results don't seem stable

GenerateData()
MultiplyImpute()

datMiss = doDataSplitOutOuter[[1]]$inDat$missing
datImput = doDataSplitOutOuter[[1]]$inDat$imputed



polyListPC = DoPolyList(missDat=datMiss, imputDatList=datImput$imputDatList, medianOrientedOrNonOrientedImputDat=datImput$medianOrientedBoxImputDat, 
	quantOrSdProp=0, scaleData=T, maxUncertainDims=maxUncertainDims, doMedian=F, doNoMiss=F, doSquarebbSd=F, doSquarebbQuant=F) 
	
polyListMed = DoPolyList(missDat=datMiss, imputDatList=datImput$imputDatList, medianOrientedOrNonOrientedImputDat=datImput$medianOrientedBoxImputDat, 
	quantOrSdProp=0, scaleData=T, maxUncertainDims=maxUncertainDims, doMedian=T, doNoMiss=F, doSquarebbSd=F, doSquarebbQuant=F) 
	
parList = list(Ccertain=1, Cuncertain=1, epsilonCertain=0, extraEpsilonUncertain=0, uncertaintySpecialTreatment=F, twoSlacks = F, linear=T) 



polyListPCSmall = polyListMedSmall = list()
polyListPCSmall[[1]] = polyListPC[[1]]
polyListMedSmall[[1]] = polyListMed[[1]]

for(i in 2:6){
	polyListPCSmall[[i]] = polyListPC[[i]]
	polyListMedSmall[[i]] = polyListMed[[i]]
	modelPC = DoTrainModel(polyListPCSmall, parList)
	modelMed = DoTrainModel(polyListMedSmall, parList)
	cat("Different with", i, "points:", sum(abs(c(modelPC$w, modelPC$w0)- c(modelMed$w, modelMed$w0))), "\n")
	
}




datStructure = doDataSplitOutOuter[[1]]
errorMedian = NULL
errorOrientedMedian = NULL
imputDatList = datStructure$inDat$imputed$imputDatList

corMatList = list()
for(index in 1:nrow(datStructure$inDat$original)){
	truePoint = datStructure$inDat$original[index,]
	missingPoint = datStructure$inDat$missing[index,]
	imputMat = GetMultipleImputSinglePoint(index, imputDatList)
	corMatList[[index]] = cor(imputMat)
	medianImput =  datStructure$inDat$imputed$medianImputDat[index,]
	medianOrientedImput =  datStructure$inDat$imputed$medianOrientedBoxImputDat[index,]
	errMedian= abs(truePoint-medianImput)
	errorMedian[index] = mean(errMedian[is.na(missingPoint)])
	errOrientedMedian= abs(truePoint-medianOrientedImput)
	errorOrientedMedian[index] = mean(errOrientedMedian[is.na(missingPoint)])
}


pp= PerformanceByParameterValue(doErrorFoldOutInnerList)

modelPC = DoTrainModel(polyListPCSmall, parList)
modelMed = DoTrainModel(polyListMedSmall, parList)

DoTrainModel(polyList, parListNonLinear)

polyListBak = polyList
polyList=list()
polyList[[1]]=polyListBak[[1]]
res = list()

for(i in 2:length(polyListBak)){
polyList[[i]]=polyListBak[[i]]

	

parListNonLinear = list(Ccertain=1, Cuncertain=1, epsilonCertain=0, extraEpsilonUncertain=0, uncertaintySpecialTreatment=F, twoSlacks = F, linear=F) 
resNonLinear = DoTrainModel(polyList, parListNonLinear); #print(sum(c(resTwoSlacks$csiPlus, resTwoSlacks$csiMinus)))

modelTwo=model
resGurobiTwo = resGurobi

parListLinear = list(Ccertain=1, Cuncertain=1, epsilonCertain=0, extraEpsilonUncertain=0, uncertaintySpecialTreatment=F, twoSlacks = F, linear=T) 
resLinear= DoTrainModel(polyList, parListLinear); #print(sum(resOneSlack$csi))


modelOne=model
resGurobiOne = resGurobi

#res[[i]] =list()	

#res[[i]]$resTwoSlacks=resTwoSlacks
#res[[i]]$resOneSlack=resOneSlack
#res[[i]]$ErrorDiff=sum(c(resTwoSlacks$csiPlus, resTwoSlacks$csiMinus))-sum(resOneSlack$csi) # expecting >=0
#cat(i, " has ", res[[i]]$ErrorDiff, "\n")

cat(i, "NonLinear, w:", resNonLinear$w, "\n")
cat(i, "Linear, w:", resLinear$w, "\n")
}

if(F){
	i=1

polyListOutIndiv = polyList[[i]]
doTrainModelOut = resOneSlack
DoMinMaxPrediction(polyListOutIndiv,doTrainModelOut)
GetSolution(modelOne,resGurobiOne, "u", i)%*%polyList[[i]]$a-GetSolution(modelOne,resGurobiOne, "w0")
GetSolution(modelOne,resGurobiOne, "v", i)%*%polyList[[i]]$a+GetSolution(modelOne,resGurobiOne, "w0")

GetSolution(modelOne,resGurobiOne, "w")%*%DoMinMaxPrediction(polyListOutIndiv,doTrainModelOut)$u$x+GetSolution(modelOne,resGurobiOne, "w0")
GetSolution(modelOne,resGurobiOne, "w")%*%DoMinMaxPrediction(polyListOutIndiv,doTrainModelOut)$v$x+GetSolution(modelOne,resGurobiOne, "w0")
GetSolution(modelOne,resGurobiOne,"csi", i)


doTrainModelOut = resTwoSlacks
DoMinMaxPrediction(polyListOutIndiv,doTrainModelOut)
GetSolution(modelTwo,resGurobiTwo, "u", i)%*%polyList[[i]]$a-GetSolution(modelTwo,resGurobiTwo, "w0")
GetSolution(modelTwo,resGurobiTwo, "v", i)%*%polyList[[i]]$a+GetSolution(modelTwo,resGurobiTwo, "w0")

GetSolution(modelTwo,resGurobiTwo, "w")%*%DoMinMaxPrediction(polyListOutIndiv,doTrainModelOut)$u$x+GetSolution(modelTwo,resGurobiTwo, "w0")
GetSolution(modelTwo,resGurobiTwo, "w")%*%DoMinMaxPrediction(polyListOutIndiv,doTrainModelOut)$v$x+GetSolution(modelTwo,resGurobiTwo, "w0")
GetSolution(modelTwo,resGurobiTwo,"csiPlus", i)
GetSolution(modelTwo,resGurobiTwo,"csiMinus", i)
}



azzOne=lapply(polyList, DoMinMaxPrediction, doTrainModelOut=NULL, w=GetSolution(modelOne,resGurobiOne, "w"), w0=GetSolution(modelOne,resGurobiOne, "w0"))
azzTwo=lapply(polyList, DoMinMaxPrediction, doTrainModelOut=NULL, w=GetSolution(modelTwo,resGurobiTwo, "w"), w0=GetSolution(modelTwo,resGurobiTwo, "w0"))

unlist(lapply(azzOne, function(x) max(x$u$err, x$v$err))) - GetSolution(modelOne,resGurobiOne, "csi")
unlist(lapply(azzTwo, function(x) max(x$u$err,0)+max(x$v$err,0))) - (GetSolution(modelTwo,resGurobiTwo, "csiPlus")+GetSolution(modelTwo,resGurobiTwo, "csiMinus"))

sum(unlist(lapply(azzOne, function(x) max(x$u$err, x$v$err))))+GetSolution(modelOne,resGurobiOne, "w")%*%GetSolution(modelOne,resGurobiOne, "w")

sum(unlist(lapply(azzTwo, function(x) max(x$u$err, x$v$err))))+GetSolution(modelTwo,resGurobiTwo, "w")%*%GetSolution(modelTwo,resGurobiTwo, "w") # disadvantageous to use the w,w0 solution of two-slack in one-slack formulation

