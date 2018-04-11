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

set.seed(seed)

numFolds = 2
scaleData = T
realData = F
injectMissingness = T
doMCAR = T ###
rejectSilently = F
maxGenerateDataAttempts=10

method = "pmm"
maxIter = 5
numImput = 11

n = 100
p = 4
meanVect = rep(0,p) 
stdVect = rep(1, p)
corVal = .5 
theoRsq = 0.9

missingVarProp = 0.8
missingObsProp = 0.8

trueW = -(1:p)
trueW0 = p/2
noiseLevel = 1

parValuesList = list(
	Ccertain=0, 
	Cuncertain=10^(-1:1),
	epsilonCertain=10^(-1:1),
	extraEpsilonUncertain =10^(-1:1), 
	uncertaintySpecialTreatment = T
)	


errMeasureVect=c("mae", "rmse", "Maxae") #maeCert #maeUncert
approachVect = c("doPCbb", "doMedian") #doSquarebb
AggregateTestError = mean
replaceImputedWithTrueY = F

maxUncertainDims = "all" # NULL #("all" considers the p+1 dims; NULL considers the actual max number of missing dims in all the data )


nRep=1

if(T){
	missingY = F # Do not modify this
}


GenerateData()
MultiplyImpute()

datMiss = doDataSplitOutOuter[[1]]$inDat$missing
datImput = doDataSplitOutOuter[[1]]$inDat$imputed


#imputDatList, medianImputDat, quantOrSdProp, scaleData, maxUncertainDims, doMedian

polyList = DoPolyList(missDat=datMiss, imputDatList=datImput$imputDatList, medianImputDat=datImput$medianImputDat, 
	quantOrSdProp=0.1, scaleData=T, maxUncertainDims=maxUncertainDims, doMedian=F, doNoMiss=F, doSquarebbSd=F, doSquarebbQuant=F) 
	
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

