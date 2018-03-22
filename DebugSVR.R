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

numFolds = 2
scaleData = T


method = "pmm"
maxIter = 5
numImput = 6

n = 50
p = 4
meanVect = rep(0,p) 
stdVect = rep(1, p)
corVal = .5 
theoRsq = 0.9

missingVarProp = 0.25
missingObsProp = 0.2

trueW = 1:p
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

datOrig =  doDataSplitOutOuter[[1]]$inDat$original 
datMiss = doDataSplitOutOuter[[1]]$inDat$missing
datImput = doDataSplitOutOuter[[1]]$inDat$imputed


#imputDatList, medianImputDat, pcProp, scaleData, maxUncertainDims, doMedian

polyList = DoPolyList(missDat=datMiss, imputDatList=datImput$imputDatList, medianImputDat=datImput$medianImputDat, 
	pcProp=1, scaleData=T, maxUncertainDims, doMedian=F, doNoMiss=F) 

parList = list(Ccertain=1, Cuncertain=1, epsilonCertain=0, extraEpsilonUncertain=0, uncertaintySpecialTreatment=T)

res = DoTrainModel(polyList, parList)