# simple, p=1 n=3 example with:
# x_1 = 2, y_1 = 4
# 2<=x_2<=4, y_2 = 6
# x_3+y_3>=11, x_3+y_3 <=13, x_3-y_3>=-5, x_3-y_3<=-3
# the uncertainty of x_2 is centered on the "true" value 3
# the uncertainty of (x_3, y_3) is centered on the "true" values (4, 8)
# epsilon is assumed 0
# So, if the uncertainty were replaced by the "true" values, the "true" values of (w, w0) would be (2, 0)

# Commented out parts below correspond to a first version where lower bounds on variables were not implemented correctly

C = 10


#lhs = matrix(0, nrow=27, ncol = 17)


lhs = matrix(0, nrow=12, ncol = 17)

lhs[1,c(1,2, 15)]= c(2,1,-1)
lhs[2,c(1,2, 15)]= c(-2,-1,-1)

lhs[3,c(2, 3,4,16)]= c(1,-2,4, -1)
lhs[4,c(1, 3,4)]= c(-1,-1,1)
lhs[5,c(2, 5,6, 16)]= c(-1,-2,4, -1)
lhs[6,c(1, 5,6)]= c(1,-1,1)

lhs[7,c(2, 7:10,17)]= c(1,-11,5,13,-3,-1)
lhs[8:9,c(1, 7:10)]= matrix(c(-1, -1, -1, 1,1, 0, -1, 1, 1, -1), byrow=T, nrow=2)
lhs[10,c(2, 11:14,17)]= c(-1,-11,5,13,-3,-1)

lhs[11:12,c(1, 11:14)]= matrix(c(1, -1, -1, 1,1, 0, -1, 1, 1, -1), byrow=T, nrow=2)

#lhs[13:27, 3:17] = diag(nrow=15, ncol=15)
#rhs = c(4, -4, 6, 0, -6, 0, 0, 0, -1, 0, 0, +1, rep(0, 15))
#direction = c("<=", "<=", "<=", "=", "<=", "=", "<=", "=", "=", "<=", "=", "=", rep(">=", 15))

rhs = c(4, -4, 6, 0, -6, 0, 0, 0, -1, 0, 0, +1)
direction = c("<=", "<=", "<=", "=", "<=", "=", "<=", "=", "=", "<=", "=", "=")
lowerBound = c(-Inf, -Inf, rep(0, 15))
startSol=c(2.0, 0.0, 0.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.5, 1.5, 0.5, 1.5, 0.0, 0.0, 0.0, 2.0, 2.0)
vBasis = rep(0, length(startSol))
vBasis[startSol==0]=-1



Q = matrix(0, nrow=ncol(lhs), ncol=ncol(lhs))
Q[1,1] = 0
coeff = rep(0, ncol(lhs))
coeff[14:17] = C
#coeff[1]=-1
model = list()
model$A = lhs
model$rhs = rhs
model$sense = direction
model$Q = Q
model$obj = coeff
model$vbasis = vBasis
model$lb = lowerBound
result <- gurobi(model)
