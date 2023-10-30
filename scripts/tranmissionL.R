Ctr=1.5
atr = 2
btr = 10
htr = function(tau) Ctr*dgamma(tau, shape = atr, rate = btr) # instantaneous rate 
Str = function(tau) exp(-Ctr*pgamma(tau, shape=atr, rate=btr)) # prob of no event from 0 to tau
loghtr = function(tau) log(Ctr)+dgamma(tau, shape = atr, rate = btr, log = TRUE) # log version
logStr = function(tau) (-Ctr*pgamma(tau,  shape=atr, rate=btr)) # log version


Cs=0.9 # Cs must now be strictly less than 1
aS = 2.5
bS= 10
hs = function(tau) Cs*dgamma(tau, shape = aS, rate = bS) # instantaneous rate 
Ss = function(tau) exp(-Cs*pgamma(tau, shape=aS, rate=bS)) # prob of no event from 0 to tau
loghs = function(tau) log(Cs)+ dgamma(tau, shape = aS, rate = bS, log = TRUE) # log version
logSs = function(tau) (-Cs*pgamma(tau,  shape=aS, rate=bS)) # log version


# recall newton's method: x_n+1 = x_n - f(x_n)/f'(x_n) to find a root of f 
# here f = x- (1-C)*exp(lambda(x-1)) and f' is 1- (1-C)*lambda*exp(lambda(x-1))
maxsteps <- 1000; n=0; 

getp0 = function(Cs, Ctr, tol=1e-6, maxsteps=1000, x0 = 0.1,output=TRUEa) { 

    f = function(x) {x - (1-Cs)*exp(Ctr*(x-1))} 
    fprime = function(x) {1 - (1-Cs)*Ctr*exp(Ctr*(x-1))}
    n=0

while(abs(f(x0)) > tol && n < maxsteps) {
  x0 = x0 - f(x0)/fprime(x0)  # Newton's method formula
  n=n+1 
}
if(n < maxsteps) {
    p0=x0
 if (output == TRUE) { 
 cat("Root found:", x0, "\n")
  cat("Function value at root:", f(x0), "\n")
  cat("Number of iterations:", n, "\n")}
} else {
  cat("The p0 algorithm did not converge after", maxstep, "iterations\n")
}
return(p0)
}

p0 = getp0(Cs, Ctr,output=TRUE)


getIndivCondition = function(p0, yri) {
    TT = 1 - exp(-(1-p0))*Str(yri)*Ss(yri)
    return(TT)
}


logLblue = logSs(3)+logStr(3)+loghtr(0.65)+loghtr(0.5) -log(getIndivCondition(p0,yri=3))
logLbob = logSs(1.2)+loghs(1.2)+logStr(3-0.5)- log(getIndivCondition(p0, yri=3-0.5))
logLeve =logSs(1.45)+loghs(1.45)+logStr(3-0.65) - log(getIndivCondition(p0, yri=3-0.65))

logEx1 = logLbob + logLeve+logLblue
logEx1



phi = 1 - p0*(1+ Ctr*(1-p0)/(1-Cs)) # ignoring that lambda might not be quite Ctr for now 
rho = 1 - exp(-phi)*Str(100)*Ss(100) # ignoring the difference between 100 and infinity for now
# i note that Str(100) is the same (to 8 digits) as Str(200), so 100 is "close to infinity". lol.
Yright = 5 # just for this example, suppose the block started 5 time units before T, so the RT time is 5. 

# I am going to use the a (a=2) and b (b=10) pars from the transmission hazard
a=atr
b=btr

 Z=0 # init for the sum 
 n=1 
 tol=1e-8
 maxn = 1e6
 term = 1 # initialize the term at something > tol 
 while ((term > tol) & (n < maxn)) {
     term =   ((1-rho)^n)*pgamma(Yright,n*a,b)
     Z = Z+term
     n=n+1
 }

cat("Approximate sum:", Z, "\n")



# recall p0 can be obtained with getp0(Cs,Ctr) 

getPhi = function( Cs, Ctr, p0) { 
    return(1 - p0*(1+ Ctr*(1-p0)/(1-Cs)))
}

getRho = function(phi) { 
    return (1 - exp(-phi)*Str(100)*Ss(100))
}

getBlockCondition = function(p0, rho, atr,btr, Yr, tol=1e-7, maxn = 1e6, output=FALSE) {
    Z =0 # init for the sum 
    n=1  
    term = 1 # initialize the term at something > tol 
    while ((term > tol) & (n < maxn)) {
        term =   ((1-rho)^n)*pgamma(Yr,n*atr,btr)
        Z = Z+term
        n=n+1
    }
    if (output==TRUE) {
        cat("Approximate sum:", Z, "\n")
        cat("Number of terms used:", n, "\n") 
        } 
    return(Z) 
    } 


getBlockLike =function(tblock, bc, Yr, rho,atr,btr, p0) {
    (1-rho^bc)*dgamma(tblock, shape=bc*atr, rate=btr) / getBlockCondition(p0,rho, atr, btr, Yr)
}






# Test 3

c1 = htr(0.85 - 0.5)*Str(1.7-0.5)*Ss(1.445-0.5)*hs(1.445-0.5) /getIndivCondition(p0,1.7-0.5)
c2 = Str(1.7-0.85)*Ss(1.006-0.85)*hs(1.006-0.85)/getIndivCondition(p0, 1.7-0.85)
c3 = Str(1.7-0.15)*htr(0.35-0.15)*Ss(0.53-0.15)*hs(0.53-0.15) /getIndivCondition(p0,1.7-0.15)
c4 = Str(1.7-0.72)*htr(0.96-0.72)*Ss(1.53-0.72)*hs(1.53-0.72)/getIndivCondition(p0, 1.7-0.72)
c5 = Str(1.7-1.33)*Ss(1.459- 1.33)*hs(1.459-1.33) /getIndivCondition(p0, 1.7-1.33)
c0n6 = Ss(1.7-0)*Str(1.7-0)*htr(0.22)*htr(0.15) /getIndivCondition(p0, 1.7)
b1 = getBlockLike(0.5-0.22, bc = 4,Yr = 1.7-0.22, rho, atr, btr, p0)
b2 = getBlockLike(0.72-0.35, bc = 3,Yr = 1.7-0.35, rho, atr, btr, p0)
b3 = getBlockLike(1.33-0.96, bc = 2,Yr = 1.7-0.96, rho, atr, btr, p0)

L = c1*c2*c3*c4*c5*c0n6*b1*b2*b3
L

## [1] 7.643214e-05

log(L)


log(c1)
log(c2)
log(c3)
log(c4)
log(c5)
log(c0n6)
log(b1)
log(b2)
log(b3)


# Test 4

c3 = Str(1.5-0.12)*Ss(0.3308-0.12)*hs(0.3308-0.12)/getIndivCondition(p0, 1.5-0.12)
c5 = Str(1.5-1.05)*Ss(1.2592-1.05)*hs(1.2592-1.05)/getIndivCondition(p0, 1.5-1.05)
c4 = Str(1.5)*htr(0.12)*htr(0.73)*Ss(1.3347)*hs(1.3347)/getIndivCondition(p0, 1.5)
b=getBlockLike(1.05-0.73, bc = 2,Yr = 1.5-0.73, rho, atr, btr, p0)

logL4 = log(c3)+log(c5)+log(c4)+log(b) 
logL4