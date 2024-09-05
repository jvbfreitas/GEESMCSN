library(qrLMM)
library(ggplot2)
library(cubature)
library(Matrix)
library(pracma)
library(MASS)
library(joineR)
library(ssym)
library(moments)
library(aod)
library(qqplotr)

data("Cholesterol")

# Organizing the covariate sex
levels(Cholesterol$sex) = c("Female", "Male")

# Exploratory plots
ggplot(Cholesterol,aes(x = year, y = log(cholst), group = newid)) +
  geom_line(alpha = 0.2)  +
  geom_point() +
  stat_summary(fun.y=mean, geom="line",lwd=1.0, aes(group= "blank")) +
  xlab("Year") + ylab("Cholesterol level") +
  scale_linetype(name="Sex")+
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  ) + facet_wrap(~sex)

t = as.vector(table(Cholesterol$newid))
n = length(t)
id = Cholesterol$newid

aux = variogram(id, Cholesterol$year,
                (Cholesterol$cholst-mean(Cholesterol$cholst))/sd(Cholesterol$cholst))
sp = spline(aux$svar)
ggplot() + geom_line(aes(x=sp$x,y=sp$y, linetype = "spline")) +
  geom_line(aes(x=aux$svar[,1],y=aux$sigma2, linetype = "Variance")) +
  theme_bw() + xlab("Lag") + ylab("Sample Variogram") +
  scale_linetype(name = "") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  )

ggplot(Cholesterol, aes(x = age, y = cholst)) + geom_point(alpha = 0.2)+
  geom_smooth(se = FALSE, colour = "Black", method = "loess")+
  xlab("Age") + ylab("Cholesterol level") +
  scale_linetype(name = "Sex") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  )


dad = Cholesterol
N = sum(t)
dad$year2 = (dad$year-5)/10
x = model.matrix(cholst~sex+year2-1, data = dad) # specification matrix
tt = dad$age
y = dad$cholst/100 # response variable
formula = y~x
all.knots = T # using all knots to obtain initial values
#nknot = 20
#family = "CST"
alpha = 100 # inital value of alpha
nformula = all.vars(formula)
p = ncol(x) # number of regression coefficients

###Get initial values

if(all.knots == F){
  t2  <- ncs(tt, nknots = nknot, lambda=1)
  tr  <- as.vector(attr(t2, "knots"))
  K   <- attr(t2, "K")
  Na   <- attr(t2, "N")
}else{
  tr  <- sort(unique(tt))
  t2  <- ncs(tt, all.knots = T, lambda = 1)
  K   <- attr(t2, "K")
  Na   <- attr(t2, "N")
}
beta = solve(t(x)%*%x)%*%t(x)%*%y
theta = solve(t(Na)%*%Na+alpha*K)%*%t(Na)%*%(y-x%*%beta)
sigma2 = var(y)
sigma = sqrt(sigma2)

gama = sign(skewness(y))*0.5
mu = x%*%beta + Na%*%theta

# Creating a new time vector
aux10 = dad$year
aux10[which(aux10==0)] = 1
aux10[which(aux10==2)] = 2
aux10[which(aux10==4)] = 3
aux10[which(aux10==6)] = 4
aux10[which(aux10==8)] = 5
aux10[which(aux10==10)] = 6
idvec=aux10

### Gamma GEE model
source("geegamma.R")
mod1 = geegama(y, mu, sigma2=1, n, p, alpha = 5000, K, Na, x,
            iter.max =65, alpha.fix= T, t, id, corstr = "exchangeable",
            alphalim = c(5000,10000), linkmu = "log", idvec)

round(mod1$beta,4) # regression coefficients
round(sqrt(diag(mod1$varnaive)[1:p]),4) # standard errors (naive)
# Wald test
wald.test(diag(mod1$varnaive)[1],mod1$beta[1], Terms = 1)
wald.test(diag(mod1$varnaive)[2],mod1$beta[2], Terms = 1)

round(sqrt(diag(mod1$vareps)[1:p]),4) # standard errors (robust)
# Wald test
wald.test(diag(mod1$vareps)[1],mod1$beta[1], Terms = 1)
wald.test(diag(mod1$vareps)[2],mod1$beta[2], Terms = 1)

mod1$sigma2 # scale parameter
mod1$df # Effective degress of freedom

# QAIC and QBIC
phig = 1/mod1$sigma2
mug = mod1$mu
QL = (1/mod1$sigma2)*(-y/mug-log(mug)+1+log(y))
-2*sum(QL)+2*(p+mod1$df+1) # QAIC
-2*sum(QL)+log(N)*(p+mod1$df+1) #QBIC

# QQ plot of the quantile residuals
dgammar = function(y,mu,phi){
  value = exp(phi*(-y/mu-log(mu))-log(gamma(phi))+phi*log(phi*y)-log(y))
}
resv = 0
for(i in 1:N){
  print(i)
  resv[i] = qnorm(integrate(function(x) dgammar(x,mug[i],phig),0,y[i])$value)
}
auxr = data.frame(resv)

ggplot(auxr,aes(sample = resv)) +
  stat_qq_band(bandType = "boot", fill = "white", colour = "black", size = 0.5) +
  stat_qq_line() +
  stat_qq_point() +
  xlab("Normal quantiles") + ylab("Quantile residuals") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  )

# Nonparametrix curves
varf = Na%*%mod1$vareps[-c(1:p),-c(1:p)]%*%t(Na)
var2 = Na%*%mod1$varnaive[-c(1:p),-c(1:p)]%*%t(Na)
ggplot() + geom_line(aes(x=tt,y=Na%*%mod1$theta), size = 1) +
  geom_line(aes(x = tt, y = Na%*%mod1$theta+1.96*sqrt(diag(var2))),
            linetype = "dotted") +
  geom_line(aes(x = tt, y = Na%*%mod1$theta-1.96*sqrt(diag(var2))),
            linetype = "dotted") +
  geom_line(aes(x = tt, y = Na%*%mod1$theta+1.96*sqrt(diag(varf))),
            linetype = "dashed") +
  geom_line(aes(x = tt, y = Na%*%mod1$theta-1.96*sqrt(diag(varf))),
            linetype = "dashed") +
  xlab("Age") + ylab("f(age)") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  )

### Inverse gaussian GEE model
source("geeig.R")
mod11 = geeig(y, mu, sigma2=1, n, p, alpha = 5000, K, Na, x,
               iter.max =65, alpha.fix= T, t, id, corstr = "exchangeable",
               alphalim = c(0.0001,0.01), linkmu = "log",idvec)

round(mod11$beta,4) # regression coefficients
round(sqrt(diag(mod11$varnaive)[1:p]),4) # standard errors (naive)
# Wald test
wald.test(diag(mod11$varnaive)[1],mod11$beta[1], Terms = 1)
wald.test(diag(mod11$varnaive)[2],mod11$beta[2], Terms = 1)

round(sqrt(diag(mod11$vareps)[1:p]),4) # standard errors (robust)
# Wald test
wald.test(diag(mod11$vareps)[1],mod11$beta[1], Terms = 1)
wald.test(diag(mod11$vareps)[2],mod11$beta[2], Terms = 1)

mod11$sigma2 # scale parameter
mod11$df # Effective degress of freedom

# QAIC and QBIC
phig = 1/mod11$sigma2
mug = mod11$mu
QL = (1/mod11$sigma2)*(1/mug - y/(2*mug^2) - 1/y + 1/(2*y))
-2*sum(QL)+2*(p+mod11$df+1) # QAIC
-2*sum(QL)+log(N)*(p+mod11$df+1) # QBIC

# QQ plot of the quantile residuals
digr = function(y,mu,phi){
  value = exp(phi*(-y/(2*mu^2)+1/mu)-.5*(log(2*pi*y^3/phi)+phi/y))
}
resv = 0
for(i in 1:N){
  print(i)
  resv[i] = qnorm(integrate(function(x) digr(x,mug[i],phig),0,y[i])$value)
}
auxr = data.frame(resv)
ggplot(auxr,aes(sample = resv)) +
  stat_qq_band(bandType = "boot", fill = "white", colour = "black", size = 0.5) +
  stat_qq_line() +
  stat_qq_point() +
  xlab("Normal quantiles") + ylab("Quantile residuals") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  )

# Nonparametrix curves
varf = Na%*%mod11$vareps[-c(1:p),-c(1:p)]%*%t(Na)
var2 = Na%*%mod11$varnaive[-c(1:p),-c(1:p)]%*%t(Na)
ggplot() + geom_line(aes(x=tt,y=Na%*%mod11$theta), size = 1) +
  geom_line(aes(x = tt, y = Na%*%mod11$theta+1.96*sqrt(diag(var2))),
            linetype = "dotted") +
  geom_line(aes(x = tt, y = Na%*%mod11$theta-1.96*sqrt(diag(var2))),
            linetype = "dotted") +
  geom_line(aes(x = tt, y = Na%*%mod11$theta+1.96*sqrt(diag(varf))),
            linetype = "dashed") +
  geom_line(aes(x = tt, y = Na%*%mod11$theta-1.96*sqrt(diag(varf))),
            linetype = "dashed") +
  xlab("age") + ylab("f(age)") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  )

### CSN GEE model
source("utils.r")
source("geecsn.R")
mod2= geecsn(y, mu, sigma2, gama, n, p, alpha = 5000, K, Na, x,
             iter.max = 50, alpha.fix= T, t, id, corstr = "exchangeable",
             alphalim = c(1,10000), link = "log",idvec)

round(mod2$beta,4) # regression coefficients
round(sqrt(diag(mod2$varnaive)[1:p]),4) # standard errors (naive)
# Wald test
wald.test(diag(mod2$varnaive)[1],mod2$beta[1], Terms = 1)
wald.test(diag(mod2$varnaive)[2],mod2$beta[2], Terms = 1)

round(sqrt(diag(mod2$vareps)[1:p]),4) # standard errors (robust)
# Wald test
wald.test(diag(mod2$vareps)[1],mod2$beta[1], Terms = 1)
wald.test(diag(mod2$vareps)[2],mod2$beta[2], Terms = 1)

mod2$sigma2 # scale parameter
mod2$gama # skewness parameter
mod2$df # Effective degress of freedom

# QAIC and QBIC
QL = -(y-mod2$mu)^2/(2*mod2$sigma2)
-2*sum(QL)+2*(p+mod2$df+2) # QAIC
-2*sum(QL)+log(N)*(p+mod2$df+2) # QBIC

# QQ plot of the quantile residuals
auxr = data.frame(resv = mod2$resq)
ggplot(auxr,aes(sample = resv)) +
  stat_qq_band(bandType = "boot", fill = "white", colour = "black", size = 0.5) +
  stat_qq_line() +
  stat_qq_point() +
  xlab("Normal quantiles") + ylab("Quantile residuals") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  )

# Nonparametrix curves
varf = Na%*%mod2$vareps[-c(1:p),-c(1:p)]%*%t(Na)
var2 = Na%*%mod2$varnaive[-c(1:p),-c(1:p)]%*%t(Na)
ggplot() + geom_line(aes(x=tt,y=Na%*%mod2$theta), size = 1) +
  geom_line(aes(x = tt, y = Na%*%mod2$theta+1.96*sqrt(diag(var2))),
            linetype = "dotted") +
  geom_line(aes(x = tt, y = Na%*%mod2$theta-1.96*sqrt(diag(var2))),
            linetype = "dotted") +
  geom_line(aes(x = tt, y = Na%*%mod2$theta+1.96*sqrt(diag(varf))),
            linetype = "dashed") +
  geom_line(aes(x = tt, y = Na%*%mod2$theta-1.96*sqrt(diag(varf))),
            linetype = "dashed") +
  xlab("age") + ylab("f(age)") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  )

### CST GEE model
source("geecst.R")
mod3 = geecst(y, mu, sigma2, gama, n, p, alpha = 5000, K, Na, x,
              iter.max = 25, alpha.fix= T, t, id, corstr = "exchangeable",
              alphalim = c(0.001,5000), linkmu = "log",idvec)

round(mod3$beta,4) # regression coefficients
round(sqrt(diag(mod3$varnaive)[1:p]),4) # standard errors (naive)
# Wald test
wald.test(diag(mod3$varnaive)[1],mod3$beta[1], Terms = 1)
wald.test(diag(mod3$varnaive)[2],mod3$beta[2], Terms = 1)

round(sqrt(diag(mod3$vareps)[1:p]),4) # standard errors (robust)
# Wald test
wald.test(diag(mod3$vareps)[1],mod3$beta[1], Terms = 1)
wald.test(diag(mod3$vareps)[2],mod3$beta[2], Terms = 1)

mod3$sigma2 # scale parameter
mod3$gama # skewness parameter
mod3$nu # shape parameter
mod3$df # Effective degress of freedom

# QAIC and QBIC
QL = -(y-mod3$mu)^2/(2*mod3$vary)
-2*sum(QL)+2*(p+mod3$df+3) # QAIC
-2*sum(QL)+log(N)*(p+mod3$df+3) # QBIC

# QQ plot of the quantile residuals
auxr = data.frame(resv = mod3$resq)
ggplot(auxr,aes(sample = resv)) +
  stat_qq_band(bandType = "boot", fill = "white", colour = "black", size = 0.5) +
  stat_qq_line() +
  stat_qq_point() +
  xlab("Normal quantiles") + ylab("Quantile residuals") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  )

# Nonparametric curves
varf = Na%*%mod3$vareps[-c(1:p),-c(1:p)]%*%t(Na)
var2 = Na%*%mod3$varnaive[-c(1:p),-c(1:p)]%*%t(Na)
ggplot() + geom_line(aes(x=tt,y=Na%*%mod3$theta), size = 1) +
  geom_line(aes(x = tt, y = Na%*%mod3$theta+1.96*sqrt(diag(var2))),
            linetype = "dotted") +
  geom_line(aes(x = tt, y = Na%*%mod3$theta-1.96*sqrt(diag(var2))),
            linetype = "dotted") +
  geom_line(aes(x = tt, y = Na%*%mod3$theta+1.96*sqrt(diag(varf))),
            linetype = "dashed") +
  geom_line(aes(x = tt, y = Na%*%mod3$theta-1.96*sqrt(diag(varf))),
            linetype = "dashed") +
  xlab("age") + ylab("f(age)") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  )

### CSS GEE model
source("geecss.R")
mod4 = geecss(y, mu, sigma2=1, gama, n, p, alpha = 5000, K, Na, x,
              iter.max = 100, alpha.fix= T, t, id,
              corstr = "exchangeable",
              linkmu = "log",idvec, alphalim = c(0.01,600000))

round(mod4$beta,4) # regression coefficients
round(sqrt(diag(mod4$varnaive)[1:p]),4) # standard errors (naive)
# Wald test
wald.test(diag(mod4$varnaive)[1],mod4$beta[1], Terms = 1)
wald.test(diag(mod4$varnaive)[2],mod4$beta[2], Terms = 1)

round(sqrt(diag(mod4$vareps)[1:p]),4) # standard errors (robust)
# Wald test
wald.test(diag(mod4$vareps)[1],mod4$beta[1], Terms = 1)
wald.test(diag(mod4$vareps)[2],mod4$beta[2], Terms = 1)

mod4$sigma2 # scale parameter
mod4$gama # skewness parameter
mod4$nu # shape parameter
mod4$df # Effective degress of freedom

# QAIC and QBIC
QL = -(y-mod4$mu)^2/(2*mod4$vary)
-2*sum(QL)+2*(p+mod4$df+3) # QAIC
-2*sum(QL)+log(N)*(p+mod4$df+3) # QBIC

# QQ plot of the quantile residuals
resq = 0
for(i in 1:N){
  resq[i] = qnorm(pCSS(y[i],mod4$mu[i],mod4$sigma2, mod4$gama, mod4$nu))
}
auxr = data.frame(resq)
ggplot(auxr,aes(sample = resq)) +
  stat_qq_band(bandType = "boot", fill = "white", colour = "black", size = 0.5) +
  stat_qq_line() +
  stat_qq_point() +
  xlab("Normal quantiles") + ylab("Quantile residuals") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  )

# Nonparametric curves
varf = Na%*%mod4$vareps[-c(1:p),-c(1:p)]%*%t(Na)
var2 = Na%*%mod4$varnaive[-c(1:p),-c(1:p)]%*%t(Na)
ggplot() + geom_line(aes(x=tt,y=Na%*%mod4$theta), size = 1) +
  geom_line(aes(x = tt, y = Na%*%mod4$theta+1.96*sqrt(diag(var2))),
            linetype = "dotted") +
  geom_line(aes(x = tt, y = Na%*%mod4$theta-1.96*sqrt(diag(var2))),
            linetype = "dotted") +
  geom_line(aes(x = tt, y = Na%*%mod4$theta+1.96*sqrt(diag(varf))),
            linetype = "dashed") +
  geom_line(aes(x = tt, y = Na%*%mod4$theta-1.96*sqrt(diag(varf))),
            linetype = "dashed") +
  xlab("age") + ylab("f(age)") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  )

### CSGT GEE model
source("geecsgt.R")
mod5 = geecsgt(y, mu, gama, n, p, alpha = 5000, K, Na, x,
              iter.max = 100, alpha.fix= T, t, id,
              corstr = "exchangeable",
              linkmu = "log",idvec)

round(mod5$beta,4) # regression coefficients
round(sqrt(diag(mod5$varnaive)[1:p]),4) # standard errors (naive)
# Wald test
wald.test(diag(mod5$varnaive)[1],mod5$beta[1], Terms = 1)
wald.test(diag(mod5$varnaive)[2],mod5$beta[2], Terms = 1)

round(sqrt(diag(mod5$vareps)[1:p]),4) # standard errors (robust)
# Wald test
wald.test(diag(mod5$vareps)[1],mod5$beta[1], Terms = 1)
wald.test(diag(mod5$vareps)[2],mod5$beta[2], Terms = 1)

mod5$sigma2 # scale parameter
mod5$gama # skewness parameter
mod5$nu # shape parameter
mod5$df # Effective degress of freedom

# QAIC and QBIC
QL = -(y-mod5$mu)^2/(2*mod5$vary)
-2*sum(QL)+2*(p+mod5$df+3) # QAIC
-2*sum(QL)+log(N)*(p+mod5$df+3) # QBIC

# QQ plot of the quantile residuals
auxr = data.frame(resq = mod5$resq)
ggplot(auxr,aes(sample = resq)) +
  stat_qq_band(bandType = "boot", fill = "white", colour = "black", size = 0.5) +
  stat_qq_line() +
  stat_qq_point() +
  xlab("Normal quantiles") + ylab("Quantile residuals") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  )

# Nonparametric curves
varf = Na%*%mod5$vareps[-c(1:p),-c(1:p)]%*%t(Na)
var2 = Na%*%mod5$varnaive[-c(1:p),-c(1:p)]%*%t(Na)
ggplot() + geom_line(aes(x=tt,y=Na%*%mod5$theta), size = 1) +
  geom_line(aes(x = tt, y = Na%*%mod5$theta+1.96*sqrt(diag(varf))),
            linetype = "dashed") +
  geom_line(aes(x = tt, y = Na%*%mod5$theta-1.96*sqrt(diag(varf))),
            linetype = "dashed") +
  geom_line(aes(x = tt, y = Na%*%mod5$theta+1.96*sqrt(diag(var2))),
            linetype = "dotted") +
  geom_line(aes(x = tt, y = Na%*%mod5$theta-1.96*sqrt(diag(var2))),
            linetype = "dotted") +
  xlab("age") + ylab("f(age)") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  )

# Local influence analysis
g13 = nthroot(mod5$gama,3)
g23 = g13^2
lambda = (s*g13)/sqrt(2/pi+s^2*g23*(2/pi - 1))
delta = lambda/sqrt(1+lambda^2)
Delta = delta/sqrt(1-b^2*delta^2)
tau = as.numeric((1-delta^2)/(1-b^2*delta^2))
Lambda = diag(as.vector(mod5$mu),N)
Lampon = Lambda
Malpha = bdiag(matrix(0,p,p),mod5$alpha*K)
M = cbind(x,Na)
Sigma = mod5$Omega
u = (y-mod5$mu)
Se = as.matrix(t(M)%*%(Lampon%*%solve(Sigma)%*%diag(u)-Lambda%*%solve(Sigma))%*%M-Malpha)

##### Case-weight perturbation
Delt = as.matrix(t(M)%*%mod5$W%*%solve(Lambda)%*%diag(u))
BG = -t(Delt)%*%solve(Se)%*%Delt
dmax = Re(eigen(BG)$vectors[,1])
qplot(id,dmax) + xlab("Subject") + ylab(expression(d[max])) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  )

##### Response perturbation
vary = mod5$vary
Delt = as.matrix(t(M)%*%mod5$W%*%solve(Lambda))*sqrt(vary)
BG = -t(Delt)%*%solve(Se)%*%Delt
dmax = Re(eigen(BG)$vectors[,1])
qplot(id,dmax) + xlab("Subject") + ylab(expression(d[max])) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  )

##### Shape1 perturbation
nu = mod5$nu
dvdo = -nu[1]*(1/(2*sqrt(vary)))*(-nu[2]/((nu[1]-2)^2))
Rrho = mod5$Rm
vary = mod5$vary
dSdo = diag(rep(dvdo,N))*Rrho*diag(rep(sqrt(vary),N)) + diag(rep(dvdo,N))*Rrho*diag(rep(sqrt(vary),N))
Sigma = diag(rep(sqrt(vary),N))*Rrho*diag(rep(sqrt(vary),N))
dSm1 = -solve(Sigma)%*%dSdo%*%solve(Sigma)
Delt = t(M)%*%Lambda%*%dSm1%*%diag(u)
BG = -t(Delt)%*%solve(Se)%*%Delt
dmax = eigen(BG)$vectors[,1]
qplot(id,Re(dmax)) + xlab("Subject") + ylab(expression(d[max])) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  )

##### Shape2 perturbation
nu = mod5$nu
dvdo = -nu[2]*(1/(2*sqrt(vary)))*(1/(nu[1]-2))
dSdo = diag(rep(dvdo,N))*Rrho*diag(rep(sqrt(vary),N)) + diag(rep(dvdo,N))*Rrho*diag(rep(sqrt(vary),N))
Sigma = diag(rep(sqrt(vary),N))*Rrho*diag(rep(sqrt(vary),N))
dSm1 = -solve(Sigma)%*%dSdo%*%solve(Sigma)
Delt = t(M)%*%Lambda%*%dSm1%*%diag(u)
BG = -t(Delt)%*%solve(Se)%*%Delt
dmax = Re(eigen(BG)$vectors[,1])
qplot(id,dmax) + xlab("Subject") + ylab(expression(d[max])) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  )

##### working correlation
R = as.matrix(mod5$Rm)
dR = -R + diag(N)
varu = vary
dSl = diag(sqrt(varu),N)%*%dR%*%diag(sqrt(varu),N)
dSigma = dSl
Sigma = mod5$Omega
invSigma = solve(Sigma)
dinvSigma <- -invSigma%*%dSigma%*%invSigma
Delt = as.matrix(t(M)%*%(Lambda)%*%dinvSigma%*%diag(u))
BG = -t(Delt)%*%solve(Se)%*%Delt
dmax = Re(eigen(BG)$vectors[,1])
qplot(id,dmax) + xlab("Subject") + ylab(expression(d[max])) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  )

### CSCN GEE model
source("geecscn.R")
mod6 = geecscn(y, mu, sigma2, gama, n, p, alpha = 5000, K, Na, x,
              iter.max = 50, alpha.fix= T, t, id,
              corstr = "exchangeable",
              linkmu = "log",idvec)

round(mod6$beta,4) # regression coefficients
round(sqrt(diag(mod6$varnaive)[1:p]),4) # standard errors (naive)
# Wald test
wald.test(diag(mod6$varnaive)[1],mod6$beta[1], Terms = 1)
wald.test(diag(mod6$varnaive)[2],mod6$beta[2], Terms = 1)

round(sqrt(diag(mod6$vareps)[1:p]),4) # standard errors (robust)
# Wald test
wald.test(diag(mod6$vareps)[1],mod6$beta[1], Terms = 1)
wald.test(diag(mod6$vareps)[2],mod6$beta[2], Terms = 1)

mod6$sigma2 # scale parameter
mod6$gama # skewness parameter
mod6$nu # shape parameter
mod6$df # Effective degress of freedom

# QAIC and QBIC
QL = -(y-mod6$mu)^2/(2*mod6$vary)
-2*sum(QL)+2*(p+mod6$df+4)# QAIC
-2*sum(QL)+log(N)*(p+mod6$df+4)# QBIC

# QQ plot of the quantile residuals
resq = 0
for(i in 1:N){
  resq[i] = qnorm(pCSCN(y[i],mod6$mu[i],mod6$sigma2, mod6$gama, mod6$nu))
}
auxr = data.frame(resq)
ggplot(auxr,aes(sample = resq)) +
  stat_qq_band(bandType = "boot", fill = "white", colour = "black", size = 0.5) +
  stat_qq_line() +
  stat_qq_point() +
  xlab("Normal quantiles") + ylab("Quantile residuals") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  )

# Nonparametric curves
varf = Na%*%mod6$vareps[-c(1:p),-c(1:p)]%*%t(Na)
var2 = Na%*%mod6$varnaive[-c(1:p),-c(1:p)]%*%t(Na)
ggplot() + geom_line(aes(x=tt,y=Na%*%mod6$theta), size = 1) +
  geom_line(aes(x = tt, y = Na%*%mod6$theta+1.96*sqrt(diag(varf))),
            linetype = "dashed") +
  geom_line(aes(x = tt, y = Na%*%mod6$theta-1.96*sqrt(diag(varf))),
            linetype = "dashed") +
  geom_line(aes(x = tt, y = Na%*%mod6$theta+1.96*sqrt(diag(var2))),
            linetype = "dotted") +
  geom_line(aes(x = tt, y = Na%*%mod6$theta-1.96*sqrt(diag(var2))),
            linetype = "dotted") +
  xlab("age") + ylab("f(age)") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  )

### CSBPN GEE model
source("geecsbpn.R")
mod7 = geecsbpn(y, mu, sigma2, gama, n, p, alpha = 5000, K, Na,
                as.matrix(x),
               iter.max = 50, alpha.fix= T, t, id,
               corstr = "unstructured",
               linkmu = "log",idvec)

round(mod7$beta,4) # regression coefficients
round(sqrt(diag(mod7$varnaive)[1:p]),4) # standard errors (naive)
# Wald test
wald.test(diag(mod7$varnaive)[1],mod7$beta[1], Terms = 1)
wald.test(diag(mod7$varnaive)[2],mod7$beta[2], Terms = 1)

round(sqrt(diag(mod7$vareps)[1:p]),4) # standard errors (robust)
# Wald test
wald.test(diag(mod7$vareps)[1],mod7$beta[1], Terms = 1)
wald.test(diag(mod7$vareps)[2],mod7$beta[2], Terms = 1)

mod7$sigma2 # scale parameter
mod7$gama # skewness parameter
mod7$nu # shape parameter
mod7$df # Effective degress of freedom

# QAIC and QBIC
QL = -(y-mod7$mu)^2/(2*mod7$vary)
-2*sum(QL)+2*(p+mod6$df+4)# QAIC
-2*sum(QL)+log(N)*(p+mod6$df+4)# QBIC

# QQ plot of the quantile residuals
resq = 0
for(i in 1:N){
  resq[i] = qnorm(pCSBPN(y[i],mod7$mu[i],mod7$sigma2, mod7$gama, mod7$nu))
}
auxr = data.frame(resq)
ggplot(auxr,aes(sample = resq)) +
  stat_qq_band(bandType = "boot", fill = "white", colour = "black", size = 0.5) +
  stat_qq_line() +
  stat_qq_point() +
  xlab("Normal quantiles") + ylab("Quantile residuals") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  )

# Nonparametric curves
varf = Na%*%mod7$vareps[-c(1:p),-c(1:p)]%*%t(Na)
var2 = Na%*%mod7$varnaive[-c(1:p),-c(1:p)]%*%t(Na)
ggplot() + geom_line(aes(x=tt,y=Na%*%mod7$theta), size = 1) +
  geom_line(aes(x = tt, y = Na%*%mod7$theta+1.96*sqrt(diag(varf))),
            linetype = "dashed") +
  geom_line(aes(x = tt, y = Na%*%mod7$theta-1.96*sqrt(diag(varf))),
            linetype = "dashed") +
  geom_line(aes(x = tt, y = Na%*%mod7$theta+1.96*sqrt(diag(var2))),
            linetype = "dotted") +
  geom_line(aes(x = tt, y = Na%*%mod7$theta-1.96*sqrt(diag(var2))),
            linetype = "dotted") +
  xlab("age") + ylab("f(age)") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  )


### CSBSN GEE model
source("geecsbsn.R")
mod8 = geecsbsn(y, mu, sigma2, gama, n, p, alpha = 5000, K, Na,
                as.matrix(x),
                iter.max = 50, alpha.fix= T, t, id,
                corstr = "exchangeable",
                linkmu = "log",idvec)

round(mod8$beta,4) # regression coefficients
round(sqrt(diag(mod8$varnaive)[1:p]),4) # standard errors (naive)
# Wald test
wald.test(diag(mod8$varnaive)[1],mod8$beta[1], Terms = 1)
wald.test(diag(mod8$varnaive)[2],mod8$beta[2], Terms = 1)

round(sqrt(diag(mod8$vareps)[1:p]),4) # standard errors (robust)
# Wald test
wald.test(diag(mod8$vareps)[1],mod8$beta[1], Terms = 1)
wald.test(diag(mod8$vareps)[2],mod8$beta[2], Terms = 1)

mod8$sigma2 # scale parameter
mod8$gama # skewness parameter
mod8$nu # shape parameter
mod8$df # Effective degress of freedom

# QAIC and QBIC
QL = -(y-mod8$mu)^2/(2*mod8$vary)
-2*sum(QL)+2*(p+mod8$df+4)# QAIC
-2*sum(QL)+log(N)*(p+mod8$df+4)# QBIC

# QQ plot of the quantile residuals
resq = 0
for(i in 1:N){
  print(i)
  resq[i] = qnorm(pCSBSN(y[i],mod8$mu[i],mod8$sigma2, mod8$gama, mod8$nu))
}
auxr = data.frame(resq)
ggplot(auxr,aes(sample = resq)) +
  stat_qq_band(fill = "white", colour = "black", size = 0.5) +
  stat_qq_line() +
  stat_qq_point() +
  xlab("Normal quantiles") + ylab("Quantile residuals") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  )

# Nonparametric curves
varf = Na%*%mod8$vareps[-c(1:p),-c(1:p)]%*%t(Na)
var2 = Na%*%mod8$varnaive[-c(1:p),-c(1:p)]%*%t(Na)
ggplot() + geom_line(aes(x=tt,y=Na%*%mod8$theta), size = 1) +
  geom_line(aes(x = tt, y = Na%*%mod8$theta+1.96*sqrt(diag(varf))),
            linetype = "dashed") +
  geom_line(aes(x = tt, y = Na%*%mod8$theta-1.96*sqrt(diag(varf))),
            linetype = "dashed") +
  geom_line(aes(x = tt, y = Na%*%mod8$theta+1.96*sqrt(diag(var2))),
            linetype = "dotted") +
  geom_line(aes(x = tt, y = Na%*%mod8$theta-1.96*sqrt(diag(var2))),
            linetype = "dotted") +
  xlab("age") + ylab("f(age)") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  )

### CSGGN GEE model
source("geecsggn.R")
mod9 = geecsggn(y, mu, gama, n, p, alpha = 10000, K, Na,
                as.matrix(x),
                iter.max = 50, alpha.fix= T, t, id,
                corstr = "exchangeable",
                linkmu = "log",idvec)

round(mod9$beta,4) # regression coefficients
round(sqrt(diag(mod9$varnaive)[1:p]),4) # standard errors (naive)
# Wald test
wald.test(diag(mod9$varnaive)[1],mod9$beta[1], Terms = 1)
wald.test(diag(mod9$varnaive)[2],mod9$beta[2], Terms = 1)

round(sqrt(diag(mod9$vareps)[1:p]),4) # standard errors (robust)
# Wald test
wald.test(diag(mod9$vareps)[1],mod9$beta[1], Terms = 1)
wald.test(diag(mod9$vareps)[2],mod9$beta[2], Terms = 1)

mod9$sigma2 # scale parameter
mod9$gama # skewness parameter
mod9$nu # shape parameter
mod9$df # Effective degress of freedom

# QAIC and QBIC
QL = -(y-mod9$mu)^2/(2*mod9$vary)
-2*sum(QL)+2*(p+mod9$df+4)# QAIC
-2*sum(QL)+log(N)*(p+mod9$df+4)# QBIC

# QQ plot of the quantile residuals
resq = 0
for(i in 1:N){
  print(i)
  resq[i] = qnorm(pCSGGN(y[i],mod9$mu[i], mod9$gama, mod9$nu))
}
auxr = data.frame(resq)
ggplot(auxr,aes(sample = resq)) +
  stat_qq_band(fill = "white", colour = "black", size = 0.5) +
  stat_qq_line() +
  stat_qq_point() +
  xlab("Normal quantiles") + ylab("Quantile residuals") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  )

# Nonparametric curves
varf = Na%*%mod9$vareps[-c(1:p),-c(1:p)]%*%t(Na)
var2 = Na%*%mod9$varnaive[-c(1:p),-c(1:p)]%*%t(Na)
ggplot() + geom_line(aes(x=tt,y=Na%*%mod9$theta), size = 1) +
  geom_line(aes(x = tt, y = Na%*%mod9$theta+1.96*sqrt(diag(varf))),
            linetype = "dashed") +
  geom_line(aes(x = tt, y = Na%*%mod9$theta-1.96*sqrt(diag(varf))),
            linetype = "dashed") +
  geom_line(aes(x = tt, y = Na%*%mod9$theta+1.96*sqrt(diag(var2))),
            linetype = "dotted") +
  geom_line(aes(x = tt, y = Na%*%mod9$theta-1.96*sqrt(diag(var2))),
            linetype = "dotted")+
  xlab("age") + ylab("f(age)") +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    text=element_text(size=15,family="serif")
  )