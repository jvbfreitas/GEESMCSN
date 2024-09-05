library(pracma)
library(rmutil)
#######################################################################################################################
## Quantities ##
#######################################################################################################################

s = nthroot(2/(4-pi),3)
b = sqrt(2/pi)

#######################################################################################################################
## AIC and BIC for CSMSN ##
#######################################################################################################################

QBICsmn=function(alpha,gama,y,mu,vary,x,Na,W,K,type){
  df = sum(diag(Na%*%solve(t(Na)%*%W%*%Na + alpha*K)%*%t(Na)%*%W))
  g13 = nthroot(gama,3)
  g23 = g13^2
  lambda = (s*g13)/sqrt(2/pi+s^2*g23*(2/pi - 1))
  delta = lambda/sqrt(1+lambda^2)
  #Delta = sigma*delta/sqrt(1-b^2*delta^2)
  #tau = as.numeric(sigma2*(1-delta^2)/(1-b^2*delta^2))
  u = y-mu
  N = length(u)
  M = cbind(x,Na)
  if(linkmu == "identity"){
    Lambda = diag(rep(1,N),N)
  }
  if(linkmu == "log"){
    Lambda = diag(as.vector(mu),N)
  }
  p = ncol(x)
  qlik = -sum((y-mu)^2/(2*vary))
  if(type=="CSS"||type=="CST"){AICm= log(N)*(p+df+3)-2*qlik}
  if(type=="CSCN"||type=="CSBPN"||type=="CSBSN"||type=="CSGT"){AICm= log(N)*(p+df+4)-2*qlik}
  if(type=="CSGGN"){AICm= log(N)*(p+df+5)-2*qlik}
  if(type=="CSN"){AICm= log(N)*(p+df+2)-2*qlik}
  return(AICm)
}

AIC.ASMCSN=function(obj){
  alpha = obj$coefficients$alpha
  lp = obj$logver
  ff = obj$coefficients$theta
  K = obj$comp$K
  N = obj$comp$N
  tau = obj$coefficients$tau
  p = ncol(obj$comp$X)
  e10 = diag(obj$expectations$S1n)
  n = nrow(e10)
  type = obj$type
  df = rep(0,length(K))
  for(i in 1:length(alpha)){
    q=ncol(K[[i]])
    L=tau*alpha[i]*K[[i]]
    aux = N[[i]]%*%solve(t(N[[i]])%*%e10%*%N[[i]] +L)%*%t(N[[i]])
    df[i]=sum(diag(aux))
    lp = lp - 0.5*alpha[i]*(t(ff[[i]])%*%K[[i]]%*%ff[[i]])
  }
  if(type=="CSS"||type=="CST"){AICm= -2*lp+2*(p+3+sum(df))}
  if(type=="CSCN"||type=="CSBPN"||type=="CSBSN"||type=="CSGT"){AICm= -2*lp+2*(p+4+sum(df))}
  if(type=="CSGGN"){AICm= -2*lp+2*(p+5+sum(df))}
  if(type=="CSN"){AICm= -2*lp+2*(p+2+sum(df))}
  return(AICm)
}

BIC.ASMCSN=function(obj){
  alpha = obj$coefficients$alpha
  lp = obj$logver
  ff = obj$coefficients$theta
  K = obj$comp$K
  N = obj$comp$N
  tau = obj$coefficients$tau
  p = ncol(obj$comp$X)
  e10 = diag(obj$expectations$S1n)
  n = nrow(e10)
  type = obj$type
  df = rep(0,length(K))
  for(i in 1:length(alpha)){
    q=ncol(K[[i]])
    L=tau*alpha[i]*K[[i]]
    aux = N[[i]]%*%solve(t(N[[i]])%*%e10%*%N[[i]] +L)%*%t(N[[i]])
    df[i]=sum(diag(aux))
    lp = lp - 0.5*alpha[i]*(t(ff[[i]])%*%K[[i]]%*%ff[[i]])
  }
  if(type=="CSS"||type=="CST"){AICm= -2*lp+log(n)*(p+3+sum(df))}
  if(type=="CSCN"||type=="CSBPN"||type=="CSBSN"||type=="CSGT"){AICm= -2*lp+log(n)*(p+4+sum(df))}
  if(type=="CSGGN"){AICm= -2*lp+log(n)*(p+5+sum(df))}
  if(type=="CSN"){AICm= -2*lp+log(n)*(p+2+sum(df))}
  return(AICm)
}

AICc.ASMCSN=function(obj){
  alpha = obj$coefficients$alpha
  lp = obj$logver
  ff = obj$coefficients$theta
  K = obj$comp$K
  N = obj$comp$N
  tau = obj$coefficients$tau
  p = ncol(obj$comp$X)
  e10 = diag(obj$expectations$S1n)
  n = nrow(e10)
  type = obj$type
  df = rep(0,length(K))
  for(i in 1:length(alpha)){
    q=ncol(K[[i]])
    L=tau*alpha[i]*K[[i]]
    aux = N[[i]]%*%solve(t(N[[i]])%*%e10%*%N[[i]] +L)%*%t(N[[i]])
    df[i]=sum(diag(aux))
    lp = lp - 0.5*alpha[i]*(t(ff[[i]])%*%K[[i]]%*%ff[[i]])
  }
  if(type=="CSS"||type=="CST"){
    kp = p+3+sum(df)
    AICm = -2*lp+2*kp*(kp+1)/(n-kp-1)
  }
  if(type=="CSCN"||type=="CSBPN"||type=="CSBSN"||type=="CSGT"){
    kp = p+4+sum(df)
    AICm = -2*lp+2*kp*(kp+1)/(n-kp-1)
  }
  if(type=="CSGGN"){
    kp = p+5+sum(df)
    AICm = -2*lp+2*kp*(kp+1)/(n-kp-1)
  }
  if(type=="CSN"){
    kp = p+2+sum(df)
    AICm = -2*lp+2*kp*(kp+1)/(n-kp-1)
  }
  return(AICm)
}

HQIC.ASMCSN=function(obj){
  alpha = obj$coefficients$alpha
  lp = obj$logver
  ff = obj$coefficients$theta
  K = obj$comp$K
  N = obj$comp$N
  tau = obj$coefficients$tau
  p = ncol(obj$comp$X)
  e10 = diag(obj$expectations$S1n)
  n = nrow(e10)
  type = obj$type
  df = rep(0,length(K))
  for(i in 1:length(alpha)){
    q=ncol(K[[i]])
    L=tau*alpha[i]*K[[i]]
    aux = N[[i]]%*%solve(t(N[[i]])%*%e10%*%N[[i]] +L)%*%t(N[[i]])
    df[i]=sum(diag(aux))
    lp = lp - 0.5*alpha[i]*(t(ff[[i]])%*%K[[i]]%*%ff[[i]])
  }
  if(type=="CSS"||type=="CST"){AICm= -2*lp+2*log(log(n))*(p+3+sum(df))}
  if(type=="CSCN"||type=="CSBPN"||type=="CSBSN"||type=="CSGT"){AICm= -2*lp+2*log(log(n))*(p+4+sum(df))}
  if(type=="CSGGN"){AICm= -2*lp+2*log(log(n))*(p+5+sum(df))}
  if(type=="CSN"){AICm= -2*lp+2*log(log(n))*(p+2+sum(df))}
  return(AICm)
}

CAIC.ASMCSN=function(obj){
  alpha = obj$coefficients$alpha
  lp = obj$logver
  ff = obj$coefficients$theta
  K = obj$comp$K
  N = obj$comp$N
  tau = obj$coefficients$tau
  p = ncol(obj$comp$X)
  e10 = diag(obj$expectations$S1n)
  n = nrow(e10)
  type = obj$type
  df = rep(0,length(K))
  for(i in 1:length(alpha)){
    q=ncol(K[[i]])
    L=tau*alpha[i]*K[[i]]
    aux = N[[i]]%*%solve(t(N[[i]])%*%e10%*%N[[i]] +L)%*%t(N[[i]])
    df[i]=sum(diag(aux))
    lp = lp - 0.5*alpha[i]*(t(ff[[i]])%*%K[[i]]%*%ff[[i]])
  }
  if(type=="CSS"||type=="CST"){AICm= -2*lp+(1+log(n))*(p+3+sum(df))}
  if(type=="CSCN"||type=="CSBPN"||type=="CSBSN"||type=="CSGT"){AICm= -2*lp+(1+log(n))*(p+4+sum(df))}
  if(type=="CSGGN"){AICm= -2*lp+(1+log(n))*(p+5+sum(df))}
  if(type=="CSN"){AICm= -2*lp+(1+log(n))*(p+2+sum(df))}
  return(AICm)
}

SABIC.ASMCSN=function(obj){
  alpha = obj$coefficients$alpha
  lp = obj$logver
  ff = obj$coefficients$theta
  K = obj$comp$K
  N = obj$comp$N
  tau = obj$coefficients$tau
  p = ncol(obj$comp$X)
  e10 = diag(obj$expectations$S1n)
  n = nrow(e10)
  type = obj$type
  df = rep(0,length(K))
  for(i in 1:length(alpha)){
    q=ncol(K[[i]])
    L=tau*alpha[i]*K[[i]]
    aux = N[[i]]%*%solve(t(N[[i]])%*%e10%*%N[[i]] +L)%*%t(N[[i]])
    df[i]=sum(diag(aux))
    lp = lp - 0.5*alpha[i]*(t(ff[[i]])%*%K[[i]]%*%ff[[i]])
  }
  if(type=="CSS"||type=="CST"){AICm= -2*lp+log((2+n)/24)*(p+3+sum(df))}
  if(type=="CSCN"||type=="CSBPN"||type=="CSBSN"||type=="CSGT"){AICm= -2*lp+log((2+n)/24)*(p+4+sum(df))}
  if(type=="CSGGN"){AICm= -2*lp+log((2+n)/24)*(p+5+sum(df))}
  if(type=="CSN"){AICm= -2*lp+log((2+n)/24)*(p+2+sum(df))}
  return(AICm)
}

#######################################################################################################################
## Densities for CSMSN ##
#######################################################################################################################

dCSN = function(y,mu,sigma2,gama){
  g13 = nthroot(gama,3)
  g23 = g13^2
  s = nthroot(2/(4-pi),3)
  sigma = sqrt(sigma2)
  xi = mu-sigma*g13*s
  om = sigma*sqrt(1+s^2*g23)
  lamb = (s*g13)/sqrt(2/pi+s^2*g23*(2/pi-1))
  dens = (2/om)*dnorm((y-xi)/om)*pnorm(lamb*(y-xi)/om)
  return(dens)
}
pCSN = function(q,mu,sigma2,gama,nu){
  value = cubintegrate(function(x) dCSN(x,mu,sigma2,gama),-Inf,q)$integral
  return(value)
}
qCSN = function(p,mu,sigma2,gama){
  value = 0
  for(i in 1:length(p)){
    value[i] = multiroot(function(x) pCSN(x,mu,sigma2,gama)-p[i], start = 0)$root
  }
  return(value)
}
dCSS = function(y,mu,sigma2,gama,nu){
  dens = 0
  for(i in 1:length(y)){
    dens[i] = integrate(function(x) dCSN(y[i],mu[i],sigma2/x,gama)*dbeta(x,nu,1),0,1)$value
  }
  return(dens)
}
pCSS = function(q,mu,sigma2,gama,nu){
  value = cubintegrate(function(x) dCSS(x,mu,sigma2,gama,nu),-Inf,q)$integral
  return(value)
}
qCSS = function(p,mu,sigma2,gama,nu){
  value = 0
  for(i in 1:length(p)){
    value[i] = multiroot(function(x) pCSS(x,mu,sigma2,gama,nu)-p[i], start = 0)$root
  }
  return(value)
}
dCSGT = function(y,mu,gama,nu){
  dens = 0
  for(i in 1:length(y)){
    dens[i] = integrate(function(x) dCSN(y[i],mu[i],1/x,gama)*dgamma(x,nu[1]/2,nu[2]/2),0,Inf)$value
  }
  return(dens)
}
pCSGT = function(q,mu,gama,nu){
  value = cubintegrate(function(x) dCSGT(x,mu,gama,nu),-Inf,q)$integral
  return(value)
}
qCSGT = function(p,mu,gama,nu){
  value = 0
  for(i in 1:length(p)){
    value[i] = multiroot(function(x) pCSGT(x,mu,gama,nu)-p[i], start = 0)$root
  }
  return(value)
}
dCST = function(y,mu,sigma2,gama,nu){
  dens = 0
  for(i in 1:length(y)){
    dens[i] = integrate(function(x) dCSN(y[i],mu[i],sigma2/x,gama)*dgamma(x,nu/2,nu/2),0,Inf,subdivisions = 10000)$value
  }
  return(dens)
}
pCST = function(q,mu,sigma2,gama,nu){
  value = cubintegrate(function(x) dCST(x,mu,sigma2,gama,nu),-Inf,q)$integral
  return(value)
}
qCST = function(p,mu,sigma2,gama,nu){
  value = 0
  for(i in 1:length(p)){
    value[i] = multiroot(function(x) pCST(x,mu,sigma2,gama,nu)-p[i], start = 0)$root
  }
  return(value)
}
dCSCN = function(y,mu,sigma2,gama,nu){
  dens = 0
  for(i in 1:length(y)){
    dens[i] = nu[1]*dCSN(y[i],mu[i],sigma2/nu[2],gama)+(1-nu[1])*dCSN(y[i],mu[i],sigma2,gama)
  }
  return(dens)
}
pCSCN = function(q,mu,sigma2,gama,nu){
  value = cubintegrate(function(x) dCSCN(x,mu,sigma2,gama,nu),-Inf,q)$integral
  return(value)
}
qCSCN = function(p,mu,sigma2,gama,nu){
  value = 0
  for(i in 1:length(p)){
    value[i] = multiroot(function(x) pCSCN(x,mu,sigma2,gama,nu)-p[i], start = 0)$root
  }
  return(value)
}
dCSBPN = function(y,mu,sigma2,gama,nu){
  dens = 0
  for(i in 1:length(y)){
    dens[i] = integrate(function(x) dCSN(y[i],mu[i],sigma2/x,gama)*dbetapr(x,nu[1],nu[2]),0,Inf)$value
  }
  return(dens)
}
pCSBPN = function(q,mu,sigma2,gama,nu){
  value = cubintegrate(function(x) dCSBPN(x,mu,sigma2,gama,nu),-Inf,q)$integral
  return(value)
}
qCSBPN = function(p,mu,sigma2,gama,nu){
  value = 0
  for(i in 1:length(p)){
    value[i] = multiroot(function(x) pCSBPN(x,mu,sigma2,gama,nu)-p[i], start = 0)$root
  }
  return(value)
}
dCSBSN = function(y,mu,sigma2,gama,nu){
  dens = 0
  for(i in 1:length(y)){
    dens[i] = integrate(function(x) dCSN(y[i],mu[i],sigma2/x,gama)*dfatigue(x,nu[1],nu[2]),0,Inf)$value
  }
  return(dens)
}
pCSBSN = function(q,mu,sigma2,gama,nu){
  value = cubintegrate(function(x) dCSBSN(x,mu,sigma2,gama,nu),-Inf,q)$integral
  return(value)
}
qCSBSN = function(p,mu,sigma2,gama,nu){
  value = 0
  for(i in 1:length(p)){
    value[i] = multiroot(function(x) pCSBSN(x,mu,sigma2,gama,nu)-p[i], start = 0)$root
  }
  return(value)
}
dCSGGN = function(y,mu,gama,nu){
  dens = 0
  for(i in 1:length(y)){
    dens[i] = integrate(function(x) dCSN(y[i],mu[i],1/x,gama)*dggamma(x,nu[3],nu[1],nu[2]),
                        0,Inf, rel.tol = 1e-6)$value
  }
  return(dens)
}
pCSGGN = function(q,mu,gama,nu){
  value = cubintegrate(function(x) dCSGGN(x,mu,gama,nu),-Inf,q)$integral
  return(value)
}
qCSGGN = function(p,mu,gama,nu){
  value = 0
  for(i in 1:length(p)){
    value[i] = multiroot(function(x) pCSGGN(x,mu,gama,nu)-p[i], start = 0)$root
  }
  return(value)
}
deltak = function(k,c,W){
  if(k <= c*W){
    return(1)
  }else{
    return(min(1/(k-c*W),1))
  }
}
daux = function(x,nu1,nu2){
  if(x==nu2){
    return(nu1)
  }
  if(x==1){
    return(1-nu1)
  }
  return(0)
}
dBSp = function(x,nu1,nu2){
  dens = 1/(2*sqrt(2*pi)*nu1*nu2)*(sqrt(nu2/x)+(nu2/x)^(3/2))*exp(-(1/(2*nu1^2))*(x/nu2+nu2/x-2))
  return(dens)
}
dGGp = function(x,nu1,nu2,nu3){
  dens = (nu2/(nu3*gamma(nu1)))*((x/nu3)^(nu1*nu2-1))*exp(-((x/nu3)^nu2))
  return(dens)
}
ghudyCST = function(y,bu,mu,sigma2,gama,nu,m){
  sigma = sqrt(sigma2)
  g13 = nthroot(gama,3)
  g23 = nthroot(gama^2,3)
  lambda = (s*g13)/sqrt(2/pi+s^2*g23*(2/pi - 1))
  delta = lambda/sqrt(1+lambda^2)
  Delta = sigma*delta/sqrt(1-b^2*delta^2)
  tau = sigma2*(1-delta^2)/(1-b^2*delta^2)
  mt = tau/(Delta^2+tau)
  un = 0
  un[1] = bu
  count = 2
  hn = 0
  while(count <= m){
    #cand = rtruncnorm(1,a=0,b=Inf,mean=bu,sd = 1)
    #c1 = propCST(y,mu,sigma2,gama,nu,cand)*dtruncnorm(un[count-1], mean = bu, sd = 1)
    #c2 = propCST(y,mu,sigma2,gama,nu,un[count-1])*dtruncnorm(cand, mean = un[count-1], sd = 1)
    cand = rgamma(1,nu/2,nu/2)
    c1 = dCSN(y,mu,sigma2/cand,gama)
    c2 = dCSN(y,mu,sigma2/un[count-1],gama)
    alfa = c1/c2

    if(is.nan(alfa)==T||is.na(alfa)==T) {alfa=0.0001}
    if (runif(1) < min(alfa, 1)){
      un[count] = cand
    }else{
      un[count] = un[count-1]
    }
    count = count + 1
  }
  mit = (Delta^2*b+Delta*sqrt(un)*(y-mu))/(Delta^2+tau)
  hn = rtruncnorm(1,a=0,b=Inf,mean = mit, sd = sqrt(mt))
  return(cbind(hn,un))
}
# JV AQUI
ghudyCSGT = function(y,bu,mu,gama,nu,m){
  sigma = sigma2 = 1
  g13 = nthroot(gama,3)
  g23 = nthroot(gama^2,3)
  lambda = (s*g13)/sqrt(2/pi+s^2*g23*(2/pi - 1))
  delta = lambda/sqrt(1+lambda^2)
  Delta = sigma*delta/sqrt(1-b^2*delta^2)
  tau = sigma2*(1-delta^2)/(1-b^2*delta^2)
  mt = tau/(Delta^2+tau)
  un = 0
  un[1] = bu
  count = 2
  hn = 0
  while(count <= m){
    #cand = rtruncnorm(1,a=0,b=Inf,mean=bu,sd = 1)
    #c1 = propCST(y,mu,sigma2,gama,nu,cand)*dtruncnorm(un[count-1], mean = bu, sd = 1)
    #c2 = propCST(y,mu,sigma2,gama,nu,un[count-1])*dtruncnorm(cand, mean = un[count-1], sd = 1)
    cand = rgamma(1,nu[1]/2,nu[2]/2)
    c1 = dCSN(y,mu,sigma2/cand,gama)
    c2 = dCSN(y,mu,sigma2/un[count-1],gama)
    alfa = c1/c2

    if(is.nan(alfa)==T||is.na(alfa)==T) {alfa=0.0001}
    if (runif(1) < min(alfa, 1)){
      un[count] = cand
    }else{
      un[count] = un[count-1]
    }
    count = count + 1
  }
  mit = (Delta^2*b+Delta*sqrt(un)*(y-mu))/(Delta^2+tau)
  hn = rtruncnorm(1,a=0,b=Inf,mean = mit, sd = sqrt(mt))
  return(cbind(hn,un))
}
ghudyCSS = function(y,bu,mu,sigma2,gama,nu,m){
  sigma = sqrt(sigma2)
  g13 = nthroot(gama,3)
  g23 = nthroot(gama^2,3)
  lambda = (s*g13)/sqrt(2/pi+s^2*g23*(2/pi - 1))
  delta = lambda/sqrt(1+lambda^2)
  Delta = sigma*delta/sqrt(1-b^2*delta^2)
  tau = sigma2*(1-delta^2)/(1-b^2*delta^2)
  mt = tau/(Delta^2+tau)
  un = 0
  un[1] = bu
  count = 2
  hn = 0
  while(count <= m){
    #cand = rtruncnorm(1,a=0,b=Inf,mean=bu,sd = 1)
    #c1 = propCST(y,mu,sigma2,gama,nu,cand)*dtruncnorm(un[count-1], mean = bu, sd = 1)
    #c2 = propCST(y,mu,sigma2,gama,nu,un[count-1])*dtruncnorm(cand, mean = un[count-1], sd = 1)
    cand = rbeta(1,nu,1)
    c1 = dCSN(y,mu,sigma2/cand,gama)
    c2 = dCSN(y,mu,sigma2/un[count-1],gama)
    alfa = c1/c2

    if(is.nan(alfa)==T||is.na(alfa)==T) {alfa=0.0001}
    if (runif(1) < min(alfa, 1)){
      un[count] = cand
    }else{
      un[count] = un[count-1]
    }
    count = count + 1
  }
  mit = (Delta^2*b+Delta*sqrt(un)*(y-mu))/(Delta^2+tau)
  hn = rtruncnorm(1,a=0,b=Inf,mean = mit, sd = sqrt(mt))
  return(cbind(hn,un))
}
ghudyCSCN = function(y,bu,mu,sigma2,gama,nu,m){
  sigma = sqrt(sigma2)
  g13 = nthroot(gama,3)
  g23 = nthroot(gama^2,3)
  lambda = (s*g13)/sqrt(2/pi+s^2*g23*(2/pi - 1))
  delta = lambda/sqrt(1+lambda^2)
  Delta = sigma*delta/sqrt(1-b^2*delta^2)
  tau = sigma2*(1-delta^2)/(1-b^2*delta^2)
  mt = tau/(Delta^2+tau)
  un = 0
  un[1] = bu
  count = 2
  hn = 0
  while(count <= m){
    #cand = rtruncnorm(1,a=0,b=Inf,mean=bu,sd = 1)
    #c1 = propCST(y,mu,sigma2,gama,nu,cand)*dtruncnorm(un[count-1], mean = bu, sd = 1)
    #c2 = propCST(y,mu,sigma2,gama,nu,un[count-1])*dtruncnorm(cand, mean = un[count-1], sd = 1)
    cand = rC(1,nu)
    c1 = dCSN(y,mu,sigma2/cand,gama)
    c2 = dCSN(y,mu,sigma2/un[count-1],gama)
    alfa = c1/c2

    if(is.nan(alfa)==T||is.na(alfa)==T) {alfa=0.0001}
    if (runif(1) < min(alfa, 1)){
      un[count] = cand
    }else{
      un[count] = un[count-1]
    }
    count = count + 1
  }
  mit = (Delta^2*b+Delta*sqrt(un)*(y-mu))/(Delta^2+tau)
  hn = rtruncnorm(m,a=0,b=Inf,mean = mit, sd = sqrt(mt))
  return(cbind(hn,un))
}
ghudyCSN = function(y,bu,mu,sigma2,gama,m){
  sigma = sqrt(sigma2)
  g13 = nthroot(gama,3)
  g23 = nthroot(gama^2,3)
  lambda = (s*g13)/sqrt(2/pi+s^2*g23*(2/pi - 1))
  delta = lambda/sqrt(1+lambda^2)
  Delta = sigma*delta/sqrt(1-b^2*delta^2)
  tau = sigma2*(1-delta^2)/(1-b^2*delta^2)
  mt = tau/(Delta^2+tau)
  un = rep(1,m)
  hn = 0
  mit = (Delta^2*b+Delta*sqrt(un)*(y-mu))/(Delta^2+tau)
  hn = rtruncnorm(1,a=0,b=Inf,mean = mit, sd = sqrt(mt))
  return(cbind(hn,un))
}
ghudyCSBPN = function(y,bu,mu,sigma2,gama,nu,m){
  sigma = sqrt(sigma2)
  g13 = nthroot(gama,3)
  g23 = nthroot(gama^2,3)
  lambda = (s*g13)/sqrt(2/pi+s^2*g23*(2/pi - 1))
  delta = lambda/sqrt(1+lambda^2)
  Delta = sigma*delta/sqrt(1-b^2*delta^2)
  tau = sigma2*(1-delta^2)/(1-b^2*delta^2)
  mt = tau/(Delta^2+tau)
  un = 0
  un[1] = bu
  count = 2
  hn = 0
  while(count <= m){
    #cand = rtruncnorm(1,a=0,b=Inf,mean=bu,sd = 1)
    #c1 = propCST(y,mu,sigma2,gama,nu,cand)*dtruncnorm(un[count-1], mean = bu, sd = 1)
    #c2 = propCST(y,mu,sigma2,gama,nu,un[count-1])*dtruncnorm(cand, mean = un[count-1], sd = 1)
    cand = rbetapr(1,nu[1],nu[2])
    c1 = dCSN(y,mu,sigma2/cand,gama)
    c2 = dCSN(y,mu,sigma2/un[count-1],gama)
    alfa = c1/c2

    if(is.nan(alfa)==T||is.na(alfa)==T) {alfa=0.0001}
    if (runif(1) < min(alfa, 1)){
      un[count] = cand
    }else{
      un[count] = un[count-1]
    }
    count = count + 1
  }
  mit = (Delta^2*b+Delta*sqrt(un)*(y-mu))/(Delta^2+tau)
  hn = rtruncnorm(1,a=0,b=Inf,mean = mit, sd = sqrt(mt))
  return(cbind(hn,un))
}
ghudyCSBSN = function(y,bu,mu,sigma2,gama,nu,m){
  sigma = sqrt(sigma2)
  g13 = nthroot(gama,3)
  g23 = nthroot(gama^2,3)
  lambda = (s*g13)/sqrt(2/pi+s^2*g23*(2/pi - 1))
  delta = lambda/sqrt(1+lambda^2)
  Delta = sigma*delta/sqrt(1-b^2*delta^2)
  tau = sigma2*(1-delta^2)/(1-b^2*delta^2)
  mt = tau/(Delta^2+tau)
  un = 0
  un[1] = bu
  count = 2
  hn = 0
  while(count <= m){
    #cand = rtruncnorm(1,a=0,b=Inf,mean=bu,sd = 1)
    #c1 = propCST(y,mu,sigma2,gama,nu,cand)*dtruncnorm(un[count-1], mean = bu, sd = 1)
    #c2 = propCST(y,mu,sigma2,gama,nu,un[count-1])*dtruncnorm(cand, mean = un[count-1], sd = 1)
    cand = rfatigue(1,nu[1],nu[2])
    c1 = dCSN(y,mu,sigma2/cand,gama)
    c2 = dCSN(y,mu,sigma2/un[count-1],gama)
    alfa = c1/c2

    if(is.nan(alfa)==T||is.na(alfa)==T) {alfa=0.0001}
    if (runif(1) < min(alfa, 1)){
      un[count] = cand
    }else{
      un[count] = un[count-1]
    }
    count = count + 1
  }
  mit = (Delta^2*b+Delta*sqrt(un)*(y-mu))/(Delta^2+tau)
  hn = rtruncnorm(1,a=0,b=Inf,mean = mit, sd = sqrt(mt))
  return(cbind(hn,un))
}
ghudyCSGGN = function(y,bu,mu,gama,nu,m){
  sigma2 = sigma = 1
  g13 = nthroot(gama,3)
  g23 = nthroot(gama^2,3)
  lambda = (s*g13)/sqrt(2/pi+s^2*g23*(2/pi - 1))
  delta = lambda/sqrt(1+lambda^2)
  Delta = sigma*delta/sqrt(1-b^2*delta^2)
  tau = sigma2*(1-delta^2)/(1-b^2*delta^2)
  mt = tau/(Delta^2+tau)
  un = 0
  un[1] = bu
  count = 2
  hn = 0
  while(count <= m){
    #cand = rtruncnorm(1,a=0,b=Inf,mean=bu,sd = 1)
    #c1 = propCST(y,mu,sigma2,gama,nu,cand)*dtruncnorm(un[count-1], mean = bu, sd = 1)
    #c2 = propCST(y,mu,sigma2,gama,nu,un[count-1])*dtruncnorm(cand, mean = un[count-1], sd = 1)
    cand = rggamma(1,nu[3],nu[1],nu[2])
    c1 = dCSN(y,mu,sigma2/cand,gama)
    c2 = dCSN(y,mu,sigma2/un[count-1],gama)
    alfa = c1/c2

    if(is.nan(alfa)==T||is.na(alfa)==T) {alfa=0.0001}
    if (runif(1) < min(alfa, 1)){
      un[count] = cand
    }else{
      un[count] = un[count-1]
    }
    count = count + 1
  }
  mit = (Delta^2*b+Delta*sqrt(un)*(y-mu))/(Delta^2+tau)
  hn = rtruncnorm(1,a=0,b=Inf,mean = mit, sd = sqrt(mt))
  return(cbind(hn,un))
}
ghudy = function(y,bu,mu,sigma2,gama,nu,family,m){
  sigma = sqrt(sigma2)
  g13 = nthroot(gama,3)
  g23 = nthroot(gama^2,3)
  lambda = (s*g13)/sqrt(2/pi+s^2*g23*(2/pi - 1))
  delta = lambda/sqrt(1+lambda^2)
  Delta = sigma*delta/sqrt(1-b^2*delta^2)
  tau = sigma2*(1-delta^2)/(1-b^2*delta^2)
  mt = tau/(Delta^2+tau)
  xi = mu-sqrt(sigma2)*g13*s
  om = sqrt(sigma2)*sqrt(1+s^2*g23)
  un = 0
  un[1] = bu
  hn = 0
  for(i in 2:m){
    if(family == "CSS"){
      us = rtruncnorm(1,a=0,b=Inf,mean=un[i-1],sd = 10)
      uu = runif(1)
      aux = (dCSN(y,mu,sigma2/us,gama)*dbeta(us,nu,1)/dtruncnorm(us,a=0,b=Inf,mean=un[i-1],sd = 10))/(dCSN(y,mu,sigma2/un[i-1],gama)*dbeta(un[i-1],nu,1)/dtruncnorm(un[i-1],a=0,b=Inf,mean=un[i-1],sd = 10))
    }
    if(family == "CSBPN"){
      us = rtruncnorm(1,a=0,b=Inf,mean=un[i-1],sd = 10)
      uu = runif(1)
      aux = (dCSN(y,mu,sigma2/us,gama)*dbetapr(us,nu[1],nu[2])/dtruncnorm(us,a=0,b=Inf,mean=un[i-1],sd = 10))/(dCSN(y,mu,sigma2/un[i-1],gama)*dbetapr(un[i-1],nu[1],nu[2])/dtruncnorm(un[i-1],a=0,b=Inf,mean=un[i-1],sd = 10))
    }
    if(family == "CSBSN"){
      us = rtruncnorm(1,a=0,b=Inf,mean=un[i-1],sd = 10)
      uu = runif(1)
      aux = (dCSN(y,mu,sigma2/us,gama)*dBSp(us,nu[1],nu[2])/dtruncnorm(us,a=0,b=Inf,mean=un[i-1],sd = 10))/(dCSN(y,mu,sigma2/un[i-1],gama)*dBSp(un[i-1],nu[1],nu[2])/dtruncnorm(un[i-1],a=0,b=Inf,mean=un[i-1],sd = 10))
    }
    if(family == "CSGGN"){
      us = rtruncnorm(1,a=0,b=Inf,mean=un[i-1],sd = 10)
      uu = runif(1)
      aux = (dCSN(y,mu,sigma2/us,gama)*dGGp(us,nu[1],nu[2],nu[3])/dtruncnorm(us,a=0,b=Inf,mean=un[i-1],sd = 10))/(dCSN(y,mu,sigma2/un[i-1],gama)*dGGp(un[i-1],nu[1],nu[2],nu[3])/dtruncnorm(un[i-1],a=0,b=Inf,mean=un[i-1],sd = 10))
    }
    if(family == "CSN"){
      un[i] = 1
    }
    if(family == "CSCN"){
      un[i] = 0.5
    }
    if(family != "CSN" && family!= "CSCN" && !is.nan(aux)){
      if(uu <= min(1,aux)){
        un[i] = us
      }else{
        un[i] = un[i-1]
      }
    }
    if(family != "CSN" && family!= "CSCN" && !is.nan(aux)){
      un[i] = un[i-1]
    }
    mit = (Delta^2*b+Delta*sqrt(un[i])*(y-mu))/(Delta^2+tau)
    hn[i] = rtruncnorm(1,a=0,b=Inf,mean = mit, sd = sqrt(mt))
  }
  return(cbind(hn,un))
}
rCSN = function(n, mu, sigma2, shape, nu=NULL){
  g13 = nthroot(shape,3)
  g23 = g13^2
  s = nthroot(2/(4-pi),3)
  lamb = (s*g13)/sqrt(2/pi+(s^2)*g23*(2/pi-1))
  xi = mu-sqrt(sigma2)*g13*s
  om = sqrt(sigma2)*sqrt(1+s^2*g23)
  delta = lamb/sqrt(1+lamb^2)
  y = xi*rep(1,n) + om*(delta*abs(rnorm(n)) + (1 - delta^2)^(1/2)*rnorm(n))
  return(y)
}
rCST = function(n, mu, sigma2, shape, nu ){
  u = rgamma(n, nu/2, nu/2)
  y = mu + rCSN(n,0,sigma2/u,shape)
  return(y)
}
rCSS = function(n, mu, sigma2, shape, nu ){
  u = rbeta(n,nu,1)
  y = mu + (1/sqrt(u))*rCSN(n,0,sigma2,shape)
  return(y)
}
rCSGT = function(n, mu, shape, nu ){
  u = rgamma(n,nu[1]/2,nu[2]/2)
  y = mu + rCSN(n,0,1/u,shape)
  return(y)
}
rCSCN = function(n, mu, sigma2, shape, nu ){
  us = runif(n)
  u = rep(0,n)
  for(i in 1:n){
    if(us[i]<nu[1]){
      u[i] = nu[2]
    }else{
      u[i] = 1
    }
  }
  y = mu + rCSN(n,0,sigma2/u,shape)
  return(y)
}
rC = function(n,nu ){
  us = runif(n)
  u = rep(0,n)
  for(i in 1:n){
    if(us[i]<nu[1]){
      u[i] = nu[2]
    }else{
      u[i] = 1
    }
  }
  return(u)
}
rCSBPN = function(n, mu, sigma2, shape, nu ){
  u = rbetapr(n,nu[1],nu[2])
  y = mu + rCSN(n,0,sigma2/u,shape)
  return(y)
}
rCSBSN = function(n, mu, sigma2, shape, nu ){
  u = rfatigue(n,nu[1],nu[2])
  y = mu + rCSN(n,0,sigma2/u,shape)
  return(y)
}
rCSGGN = function(n, mu, shape, nu ){
  u = rggamma(n,nu[3],nu[1],nu[2])
  y = mu + rCSN(n,0,1/u,shape)
  return(y)
}
csnm = function(y, X, N, K, alpha0 = 1, nknot, tol = 0.001,
                iter.max = 300, k = 1/3, alpha.fix = F, a0, a1){
  p = ncol(X) # Número de parâmetros
  n = nrow(X) # Número de unidades experimentais
  reg = lm(y~X)

  #=== Initial values

  beta   = solve(t(X)%*%X)%*%t(X)%*%y
  sigma2 = as.numeric(t((y-X%*%beta))%*%(y-X%*%beta)/(n-2))
  sigma = sqrt(sigma2)
  lambi = (sum((y - X%*%beta)^3)/n)/(sum((y - X%*%beta)^2)/n)^(3/2)
  delti = lambi/sqrt(1+lambi^2)
  gama = ((4-pi)/2)*((sqrt(2/pi)*delti)^3)/((1-(2/pi)*(delti^2))^(3/2))
  print(gama)
  theta = list(NULL)
  prednp = list(NULL)
  for(i in 1:length(N)){
    theta[[i]]  = solve(t(N[[i]])%*%N[[i]]+alpha0[i]*K[[i]])%*%t(N[[i]])%*%(y-X%*%beta)
    prednp[[i]] = N[[i]]%*%theta[[i]]
  }
  mu = X%*%beta+Reduce("+",prednp)

  S2o = 0
  S3o = 0
  criterio = 1
  count    = 1
  m = 20
  betao = matrix(0,iter.max,p)
  sigma2o = rep(0,iter.max)
  gamao = rep(0,iter.max)
  theta0 = list(NULL)
  alpha = alpha0
  pb = txtProgressBar(min = 0, max = iter.max, initial = 1)
  while((count <= iter.max)){
    g13 = nthroot(gama,3)
    g23 = nthroot(gama^2,3)
    lambda = (s*g13)/sqrt(2/pi+s^2*g23*(2/pi - 1))
    delta = lambda/sqrt(1+lambda^2)
    Delta = sigma*delta/sqrt(1-b^2*delta^2)
    tau = sigma2*(1-delta^2)/(1-b^2*delta^2)
    mt = tau/(Delta^2+tau)
    esq = matrix(0,n,m)
    e10 = matrix(0,n,m)
    e51 = matrix(0,n,m)
    e02 = matrix(0,n,m)
    for(i in 1:n){
      mit = (Delta^2*b+Delta*(y[i]-mu[i]))/(Delta^2+tau)
      eh1 = rtruncnorm(1,a=0,b=Inf,mean = mit, sd = sqrt(mt))
      e02[i,] = (eh1-b)^2
      e51[i,] = (eh1-b)
    }
    ss2n = 0
    ss3n = 0
    for(i in 1:m){
      ss2n = ss2n + e02[,m]
      ss3n = ss3n + e51[,m]
    }
    S2n = S2o + deltak(count,k,iter.max)*((1/m)*ss2n-S2o)
    S3n = S3o + deltak(count,k,iter.max)*((1/m)*ss3n-S3o)
    Delta0 = (1/sum(S2n))*sum((y-mu)*S3n)
    tau0 = (1/n)*sum((y-mu)^2-2*(y-mu)*Delta0*S3n+Delta0^2*S2n)
    beta0 = solve(t(X)%*%X)%*%t(X)%*%(y-Reduce("+",prednp)-S3n*Delta0)
    if(length(N)>1){
      for(i in 1:length(N)){
        theta0[[i]] = solve(as.matrix(t(N[[i]])%*%N[[i]]+alpha[i]*tau0*K[[i]]))%*%t(N[[i]])%*%(y-X%*%beta0-Reduce("+",prednp[-i])-S3n*Delta0)
        prednp[[i]] = N[[i]]%*%theta0[[i]]
      }
    }else{
      theta0[[1]] = solve(as.matrix(t(N[[1]])%*%N[[1]]+alpha[1]*tau0*K[[1]]))%*%t(N[[1]])%*%(y-X%*%beta0-S3n*Delta0)
    }
    s20 = tau0+Delta0^2*(1-b^2)
    delta0 = Delta0/sqrt(tau0+Delta0^2)
    gama0 = ((4-pi)/2)*b^3*(delta0^3)/(sqrt((1-(2/pi)*delta0^2)^3))
    for(i in 1:length(N)){
      prednp[[i]] = N[[i]]%*%theta0[[i]]
    }
    mu0 = X%*%beta0+Reduce("+",prednp)
    ver_novo = sum(log(dCSN(y, mu0, s20, gama0)))
    if(alpha.fix == F){
      ft = function(alpha0) BICsmn(alpha = alpha0,logver = ver_novo,
                                    ff = theta0, K = K, N = N, tau = tau0,
                                    p = p, e10 = diag(1,n),
                                    type = "CSN")
      alpha = optim(alpha, ft, lower = a0, upper = a1,
                     method = "L-BFGS-B")$par
    }
    betao[count,] = beta0
    sigma2o[count] = s20
    gamao[count] = gama0
    count = (count+1)
    beta = beta0
    sigma2 = s20
    sigma = sqrt(sigma2)
    theta = theta0
    gama = gama0
    mu  = mu0
    S2o = S2n
    S3o = S3n
    progress(count,iter.max, progress.bar = T)
  }

  fit = list()
  attr(fit, "class") = c("csnm")
  fit$title = "csnm:  CENTERED SKEW-NORMAL MODEL"
  fit$model = list()
  fit$call = call
  fit$n = n
  fit$beta = beta
  mu = as.vector(mu)
  fit$fitted.values = mu
  fit$family = "Centered Skew-Normal"
  fit$y = as.vector(y)
  fit$sigma2 = sigma2
  fit$gama = gama
  fit$theta = theta
  fit$log_vero = ver_novo
  fit$comp$Delta = Delta
  fit$comp$tau = tau
  fit$comp$S2n = S2n
  fit$comp$S3n = S3n
  fit$comp$K = K
  fit$alpha = alpha
  fit$X = X
  fit$N = N
  fit$alpha = alpha
  fit$betao = betao
  fit$sigma2o = sigma2o
  fit$gamao = gamao
  return(fit)
}

SAEM.CSS = function(y,mu,theta,sigma2,gama,n,p,alpha0,K,N,X,k,iter.max,
                    alpha.fix,a0,a1){
  nllg  = function(nu,y,mu,sigma2,gama){sum(log(dCSS(y, mu, sigma2, gama, nu)))}
  nu = optim(par = 2, fn = nllg, y = y, mu = mu, sigma2 = sigma2, gama = gama,
             method = "L-BFGS-B", lower = 2, upper = 15, control = list(fnscale = -1))$par
  S1o = 0
  S2o = 0
  S3o = 0
  S4o = 0
  criterio = 1
  count    = 1
  bu = rbeta(n,nu,1)
  m = 20
  betao = matrix(0,iter.max,p)
  sigma2o = rep(0,iter.max)
  gamao = rep(0,iter.max)
  nuo = matrix(0,iter.max,1)
  prednp = list(NULL)
  for(i in 1:length(N)){
    prednp[[i]] = N[[i]]%*%theta[[i]]
  }
  sigma = sqrt(sigma2)
  theta0 = list(NULL)
  alpha = alpha0
  while((count <= iter.max)){
    g13 = nthroot(gama,3)
    g23 = g13^2
    lambda = (s*g13)/sqrt(2/pi+s^2*g23*(2/pi - 1))
    delta = lambda/sqrt(1+lambda^2)
    Delta = sigma*delta/sqrt(1-b^2*delta^2)
    tau = sigma2*(1-delta^2)/(1-b^2*delta^2)
    esq = matrix(0,n,m)
    e10 = matrix(0,n,m)
    e51 = matrix(0,n,m)
    e02 = matrix(0,n,m)
    for(i in 1:n){
      aux2 = ghudyCSS(y[i],bu[i],mu[i],sigma2,gama,nu,m)
      eh1 = aux2[,1]
      e02[i,] = (eh1-b)^2
      e10[i,] = aux2[,2]
      esq[i,] = sqrt(e10[i,])
      e51[i,] = esq[i,]*(eh1-b)
    }
    ss1n = 0
    ss2n = 0
    ss3n = 0
    for(i in 1:m){
      ss1n = ss1n + e10[,m]
      ss2n = ss2n + e02[,m]
      ss3n = ss3n + e51[,m]
    }
    S1n = S1o + deltak(count,k,iter.max)*((1/m)*ss1n-S1o)
    S2n = S2o + deltak(count,k,iter.max)*((1/m)*ss2n-S2o)
    S3n = S3o + deltak(count,k,iter.max)*((1/m)*ss3n-S3o)
    Delta0 = (1/sum(S2n))*sum((y-mu)*S3n)
    tau0 = (1/n)*sum((y-mu)^2*S1n-2*(y-mu)*Delta0*S3n+Delta0^2*S2n)
    beta0 = solve(t(X)%*%diag(S1n)%*%X)%*%t(X)%*%(diag(S1n)%*%y-diag(S1n)%*%Reduce("+",prednp)-S3n*Delta0)
    if(length(N)>1){
      for(i in 1:length(N)){
        theta0[[i]] = solve(t(N[[i]])%*%diag(S1n)%*%N[[i]]+tau0*alpha[[i]]*K[[i]])%*%t(N[[i]])%*%(diag(S1n)%*%y-diag(S1n)%*%X%*%beta0-S3n*Delta0-diag(S1n)%*%Reduce("+",prednp[-i]))
        prednp[[i]] = N[[i]]%*%theta0[[i]]
      }
    }else{
      theta0[[1]] = solve(t(N[[1]])%*%diag(S1n)%*%N[[1]]+tau0*alpha[[1]]*K[[1]])%*%t(N[[1]])%*%(diag(S1n)%*%y-diag(S1n)%*%X%*%beta0-S3n*Delta0)
    }

    s20 = tau0+Delta0^2*(1-b^2)
    delta0 = Delta0/sqrt(tau0+Delta0^2)
    gama0 = ((4-pi)/2)*b^3*(delta0^3)/((1-(2/pi)*delta0^2)^(3/2))
    for(i in 1:length(N)){
      prednp[[i]] = N[[i]]%*%theta0[[i]]
    }
    mu0 = X%*%beta0+Reduce("+",prednp)
    nu0 = optim(par = nu, fn = nllg, y = y, mu = mu0, sigma2 = s20,
                gama = gama0, method = "L-BFGS-B", lower = 2, upper = 15,
                control = list(fnscale = -1))$par
    nu = nu0
    ver_novo = nllg(nu,y,mu0,s20,gama0)

    if(alpha.fix == F){
      ft = function(alpha0) BICsmn(alpha0,ver_novo,theta0,K,N,tau0,p,
                                    diag(S1n),"CSS")
      alpha = optim(alpha, ft, method = "L-BFGS-B", lower = a0,
                     upper = a1)$par
    }

    betao[count,] = beta0
    sigma2o[count] = s20
    gamao[count] = gama0
    nuo[count,] = nu
    count = (count+1)
    beta = beta0
    theta = theta0
    sigma2 = s20
    sigma = sqrt(sigma2)
    gama = gama0
    mu  = mu0
    bu = S1n
    S1o = S1n
    S2o = S2n
    S3o = S3n
    progress(count,iter.max, progress.bar = T)
  }
  return(list(beta = beta, theta = theta, sigma2 = sigma2, gama = gama, mu = mu, betao = betao,
              sigma2o = sigma2o, gamao = gamao, nuo = nuo, nu = nu, alpha = alpha, S1n = S1n,
              S2n = S2n, S3n = S3n))
}
SAEM.CSGT = function(y,mu,theta,gama,n,p,alpha0,K,N,X,k,iter.max, alpha.fix,
                     a0, a1){
  nllg  = function(nu,y,mu,gama){
    nu1=nu[1]
    nu2=nu[2]
    return(sum(log(dCSGT(y, mu,gama,c(nu1,nu2)))))
  }
  nu = optim(par = c(5,5), fn = nllg, y = y, mu = mu, gama = gama,
             method = "L-BFGS-B", lower = c(3,1), upper = c(15,15),
             control = list(fnscale = -1))$par

  S1o = 0
  S2o = 0
  S3o = 0
  S4o = 0
  criterio = 1
  count    = 1
  bu = rgamma(n,nu[1]/2,nu[2]/2)
  m = 30
  betao = matrix(0,iter.max,p)
  gamao = rep(0,iter.max)
  nuo = matrix(0,iter.max,2)
  sigma2 = sigma = 1
  prednp = list(NULL)
  for(i in 1:length(N)){
    prednp[[i]] = N[[i]]%*%theta[[i]]
  }
  theta0 = list(NULL)
  alpha = alpha0
  while((count <= iter.max)){
    g13 = nthroot(gama,3)
    g23 = g13^2
    lambda = (s*g13)/sqrt(2/pi+s^2*g23*(2/pi - 1))
    delta = lambda/sqrt(1+lambda^2)
    Delta = sigma*delta/sqrt(1-b^2*delta^2)
    tau = sigma2*(1-delta^2)/(1-b^2*delta^2)
    esq = matrix(0,n,m)
    e10 = matrix(0,n,m)
    e51 = matrix(0,n,m)
    e02 = matrix(0,n,m)
    for(i in 1:n){
      aux2 = ghudyCSGT(y[i],bu[i],mu[i],gama,nu,m)
      eh1 = aux2[,1]
      e02[i,] = (eh1-b)^2
      e10[i,] = aux2[,2]
      esq[i,] = sqrt(e10[i,])
      e51[i,] = esq[i,]*(eh1-b)
    }
    ss1n = 0
    ss2n = 0
    ss3n = 0
    for(i in 1:m){
      ss1n = ss1n + e10[,m]
      ss2n = ss2n + e02[,m]
      ss3n = ss3n + e51[,m]
    }
    S1n = S1o + deltak(count,k,iter.max)*((1/m)*ss1n-S1o)
    S2n = S2o + deltak(count,k,iter.max)*((1/m)*ss2n-S2o)
    S3n = S3o + deltak(count,k,iter.max)*((1/m)*ss3n-S3o)
    Delta0 = (1/sum(S2n))*sum((y-mu)*S3n)
    tau0 = (1/n)*sum((y-mu)^2*S1n-2*(y-mu)*Delta0*S3n+Delta0^2*S2n)
    beta0 = solve(t(X)%*%diag(S1n)%*%X)%*%t(X)%*%(diag(S1n)%*%y-diag(S1n)%*%Reduce("+",prednp)-S3n*Delta0)
    if(length(N)>1){
      for(i in 1:length(N)){
        theta0[[i]] = solve(t(N[[i]])%*%diag(S1n)%*%N[[i]]+tau0*alpha[[i]]*K[[i]])%*%t(N[[i]])%*%(diag(S1n)%*%y-diag(S1n)%*%X%*%beta0-S3n*Delta0-diag(S1n)%*%Reduce("+",prednp[-i]))
        prednp[[i]] = N[[i]]%*%theta0[[i]]
      }
    }else{
      theta0[[1]] = solve(t(N[[1]])%*%diag(S1n)%*%N[[1]]+tau0*alpha[[1]]*K[[1]])%*%t(N[[1]])%*%(diag(S1n)%*%y-diag(S1n)%*%X%*%beta0-S3n*Delta0)
    }
    delta0 = Delta0/sqrt(tau0+Delta0^2)
    gama0 = ((4-pi)/2)*b^3*(delta0^3)/((1-(2/pi)*delta0^2)^(3/2))
    for(i in 1:length(N)){
      prednp[[i]] = N[[i]]%*%theta0[[i]]
    }
    mu0 = X%*%beta0+Reduce("+",prednp)

    nu0 = optim(par = nu, fn = nllg, y = y, mu = mu0, gama = gama0,
                method = "L-BFGS-B", lower = c(3,1), upper = c(15,15),
                control = list(fnscale = -1))$par
    nu = nu0
    ver_novo = nllg(nu,y,mu0,gama0)

    if(alpha.fix == F){
      ft = function(alpha0) BICsmn(alpha0,ver_novo,theta0,K,N,tau0,p,
                                    diag(S1n),"CSGT")
      alpha = optim(alpha, ft, method = "L-BFGS-B", lower = a0,
                     upper = a1)$par
    }

    betao[count,] = beta0
    gamao[count] = gama0
    nuo[count,] = nu
    count = (count+1)
    beta = beta0
    theta = theta0
    gama = gama0
    mu  = mu0
    bu = S1n
    S1o = S1n
    S2o = S2n
    S3o = S3n
    progress(count,iter.max, progress.bar = T)
  }
  return(list(beta = beta, theta = theta, sigma2 = sigma2, gama = gama, mu = mu, betao = betao,
              gamao = gamao, nuo = nuo, nu = nu, alpha = alpha, delta = delta,
              S1n = S1n, S2n = S2n, S3n = S3n))
}

SAEM.CST = function(y,mu,theta,sigma2,gama,n,p,alpha0,K,N,X,k,iter.max,
                    alpha.fix,a0,a1){
  nllg  = function(nu,y,mu,sigma2,gama){sum(log(dCST(y, mu, sigma2, gama, nu)))}
  nu = optim(par = 5, fn = nllg, y = y, mu = mu, sigma2 = sigma2, gama = gama,
             method = "L-BFGS-B", lower = 2, upper = 15,
             control = list(fnscale = -1))$par
  S1o = 0
  S2o = 0
  S3o = 0
  S4o = 0
  criterio = 1
  count    = 1
  bu = rgamma(n,nu/2,nu/2)
  m = 30
  betao = matrix(0,iter.max,p)
  sigma2o = rep(0,iter.max)
  gamao = rep(0,iter.max)
  nuo = matrix(0,iter.max,1)
  prednp = list(NULL)
  for(i in 1:length(N)){
    prednp[[i]] = N[[i]]%*%theta[[i]]
  }
  sigma = sqrt(sigma2)
  theta0 = list(NULL)
  alpha = alpha0
  while((count <= iter.max)){
    g13 = nthroot(gama,3)
    g23 = g13^2
    lambda = (s*g13)/sqrt(2/pi+s^2*g23*(2/pi - 1))
    delta = lambda/sqrt(1+lambda^2)
    Delta = sigma*delta/sqrt(1-b^2*delta^2)
    tau = sigma2*(1-delta^2)/(1-b^2*delta^2)
    esq = matrix(0,n,m)
    e10 = matrix(0,n,m)
    e51 = matrix(0,n,m)
    e02 = matrix(0,n,m)
    for(i in 1:n){
      aux2 = ghudyCST(y[i],bu[i],mu[i],sigma2,gama,nu,m)
      eh1 = aux2[,1]
      e02[i,] = (eh1-b)^2
      e10[i,] = aux2[,2]
      esq[i,] = sqrt(e10[i,])
      e51[i,] = esq[i,]*(eh1-b)
    }
    ss1n = 0
    ss2n = 0
    ss3n = 0
    for(i in 1:m){
      ss1n = ss1n + e10[,m]
      ss2n = ss2n + e02[,m]
      ss3n = ss3n + e51[,m]
    }
    S1n = S1o + deltak(count,k,iter.max)*((1/m)*ss1n-S1o)
    S2n = S2o + deltak(count,k,iter.max)*((1/m)*ss2n-S2o)
    S3n = S3o + deltak(count,k,iter.max)*((1/m)*ss3n-S3o)
    Delta0 = (1/sum(S2n))*sum((y-mu)*S3n)
    tau0 = (1/n)*sum((y-mu)^2*S1n-2*(y-mu)*Delta0*S3n+Delta0^2*S2n)
    beta0 = solve(t(X)%*%diag(S1n)%*%X)%*%t(X)%*%(diag(S1n)%*%y-diag(S1n)%*%Reduce("+",prednp)-S3n*Delta0)
    if(length(N)>1){
      for(i in 1:length(N)){
        theta0[[i]] = solve(t(N[[i]])%*%diag(S1n)%*%N[[i]]+tau0*alpha[[i]]*K[[i]])%*%t(N[[i]])%*%(diag(S1n)%*%y-diag(S1n)%*%X%*%beta0-S3n*Delta0-diag(S1n)%*%Reduce("+",prednp[-i]))
        prednp[[i]] = N[[i]]%*%theta0[[i]]
      }
    }else{
      theta0[[1]] = solve(t(N[[1]])%*%diag(S1n)%*%N[[1]]+tau0*alpha[[1]]*K[[1]])%*%t(N[[1]])%*%(diag(S1n)%*%y-diag(S1n)%*%X%*%beta0-S3n*Delta0)
    }
    s20 = tau0+Delta0^2*(1-b^2)
    delta0 = Delta0/sqrt(tau0+Delta0^2)
    gama0 = ((4-pi)/2)*b^3*(delta0^3)/((1-(2/pi)*delta0^2)^(3/2))
    for(i in 1:length(N)){
      prednp[[i]] = N[[i]]%*%theta0[[i]]
    }
    mu0 = X%*%beta0+Reduce("+",prednp)

    nu0 = optim(par = nu, fn = nllg, y = y, mu = mu0, sigma2 = s20,
                gama = gama0, method = "L-BFGS-B", lower = 2, upper = 15,
                control = list(fnscale = -1))$par
    nu = nu0
    ver_novo = nllg(nu,y,mu0,s20,gama0)

    if(alpha.fix == F){
      ft = function(alpha0) BICsmn(alpha0,ver_novo,theta0,K,N,tau0,p,
                                    diag(S1n),"CST")
      alpha = optim(alpha, ft, lower = a0, upper = a1,
                     method = "L-BFGS-B")$par
    }

    betao[count,] = beta0
    sigma2o[count] = s20
    gamao[count] = gama0
    nuo[count] = nu
    count = (count+1)
    beta = beta0
    theta = theta0
    sigma2 = s20
    sigma = sqrt(sigma2)
    gama = gama0
    mu  = mu0
    bu = S1n
    S1o = S1n
    S2o = S2n
    S3o = S3n
    progress(count,iter.max, progress.bar = T)
  }
  return(list(beta = beta, theta = theta, sigma2 = sigma2, gama = gama, mu = mu, betao = betao,
              sigma2o = sigma2o, gamao = gamao, nuo = nuo, nu = nu, alpha = alpha,
              delta = delta0,S1n = S1n, S2n = S2n, S3n = S3n))
}

SAEM.CSCN = function(y,mu,theta,sigma2,gama,n,p,alpha0,K,N,X,k,iter.max,
                      alpha.fix,a0,a1){
  nllg  = function(nu,y,mu,sigma2,gama){
    nu1=nu[1]
    nu2=nu[2]
    return(sum(log(dCSCN(y, mu, sigma2, gama,c(nu1,nu2)))))
  }
  nu = optim(par = c(0.5,0.5), fn = nllg, y = y, mu = mu, sigma2 = sigma2,
             gama = gama, method = "L-BFGS-B", lower = c(0.1,0.1),
             upper = c(0.9,0.9), control = list(fnscale = -1))$par
  S1o = 0
  S2o = 0
  S3o = 0
  S4o = 0
  criterio = 1
  count    = 1
  bu = rep(0.5,n)
  m = 30
  betao = matrix(0,iter.max,p)
  sigma2o = rep(0,iter.max)
  gamao = rep(0,iter.max)
  nuo = matrix(0,iter.max,2)
  prednp = list(NULL)
  for(i in 1:length(N)){
    prednp[[i]] = N[[i]]%*%theta[[i]]
  }
  sigma = sqrt(sigma2)
  theta0 = list(NULL)
  alpha = alpha0
  while((count <= iter.max)){
    g13 = nthroot(gama,3)
    g23 = g13^2
    lambda = (s*g13)/sqrt(2/pi+s^2*g23*(2/pi - 1))
    delta = lambda/sqrt(1+lambda^2)
    Delta = sigma*delta/sqrt(1-b^2*delta^2)
    tau = sigma2*(1-delta^2)/(1-b^2*delta^2)
    esq = matrix(0,n,m)
    e10 = matrix(0,n,m)
    e51 = matrix(0,n,m)
    e02 = matrix(0,n,m)
    for(i in 1:n){
      aux2 = ghudyCSCN(y[i],bu[i],mu[i],sigma2,gama,nu,m)
      eh1 = aux2[,1]
      e02[i,] = (eh1-b)^2
      e10[i,] = aux2[,2]
      esq[i,] = sqrt(e10[i,])
      e51[i,] = esq[i,]*(eh1-b)
    }
    ss1n = 0
    ss2n = 0
    ss3n = 0
    for(i in 1:m){
      ss1n = ss1n + e10[,m]
      ss2n = ss2n + e02[,m]
      ss3n = ss3n + e51[,m]
    }
    S1n = S1o + deltak(count,k,iter.max)*((1/m)*ss1n-S1o)
    S2n = S2o + deltak(count,k,iter.max)*((1/m)*ss2n-S2o)
    S3n = S3o + deltak(count,k,iter.max)*((1/m)*ss3n-S3o)
    Delta0 = (1/sum(S2n))*sum((y-mu)*S3n)
    tau0 = (1/n)*sum((y-mu)^2*S1n-2*(y-mu)*Delta0*S3n+Delta0^2*S2n)
    beta0 = solve(t(X)%*%diag(S1n)%*%X)%*%t(X)%*%(diag(S1n)%*%y-diag(S1n)%*%Reduce("+",prednp)-S3n*Delta0)
    if(length(N)>1){
      for(i in 1:length(N)){
        theta0[[i]] = solve(t(N[[i]])%*%diag(S1n)%*%N[[i]]+tau0*alpha[[i]]*K[[i]])%*%t(N[[i]])%*%(diag(S1n)%*%y-diag(S1n)%*%X%*%beta0-S3n*Delta0-diag(S1n)%*%Reduce("+",prednp[-i]))
        prednp[[i]] = N[[i]]%*%theta0[[i]]
      }
    }else{
      theta0[[1]] = solve(t(N[[1]])%*%diag(S1n)%*%N[[1]]+tau0*alpha[[1]]*K[[1]])%*%t(N[[1]])%*%(diag(S1n)%*%y-diag(S1n)%*%X%*%beta0-S3n*Delta0)
    }

    s20 = tau0+Delta0^2*(1-b^2)
    delta0 = Delta0/sqrt(tau0+Delta0^2)
    gama0 = ((4-pi)/2)*b^3*(delta0^3)/((1-(2/pi)*delta0^2)^(3/2))
    for(i in 1:length(N)){
      prednp[[i]] = N[[i]]%*%theta0[[i]]
    }
    mu0 = X%*%beta0+Reduce("+",prednp)

    nu0 = optim(par = nu, fn = nllg, y = y, mu = mu0, sigma2 = s20,
                gama = gama0, method = "L-BFGS-B", lower = c(0.1,0.1),
                upper = c(0.9,0.9), control = list(fnscale = -1))$par
    nu = nu0
    ver_novo = nllg(nu,y,mu0,s20,gama0)

    if(alpha.fix == F){
      ft = function(alpha0) BICsmn(alpha0,ver_novo,theta0,K,N,tau0,p,
                                    diag(S1n),"CSCN")
      alpha = optim(alpha, ft, method = "L-BFGS-B", lower = a0,
                     upper = a1)$par
    }

    betao[count,] = beta0
    sigma2o[count] = s20
    gamao[count] = gama0
    nuo[count,] = nu
    count = (count+1)
    beta = beta0
    theta = theta0
    sigma2 = s20
    sigma = sqrt(sigma2)
    gama = gama0
    mu  = mu0
    bu = S1n
    S1o = S1n
    S2o = S2n
    S3o = S3n
    progress(count,iter.max, progress.bar = T)
  }
  return(list(beta = beta, theta = theta, sigma2 = sigma2, gama = gama, mu = mu, betao = betao,
              sigma2o = sigma2o, gamao = gamao, nuo = nuo, nu = nu, alpha = alpha,
              delta = delta0,S1n = S1n, S2n = S2n, S3n = S3n))
}

SAEM.CSBPN = function(y,mu,theta,sigma2,gama,n,p,alpha0,K,N,X,k,iter.max,
                     alpha.fix,a0,a1){
  nllg  = function(nu,y,mu,sigma2,gama){
    nu1=nu[1]
    nu2=nu[2]
    return(sum(log(dCSBPN(y, mu, sigma2,gama,c(nu1,nu2)))))
  }
  nu = optim(par = c(5,5), fn = nllg, y = y, mu = mu, sigma2 = sigma2, gama = gama,
             method = "L-BFGS-B", lower = c(2,1), upper = c(5,5),
             control = list(fnscale = -1))$par

  S1o = 0
  S2o = 0
  S3o = 0
  S4o = 0
  criterio = 1
  count    = 1
  bu = rbetapr(n,nu[1],nu[2])
  m = 20
  betao = matrix(0,iter.max,p)
  sigma2o = rep(0,iter.max)
  gamao = rep(0,iter.max)
  nuo = matrix(0,iter.max,2)
  prednp = list(NULL)
  for(i in 1:length(N)){
    prednp[[i]] = N[[i]]%*%theta[[i]]
  }
  sigma = sqrt(sigma2)
  theta0 = list(NULL)
  alpha = alpha0
  while((count <= iter.max)){
    g13 = nthroot(gama,3)
    g23 = g13^2
    lambda = (s*g13)/sqrt(2/pi+s^2*g23*(2/pi - 1))
    delta = lambda/sqrt(1+lambda^2)
    Delta = sigma*delta/sqrt(1-b^2*delta^2)
    tau = sigma2*(1-delta^2)/(1-b^2*delta^2)
    esq = matrix(0,n,m)
    e10 = matrix(0,n,m)
    e51 = matrix(0,n,m)
    e02 = matrix(0,n,m)
    for(i in 1:n){
      aux2 = ghudyCSBPN(y[i],bu[i],mu[i],sigma2,gama,nu,m)
      eh1 = aux2[,1]
      e02[i,] = (eh1-b)^2
      e10[i,] = aux2[,2]
      esq[i,] = sqrt(e10[i,])
      e51[i,] = esq[i,]*(eh1-b)
    }
    ss1n = 0
    ss2n = 0
    ss3n = 0
    for(i in 1:m){
      ss1n = ss1n + e10[,m]
      ss2n = ss2n + e02[,m]
      ss3n = ss3n + e51[,m]
    }
    S1n = S1o + deltak(count,k,iter.max)*((1/m)*ss1n-S1o)
    S2n = S2o + deltak(count,k,iter.max)*((1/m)*ss2n-S2o)
    S3n = S3o + deltak(count,k,iter.max)*((1/m)*ss3n-S3o)
    Delta0 = (1/sum(S2n))*sum((y-mu)*S3n)
    tau0 = (1/n)*sum((y-mu)^2*S1n-2*(y-mu)*Delta0*S3n+Delta0^2*S2n)
    beta0 = solve(t(X)%*%diag(S1n)%*%X)%*%t(X)%*%(diag(S1n)%*%y-diag(S1n)%*%Reduce("+",prednp)-S3n*Delta0)
    if(length(N)>1){
      for(i in 1:length(N)){
        theta0[[i]] = solve(t(N[[i]])%*%diag(S1n)%*%N[[i]]+tau0*alpha[[i]]*K[[i]])%*%t(N[[i]])%*%(diag(S1n)%*%y-diag(S1n)%*%X%*%beta0-S3n*Delta0-diag(S1n)%*%Reduce("+",prednp[-i]))
        prednp[[i]] = N[[i]]%*%theta0[[i]]
      }
    }else{
      theta0[[1]] = solve(t(N[[1]])%*%diag(S1n)%*%N[[1]]+tau0*alpha[[1]]*K[[1]])%*%t(N[[1]])%*%(diag(S1n)%*%y-diag(S1n)%*%X%*%beta0-S3n*Delta0)
    }

    s20 = tau0+Delta0^2*(1-b^2)
    delta0 = Delta0/sqrt(tau0+Delta0^2)
    gama0 = ((4-pi)/2)*b^3*(delta0^3)/((1-(2/pi)*delta0^2)^(3/2))
    for(i in 1:length(N)){
      prednp[[i]] = N[[i]]%*%theta0[[i]]
    }
    mu0 = X%*%beta0+Reduce("+",prednp)

    nu0 = optim(par = nu, fn = nllg, y = y, mu = mu0, sigma2 = s20,
                gama = gama0,
                method = "L-BFGS-B", lower = c(2,1), upper = c(5,5),
                control = list(fnscale = -1))$par
    nu = nu0
    ver_novo = nllg(nu,y,mu0,s20,gama0)

    if(alpha.fix == F){
      ft = function(alpha0) BICsmn(alpha0,ver_novo,theta0,K,N,tau0,p,
                                    diag(S1n),"CSBPN")
      alpha = optim(alpha, ft, method = "L-BFGS-B", lower = a0,
                     upper = a1)$par
    }

    betao[count,] = beta0
    sigma2o[count] = s20
    gamao[count] = gama0
    nuo[count,] = nu
    count = (count+1)
    beta = beta0
    theta = theta0
    sigma2 = s20
    sigma = sqrt(sigma2)
    gama = gama0
    mu  = mu0
    bu = S1n
    S1o = S1n
    S2o = S2n
    S3o = S3n
    progress(count,iter.max, progress.bar = T)
  }
  return(list(beta = beta, theta = theta, sigma2 = sigma2, gama = gama, mu = mu, betao = betao,
              sigma2o = sigma2o, gamao = gamao, nuo = nuo, nu = nu,
              alpha = alpha, delta = delta0, S1n = S1n, S2n = S2n, S3n = S3n))
}


SAEM.CSBSN = function(y,mu,theta,sigma2,gama,n,p,alpha0,K,N,X,k,iter.max,
                     alpha.fix, a0, a1){
  nllg  = function(nu,y,mu,sigma2,gama){
    nu1=nu[1]
    nu2=nu[2]
    return(sum(log(dCSBSN(y, mu, sigma2,gama,c(nu1,nu2)))))
  }
  nu = optim(par = c(5,5), fn = nllg, y = y, mu = mu, sigma2 = sigma2, gama = gama,
             method = "L-BFGS-B", lower = c(2,1), upper = c(10,10),
             control = list(fnscale = -1, maxit = 50))$par

  S1o = 0
  S2o = 0
  S3o = 0
  S4o = 0
  criterio = 1
  count    = 1
  bu = rfatigue(n,nu[1],nu[2])
  m = 20
  betao = matrix(0,iter.max,p)
  sigma2o = rep(0,iter.max)
  gamao = rep(0,iter.max)
  nuo = matrix(0,iter.max,2)
  prednp = list(NULL)
  for(i in 1:length(N)){
    prednp[[i]] = N[[i]]%*%theta[[i]]
  }
  sigma = sqrt(sigma2)
  theta0 = list(NULL)
  alpha = alpha0
  while((count <= iter.max)){
    g13 = nthroot(gama,3)
    g23 = g13^2
    lambda = (s*g13)/sqrt(2/pi+s^2*g23*(2/pi - 1))
    delta = lambda/sqrt(1+lambda^2)
    Delta = sigma*delta/sqrt(1-b^2*delta^2)
    tau = sigma2*(1-delta^2)/(1-b^2*delta^2)
    esq = matrix(0,n,m)
    e10 = matrix(0,n,m)
    e51 = matrix(0,n,m)
    e02 = matrix(0,n,m)
    for(i in 1:n){
      aux2 = ghudyCSBSN(y[i],bu[i],mu[i],sigma2,gama,nu,m)
      eh1 = aux2[,1]
      e02[i,] = (eh1-b)^2
      e10[i,] = aux2[,2]
      esq[i,] = sqrt(e10[i,])
      e51[i,] = esq[i,]*(eh1-b)
    }
    ss1n = 0
    ss2n = 0
    ss3n = 0
    for(i in 1:m){
      ss1n = ss1n + e10[,m]
      ss2n = ss2n + e02[,m]
      ss3n = ss3n + e51[,m]
    }
    S1n = S1o + deltak(count,k,iter.max)*((1/m)*ss1n-S1o)
    S2n = S2o + deltak(count,k,iter.max)*((1/m)*ss2n-S2o)
    S3n = S3o + deltak(count,k,iter.max)*((1/m)*ss3n-S3o)
    Delta0 = (1/sum(S2n))*sum((y-mu)*S3n)
    tau0 = (1/n)*sum((y-mu)^2*S1n-2*(y-mu)*Delta0*S3n+Delta0^2*S2n)
    beta0 = solve(t(X)%*%diag(S1n)%*%X)%*%t(X)%*%(diag(S1n)%*%y-diag(S1n)%*%Reduce("+",prednp)-S3n*Delta0)
    if(length(N)>1){
      for(i in 1:length(N)){
        theta0[[i]] = solve(t(N[[i]])%*%diag(S1n)%*%N[[i]]+tau0*alpha[[i]]*K[[i]])%*%t(N[[i]])%*%(diag(S1n)%*%y-diag(S1n)%*%X%*%beta0-S3n*Delta0-diag(S1n)%*%Reduce("+",prednp[-i]))
        prednp[[i]] = N[[i]]%*%theta0[[i]]
      }
    }else{
      theta0[[1]] = solve(t(N[[1]])%*%diag(S1n)%*%N[[1]]+tau0*alpha[[1]]*K[[1]])%*%t(N[[1]])%*%(diag(S1n)%*%y-diag(S1n)%*%X%*%beta0-S3n*Delta0)
    }

    s20 = tau0+Delta0^2*(1-b^2)
    delta0 = Delta0/sqrt(tau0+Delta0^2)
    gama0 = ((4-pi)/2)*b^3*(delta0^3)/((1-(2/pi)*delta0^2)^(3/2))
    for(i in 1:length(N)){
      prednp[[i]] = N[[i]]%*%theta0[[i]]
    }
    mu0 = X%*%beta0+Reduce("+",prednp)

    nu0 = optim(par = nu, fn = nllg, y = y, mu = mu0, sigma2 = s20,
                gama = gama0, method = "L-BFGS-B", lower = c(2,1),
                upper = c(10,10), control = list(fnscale = -1,maxit = 50))$par
    nu = nu0
    ver_novo = nllg(nu,y,mu0,s20,gama0)

    if(alpha.fix == F){
      ft = function(alpha0) BICsmn(alpha0,ver_novo,theta0,K,N,tau0,p,diag(S1n),"CSBSN")
      alpha = optim(alpha, ft, method = "L-BFGS-B", lower = a0,
                     upper = a1)$par
    }

    betao[count,] = beta0
    sigma2o[count] = s20
    gamao[count] = gama0
    nuo[count,] = nu
    count = (count+1)
    beta = beta0
    theta = theta0
    sigma2 = s20
    sigma = sqrt(sigma2)
    gama = gama0
    mu  = mu0
    bu = S1n
    S1o = S1n
    S2o = S2n
    S3o = S3n
    progress(count,iter.max, progress.bar = T)
  }
  return(list(beta = beta, theta = theta, sigma2 = sigma2, gama = gama, mu = mu, betao = betao,
              sigma2o = sigma2o, gamao = gamao, nuo = nuo, nu = nu,
              alpha = alpha, S1n = S1n, S2n = S2n, S3n = S3n))
}

SAEM.CSGGN = function(y,mu,theta,gama,n,p,alpha0,K,N,X,k,iter.max, alpha.fix,
                     a0, a1){
  nllg  = function(nu,y,mu,gama){
    nu1=nu[1]
    nu2=nu[2]
    nu3 = nu[3]
    return(sum(log(dCSGGN(y, mu,gama,c(nu1,nu2,nu3)))))
  }
  nu = optim(par = c(2,2,2), fn = nllg, y = y, mu = mu, gama = gama,
             method = "L-BFGS-B", lower = c(1,1,0.1), upper = c(5,5,5),
             control = list(fnscale = -1, maxit = 50))$par

  S1o = 0
  S2o = 0
  S3o = 0
  S4o = 0
  criterio = 1
  count    = 1
  bu = rggamma(n,nu[3],nu[1],nu[2])
  m = 20
  betao = matrix(0,iter.max,p)
  gamao = rep(0,iter.max)
  nuo = matrix(0,iter.max,3)
  sigma2 = sigma = 1
  prednp = list(NULL)
  for(i in 1:length(N)){
    prednp[[i]] = N[[i]]%*%theta[[i]]
  }
  theta0 = list(NULL)
  alpha = alpha0
  while((count <= iter.max)){
    g13 = nthroot(gama,3)
    g23 = g13^2
    lambda = (s*g13)/sqrt(2/pi+s^2*g23*(2/pi - 1))
    delta = lambda/sqrt(1+lambda^2)
    Delta = sigma*delta/sqrt(1-b^2*delta^2)
    tau = sigma2*(1-delta^2)/(1-b^2*delta^2)
    esq = matrix(0,n,m)
    e10 = matrix(0,n,m)
    e51 = matrix(0,n,m)
    e02 = matrix(0,n,m)
    for(i in 1:n){
      aux2 = ghudyCSGGN(y[i],bu[i],mu[i],gama,nu,m)
      eh1 = aux2[,1]
      e02[i,] = (eh1-b)^2
      e10[i,] = aux2[,2]
      esq[i,] = sqrt(e10[i,])
      e51[i,] = esq[i,]*(eh1-b)
    }
    ss1n = 0
    ss2n = 0
    ss3n = 0
    for(i in 1:m){
      ss1n = ss1n + e10[,m]
      ss2n = ss2n + e02[,m]
      ss3n = ss3n + e51[,m]
    }
    S1n = S1o + deltak(count,k,iter.max)*((1/m)*ss1n-S1o)
    S2n = S2o + deltak(count,k,iter.max)*((1/m)*ss2n-S2o)
    S3n = S3o + deltak(count,k,iter.max)*((1/m)*ss3n-S3o)
    Delta0 = (1/sum(S2n))*sum((y-mu)*S3n)
    tau0 = (1/n)*sum((y-mu)^2*S1n-2*(y-mu)*Delta0*S3n+Delta0^2*S2n)
    beta0 = solve(t(X)%*%diag(S1n)%*%X)%*%t(X)%*%(diag(S1n)%*%y-diag(S1n)%*%Reduce("+",prednp)-S3n*Delta0)
    if(length(N)>1){
      for(i in 1:length(N)){
        theta0[[i]] = solve(t(N[[i]])%*%diag(S1n)%*%N[[i]]+tau0*alpha[[i]]*K[[i]])%*%t(N[[i]])%*%(diag(S1n)%*%y-diag(S1n)%*%X%*%beta0-S3n*Delta0-diag(S1n)%*%Reduce("+",prednp[-i]))
        prednp[[i]] = N[[i]]%*%theta0[[i]]
      }
    }else{
      theta0[[1]] = solve(t(N[[1]])%*%diag(S1n)%*%N[[1]]+tau0*alpha[[1]]*K[[1]])%*%t(N[[1]])%*%(diag(S1n)%*%y-diag(S1n)%*%X%*%beta0-S3n*Delta0)
    }

    delta0 = Delta0/sqrt(tau0+Delta0^2)
    gama0 = ((4-pi)/2)*b^3*(delta0^3)/((1-(2/pi)*delta0^2)^(3/2))
    for(i in 1:length(N)){
      prednp[[i]] = N[[i]]%*%theta0[[i]]
    }
    mu0 = X%*%beta0+Reduce("+",prednp)

    nu0 = optim(par = nu, fn = nllg, y = y, mu = mu0, gama = gama0,
                method = "L-BFGS-B", lower = c(1,1,0.1), upper = c(5,5,5),
                control = list(fnscale = -1,maxit = 50))$par
    nu = nu0
    ver_novo = nllg(nu,y,mu0,gama0)

    if(alpha.fix == F){
      ft = function(alpha0) BICsmn(alpha0,ver_novo,theta0,K,N,tau0,p,diag(S1n),"CSGGN")
      alpha = optim(alpha, ft, method = "L-BFGS-B", lower = a0,
                     upper = a1)$par
    }

    betao[count,] = beta0
    gamao[count] = gama0
    nuo[count,] = nu
    count = (count+1)
    beta = beta0
    theta = theta0
    gama = gama0
    mu  = mu0
    bu = S1n
    S1o = S1n
    S2o = S2n
    S3o = S3n
    progress(count,iter.max, progress.bar = T)
  }
  return(list(beta = beta, theta = theta, sigma2 = sigma2, gama = gama, mu = mu, betao = betao,
              gamao = gamao, nuo = nuo, nu = nu, alpha = alpha, S1n = S1n, S2n = S2n, S3n = S3n))
}

dbetapr = function(x,nu1,nu2){
  res = x^(nu1 - 1)/((1+x)^(nu1+nu2)*beta(nu1,nu2))
  return(res)
}
rbetapr = function(n, nu1, nu2){
  x = rbeta(n, nu1, nu2)
  return(x/(1-x))
}
dfatigue = function(x,alpha,beta){
  res = ((sqrt((x)/beta) + sqrt(beta/(x)))/(2*alpha*(x))) * dnorm((sqrt((x)/beta) - sqrt(beta/(x)))/alpha)
  return(res)
}
rfatigue = function(n, alpha, beta){
  z = rnorm(n)
  res = (alpha/2*z + sqrt((alpha/2*z)^2 + 1))^2 * beta
  return(res)
}

