geecsbsn = function(y,mu,sigma2,gama,n,p,alpha,K,Na,x,iter.max, alpha.fix,
                    t, id, corstr = "independence",
                    linkmu = "identity",idvec){
  nllg  <- function(nut,y,mu,sigma2,gama){
    nut1 = nut[1]
    nut2 = nut[2]
    return(-sum(log(dCSBSN(y, mu, sigma2, gama, c(1/nut1,1/nut2)))))
  }
  nllg3 <- function(parms,y,mu){
    nut1 = parms[1]
    nut2 = parms[2]
    sigma2 = parms[3]
    gama = parms[4]
    return(-sum(log(dCSBSN(y, mu, sigma2, gama, c(1/nut1,1/nut2)))))
  } #####

  nut = optim(par = c(0.15,0.15), fn = nllg, y = y, mu = mu, sigma2 = sigma2,
              gama = gama, method = "L-BFGS-B", lower = c(0.05,0.05),
              upper = c(1,1))$par #####
  nu = 1/nut
  if(linkmu == "identity"){
    yvn = y
  }
  if(linkmu == "log"){
    yvn = log(y)
  }
  beta = solve(t(x)%*%x)%*%t(x)%*%yvn
  theta = solve(t(Na)%*%Na+alpha*K)%*%t(Na)%*%(yvn-x%*%beta)
  sigma = sqrt(sigma2)
  N = sum(t)
  p = ncol(x)
  s = nthroot(2/(4-pi),3)
  b = sqrt(2/pi)
  criterio <- 1
  count    <- 1
  alphac = 0
  df = 1
  cat("Starting algorithm")
  while((count <= iter.max)){
    print(count)
    print(nu)
    eta = x%*%beta + Na%&%theta
    if(linkmu == "identity"){
      mu = as.vector(eta) # mi para a ligação logarítmica
    }
    if(linkmu == "log"){
      mu = exp(eta)
    }
    g13 = nthroot(gama,3)
    g23 = g13^2
    lambda = (s*g13)/sqrt(2/pi+s^2*g23*(2/pi - 1))
    delta = lambda/sqrt(1+lambda^2)
    Delta = sigma*delta/sqrt(1-b^2*delta^2)
    tau = as.numeric(sigma2*(1-delta^2)/(1-b^2*delta^2))
    u = y-mu
    M = cbind(x,Na)
    if(linkmu == "identity"){
      Lambda = diag(rep(1,N),N)
    }
    if(linkmu == "log"){
      Lambda = diag(as.vector(mu),N)
    }
    vary = sigma2*(nu[1]^2+2)/(2*nu[2]) #####
    varu = vary
    uc2 = split(u,id)
    uc = list(NULL)
    for(i in 1:n){
      idaux = idvec[id==i]
      id2 = rep(0,max(t))
      id2[idaux] = 1
      cid = 0
      uc[[i]] = rep(0,max(t))
      for(j in 1:max(t)){
        if(id2[j]==1){
          cid = cid +1
          uc[[i]][j] = uc2[[i]][cid]
        }else{
          uc[[i]][j] = 0
        }
      }
    }
    if(corstr == "unstructured"){
      Rg = matrix(0,max(t),max(t))
      cnum = den1 = den2 = 0
      for(j in 1:(max(t))){
        for(k in j:(max(t))){
          for(i in 1:n){
            if(is.na(uc[[i]][j])||is.na(uc[[i]][k])){
              cnum = cnum
            }
            else{
              cnum = cnum + (uc[[i]][j])*(uc[[i]][k])
              den1 = den1 + (uc[[i]][j])^2
              den2 = den2 + (uc[[i]][k])^2
            }
          }
          Rg[j,k] = cnum/(sqrt(den1)*sqrt(den2))
          Rg[k,j] = Rg[j,k]
        }
      }
      diag(Rg) = 1
      R = list(NULL)
      for(i in 1:n){
        R[[i]] = Rg[idvec[id==i],idvec[id==i]]
      }
      Rm = bdiag(R)
    }
    if(corstr == "AR-1"){
      cnum = den1 = den2 = 0
      for(i in 1:n){
        if(t[i]==1){
          cnum = cnum
          den1 = den1
          den2 = den2
        }else{
          for(j in 1:(max(t)-1)){
            cnum = cnum + uc[[i]][j]*uc[[i]][j+1]
            den1 = den1 + uc[[i]][j]^2
            den2 = den2 + uc[[i]][j+1]^2
          }
        }
      }
      alphac = cnum/sqrt(den1*den2)
      Rm = matrix(0,N,N)
      diag(Rm) = 1
      R = list(NULL)
      for(i in 1:n){
        R[[i]] = matrix(0,t[i],t[i])
        for(j in 1:t[i]){
          for(l in 1:t[i]){
            R[[i]][j,l] = alphac^(abs(idvec[id==i][j]-idvec[id==i][l]))
          }
        }
      }
      # Matriz de correlação AR-1
      Rm = as.matrix(bdiag(R))
      R=R[[1]]
    }
    if(corstr == "exchangeable"){
      cnum = den =  0
      for(i in 1:n){
        for(j in 1:max(t)){
          for(k in 1:max(t)){
            if(is.na(uc[[i]][j])||is.na(uc[[i]][k])){
              cnum = cnum
              den = den
            }else{
              if(j>k){
                cnum = cnum + uc[[i]][j]*uc[[i]][k]
              }else{
                cnum = cnum
              }
            }
          }
          den = den + uc[[i]][j]^2
        }
      }
      alphac = (cnum/sum(.5*t*(t-1)))/(den/sum(t))
      Rm = matrix(0,N,N)
      diag(Rm) = 1
      R = list(NULL)
      for(i in 1:n){
        R[[i]] = matrix(0,t[i],t[i])
        for(j in 1:t[i]){
          for(l in 1:t[i]){
            R[[i]][j,l] = alphac
          }
        }
        diag(R[[i]]) = 1
      }
      # Matriz de correlação AR-1
      Rm = as.matrix(bdiag(R))
      R=R[[1]]
      diag(Rm) = 1
    }
    Sigma = diag(rep(sqrt(varu),N),N)%*%Rm%*%diag(rep(sqrt(varu),N),N)
    W = t(Lambda)%*%solve(Sigma)%*%Lambda
    z = eta + solve(Lambda)%*%u
    beta0 = as.vector(solve(t(x)%*%W%*%x)%*%t(x)%*%W%*%(z-Na%*%theta))
    theta0 = as.vector(solve(t(Na)%*%W%*%Na+alpha*K)%*%t(Na)%*%W%*%(z-x%*%beta0))
    df = sum(diag(Na%*%solve(t(Na)%*%W%*%Na + alpha*K)%*%t(Na)%*%W))
    if(linkmu == "identity"){
      mu0 = as.vector(x%*%beta0+Na%*%theta0)
    }
    if(linkmu == "log"){
      mu0 = exp(as.vector(x%*%beta0+Na%*%theta0))
    }
    parmse = nlminb(c(nut[1],nut[2],sigma2,gama),nllg3, y = y, mu = mu0,
                    lower = c(0.03,0.03,0.1,-0.9),
                    upper = c(1,1,10,0.9))$par #####
    nut  = c(parmse[1],parmse[2])
    s20 = parmse[3]
    gama0 = parmse[4]
    nu = 1/nut

    count <- (count+1)
    dife = sum(abs(c(beta,theta,sigma2,gama)-c(beta0,theta0,s20,gama0)))
    print(dife)
    beta = beta0
    theta = theta0
    sigma2 = s20
    sigma = sqrt(sigma2)
    gama = gama0
    mu  <- mu0
    if(dife<0.0001){
      break
    }
  }
  auxt = c(0,cumsum(t))
  ploglik = 0
  for(i in 1:n){
    ploglik = ploglik -(t[i]/2)*log(2*pi)-.5*determinant(as.matrix(Sigma[(auxt[i]+1):auxt[i+1],(auxt[i]+1):auxt[i+1]]))$modulus-t(u[(auxt[i]+1):auxt[i+1]])%*%solve(Sigma[(auxt[i]+1):auxt[i+1],(auxt[i]+1):auxt[i+1]])%*%u[(auxt[i]+1):auxt[i+1]]
  }
  pAIC = 2*(p+df+3)-2*ploglik
  pBIC = log(N)*(p+df+3)-2*ploglik
  print(pAIC)
  print(pBIC)
  g13 = nthroot(gama,3)
  g23 = g13^2
  lambda = (s*g13)/sqrt(2/pi+s^2*g23*(2/pi - 1))
  delta = lambda/sqrt(1+lambda^2)
  wm = sqrtm(as.matrix(W))$B
  eta = x%*%beta+Na%*%theta
  Ha = Na%*%solve(t(Na)%*%W%*%Na + alpha*K)%*%t(Na)%*%W
  resido = wm%*%(z-eta)
  g13 = nthroot(gama,3)
  g23 = g13^2
  lambda = (s*g13)/sqrt(2/pi+s^2*g23*(2/pi - 1))
  delta = lambda/sqrt(1+lambda^2)
  Delta = sigma*delta/sqrt(1-b^2*delta^2)
  tau = as.numeric(sigma2*(1-delta^2)/(1-b^2*delta^2))
  vary = sigma2*(nu[1]^2+2)/(2*nu[2])
  varu = vary
  Omega = as.matrix(diag(rep(sqrt(vary),N),N)%*%Rm%*%diag(rep(sqrt(vary),N),N))
  varep = (diag(1,N)-Ha)%*%wm%*%solve(Lambda)%*%Omega%*%solve(Lambda)%*%wm%*%(diag(1,N)-Ha)
  resid = resido/sqrt(diag(varep))
  resnai = resido/sqrt(1-diag(Ha))
  Malpha = bdiag(matrix(0,p,p),alpha*K)
  Seps = -t(M)%*%W%*%M-as.matrix(Malpha)
  u = (y-mu)
  vareps = solve(Seps)%*%t(M)%*%W%*%solve(Lambda)%*%u%*%t(u)%*%solve(Lambda)%*%t(W)%*%M%*%solve(Seps)
  varnaive = solve(Seps)%*%t(M)%*%W%*%solve(Lambda)%*%Omega%*%solve(Lambda)%*%t(W)%*%M%*%solve(Seps)
  return(list(beta = beta, theta = theta, sigma2 = sigma2, gama = gama, mu = mu,
              nu = nu, alpha = alpha, delta = delta, alphac = alphac,pAIC,
              pBIC,Rm = Rm, nobs = N, nclusters = n, Omega = Sigma,
              clusters = t,
              resnai = resnai, sdeps = sqrt(diag(vareps)),
              vareps = vareps, resrob = resid, df = df, vary = vary,
              varnaive = varnaive))
}
