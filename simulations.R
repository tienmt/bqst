

rwMH = rwMH.eig =c()
for(rr in 1:50){
  rho.old = matrix(0,nr=d,nc=d)
  Te = rexp(d)
  Id = diag(d)
  U <- sapply(1:d, function(i)norm.complex(rcmvnorm(1,sigma=diag(d))))
  Lamb = c(Te/sum(Te))
  ro = 1
  be = 1
  gamm = n.size/2
  entry = c()
  #MCMC loop
  Iter = 2000
  burnin = 100
  
  for(t in 1:(Iter+burnin)){
    for(j in 1:d){
      #update Lambda
      Te.can = Te
      Te.can[j] = Te[j]*exp(be*runif(1,min=-0.5,.5))
      L.can = Te.can/sum(Te.can)
      tem.can = c(U%*%diag(L.can)%*%Conj(t(U)))
      tem = c(U%*%diag(Lamb)%*%Conj(t(U)))
      ss = sum((sapply(Pra, function(x)crossprod(tem.can,x)) -p_ra1)^2-
                 (sapply(Pra, function(x)crossprod(tem,x)) -p_ra1)^2)
      r.prior = (ro-1)*log(Te.can[j]/Te[j]) -Te.can[j] +Te[j]
      ap = -gamm*Re(ss)
      if(log(runif(1)) <= ap + r.prior){
        Te<- Te.can
      }
      Lamb = c(Te/sum(Te))
    }
    #update U
    for(j in 1:d){
      U.can = U
      U.can[,j] = norm.complex(U[,j]+rcmvnorm(1,sigma=diag(d))/100)
      tem.can =c(U.can%*%diag(Lamb)%*%Conj(t(U.can)))
      tem = c(U%*%diag(Lamb)%*%Conj(t(U)))
      ss = sum((sapply(Pra, function(x)crossprod(tem.can,x)) -p_ra1)^2-
                 (sapply(Pra, function(x)crossprod(tem,x)) -p_ra1)^2)
      ap = Re(-gamm*ss)
      if(log(runif(1)) <= ap){U<- U.can}
    }
    #approx rho
    if(t>burnin){
      rho.old = U%*%diag(Lamb)%*%Conj(t(U))/(t-burnin) + rho.old*(1-1/(t-burnin))
      #mse.r[t] = Re(mean((dens.ma-rho.old)%*%Conj(t((dens.ma-rho.old))) ) )
    }
  }
  rwMH[rr] = Re(mean((dens.ma-rho.old)%*%Conj(t((dens.ma-rho.old))) ) )
  rwMH.eig[rr] = mean( abs(eigen(rho.old,only.values = T)$values -eigen(dens.ma,only.values = T)$values) )
}




adMH = adMH.eig = c()
for(rr in 1:50){
  rho.fast = matrix(0,nr=d,nc=d)
  Te = rexp(d)
  U <-  sapply(1:d, function(i)norm.complex(rcmvnorm(1,sigma=diag(d))))
  Lamb = c(Te/sum(Te))
  gamm = n.size/2
  #MCMC loop
  Iter = 1000
  burnin = 100
  be.y = 0.05
  beta.z = .04
  ro = 1/4
  ac = 0
  U.can = U
  Te.can = Te
  for(t in 1:(Iter+burnin)){
    Te.can = Te*exp(be.y*runif(d,min=-0.5,.5) )
    U.can = sqrt(1-beta.z^2)*U + beta.z*rcmvnorm(d,sigma=diag(d))
    L.can = Te.can/sum(Te.can)
    tam1 = apply(U.can,2,norm.complex)
    tem.can = c(tam1%*% ( Conj(t(tam1 ))* L.can ) )
    tam2 = apply(U,2,norm.complex)
    tem =  c(tam2%*% ( Conj(t(tam2 ))* Lamb ) )
    ss = sum((sapply(Pra, function(x)crossprod(tem.can,x)) -p_ra1)^2-
               (sapply(Pra, function(x)crossprod(tem,x)) -p_ra1)^2)
    r.prior = sum( (ro)*log(Te.can/Te) -Te.can +Te,na.rm = T )
    ap = -gamm*Re(ss)
    if(log(runif(1)) <= ap + r.prior){
      Te<- Te.can
      U = U.can
      ac = ac +1
    }
    Lamb = c(Te/sum(Te))
    
    #approx rho
    if(t>burnin){
      U.n = apply(U,2,norm.complex)
      rho.fast = U.n %*% (Conj(t(U.n ))*Lamb)/(t-burnin) + rho.fast*(1-1/(t-burnin))
      #mse.a[t] = Re(mean(tcprod.Cplex(dens.ma-rho.fast)))
    }
  }
  adMH[rr]= Re(mean((dens.ma-rho.fast)%*%Conj(t((dens.ma-rho.fast)))) )
  adMH.eig[rr] = mean( abs(eigen(rho.fast,only.values = T)$values -eigen(dens.ma,only.values = T)$values) )
print(rr)
  }


