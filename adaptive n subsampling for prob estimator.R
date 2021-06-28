################################################
####### MAIN CODES: ############################
rho.fast.sub = matrix(0,nr=d,nc=d)
Te = rexp(d)
U <- sapply(1:d, function(i)norm.complex(rcmvnorm(1,sigma=diag(d))))
Lamb = c(Te/sum(Te))
gamm = n.size/2
#MCMC loop
Iter = 1000
burnin = 200
be.y = 0.01
beta.z = .03
ro = 1/2
accc = 0
U.can = U
Te.can = Te
sub_sample = 0.6
start_time <- Sys.time()
for(t in 1:(Iter+burnin)){
  Te.can = Te*exp(be.y*runif(d,min=-0.5,.5) )
  U.can = sqrt(1-beta.z^2)*U + beta.z*rcmvnorm(d,sigma=diag(d))
  L.can = Te.can/sum(Te.can)
  tam1 = apply(U.can,2,norm.complex)
  tem.can = c(tam1%*% ( Conj(t(tam1 ))* L.can ) )
  tam2 = apply(U,2,norm.complex)
  tem =  c(tam2%*% ( Conj(t(tam2 ))* Lamb ) )
  #subsampling
  rd = sample(A*R,A*R*sub_sample)
  Pra.rd = Pra[rd]
  p_ra1.rd = p_ra1[rd]
  ss = sum((sapply(Pra.rd, function(x)crossprod(tem.can,x)) -p_ra1.rd)^2-
             (sapply(Pra.rd, function(x)crossprod(tem,x)) -p_ra1.rd)^2)
  r.prior = sum( (ro)*log(Te.can/Te) -Te.can +Te,na.rm = T )
  ap = -gamm*Re(ss)
  if(log(runif(1)) <= ap + r.prior){
    Te<- Te.can
    U = U.can
    accc = accc +1
  }
  Lamb = c(Te/sum(Te))
  
  #approx rho
  if(t>burnin){
    U.n = apply(U,2,norm.complex)
    rho.fast.sub = U.n %*% (Conj(t(U.n ))*Lamb)/(t-burnin) + rho.fast.sub*(1-1/(t-burnin))
  }
}
end_time <- Sys.time()
end_time - start_time
accc/(Iter+burnin)
Re(mean(tcprod.Cplex(dens.ma-rho.old)))
Re(mean(tcprod.Cplex(dens.ma-rho.fast)))
Re(mean(tcprod.Cplex(dens.ma-rho.fast.sub)))
Re(mean(tcprod.Cplex(dens.ma-rho.hat)))


