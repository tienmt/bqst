##########################
# code for distance btw probability matrices
##########################
library(gtools)
library(cmvnorm)
library(rTensor)


######## Auxiliary functions to speed up ##############
# normalize a complex vector #
norm.complex = function(x)return(x/sqrt(sum(Mod(x)^2)))
norm.complex <- compiler::cmpfun(norm.complex)
# taking the diagonal faster #
Di.ag = function(A,B)return(sapply(1:nrow(A),function(s) A[s,]%*%B[,s]))
Di.ag <- compiler::cmpfun(Di.ag)
# return the eigen vectors
eigen.vec = function(A)return(eigen(A)$vectors)
eigen.vec <- compiler::cmpfun(eigen.vec)
# faster conjugate tcrossprod
tcprod.Cplex <- function(x)return(tcrossprod(x,Conj(x)))
tcprod.Cplex <- compiler::cmpfun(tcprod.Cplex )


### Pauli basis 2x2  ##
sx = matrix(c(0,1,1,0),nr=2)
sy = matrix(c(0,1i,-1i,0),nr=2)
sz = matrix(c(1,0,0,-1),nr=2)
basis = list(diag(2),sx,sy,sz)

### number of qubits  ###
n = 2
J = 4^n
I = 6^n
d = R = 2^n
A = 3^n
## total number of a, b and r
b = permutations(4,n,  repeats.allowed=T)
a = permutations(3,n, v=c(2,3,4), repeats.allowed=T)
r = permutations(2,n, v = c(-1,1),  repeats.allowed=T)

## return the projectors
projectors = function(a,r){
  tem1 = lapply(basis[a],eigen.vec)
  tem3 = lapply(1:length(tem1), function(s) return(tem1[[s]][r[s],]))
  lapply(tem3, function(x) tcprod.Cplex(x) )
}

### Pauli basis for n qubit
sig_b = list()
for (i in 1:J){
  sig_b[[i]] = kronecker_list(basis[b[i,]])
}


# The projectors matrices Pra
Pra = list()
count = 0
for(j in 1:A){   
  for(i in 1:R){
    count = count + 1
    Pra[[count]] = c(kronecker_list(projectors(a =a[j,],r = r[i,])))
  }
}


### simulate The "true" density matrix  ###
## pure state
#dens.ma = matrix(0,nr=d,nc=d)
#dens.ma[1,1] = 1
## mixed state
#u = sapply(1:d, function(i)norm.complex(rcmvnorm(1,sigma=diag(d)/100)))
#dens.ma = Conj(t(u))%*%u/d
## rank-2 state
v1 = t(rep(0,d));v1[1:(d/2)]=1;v1 = norm.complex(v1)
v2 = t(rep(0,d));v2[d:(d/2+1)] = 1;v2 = norm.complex(v2)
dens.ma = Conj(t(v1))%*%v1*0.5 + Conj(t(v2))%*%v2*0.5 #+(1-2*0.49)*diag(d)/d


### calculate the probabilities matrix P.ar, used to simulate data
Prob.ar = matrix(0,nr=A,nc=R)
if(n==1){
  for(i in 1:A){
    for(j in 1:R){
      Prob.ar[i,j] <- c(dens.ma)%*%unlist(projectors(a = a[i,],r = r[j,]))
    }  }
}else{
  for(i in 1:A){
    for(j in 1:R){
      Prob.ar[i,j] <- sum( dens.ma * kronecker_list(projectors(a = a[i,],r = r[j,])) )
    }}}
Prob.ar = Re(Prob.ar)

### calculate the frequency matrix after simulating sample
n.size = 1000 ## numbers of repeat the measurements
p_ra = apply(Prob.ar, 1,function(x){
  H = sample(1:R,n.size,prob = x,replace=TRUE)
  return(sapply(1:R,function(s)sum(H==s)/n.size))
})
p_ra1 = c(p_ra)  
Id = diag(d)


################################################
####### MAIN CODES: ############################
rho.fast = matrix(0,nr=d,nc=d)
Te = rexp(d)
U <- sapply(1:d, function(i)norm.complex(rcmvnorm(1,sigma=diag(d))))
Lamb = c(Te/sum(Te))
gamm = n.size/2
#MCMC loop
Iter = 1000
burnin = 100
be.y = 0.3
beta.z = .03
ro = 1
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
  }
}

Re(mean((dens.ma-rho.old)%*%Conj(t((dens.ma-rho.old))) )  )
Re(mean(tcprod.Cplex(dens.ma-rho.fast)))
Re(mean((rho.hat-dens.ma)%*%Conj(t((rho.hat-dens.ma)))) )

mean( abs(eigen(rho.fast,only.values = T)$values -eigen(dens.ma,only.values = T)$values) )
ac/(Iter+burnin)

# n =2, be.y = .33 ; beta.z = .2 ; ro = 1
# n = 3 ; be.y = .1 ; beta.z = .09
# n = 4 ; be.y = .1 ; beta.z = .04


