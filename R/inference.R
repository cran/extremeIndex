#'@export
#'@import gmm
#'@importFrom stats ecdf pbeta runif
#########################
###                   ###
### Inference via PMW ###
###                   ###
#########################

.PWM <- function(orders,kappa=NA,sigma=NA,xi=NA,empiric=FALSE,Unif=NULL,NbSamples=10^4,N=200){
  H.L <- 0
  H.U <- 1

  F.L <- 0
  F.U <- 1
  prob.LU <- F.U-F.L

  if(!empiric){
    E2 <- ((1-F.L)^(orders+1)-(1-F.U)^(orders+1))/(prob.LU*(orders+1))


    beta.coefs <- beta(c(1:(max(orders)+1))*kappa,1-xi)
    probs.beta <- pbeta(H.U,shape1=c(1:(max(orders)+1))*kappa,shape2=1-xi)-pbeta(H.L,shape1=c(1:(max(orders)+1))*kappa,shape2=1-xi)
    E1 <- c()
    for(i in 1:length(orders)){
      E1[i] <- (kappa/prob.LU)*sum((-1)^c(0:orders[i])*choose(orders[i],c(0:orders[i]))*beta.coefs[1:(orders[i]+1)]*probs.beta[1:(orders[i]+1)])
    }
    return( (sigma/xi)*(E1-E2) )

  } else{
    if(is.null(Unif)){
      Unif <- runif(NbSamples)
    }

    V <- .rG(length(Unif),kappa,Unif)

    X <- .qgpd2(V,scale=sigma,shape=xi)
    res <- c()
    for(i in 1:length(orders)){
      res[i] <- mean(X*(1-(F.L + (F.U-F.L)*Unif))^orders[i],na.rm=TRUE)
    }
    return( res )
  }
}

.EGP.fitPWM <- function(x,kappa0=NA,sigma0=NA,xi0=NA,empiric=FALSE,Unif=NULL,NbSamples=10^4){
  Fn = ecdf(x)

  mu0hat = mean(x)
  mu1hat = mean(x*(1-Fn(x)))
  mu2hat = mean(x*(1-Fn(x))^2)
  mu3hat = mean(x*(1-Fn(x))^3)
  mu4hat = mean(x*(1-Fn(x))^4)

  fct <- function(theta,x){
    pwm.theor <- .PWM(orders=c(0:2),kappa=theta[1],sigma=theta[2],xi=theta[3],empiric=empiric,Unif=Unif,NbSamples=NbSamples)
    pwm.empir <- c(mu0hat,mu1hat,mu2hat)
    return(matrix(pwm.theor - pwm.empir,ncol=3))
  }
  theta0 <- c(kappa0,sigma0,xi0)
  res <- gmm::gmm(fct, x, theta0, optfct = "nlminb", lower = c(0.0001, 0.0001, 10^(-6)),  upper = c(Inf,Inf, .99),vcov="iid")
  thetahat <- res$coefficients
  names(thetahat) <- c("kappa","sigma","xi")
  return(thetahat)

}


.EGP.fitPWM.boot <- function(data,i,kappa0=NA,sigma0=NA,xi0=NA,empiric=FALSE,Unif=NULL,NbSamples=10^4){
  return( .EGP.fitPWM(data[i],kappa0=kappa0,sigma0=sigma0,xi0=xi0,empiric=empiric,Unif=Unif,NbSamples=NbSamples) )
}




