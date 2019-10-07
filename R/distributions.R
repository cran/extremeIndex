#'@export
#'@importFrom evd pgpd qgpd
#'@importFrom stats ecdf pbeta runif
################################################################
###                                                          ###
### Distribution/Density/Quantile/Random generator functions ###
###                                                          ###
################################################################

### GPD
.qgpd2 <- function (p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE)
{
  inds <- p <=0 | p>=1
  res <- rep(NA,length(p))
  if (min(scale) < 0)
    stop("invalid scale")
  if (length(shape) != 1)
    stop("invalid shape")
  if (lower.tail)
    p <- 1 - p
  if (shape == 0)
    res[!inds] <- loc - scale * log(p[!inds])
  else res[!inds] <- loc + scale * (p[!inds]^(-shape) - 1)/shape

  return(res)
}
#'@export
.qgpd.fullrange <- function (p, loc = 0, scale = 1, shape = 0, prob.loc=0.95)
{
  res <- loc + .qgpd2((p-prob.loc)/(1-prob.loc),scale=scale,shape=shape)
  return(res)
}


#'@export
.pG <- function(u,kappa=NA){

  return(u^kappa)

}
#'@export
.dG <- function(u,kappa=NA,log=FALSE){
  if(log==FALSE){

    return(kappa*u^(kappa-1))

  } else{

    return(log(kappa) + (kappa-1)*log(u))

  }
}
#'@export
.qG <- function(u,kappa=NA){

  return(u^(1/kappa))

}
#'@export
.rG <- function(n,kappa=NA,Unif=NULL){
  if(is.null(Unif)){
    Unif <- runif(n)
  }


  return( .qG(Unif,kappa) )


}


#'@export
.pEGP <- function(x,kappa=NA,sigma=NA,xi=NA){
  return( .pG(evd::pgpd(x,scale=sigma,shape=xi),kappa) )
}
#'@export
.dEGP <- function(x,kappa=NA,sigma=NA,xi=NA,log=FALSE){
  if(log==FALSE){
    return( .dG(evd::pgpd(x,scale=sigma,shape=xi),kappa)*evd::dgpd(x,scale=sigma,shape=xi) )
  } else{
    return( .dG(evd::pgpd(x,scale=sigma,shape=xi),kappa,log=TRUE) + evd::dgpd(x,scale=sigma,shape=xi,log=TRUE) )
  }
}
#'@export
.qEGP <- function(p,kappa=NA,sigma=NA,xi=NA){
  return( .qgpd2(.qG(p,kappa),scale=sigma,shape=xi) )
}
#'@export
.rEGP <- function(n,kappa=NA,sigma=NA,xi=NA,Unif=NULL){

  return( .qgpd2(.rG(n,kappa,Unif),scale=sigma,shape=xi) )

}
