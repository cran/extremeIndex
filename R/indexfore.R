#' Function for computing the index for a forecast system vs. climatological forecast. You must provide an indexclim object.
#'
#' @param clim an indexclim object coming from indexclim
#' @param score_fore the time serie of the ensemble forecast's CRPS/MAE. Be careful that score_fore is consistent with "score" in indexclim
#'
#' @return an indexfore object with the index computed vs. climatological forecast
#' @export
#'
#' @import goftest
#' @importFrom evir gpd

indexfore=function(score_fore,clim){
  stopifnot(class(clim)=="indexclim")
  


  testdata=function(yy){

    exc=score_fore[which(score_fore>yy)]
 
    pa=c(clim$xi,evir::gpd(clim$obs,threshold=yy,method="pwm")$par.ests[2])


    cv=goftest::cvm.test(exc,function (q,xi,mu,beta)
    { if (xi>0.01){(1 - (1 + (xi * (q - mu))/beta)^(-1/xi))}
      else{1-exp(-(q-mu)/beta)}
    },xi=pa[1],mu=0,beta=pa[2])$p.value
    return(cv)

  }
  cvm=sapply(clim$quantiles,function(q) testdata(q))
  result=list(quantiles=clim$quantiles,index=1-cvm/clim$index,obs=clim$obs,clim=clim,score=score_fore)
  class(result)="indexfore"
  return(result)
}
