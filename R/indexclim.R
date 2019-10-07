#' Function which computes the index for the climatological CRPS/MAE. You must provide the observations. If you computes climatological CRPS/MAE previously, you can add the corresponding vector
#'
#' @param y The observations
#' @param thresh Vector of thresholds where you want to compute the index
#' @param score_clim If not NULL, must be the time serie of the CRPS/MAE of the climatology. It is recommended to compute CRPS/MAE out of this function
#' @param xi Shape parameter of the GP
#' @param score A character string indicating if you want to work with CRPS ("crps") or MAE ("mae"), by default "crps"

#' @return An indexclim object containing xi, y, the score time serie, the score considered, the index values, and the corresponding quantiles of the observations
#' @export
#'
#' @import goftest
#' @importFrom evir gpd
indexclim=function(y,thresh=NULL,score_clim=NULL,xi=NULL,score="crps"){

  stopifnot(is.numeric(y))
  stopifnot(!is.null(thresh))
stopifnot(!is.null(xi))
if ((score!="crps")&(score!="mae")){stop("score must be crps or mae")}

  if (is.null(score_clim)){

	if (score=="crps"){
    alpha=mean(y)-2*mean(seq(0,1,length.out=length(y))*sort(y))
    score_clim=sapply(y,function(i) mean(abs(y-i))+alpha)
    }
	if (score=="mae"){
    score_clim=sapply(y,function(i) mean(abs(y-i)))
    }

  }



  testdata=function(yy){

    exc=score_clim[which(score_clim>yy)]
    pa=c(xi,evir::gpd(y,threshold=yy,method="pwm")$par.ests[2])

    cv=goftest::cvm.test(exc,function (q,xi,mu,beta)
    { if (xi>0.01){(1 - (1 + (xi * (q - mu))/beta)^(-1/xi))}
      else{1-exp(-(q-mu)/beta)}
    },xi=pa[1],mu=0,beta=pa[2])$p.value
    return(cv)

  }
  cvm=sapply(thresh,function(q) testdata(q))


  result=list(quantiles=thresh,crps=score_clim,index=cvm,obs=y,xi=xi,score=score)
  class(result)="indexclim"
  return(result)
}
