#' Function which plots the index for differents forecasts sharing the same observations
#'
#' @param forecasts list of "indexfore" objects, all forecasts must be computed on the same climatology and thresholds
#' @param col colors of the differents forecasts for the plot
#' @param leg legend of the plot
#' @param xtypq the x-axis of the plot is quantiles values or orders (TRUE for quantiles)
#' @param ... other arguments for the plot
#'
#' @return a plot of the indices and a matrix containing the indexes for each threshold/order
#' @export
#'
#' @importFrom graphics abline arrows legend lines par plot points polygon title
#' @importFrom grDevices rainbow
#' @examples
#' data("crps")
#' y=crps[1:500,1]
#' cli=indexclim(y,thresh=seq(3,quantile(y,probs=0.995),length=2),xi=0.2)
#' frcst=crps[1:500,2]
#' idf=indexfore(frcst,cli)
#' frcst=crps[1:500,3]
#' idf2=indexfore(frcst,cli)
#'fore=list(idf,idf2)
#'idxp2=index.plot(fore,col=c("red","blue"),leg=c("forecast 1",
#'"forecast 2"),main="Index plot")
index.plot=function(forecasts,col=NULL,leg=NULL,xtypq=TRUE,...){
  stopifnot(all(sapply(forecasts,class)=="indexfore"))
  if(is.null(col)){col=rainbow(length(forecasts))}
  else{stopifnot(length(col)==length(forecasts))}
  inds=matrix(NA,ncol=length(forecasts),nrow=length(forecasts[[1]]$quantiles))

if (xtypq==FALSE){
forecasts[[1]]$quantiles=ecdf(forecasts[[1]]$clim$obs)(forecasts[[1]]$quantiles)
}

plot(forecasts[[1]]$quantiles,forecasts[[1]]$index,ylim=c(0,1),col=col[1],type="l",xlab="Threshold",ylab="Index",...)
abline(h = 0, col="black",lty=2)
inds[,1]=forecasts[[1]]$index
if(length(forecasts)>1){
for( i in 2:length(forecasts)){
 lines(forecasts[[1]]$quantiles,forecasts[[i]]$index,type="l",col=col[i],...)
  inds[,i]=forecasts[[i]]$index

}
}
if(length(leg)!=length(col)){leg=NULL}
if(is.null(leg)){leg=as.character(1:length(forecasts))}
legend("topleft",legend=leg,col=col,bty="n",lwd=par("lwd"))
return(invisible(list(thresh=forecasts[[1]]$quantiles,inds=inds)))
}
