#' Function for heuristically choosing the domain where extreme value theory can be applied
#'
#' @param data a numeric vector containing the observation used for verification
#' @param thresh vector of thresholds to try
#' @param guess starting values for GPD's sigma and xi (0<xi<1)
#' @param plots which parameter plots do you want
#' @param R number of bootstrap estimates for confidence intervals
#' @param ncpus if you want to make bootstrap on several cores
#' @return three plots summarizing the stability of the parameters to threshold. The starting threshold admits kappa=1 and its confidence interval ;
#' according Papastathopoulos & Tawn (2013)
#' @return a list with thresholds used, GP parameters and CIs, p-values of Cramer von Mises test (accordance of thresholded data with GP), optimal threshold and xi.
#' @export
#'
#' @import goftest
#' @importFrom graphics abline arrows legend lines par plot points polygon title
choosethres=function(data,thresh,guess=c(1,0.1),plots=1:3,R=200,ncpus=1){
thresh=sort(thresh)
  pe=matrix(0,ncol=3,nrow=length(thresh))
  se=array(0,dim=c(length(thresh),3,2))
  cv=rep(0,length(thresh))
  fi=suppressWarnings(.gp.pw.fit(data[data>thresh[1]]-thresh[1],init=c(1,guess),CI=T,R=R,ncpus=ncpus))

  pe[1,]=fi$fit$PWM
  se[1,,]=t(fi$CI$PWM)
  cv[1]=goftest::cvm.test(data[data>thresh[1]]-thresh[1],function (q,xi,mu,beta)
  {
    (1 - (1 + (xi * (q - mu))/beta)^(-1/xi))
  },xi=pe[1,3],mu=0,beta=pe[1,2])$p.value
  for (i in 2:length(thresh)){

    fi=suppressWarnings(.gp.pw.fit(data[data>thresh[i]]-thresh[i],init=c(1,guess),CI=T,R=R,ncpus=ncpus))
    pe[i,]=fi$fit$PWM
    se[i,,]=t(fi$CI$PWM)
    cv[i]=goftest::cvm.test(data[data>thresh[i]]-thresh[i],function (q,xi,mu,beta)
    {
      (1 - (1 + (xi * (q - mu))/beta)^(-1/xi))
    },xi=pe[i,3],mu=0,beta=pe[i,2])$p.value
  }

  iopt=which((se[,1,1]<1)&(se[,1,2]>1))[1]
  if (is.na(iopt)) {
    stop("A pareto approximation of data is not available here")
  }


  if (!(length(plots) == 1 && is.na(plots))) {
    plots <- sort(unique(plots))
    if (!isTRUE(all(plots %in% 1:3))) {
      stop("Invalid plot selection. Must be a vector of integers containing indices 1, 2 or 3.")
    }
    old.par <- par(no.readonly = TRUE)
    par(mfrow = c(length(plots), 1), mar = c(4.5, 4.5, 3.1,
                                             0.1))
	    on.exit(par(old.par))
    for (i in plots) {

      ylims = c(min(se[, i,],pe[,i]), max(se[, i,],pe[,i]))
      plot(x = thresh, y = pe[, i], pch = 20, xlab = "Threshold",
           bty = "l", ylab = switch(i, expression(kappa),
                                    expression(sigma), expression(xi)),
           ylim = ylims, type = "n")
      polygon(c(thresh, rev(thresh)), c( se[, i,1], rev(se[, i,2])),
              col = "gray57", border = FALSE)
      abline(v=thresh[iopt],lty=2,col="seagreen")
      if (i == min(plots)) {
        title(paste0("Parameter stability plots for EGP"))
      }
      if (i == 1) {
        abline(h = 1, lwd = 0.5, col = "gray20", lty = 2)
      }
      arrows(x0 = thresh, y0 = se[, i,1], y1 = se[, i,2],
             length = 0.05, angle = 90, code = 3)
      points(x = thresh, y = pe[, i], type = "b", pch = 20)
    }

  }
  return(invisible(list(par = pe, CI = se, thresh = thresh, cvm=cv, opt=c(thresh[iopt],pe[iopt,3]))))
}
