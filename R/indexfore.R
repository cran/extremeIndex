#' Function for computing the index for a forecast system vs. climatological forecast. You must provide an indexclim object.
#'
#' @param clim an indexclim object coming from indexclim
#' @param score_fore the time serie of the ensemble forecast's CRPS/MAE. Be careful that score_fore is consistent with "score" in indexclim
#'
#' @return an indexfore object with the index computed vs. climatological forecast and the statistic omega2
#' @export
#'
#' @importFrom evir gpd

indexfore=function (score_fore, clim)
{
  stopifnot(class(clim) == "indexclim")
  testdata = function(yy) {
    exc = score_fore[which(clim$obs > yy)]
    if (clim$estim_xi == TRUE) {
      pa = evir::gpd(clim$obs, threshold = yy, method = "pwm")$par.ests[c(1,2)]
    }
    else {
      pa = c(clim$xi, evir::gpd(clim$obs, threshold = yy,
                                method = "pwm")$par.ests[2])
    }


    ##### comme on conditionne les crps aux obs sup à yy, il peut arriver que crps < yy ce qui posera problème dans le calcul de pgpd : on ramene donc tous les crps < yy à yy :
    exc[exc<yy]=yy
    U=sapply(exc, function(q, xi, mu, beta) {
      if (xi > 0.01) {
        (1 - (1 + (xi * (q - mu))/beta)^(-1/xi))
      }
      else {
        1 - exp(-(q - mu)/beta)
      }
    }, xi = pa[1], mu = yy, beta = pa[2])

    U = sort(U)
    n = length(U)
    k = seq_len(n)
    omega2 = 1/(12 * n) + sum((U - (2 * k - 1)/(2 * n))^2)



    return(omega2)
  }
  cvm = sapply(clim$quantiles, function(q) testdata(q))
  result = list(quantiles = clim$quantiles, index = 1 - clim$index/cvm,
                obs = clim$obs, clim = clim, score = score_fore, estim_xi = clim$estim_xi, omega2=cvm)
  class(result) = "indexfore"
  return(result)
}
