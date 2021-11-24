#' Function which computes the index for the climatological CRPS/MAE. You must provide the observations. If you computes climatological CRPS/MAE previously, you can add the corresponding vector
#'
#' @param y The observations
#' @param thresh Vector of thresholds where you want to compute the index
#' @param score_clim If not NULL, must be the time serie of the CRPS/MAE of the climatology. It is recommended to compute CRPS/MAE out of this function
#' @param xi Shape parameter of the GP ( xi > 0)
#' @param score A character string indicating if you want to work with CRPS ("crps") or MAE ("mae"), by default "crps"
#' @param estim_xi If you want xi estimated for each threshold (for numerical reasons for instance)
#' @return An indexclim object containing xi, y, the score time serie, the score considered, the index values, and the corresponding quantiles of the observations
#' @export
#'
#' @importFrom evir gpd
indexclim=function (y, thresh = NULL, score_clim = NULL, xi = NULL, score = "crps", estim_xi = FALSE)
{
  stopifnot(is.numeric(y))
  stopifnot(!is.null(thresh))
  stopifnot(!is.null(xi))
  compscore = is.null(score_clim)
  thresh = sort(thresh)
  if ((score != "crps") & (score != "mae")) {
    stop("score must be crps or mae")
  }
  if (compscore) {
    if (score == "crps") {
      alpha = mean(y) - 2 * mean(seq(0, 1, length.out = length(y)) *
                                   sort(y))
      score_clim = sapply(y[y > thresh[1]], function(i) mean(abs(y -
                                                                   i)) + alpha)
    }
    if (score == "mae") {
      score_clim = sapply(y[y > thresh[1]], function(i) mean(abs(y -
                                                                   i)))
    }
  }
  testdata = function(yy) {
    if (compscore) {
      exc = score_clim[which(y[y > thresh[1]] >= yy)]
    }
    else {
      exc = score_clim[which(y > yy)]
    }
    #pa = c(xi, evir::gpd(y, threshold = yy, method = "ml")$par.ests[2])
    if (estim_xi == TRUE) {
      pa = evir::gpd(y, threshold = yy, method = "ml")$par.ests[c(1,
                                                                  2)]
    }
    else {
      pa = c(xi, 1+xi*yy)
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
  cvm = sapply(thresh, function(q) testdata(q))
  result = list(quantiles = thresh, crps = score_clim, index = cvm,
                obs = y, xi = xi, score = score, estim_xi = estim_xi)
  class(result) = "indexclim"
  return(result)
}
