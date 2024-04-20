#' Function to apply thresholding on warped wavelet coefficients.
#' The thresholding steps are based on the estimator in Eq. (16) of
#' Kerkyacharian and Picard (2004, Bernoulli).
#'
#' @param x a vector with observations of the covariate/feature.
#'
#' @param y a vector with observations of the response variable.
#'
#' @param wm_filter_num the filter.number variable in the wavethresh::wd
#' function. The higher it is, the smoother the fitted functions are.
#'
#' @param wm_family the wavelet family used in the wavethresh::wd function.
#'
#' @param thresh_par tuning parameter used in the threshold of wavelet
#' coefficients. The threshold value is thresh_par * sqrt(log(n)/n).
#'
#' @param n_levels number of levels used in the warped wavelet estimator.
#'
#' @param boundary_handling boundary handling used in wavethresh 'bc' argument.
#'
#' @return thresholded wavelet coefficients
#'
#' @export

apply_thresh_wavelet_coeff <- function(x=NULL, y=NULL,
                                       wm_filter_num=NULL,
                                       wm_family=NULL,
                                       threshold_type=NULL,
                                       thresh_par=NULL,
                                       n_levels=NULL,
                                       boundary_handling="periodic"){
  # number of observations
  n <- length(y)

  # wavelet transform
  wavelet_decomp <- wavethresh::wd(y, filter.number = wm_filter_num,
                                   family = wm_family, bc=boundary_handling)
  # number of levels
  J_max <- wavelet_decomp$nlevels

  # number of levels considered for the wavelet estimator
  if(is.null(n_levels)){
    J <- ceiling(log2(sqrt(n/log(n))))
  }else{
    if(n_levels > J_max){
      stop(paste("The number of levels must be lower than",J_max))
    }else{
      J <- n_levels
    }
  }

  # the threshold used for wavelet coefficients
  coef_thresh <- thresh_par * sqrt(log(n)/n)

  #1st threshold: zeroing coefficients at scales greater than J
  wav_thr_uppLevels <- wavethresh::threshold(wavelet_decomp, type = "hard",
                                             policy = "manual",
                                             value = Inf,
                                             levels = (J+1):(J_max-1))

  #2nd threshold: thresholding remaining coefficients
  wav_thr <- wavethresh::threshold(wav_thr_uppLevels,
                                   type = threshold_type,
                                   policy = "manual",
                                   value = coef_thresh,
                                   levels = 0:J)

  return(wav_thr)
}
