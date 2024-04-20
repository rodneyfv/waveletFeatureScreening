#' Computes warped wavelet coefficients and then apply thresholding
#' on them. The computed estimator is presented in Eq. (16) of
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
#' @return A list with three elemets:
#' wavelet_coeff_thr: thresholded wavelet coefficients;
#' n_first_extend_values: number of observations appended at the beginning of
#' the vector if it does not have a length that is a power of two;
#' n_last_extend_values: number of observations appended at the end of
#' the vector if it does not have a length that is a power of two.
#'
#' @export

compute_warped_wavelet_coefficients <- function(x=NULL, y=NULL,
                                                wm_filter_num=1,
                                                wm_family="DaubExPhase",
                                                thresh_type = "hard",
                                                thresh_par=NULL,
                                                n_levels=NULL,
                                                boundary_handling="periodic"){
  if(length(y)!=length(x)){
    stop("The length of x and y do not match.")
  }
  n <- length(y)
  # sorting y according to increasing order of x
  y_xOrder <- y[order(x)]

  # if the sample size is not a power of two, we extend the vector with a
  # symmetric repetition of values in the borders.
  extension_y <- extend_vector_length_topower2(y_xOrder)
  y_xOrder <- extension_y$y_extended
  n_first_extend_values <- extension_y$n_first_extend_values
  n_last_extend_values <- extension_y$n_last_extend_values

  # number of observations in the ordered vector. Notice that if n is not a
  # power of two, then 'n_xOrder' is different than 'n'
  n_xOrder <- length(y_xOrder)

  # applying threshold on wavelet coefficients
  wavelet_coeff_thr <- apply_thresh_wavelet_coeff(y=y_xOrder,
                                                  wm_filter_num = wm_filter_num,
                                                  wm_family = wm_family,
                                                  threshold_type = thresh_type,
                                                  thresh_par = thresh_par,
                                                  n_levels = n_levels,
                                                  boundary_handling = boundary_handling)

  return(list(wavelet_coeff_thr = wavelet_coeff_thr,
              n_first_extend_values = n_first_extend_values,
              n_last_extend_values = n_last_extend_values))
}
