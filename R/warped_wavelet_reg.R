#' Function to fit a warped wavelet regression using the
#' wd function of the wavethresh package.
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
#' x: the sorted values of the feature
#' y_fit: fitted (predicted) values for sorted values in x
#' residuals: residuals of the regression (y - y_fit)
#'
#' @examples
#' x <- runif(50)
#' y = sin(2*pi * x) + rnorm(50)*.1
#' plot(x,y)
#' wave_reg <- warped_wavelet_reg(x=x, y=y, wm_filter_num=8,
#'                                wm_family="DaubLeAsymm", thresh_type = "hard",
#'                                thresh_par=1, n_levels=2)
#' lines(wave_reg$x,wave_reg$y_fit)
#'
#' @export

warped_wavelet_reg <- function(x=NULL, y=NULL,
                               wm_filter_num=9, wm_family="DaubLeAsymm",
                               thresh_type = "hard",
                               thresh_par=1,
                               n_levels=1,
                               boundary_handling="periodic"){
  if(length(y)!=length(x)){
    stop("The length of x and y do not match.")
  }
  n <- length(y)
  order_x <- order(x)

  wav_thr <- compute_warped_wavelet_coefficients(x=x, y=y,
                                                 wm_filter_num = wm_filter_num,
                                                 wm_family = wm_family,
                                                 thresh_type = thresh_type,
                                                 thresh_par = thresh_par,
                                                 n_levels = n_levels,
                                                 boundary_handling = boundary_handling)

  # fitted values in the same order of input y
  y_fit <- rep(0, n)
  wavelet_fit <- wavethresh::wr(wav_thr$wavelet_coeff_thr)

  # checking if the observations had to be duplicated, i.e., the sample size is
  # not a power of two
  if(log2(n)%%1 != 0){
    id_begin <- wav_thr$n_first_extend_values + 1
    id_end <- wav_thr$n_first_extend_values + n
    y_fit <- wavelet_fit[id_begin:id_end]
  }else{
    y_fit <- wavelet_fit
  }
  # residuals
  res_fit <- y[order_x] - y_fit

  return(list(x = sort(x),
              y_fit = y_fit,
              residuals = res_fit))
}
