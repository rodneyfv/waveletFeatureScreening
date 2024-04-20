#' Feature screening via warped wavelets.
#'
#' @param x nxp matrix of covariates.
#'
#' @param y n-dimensional vector with response variable.
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
#' @return A list with two elemets:
#' rank: importance of the variable (important variables are expected to
#' have the lowest indices);
#' measurement: the screening utility used to compute the rank.
#'
#' @examples
#' set.seed(2024)
#' sample_size = 100
#' number_parameters = 500
#' matrix_features = matrix(data = runif(sample_size * number_parameters),
#'                          nrow=sample_size,
#'                          ncol=number_parameters)
#' response_variable = sin(2*pi*matrix_features[,1]) +
#'                       3*(matrix_features[,2]^2) + 0.5*rnorm(sample_size)
#' wws <- wavelet_screening(matrix_features, response_variable)
#' wws$rank[1:20]
#'
#' @export

wavelet_screening <- function(x, y, wm_filter_num=1,
                                  wm_family="DaubExPhase",
                                  thresh_type = "hard",
                                  thresh_par = 1, n_levels = 1,
                                  boundary_handling="periodic"){
  # sample size
  n <- nrow(x)
  # dimension
  p <- ncol(x)

  # vector to keep the norm of warped wavelet coefficients
  # of each feature
  Energyvector   <- rep(0, p)
  # normalized response
  y_normalized <- (y - mean(y))/sd(y)
  for (i in 1:p){
    # thresholded warped wavelet coefficients
    wav_thr <- compute_warped_wavelet_coefficients(x=x[,i],
                                                   y=y_normalized,
                                                   wm_filter_num = wm_filter_num,
                                                   wm_family = wm_family,
                                                   thresh_type = thresh_type,
                                                   thresh_par = thresh_par,
                                                   n_levels = n_levels,
                                                   boundary_handling=boundary_handling)
    # squared norm of thresholded warped wavelet coefficients
    Energyvector[i] <- sum(wav_thr$wavelet_coeff_thr$C^2) +
      sum(wav_thr$wavelet_coeff_thr$D^2)
  }

  scree_meas <- Energyvector
  names(scree_meas) <- colnames(x)
  # rank of the features (rank 1 has the highest score, etc)
  scree_rank <- p - rank(Energyvector) + 1
  names(scree_rank) <- colnames(x)

  return(list(rank = scree_rank, measurement = scree_meas))
}
