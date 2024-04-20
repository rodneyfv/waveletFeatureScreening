#' Function that extends a vector to make it have a length that is a power of
#' two. The extension is done via a symmetric repetition of values in the
#' border. This function is similar to "wextend('1D','sym',x,len)" in Matlab.
#'
#' @param y vector to be extended.
#'
#' @return A list with three elemets:
#' y_extended: the extended vector;
#' n_first_extend_values: number of repeated values in the beginning of y_extended;
#' n_last_extend_values: number of repeated values in the end of y_extended.
#'
#' @examples
#' y = c(1:9)
#' extend_vector_length_topower2(y)
#' y = c(1:7)
#' extend_vector_length_topower2(y)
#'
#' @export

extend_vector_length_topower2 <- function(y=null){

  n <- length(y)

  # if the length of y is not a power of two, we extend it with a
  # symmetric repetition of values in the border.
  if(log2(n)%%1 != 0){
    # power of 2 that will be the length of the extended vector
    pow2 <- 2^(ceiling(log2(n)))
    # number of observations duplicated at the first entries
    n_first_extend_values <- ceiling((pow2 - n)/2)
    # number of observations duplicated at the last entries
    n_last_extend_values <- pow2 - n - n_first_extend_values
    # vector with extended observations
    y_ext <- rep(0,pow2)

    if((pow2-n)==1){
      y_extended <- c(y[n_first_extend_values:1], y)
      n_last_extend_values <- 0
    }else{
      y_extended <- c(y[n_first_extend_values:1], y, y[n:(n-n_last_extend_values+1)])
    }
  }else{
    # if the length of y is already a power of two, no extension is made
    y_extended <- y
    n_first_extend_values <- 0
    n_last_extend_values <- 0
  }

  return(list(y_extended = y_extended,
              n_first_extend_values = n_first_extend_values,
              n_last_extend_values = n_last_extend_values))
}
