# waveletFeatureScreening

This package implements a feature screening method that uses a nonparametric wavelet regression. 
It can be used to identify features with a highly nonlinear relation with the response variable.
For more details, see [Fonseca, Morettin, Pinheiro (2024, JCGS)](https://doi.org/10.1080/10618600.2024.2342984).

## Installation

Install and load the _devtools_ package in R and run:
```R
# install.packages("devtools")
devtools::install_github("rodneyfv/waveletFeatureScreening")
```
## An example

```R
set.seed(2024)
sample_size = 100
number_parameters = 500
matrix_features = matrix(data = runif(sample_size * number_parameters),
                         nrow=sample_size,
                         ncol=number_parameters)
# The first two features are used to generate a response
response_variable = sin(2*pi*matrix_features[,1]) +
                        3*(matrix_features[,2]^2) + 0.5*rnorm(sample_size)
# Using the wavelet feature screening to rank the features
wws <- wavelet_screening(matrix_features, response_variable)
# most relevant features
wws$rank[1:20]

```

