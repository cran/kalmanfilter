#kalmanfilter 2.1

## Major changes

* added kalman_filter_tvp for time varying parameters

kalmanfilter 2.0.2

## Minor changes

* updated to use sentinel "_PACKAGE"

#kalmanfilter 2.0.1

## Minor changes

* Updated Rcpp Rginv function to matrix that has no positive values in SVD.

#kalmanfilter 2.0

## Minor changes

* Updated Rcpp code to handle deprecation of << assignment operator

## Major changes

* Updated Rcpp code to handle NULL values for the exogenous and weight matrices
* Updated Rcpp code to handle cases without exogenous coefficient matrices in the state space model
* Removed R wrapper functions
* Added Stock and Watson dynamic common factor model example to the vignette

# kalmanfilter 1.0

## Bug fixes

* none

## Major changes

* Initial release
