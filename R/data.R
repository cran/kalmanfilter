#' Treasuries
#' 
#' @format data.table with columns DATE, VARIABLE, VALUE, and MATURITY 
#' The data is quarterly frequency with variables DGS1MO, DGS3MO, DGS6MO, DGS1, 
#' DGS2, DGS3, DGS5, DGS7, DGS10, DGS20, and DGS30
#' @source FRED
#' @usage data(treasuries)
"treasuries"

#' Stock and Watson Dynamic Common Factor Data Set
#' 
#' @format data.table with columns DATE, VARIABLE, VALUE, and MATURITY 
#' The data is monthly frequency with variables ip (industrial production), 
#' gmyxpg (total personal income less transfer payments in 1987 dollars), 
#' mtq (total manufacturing and trade sales in 1987 dollars), 
#' lpnag (employees on non-agricultural payrolls), and
#' dcoinc (the coincident economic indicator)
#' @source Kim, Chang-Jin and Charles R. Nelson (1999) "State-Space Models with Regime Switching: Classical and Gibbs-Sampling Approaches with Applications" <doi:10.7551/mitpress/6444.001.0001><http://econ.korea.ac.kr/~cjkim/>. 
#' @usage data(sw_dcf)
"sw_dcf"

