#include <RcppArmadillo.h>
#include <Rcpp.h>

// Correctly setup the build environment
// [[Rcpp::depends(RcppArmadillo)]]

// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

//' R's implementation of the Moore-Penrose pseudo matrix inverse
//'
//' @param m matrix
//' @return matrix inverse of m
// [[Rcpp::export]]
arma::mat Rginv(const arma::mat& m){
  arma::mat U, V;
  arma::vec S;
  arma::svd(U, S, V, m, "dc");
  arma::uvec Positive = arma::find(S > 0.0);
  if(Positive.size() == 0){
    return arma::zeros(m.n_rows, m.n_cols);
  }else if(all(Positive)){
    arma::mat D = diagmat(S);
    return V * (1/D * U.t());
  }else if(!any(Positive)){
    return arma::zeros(m.n_rows, m.n_cols);
  }else{
    S.elem(Positive) = 1/S.elem(Positive);
    arma::mat D = diagmat(S);
    return V * D * U.t();
  }
}

//' Generalized matrix inverse
//'
//' @param m matrix
//' @return matrix inverse of m
// [[Rcpp::export]]
arma::mat gen_inv(arma::mat& m){
  arma::mat out(m.n_rows, m.n_cols);
  try{
    out = inv(m);
  }catch(std::exception &ex){
    out = Rginv(m);
  }
  return out;
}

//' Check if list contains a name
//' 
//' @param s a string name
//' @param L a list object
//' @return boolean
// [[Rcpp::export]]                                                                                                                                           
bool contains(std::string s, Rcpp::List L){                                                                                                                  
  Rcpp::CharacterVector nv = L.names();                                                                                                                     
  for(int i = 0; i < nv.size(); i++){                                                                                                                         
    if(std::string(nv[i]) == s){                                                                                                                        
      return true;                                                                                                                                      
    }                                                                                                                                                     
  }                                                                                                                                                         
  return false;                                                                                                                                             
}     

//' Kalman Filter
//'
//' @param ssm list describing the state space model, must include names
//' B0 - N_b x 1 matrix, initial guess for the unobserved components 
//' P0 - N_b x N_b matrix, initial guess for the covariance matrix of the unobserved components
//' Dm - N_b x 1 matrix, constant matrix for the state equation
//' Am - N_y x 1 matrix, constant matrix for the observation equation
//' Fm - N_b X p matrix, state transition matrix
//' Hm - N_y x N_b matrix, observation matrix
//' Qm - N_b x N_b matrix, state error covariance matrix
//' Rm - N_y x N_y matrix, state error covariance matrix
//' betaO - N_y x N_o matrix, coefficient matrix for the observation exogenous data
//' betaS - N_b x N_s matrix, coefficient matrix for the state exogenous data
//' @param yt N x T matrix of data
//' @param Xo N_o x T matrix of exogenous observation data
//' @param Xs N_s x T matrix of exogenous state 
//' @param weight column matrix of weights, T x 1
//' @param smooth boolean indication whether to run the backwards smoother
//' @return list of matrices and cubes output by the Kalman filter
//' @examples
//' #Nelson-Siegel dynamic factor yield curve
//' library(kalmanfilter)
//' library(data.table)
//' data(treasuries)
//' tau = unique(treasuries$maturity)
//'
//' #Set up the state space model
//' ssm = list()
//' ssm[["Fm"]] = rbind(c(0.9720, -0.0209, -0.0061), 
//'                     c(0.1009 , 0.8189, -0.1446), 
//'                     c(-0.1226, 0.0192, 0.8808))
//' ssm[["Dm"]] = matrix(c(0.1234, -0.2285, 0.2020), nrow = nrow(ssm[["Fm"]]), ncol = 1)
//' ssm[["Qm"]] = rbind(c(0.1017, 0.0937, 0.0303), 
//'                     c(0.0937, 0.2267, 0.0351), 
//'                     c(0.0303, 0.0351, 0.7964))
//' ssm[["Hm"]] = cbind(rep(1, 11),
//'                     -(1 - exp(-tau*0.0423))/(tau*0.0423), 
//'                     (1 - exp(-tau*0.0423))/(tau*0.0423) - exp(-tau*0.0423))
//' ssm[["Am"]] = matrix(0, nrow = length(tau), ncol = 1)
//' ssm[["Rm"]] = diag(c(0.0087, 0, 0.0145, 0.0233, 0.0176, 0.0073, 
//'                      0, 0.0016, 0.0035, 0.0207, 0.0210))
//' ssm[["B0"]] = matrix(c(5.9030, -0.7090, 0.8690), nrow = nrow(ssm[["Fm"]]), ncol = 1)
//' ssm[["P0"]] = diag(rep(0.0001, nrow(ssm[["Fm"]])))
//'     
//' #Convert to an NxT matrix
//' yt = dcast(treasuries, "date ~ maturity", value.var = "value")
//' yt = t(yt[, 2:ncol(yt)])
//' kf = kalman_filter(ssm, yt, smooth = TRUE)   
// [[Rcpp::export]]
Rcpp::List kalman_filter_cpp(Rcpp::List& ssm, const arma::mat& yt, 
                         Rcpp::Nullable<Rcpp::NumericMatrix> Xo = R_NilValue,
                         Rcpp::Nullable<Rcpp::NumericMatrix> Xs = R_NilValue,
                         Rcpp::Nullable<Rcpp::NumericMatrix> weight = R_NilValue,
                         bool smooth = false){
  
  //Initialize matrices
  int n_cols = yt.n_cols;
  int n_rows = yt.n_rows;
  arma::mat B0 = ssm["B0"];
  arma::mat P0 = ssm["P0"];
  arma::mat Dm = ssm["Dm"];
  arma::mat Am = ssm["Am"];
  arma::mat Fm = ssm["Fm"];
  arma::mat Hm = ssm["Hm"];
  arma::mat Qm = ssm["Qm"];
  arma::mat Rm = ssm["Rm"];
  arma::mat betaO;
  arma::mat betaS;
  arma::mat X_o;
  arma::mat X_s;
  
  //Set the exogenous matrices
  if(Xo.isNotNull()){
    X_o = Rcpp::as<arma::mat>(Xo);
  }else{
    X_o = arma::zeros(1, n_cols);
  }
  if(Xs.isNotNull()){
    X_s = Rcpp::as<arma::mat>(Xs);
  }else{
    X_s = arma::zeros(1, n_cols);
  }
  
  //Set the exogenous coefficient matrices
  if(contains("betaO", ssm)){
    betaO = Rcpp::as<arma::mat>(ssm["betaO"]);
  }else{
    betaO = arma::zeros(n_rows, 1);
  }
  if(contains("betaS", ssm)){
    betaS = Rcpp::as<arma::mat>(ssm["betaS"]);
  }else{
    betaS = arma::zeros(Fm.n_rows, 1);
  }
  
  //Initialize the filter
  arma::mat B_LL = B0;
  arma::mat P_LL = P0;
  arma::mat lnl(1, n_cols, arma::fill::zeros);
  arma::mat w(n_cols, 1);
  
  //Rescale the weights
  if(weight.isNotNull()){
    w = Rcpp::as<arma::mat>(weight);
  }else{
    w = arma::ones(n_cols, 1);
  }
  w = w * n_cols/arma::as_scalar(sum(w));  
  
  //Define the storage matrices
  arma::mat B_tt(Fm.n_rows, n_cols);
  arma::mat B_tl(B_tt.n_rows, n_cols);
  arma::cube P_tt(Fm.n_rows, Fm.n_rows, n_cols);
  arma::cube P_tl(Fm.n_rows, Fm.n_rows, n_cols);
  arma::mat N_t(n_rows, n_cols);
  arma::cube F_t = arma::ones(n_rows, n_rows, n_cols)*R_PosInf;
  arma::cube K_t(Fm.n_rows, n_rows, n_cols);
  arma::uvec non_na_idx;
  arma::uvec iv(1);
  arma::mat y_tl(n_rows, n_cols);
  arma::mat y_tt(n_rows, n_cols);
  
  //Define some matrix transforms
  arma::mat Fm_t = Fm.t();
  arma::mat Hm_t = Hm.t();
  arma::mat F_t_inv = arma::zeros(F_t.n_rows, F_t.n_cols);
  arma::mat F_t_submat;
  
  //Kalman filter routine
  for(int i = 0; i < n_cols; i++){
    //Find the non-missing values
    non_na_idx = arma::find_finite(yt.col(i));
    
    //Initial estimates conditional on t-1
    //Initial estimate of unobserved values conditional on t-1
    B_tl.col(i) = Dm + Fm * B_LL + betaS * X_s.col(i); 
    //Initial estimate of the covariance matrix conditional on t-1
    P_tl.slice(i) = Fm * P_LL * Fm_t + Qm; 
    //Prediction error conditional on t-1
    N_t.col(i) = yt.col(i) - Am - Hm * B_tl.col(i) - betaO * X_o.col(i); 
    F_t_inv = arma::zeros(F_t.slice(i).n_rows, F_t.slice(i).n_cols);

    if(!non_na_idx.is_empty()){
      //Variance of the prediction error conditional on t-1
      F_t.slice(i).submat(non_na_idx, non_na_idx) = Hm.rows(non_na_idx) * P_tl.slice(i) * Hm_t.cols(non_na_idx) + Rm.submat(non_na_idx, non_na_idx);
      F_t_submat = F_t.slice(i).submat(non_na_idx, non_na_idx);
      F_t_inv.submat(non_na_idx, non_na_idx) = gen_inv(F_t_submat);
    }
    //Kalman gain conditional on t-1
    K_t.slice(i) = P_tl.slice(i) * Hm_t * F_t_inv; 
    
    //Final estimates conditional on t
    if(!non_na_idx.is_empty()){
      iv[0] = i;
      //Final estimate of the unobserved values
      B_tt.col(i) = B_tl.col(i) + K_t.slice(i).cols(non_na_idx) * N_t.submat(non_na_idx, iv); 
      //Final estimate of the covariance matrix
      P_tt.slice(i) = P_tl.slice(i) - K_t.slice(i).cols(non_na_idx) * Hm.rows(non_na_idx) * P_tl.slice(i); 
      //Update the log likelihood
      lnl.col(i) = 0.5*arma::as_scalar((log(det(F_t.slice(i).submat(non_na_idx, non_na_idx))) + N_t.submat(non_na_idx, iv).t() * F_t_inv.submat(non_na_idx, non_na_idx) * N_t.submat(non_na_idx, iv)));
    }else{
      B_tt.col(i) = B_tl.col(i);
      P_tt.slice(i) = P_tl.slice(i);
    }
    
    //Get the finalized predictions
    y_tl.col(i) = Am + Hm * B_tl.col(i) + betaO * X_o.col(i); 
    y_tt.col(i) = Am + Hm * B_tt.col(i) + betaO * X_o.col(i);
    
    //Reinitialize for the next iteration
    B_LL = B_tt.col(i);
    P_LL = P_tt.slice(i);
  }
  
  if(smooth == true){
    int t = B_tt.n_cols - 1;
    arma::mat Ptt_x_Ft_x_PtInv = P_tt.slice(t - 1) * Fm_t * Rginv(P_tl.slice(t));
    
    //Update the estimates of B_tt and P_tt using future information
    for(int i = t - 1; i >= 0; i--){
      Ptt_x_Ft_x_PtInv = P_tt.slice(i) * Fm_t * Rginv(P_tl.slice(i + 1));
      B_tt.col(i) = B_tt.col(i) + Ptt_x_Ft_x_PtInv * (B_tt.col(i + 1) - B_tl.col(i + 1));
      P_tt.slice(i) = P_tt.slice(i) + Ptt_x_Ft_x_PtInv * (P_tt.slice(i + 1) - P_tl.slice(i + 1)) * Ptt_x_Ft_x_PtInv.t();
      
      //Get the finalized predictions
      y_tt.col(i) = Am + Hm * B_tt.col(i) + betaO * X_o.col(i); 
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("lnl") = arma::as_scalar(-(lnl * w)),
                            Rcpp::Named("y_tl") = y_tl,
                            Rcpp::Named("y_tt") = y_tt,
                            Rcpp::Named("B_tl") = B_tl,
                            Rcpp::Named("B_tt") = B_tt,
                            Rcpp::Named("P_tl") = P_tl,
                            Rcpp::Named("P_tt") = P_tt,
                            Rcpp::Named("F_t") = F_t,
                            Rcpp::Named("N_t") = N_t,
                            Rcpp::Named("K_t") = K_t);
}

 //' Kalman Filter for Time Varying Parameters
 //' 
 //' @param ssm list describing the state space model, must include names
 //' B0 - N_b x 1 matrix, initial guess for the unobserved components 
 //' P0 - N_b x N_b matrix, initial guess for the covariance matrix of the unobserved components
 //' Dm - N_b x 1 x T array or matrices, constant matrix for the state equation
 //' Am - N_y x 1 x T array of matrices, constant matrix for the observation equation
 //' Fm - N_b X p x T array of matrices, state transition matrix
 //' Hm - N_y x N_b x T array of matrices, observation matrix
 //' Qm - N_b x N_b x T array of matrices, state error covariance matrix
 //' Rm - N_y x N_y x T array of matrices, state error covariance matrix
 //' betaO - N_y x N_o x T array of matrices, coefficient matrix for the observation exogenous data
 //' betaS - N_b x N_s x T array of matrices, coefficient matrix for the state exogenous data
 //' @param yt N x T matrix of data
 //' @param Xo N_o x T matrix of exogenous observation data
 //' @param Xs N_s x T matrix of exogenous state 
 //' @param weight column matrix of weights, T x 1
 //' @param smooth boolean indication whether to run the backwards smoother
 //' @return list of cubes and matrices output by the Kim filter
 //' @examples
 // [[Rcpp::export]]
 Rcpp::List kalman_filter_tvp_cpp(Rcpp::List& ssm, const arma::mat& yt, 
                              Rcpp::Nullable<Rcpp::NumericMatrix> Xo = R_NilValue,
                              Rcpp::Nullable<Rcpp::NumericMatrix> Xs = R_NilValue,
                              Rcpp::Nullable<Rcpp::NumericMatrix> weight = R_NilValue,
                              bool smooth = false){
   //Initialize matrices
   int n_cols = yt.n_cols;
   int n_rows = yt.n_rows;
   arma::mat B0 = ssm["B0"];
   arma::mat P0 = ssm["P0"];
   arma::cube Dm = ssm["Dm"];
   arma::cube Am = ssm["Am"];
   arma::cube Fm = ssm["Fm"];
   arma::cube Hm = ssm["Hm"];
   arma::cube Qm = ssm["Qm"];
   arma::cube Rm = ssm["Rm"];
   arma::cube betaO;
   arma::cube betaS;
   arma::mat X_o;
   arma::mat X_s;
   
   //Set the exogenous matrices
   if(Xo.isNotNull()){
     X_o = Rcpp::as<arma::mat>(Xo);
   }else{
     X_o = arma::zeros(1, n_cols);
   }
   if(Xs.isNotNull()){
     X_s = Rcpp::as<arma::mat>(Xs);
   }else{
     X_s = arma::zeros(1, n_cols);
   }
   
   //Set the exogenous coefficient matrices
   if(contains("betaO", ssm)){
     betaO = Rcpp::as<arma::cube>(ssm["betaO"]);
   }else{
     betaO = arma::zeros(n_rows, 1, n_cols);
   }
   if(contains("betaS", ssm)){
     betaS = Rcpp::as<arma::cube>(ssm["betaS"]);
   }else{
     betaS = arma::zeros(Fm.slice(1).n_rows, 1, n_cols);
   }
   
   //Define the storage matrices
   arma::mat B_tt(Fm.slice(1).n_rows, n_cols);
   arma::mat B_tl(B_tt.n_rows, n_cols);
   arma::cube P_tt(Fm.n_rows, Fm.n_rows, n_cols);
   arma::cube P_tl(Fm.slice(1).n_rows, Fm.slice(1).n_rows, n_cols);
   arma::mat N_t(n_rows, n_cols);
   arma::cube F_t = arma::ones(n_rows, n_rows, n_cols)*R_PosInf;
   arma::cube K_t(Fm.slice(1).n_rows, n_rows, n_cols);
   arma::uvec non_na_idx;
   arma::uvec iv(1);
   arma::mat y_tl(n_rows, n_cols);
   arma::mat y_tt(n_rows, n_cols);
   
   //Define some matrix transforms
   arma::mat F_t_inv = arma::zeros(Fm.slice(1).n_cols, Fm.slice(1).n_rows);
   arma::mat F_t_submat;
   
   //Initialize the filter
   arma::mat B_LL = B0;
   arma::mat P_LL = P0;
   arma::mat lnl(1, n_cols, arma::fill::zeros);
   arma::mat w(n_cols, 1);
   
   //Rescale the weights
   if(weight.isNotNull()){
     w = Rcpp::as<arma::mat>(weight);
   }else{
     w = arma::ones(n_cols, 1);
   }
   w = w * n_cols/arma::as_scalar(sum(w));  
   
   //Kalman filter routine
   for(int i = 0; i < n_cols; i++){
     //Find the non-missing values
     non_na_idx = arma::find_finite(yt.col(i));
     
     //Initial estimates conditional on t-1
     //Initial estimate of unobserved values conditional on t-1
     B_tl.col(i) = Dm.slice(i) + Fm.slice(i) * B_LL + betaS.slice(i) * X_s.col(i); 
     //Initial estimate of the covariance matrix conditional on t-1
     P_tl.slice(i) = Fm.slice(i) * P_LL * Fm.slice(i).t() + Qm.slice(i); 
     //Prediction error conditional on t-1
     N_t.col(i) = yt.col(i) - Am.slice(i) - Hm.slice(i) * B_tl.col(i) - betaO.slice(i) * X_o.col(i); 
     F_t_inv = arma::zeros(F_t.slice(i).n_rows, F_t.slice(i).n_cols);
     
     if(!non_na_idx.is_empty()){
       //Variance of the prediction error conditional on t-1
       F_t.slice(i).submat(non_na_idx, non_na_idx) = Hm.slice(i).rows(non_na_idx) * P_tl.slice(i) * Hm.slice(i).rows(non_na_idx).t() + Rm.slice(i).submat(non_na_idx, non_na_idx);
       F_t_submat = F_t.slice(i).submat(non_na_idx, non_na_idx);
       F_t_inv.submat(non_na_idx, non_na_idx) = gen_inv(F_t_submat);
     }
     //Kalman gain conditional on t-1
     K_t.slice(i) = P_tl.slice(i) * Hm.slice(i).t() * F_t_inv; 
     
     //Final estimates conditional on t
     if(!non_na_idx.is_empty()){
       iv[0] = i;
       //Final estimate of the unobserved values
       B_tt.col(i) = B_tl.col(i) + K_t.slice(i).cols(non_na_idx) * N_t.submat(non_na_idx, iv); 
       //Final estimate of the covariance matrix
       P_tt.slice(i) = P_tl.slice(i) - K_t.slice(i).cols(non_na_idx) * Hm.slice(i).rows(non_na_idx) * P_tl.slice(i); 
       //Update the log likelihood
       lnl.col(i) = 0.5*arma::as_scalar((log(det(F_t.slice(i).submat(non_na_idx, non_na_idx))) + N_t.submat(non_na_idx, iv).t() * F_t_inv.submat(non_na_idx, non_na_idx) * N_t.submat(non_na_idx, iv)));
     }else{
       B_tt.col(i) = B_tl.col(i);
       P_tt.slice(i) = P_tl.slice(i);
     }
     
     //Get the finalized predictions
     y_tl.col(i) = Am.slice(i) + Hm.slice(i) * B_tl.col(i) + betaO.slice(i) * X_o.col(i); 
     y_tt.col(i) = Am.slice(i) + Hm.slice(i) * B_tt.col(i) + betaO.slice(i) * X_o.col(i);
     
     //Reinitialize for the next iteration
     B_LL = B_tt.col(i);
     P_LL = P_tt.slice(i);
   }
   
   if(smooth == true){
     int t = B_tt.n_cols - 1;
     arma::mat Ptt_x_Ft_x_PtInv = P_tt.slice(t - 1) * Fm.slice(t).t() * Rginv(P_tl.slice(t));
     
     //Update the estimates of B_tt and P_tt using future information
     for(int i = t - 1; i >= 0; i--){
       Ptt_x_Ft_x_PtInv = P_tt.slice(i) * Fm.slice(i).t() * Rginv(P_tl.slice(i + 1));
       B_tt.col(i) = B_tt.col(i) + Ptt_x_Ft_x_PtInv * (B_tt.col(i + 1) - B_tl.col(i + 1));
       P_tt.slice(i) = P_tt.slice(i) + Ptt_x_Ft_x_PtInv * (P_tt.slice(i + 1) - P_tl.slice(i + 1)) * Ptt_x_Ft_x_PtInv.t();
       
       //Get the finalized predictions
       y_tt.col(i) = Am.slice(i) + Hm.slice(i) * B_tt.col(i) + betaO.slice(i) * X_o.col(i); 
     }
   }
   
   return Rcpp::List::create(Rcpp::Named("lnl") = arma::as_scalar(-(lnl * w)),
                             Rcpp::Named("y_tl") = y_tl,
                             Rcpp::Named("y_tt") = y_tt,
                             Rcpp::Named("B_tl") = B_tl,
                             Rcpp::Named("B_tt") = B_tt,
                             Rcpp::Named("P_tl") = P_tl,
                             Rcpp::Named("P_tt") = P_tt,
                             Rcpp::Named("F_t") = F_t,
                             Rcpp::Named("N_t") = N_t,
                             Rcpp::Named("K_t") = K_t);
 }

//git config remote.origin.url https://ajhubb@bitbucket.org/ajhubb/kalmanfilter.git
