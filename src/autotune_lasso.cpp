#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List autotune_lasso(SEXP xin,
                    SEXP yin,
                    float alpha = 0.01,
                    bool standardize = true,
                    bool standardize_response = true,
                    bool intercept = true,
                    bool active = false,
                    bool trace_it = false,
                    double tolerance = 1e-4,
                    double beta_tolerance = 1e-3,
                    short int iter_max = 30,
                    short int beta_iter_max = 40,
                    bool PR_norm_l2 = false) {
  
  NumericVector y;
  
  if (!Rf_isVector(yin)) {
    stop("y must be a numeric vector.");
  } else {
    y = yin;
  }
  
  NumericMatrix x;
  
  if (!Rf_isMatrix(xin) and !Rf_isVector(xin)) {
    stop("X must be a numeric matrix");
  } else if (Rf_isMatrix(xin)) {
    x = xin;
  } else {
    stop("X must have more than 1 variables");
  }
  if (Rcpp::as<bool>(any(is_na(y))) || Rcpp::as<bool>(any(is_na(x)))) {
    stop("x or y contains NA values.");
  }
  
  int n = x.nrow();
  int p = x.ncol();
  
  if( y.size() != n ) {
    stop("Length of y doesn't match with the no of observations in x matrix.");
  }
  
  NumericVector beta(p, 0.0);
  NumericVector predmeans(p, 0.0);
  NumericVector predsds(p, 0.0);
  double mu_j, sd_j;
  
  if(standardize) {
    for (int j = 0; j < p; j++) {
      mu_j = mean(x(_, j));
      predmeans[j] = mu_j;
      sd_j = sd(x(_, j));
      predsds[j] = sd_j;
      if(sd_j > 0.0) {
        x(_,j) = (x(_,j) - mu_j) / sd_j;
      } else {
        NumericVector zeroVector(n, 0.0);
        x(_, j) = zeroVector;
      }
    }
  }
  
  double y_mean = 0.0;
  if(standardize_response) {
    y_mean = mean(y);
    y = y - rep(y_mean, y.size());
  }
  
  NumericVector r = clone(y);
  double sigma2est = var(r);
  IntegerVector active_indices = seq_len(p) - 1;
  // NumericMatrix resi_mat(iter_max+1, n);
  // NumericMatrix beta_mat(iter_max+1, p);
  NumericVector vec_sig_beta_count(iter_max);
  NumericVector sigma2_seq(iter_max);
  IntegerVector support_set, old_support_set;
  int max_no_of_preds = std::min(p / 2, n / 2);
  NumericMatrix u(n, max_no_of_preds);
  double init_lambda, lambda_value;
  short int flag = 1;
  short int s;
  bool null_support = FALSE;
  
  NumericVector temp(p);
  for (int j = 0; j < p; j++) {
    temp[j] = sum(x(_, j) * y);
  }
  init_lambda = max(abs(temp)) / n;
  lambda_value = init_lambda * (1.0 / (sigma2est));
  
  // if (lambda0.isNull()) {
  // } else {
  //   lambda_value = as<double>(lambda0);
  // }
  
  int idx;
  int iteration = 1;
  double error = R_PosInf;
  NumericMatrix temp_stor(n, p);
  double lambda_effective, beta_temp;
  NumericVector old_beta(p), partial_res_l1(p), change_in_xbeta;
  double mean_abs_old_beta, mean_abs_diff;
  
  while (error > tolerance && iteration <= iter_max) {
    old_beta = clone(beta);
    old_support_set = clone(support_set);
    s = support_set.size();
    if(s > 0) {
      support_set.erase(support_set.begin(), support_set.end());
    }
    lambda_effective = lambda_value * sigma2est / 2.0;
    
    for (int j = 0; j < p; j++) {
      idx = active_indices[j];
      beta_temp = (sum(x(_, idx) * r) / n) + beta[idx];
      
      if (std::abs(beta_temp) > lambda_effective) {
        beta[idx] = beta_temp - (beta_temp > 0 ? 1 : -1) * lambda_effective;
      } else {
        beta[idx] = 0.0;
      }
      
      change_in_xbeta = x(_, idx) * (beta[idx] - old_beta[idx]);
      r = r - change_in_xbeta;
    }
    
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < p; j++) {
        temp_stor(i, j) = x(i, j) * beta[j] + r[i];
      }
    }
    
    if(!PR_norm_l2) {
      for (int j = 0; j < p; j++) {
        partial_res_l1[j] = sum(abs(temp_stor(_, j)));
      }
    } else {
      for (int j = 0; j < p; j++) {
        partial_res_l1[j] = sd(temp_stor(_, j));
      }
    }
    
    active_indices = seq_len(p) - 1;
    std::sort(active_indices.begin(), active_indices.end(),
              [&partial_res_l1](int a, int b) {
                return partial_res_l1[a] > partial_res_l1[b];
              });
    
    NumericVector ytemp = clone(y);
    int iter = 0;
    NumericVector xk(n);
    NumericVector uk(n);
    NumericVector proj_coeffs(max_no_of_preds - 1);
    
    int idx;
    for (int k = 0; k < max_no_of_preds; k++) {
      idx = active_indices[k];
      xk = x(_, idx);
      u(_, k) = xk;
      if (k > 0){
        for (int k1 = 0; k1 < k; k1++) {
          proj_coeffs[k1] = sum(xk * u(_,k1))/sum(u(_,k1) * u(_,k1));
        }
        for (int k1 = 0; k1 < k; k1++) {
          u(_, k) = u(_, k) - proj_coeffs[k1] * u(_, k1);
        }
      }
      
      uk = u(_, k);
      NumericVector yhat = uk * (sum(uk * ytemp)) / (sum(uk * uk));
      double rss = sum(pow(ytemp - yhat, 2));
      double sst = sum(pow(ytemp - mean(ytemp), 2));
      double ss_reg = sst - rss;
      double f_stat = ((n - k - 1) * ss_reg)/rss;
      double cutoff = R::qf(1 - alpha, 1, n - k - 1, true, false);
      
      if (f_stat < cutoff) {
        break;
      } else {
        ytemp = ytemp - yhat;
        support_set.push_back(idx);
      }
      iter++;
    }
    
    sigma2est = ((n - 1) / std::max(n - iter, 2)) * var(ytemp);
    
    // resi_mat(iteration - 1, _) = r;
    // beta_mat(iteration - 1, _) = beta;
    vec_sig_beta_count[iteration - 1] = iter;
    sigma2_seq[iteration - 1] = sigma2est;
    
    // Calculate error for stopping criteria
    mean_abs_old_beta = mean(abs(old_beta));
    mean_abs_diff = mean(abs(old_beta - beta));
    error = mean_abs_diff / std::max(mean_abs_old_beta, 1e-8);
    
    
    if (trace_it) {Rcout << "\rIteration: " << iteration << std::flush;}
    iteration++;
    
    if (setdiff(support_set, old_support_set).size() == 0) {
      flag--;
    } else {
      flag = 1;
    }
    
    if (flag == 0) {
      break;
    }
  }
  
  iteration--;
  
  if(support_set.size() == 0) {
    null_support = TRUE;
    sigma2est = sigma2_seq[iteration - 2];
    active_indices = active_indices[beta[active_indices] != 0];
    p = active_indices.size();
  }
  
  if (iteration < iter_max) {
    // resi_mat = resi_mat(Range(0, iteration), _);
    // beta_mat = beta_mat(Range(0, iteration), _);
    vec_sig_beta_count = vec_sig_beta_count[Range(0, iteration - 1)];
    sigma2_seq = sigma2_seq[Range(0, iteration - 1)];
  }
  
  
  lambda_effective = lambda_value * sigma2est / 2.0;
  int beta_iteration = 1;
  error = R_PosInf;
  IntegerVector active_set = support_set;
  s = support_set.size();
  old_beta = clone(beta);
  // NumericVector act_pred_count(beta_iter_max);
  
  if(active) {
    while( beta_iteration < 2 ) {
      
      for(int& idx : active_set) {
        beta_temp = (sum(x(_, idx) * r) / n) + beta[idx];
        
        if (std::abs(beta_temp) > lambda_effective) {
          beta[idx] = beta_temp - (beta_temp > 0 ? 1 : -1) * lambda_effective;
        } else {
          beta[idx] = 0.0;
        }
        
        NumericVector change_in_xbeta = x(_, idx) * (beta[idx] - old_beta[idx]);
        r = r - change_in_xbeta;
      }
      
      for (int j = active_set.size(); j < p; j++) {
        idx = active_indices[j];
        if((sum(x(_, idx) * r) - beta[idx])/n >= lambda_effective) {
          active_set.push_back(idx);
        } else {
          r += x(_, idx) * beta[idx];
          beta[idx] = 0;
        }
      }
      
      old_beta = clone(beta);
      // act_pred_count[beta_iteration - 1] = active_set.size();
      beta_iteration++;
    }
    
    while(error > beta_tolerance && beta_iteration <= beta_iter_max) {
      NumericVector old_beta = clone(beta);
      
      for(int& idx : active_set){
        beta_temp = (sum(x(_, idx) * r) / n) + beta[idx];
        
        if (std::abs(beta_temp) > lambda_effective) {
          beta[idx] = beta_temp - (beta_temp > 0 ? 1 : -1) * lambda_effective;
        } else {
          beta[idx] = 0.0;
        }
        
        NumericVector change_in_xbeta = x(_, idx) * (beta[idx] - old_beta[idx]);
        r = r - change_in_xbeta;
      }
      // act_pred_count[beta_iteration - 1] = active_set.size();
      mean_abs_old_beta = mean(abs(old_beta));
      mean_abs_diff = mean(abs(old_beta - beta));
      error = mean_abs_diff / std::max(mean_abs_old_beta, 1e-8);
      
      beta_iteration++;
    }
  } else {
    
    while(error > beta_tolerance && beta_iteration <= beta_iter_max) {
      NumericVector old_beta = clone(beta);
      
      for (int j = 0; j < p; j++) {
        idx = active_indices[j];
        beta_temp = (sum(x(_, idx) * r) / n) + beta[idx];
        
        if (std::abs(beta_temp) > lambda_effective) {
          beta[idx] = beta_temp - (beta_temp > 0 ? 1 : -1) * lambda_effective;
        } else {
          beta[idx] = 0.0;
        }
        
        NumericVector change_in_xbeta = x(_, idx) * (beta[idx] - old_beta[idx]);
        r = r - change_in_xbeta;
      }
      
      mean_abs_old_beta = mean(abs(old_beta));
      mean_abs_diff = mean(abs(old_beta - beta));
      error = mean_abs_diff / std::max(mean_abs_old_beta, 1e-8);
      if (trace_it) {Rcout << "\rLambda converged, Iteration: " << iteration + beta_iteration << std::flush;}
      beta_iteration++;
    }
  }
  
  Rcout << std::endl;
  
  // resi_mat(iteration, _) = r;
  // beta_mat(iteration, _) = beta;
  
  if(standardize) {
    for(int i = 0; i <= iteration; i++) {
      // NumericVector col = beta_mat(i, _);
      for(int j = 0; j < p; j++) {
        if(predsds[j] > 0) {
          beta[j] = beta[j]/predsds[j];
        }
      }
    }
    // beta = beta_mat(iteration, _);
    for(int j = 0; j < p; j++) {
      x(_, j) = x(_, j) * predsds[j] + predmeans[j];
    }
  }
  
  if(standardize_response) {
    y = y + rep(y_mean, y.size());
  }
  
  beta_iteration--;
  
  if (trace_it) {
    Rcout << "\nNo of predictors significant for sigma estimation: " << vec_sig_beta_count[iteration - 1] << std::endl;
  }
  
  double intercept_estimate = 0.0;
  if(intercept) {
    intercept_estimate = y_mean - sum(beta * predmeans);
  }
  
  double final_sigma = rev(sigma2_seq)[0];
  // int total_iter;
  // total_iter = iteration + beta_iteration;
  
  List cd_path_details = List::create(
    _["sorted_predictors"] = active_indices + 1,
    _["sigma_sq_seq"] = sigma2_seq,
    _["no_of_iter_before_lambda_conv"] = iteration,
    _["no_of_iter_after_lambda_conv"] = beta_iteration,
    _["no_of_iterations"] = iteration + beta_iteration,
    _["support_set"] = support_set + 1,
    _["active_set"] = active_set + 1,
    _["count_sig_beta"] = vec_sig_beta_count,
    _["lambda0"] = lambda_value/2.0,
    _["null_support"] = null_support
  );
  
  return List::create(
    // _["residual_matrix"] = resi_mat,
    // _["beta_matrix"] = beta_mat,
    _["beta"] = beta,
    _["a0"] = intercept_estimate,
    _["lambda"] = lambda_effective,
    _["sigma_sq"] = final_sigma,
    _["CD.path.details"] = cd_path_details
  );
}


