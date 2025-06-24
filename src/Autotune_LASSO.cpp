#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List autotune_lasso_l2(NumericMatrix x,
                      NumericVector y,
                      float alpha = 0.01,
                      Nullable<double> lambda = R_NilValue,
                      double tolerance = 1e-4,
                      double beta_tolerance = 1e-3,
                      short int iter_max = 30,
                      short int beta_iter_max = 40,
                      bool verbose = false) {
  
  int n = x.nrow();
  int p = x.ncol();
  NumericVector beta(p, 0.0);
  NumericVector r = clone(y);
  double sigma2est = var(r);
  IntegerVector active_indices = seq_len(p) - 1;
  NumericMatrix resi_mat(iter_max+1, n);
  NumericMatrix beta_mat(iter_max+1, p);
  NumericVector vec_sig_beta_count(iter_max);
  NumericVector sigma2_seq(iter_max);
  IntegerVector support_set, old_support_set;
  int max_no_of_preds = std::min(p / 2, n / 2);
  NumericMatrix u(n, max_no_of_preds);
  double init_lambda, lambda_value;
  int flag = 1;
  int s;
  
  if (lambda.isNull()) {
    NumericVector temp(p);
    for (int j = 0; j < p; j++) {
      temp[j] = sum(x(_, j) * y);
    }
    init_lambda = max(abs(temp)) / n;
    lambda_value = init_lambda * (1.0 / (sigma2est));
  } else {
    lambda_value = as<double>(lambda);
  }
  
  int idx;
  int iteration = 1;
  double error = R_PosInf;
  NumericMatrix temp_stor(n, p);
  double lambda_effective, beta_temp;
  NumericVector old_beta(p), partial_res_sd(p), change_in_xbeta;
  double mean_abs_old_beta, mean_abs_diff;
  
  while (error > tolerance && iteration <= iter_max) {
    old_beta = clone(beta);
    old_support_set = clone(support_set);
    s = support_set.size();
    if(s > 0) {
      support_set.erase(0, s);
    }
    lambda_effective = lambda_value * sigma2est / 2.0;
    
    for (int j = 0; j < p; j++) {
      idx = active_indices[j];
      beta_temp = sum(x(_, idx) * r) / (n - 1) + beta[idx];
      
      if (std::abs(beta_temp) > lambda_effective) {
        beta[idx] = beta_temp - (beta_temp > 0 ? 1 : -1) * lambda_effective;
      } else {
        beta[idx] = 0.0;
      }
      
      change_in_xbeta = x(_, idx) * (beta[idx] - old_beta[idx]);
      r = r - change_in_xbeta;
    }
    
    if (iteration < 60) {
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < p; j++) {
          temp_stor(i, j) = x(i, j) * beta[j] + r[i];
        }
      }
      
      for (int j = 0; j < p; j++) {
        partial_res_sd[j] = sd(temp_stor(_, j));
      }
      
      active_indices = seq_len(p) - 1;
      std::sort(active_indices.begin(), active_indices.end(),
                [&partial_res_sd](int a, int b) {
                  return partial_res_sd[a] > partial_res_sd[b];
                });
    }
    
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
    
    resi_mat(iteration - 1, _) = r;
    beta_mat(iteration - 1, _) = beta;
    vec_sig_beta_count[iteration - 1] = iter;
    sigma2_seq[iteration - 1] = sigma2est;
    
    // Calculate error for stopping criteria
    mean_abs_old_beta = mean(abs(old_beta));
    mean_abs_diff = mean(abs(old_beta - beta));
    error = mean_abs_diff / std::max(mean_abs_old_beta, 1e-8);
    
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
  
  if (iteration < iter_max) {
    resi_mat = resi_mat(Range(0, iteration), _);
    beta_mat = beta_mat(Range(0, iteration), _);
    vec_sig_beta_count = vec_sig_beta_count[Range(0, iteration - 1)];
    sigma2_seq = sigma2_seq[Range(0, iteration - 1)];
  }
  
  
  lambda_effective = lambda_value * sigma2est / 2.0;
  int beta_iteration = 1;
  error = R_PosInf;
  
  while(error > beta_tolerance && beta_iteration <= beta_iter_max) {
    NumericVector old_beta = clone(beta);
    
    for (int j = 0; j < p; j++) {
      idx = active_indices[j];
      beta_temp = sum(x(_, idx) * r) / (n - 1) + beta[idx];
      
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
    
    beta_iteration++;
  }
  
  resi_mat(iteration, _) = r;
  beta_mat(iteration, _) = beta;
  
  if (verbose) {
    Rcout << "\nNo of iterations: " << iteration << "\nNo of Beta iterations: " << beta_iteration
          << "\nNo of outliers in last step: " << vec_sig_beta_count[iteration - 1] << std::endl;
  }
  
  return List::create(
    _["residual_matrix"] = resi_mat,
    _["beta_matrix"] = beta_mat,
    _["sigma2_seq"] = sigma2_seq,
    _["beta"] = beta,
    _["lambda"] = lambda_value,
    _["lambda_for_soft_thres"] = lambda_effective,
    _["sorted_predictors"] = active_indices,
    _["no_of_iter"] = iteration,
    _["no_of_beta_iter"] = beta_iteration,
    _["count_sig_beta"] = vec_sig_beta_count
  );
}


