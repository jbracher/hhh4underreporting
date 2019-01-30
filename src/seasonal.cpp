#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector back_seas_cpp(NumericVector vect, int to, int lgt) {
  int L = vect.size();
  NumericVector vect_long = rep(vect, lgt/L + 2);
  NumericVector vect_back;
  if(to == L){
    vect_back = rev(tail(vect_long, lgt));
  }else{
    vect_back = rev(head(tail(vect_long, L - to + lgt), lgt));
  }
  return vect_back;
}

// [[Rcpp::export]]
NumericMatrix get_weight_matrix_seas_cpp(NumericVector phi, NumericVector kappa, int max_lag){
  int L = phi.size();
  NumericVector phi_long = rep(phi, max_lag/L + 3);
  NumericVector kappa_long = rep(kappa, max_lag/L + 3);
  int shift = L*(max_lag/L + 1);
  NumericMatrix wgts(L, max_lag);
  for(int i = 0; i < L; i++){// fill in phi values
    for(int j = 0; j < max_lag; j++){
      wgts(i, j) = phi_long(phi_long.size() - shift + i - j);
    }
    for(int k = 1; k < max_lag; k++){
      for(int l = k; l < max_lag; l++){
        wgts(i, l) = wgts(i, l)*kappa_long(kappa_long.size() - shift + i - k - L + 1);
      }
    }
  }
  return wgts;
}

// [[Rcpp::export]]
double prod(NumericVector vect){
  double pr = vect(0);
  for(int i = 1; i < vect.size(); i++){
    pr = pr*vect(i);
  }
  return pr;
}

// [[Rcpp::export]]
NumericVector order_back_in_time(NumericVector vect, int from){
  NumericVector vect_reordered(vect.size());
  for(int i = 0; i < vect.size() - from; i++){
    vect_reordered(i) = vect(from + i);
  }
  for(int i = 0; i < from; i++){
    vect_reordered(vect.size() - from + i) = vect(i);
  }
  return vect_reordered;
}

// [[Rcpp::export]]
NumericVector stat_mean_seas_cpp(NumericVector nu, NumericVector phi, NumericVector kappa){
  int L = nu.size();
  NumericVector mu(L);
  for(int i = L; i > 1; i--){
    mu(L - 1) = (phi(L - i + 1) + kappa(L - i + 1))*(mu(L - 1) + nu(L - i));
  }
  mu[L - 1] = (mu(L - 1) + nu(L - 1))/(1 - prod(phi + kappa));
  mu(0) = nu(0) + (phi(0) + kappa(0))*mu(L - 1);
  for(int i = 1; i < L - 1; i++){
    mu(i) = nu(i) + (phi(i) + kappa(i))*mu(i - 1);
  }
  return mu;
}

// [[Rcpp::export]]
List stat_var_seas_cpp(NumericVector nu, NumericVector phi, NumericVector kappa, NumericVector psi){
  int L = nu.size();
  NumericVector mu = stat_mean_seas_cpp(nu, phi, kappa);
  NumericVector v = psi*pow(mu, 2) + mu;

  NumericVector h(L);
  h(0) = pow(phi(0) + kappa(0), 2) + pow(phi(0), 2)*psi(L - 1);
  for(int i = 1; i < L; i++){
    h(i) = pow(phi(i) + kappa(i), 2) + pow(phi(i), 2)*psi(i - 1);
  }

  NumericVector wgts(L);
  wgts(L - 1) = 1;
  for(int i = 2; i <= L; i++){
    wgts(L - i) = h(L - i + 1)*wgts(L - i + 1);
  }

  NumericVector summands(L);
  summands(0) = pow(phi(0), 2)*v(L - 1);
  for(int i = 1; i < L; i++){
    summands(i) = pow(phi(i), 2)*v(i - 1);
  }

  NumericVector stat_var_lambda(L);
  stat_var_lambda(L - 1) = sum(wgts*summands)/(1 - prod(h));
  stat_var_lambda(0) = pow(phi(0), 2)*v(L - 1) + (pow(phi(0) + kappa(0), 2) + pow(phi(0), 2)*psi(L - 1))*stat_var_lambda(L - 1);
  for(int i = 1; i < L - 1; i++){
    stat_var_lambda(i) = pow(phi(i), 2)*v(i - 1) + (pow(phi(i) + kappa(i), 2) + pow(phi(i), 2)*psi(i - 1))*stat_var_lambda(i - 1);
  }
  NumericVector stat_var_X = (1 + psi)*stat_var_lambda + v;
  return Rcpp::List::create(Rcpp::Named("stat_var_lambda") = stat_var_lambda,
                            Rcpp::Named("stat_var_X") = stat_var_X);
}

// [[Rcpp::export]]
NumericVector nu_to_nu_tilde_seas_cpp(NumericVector nu, NumericVector kappa, int max_lag){
  return stat_mean_seas_cpp(nu, NumericVector(nu.size()), kappa);
}

// [[Rcpp::export]]
List compute_sop_seas_cpp(NumericVector nu, NumericVector phi, NumericVector kappa, NumericVector psi, double p){
  int L = nu.size();
  NumericVector mu_X = stat_mean_seas_cpp(nu, phi, kappa);
  NumericVector sigma2_X = stat_var_seas_cpp(nu, phi, kappa, psi)(1);
  NumericVector cov1_X(L);
  cov1_X(0) = phi(0)*sigma2_X(L - 1) +
    kappa(0)*(sigma2_X(L - 1) - mu_X(L - 1) - psi(L - 1)*pow(mu_X(L - 1), 2))/
      (1 + psi(L - 1));
  for(int i = 1; i < L; i++){
    cov1_X(i) = phi(i)*sigma2_X(i - 1) +
      kappa(i)*(sigma2_X(i - 1) - mu_X(i - 1) - psi(i - 1)*pow(mu_X(i - 1), 2))/
        (1 + psi(i - 1));
  }
  NumericVector decay_X = phi + kappa;

  NumericVector mu_Y = p*mu_X;
  NumericVector sigma2_Y = pow(p, 2)*sigma2_X + p*(1 - p)*mu_X;
  NumericVector cov1_Y = pow(p, 2)*cov1_X;
  NumericVector decay_Y = decay_X;

  return Rcpp::List::create(Rcpp::Named("mu_Y") = mu_Y,
                            Rcpp::Named("v_Y") = sigma2_Y,
                            Rcpp::Named("cov1_Y") = cov1_Y,
                            Rcpp::Named("decay_cov_Y") = decay_Y);
}

// [[Rcpp::export]]
List reparam_seas_cpp(NumericVector nu, NumericVector phi, NumericVector kappa, NumericVector psi, double q){
  int L = nu.size();

  if(q == 1){
    return Rcpp::List::create(Rcpp::Named("nu_star") = nu,
                              Rcpp::Named("phi_star") = phi,
                              Rcpp::Named("kappa_star") = kappa,
                              Rcpp::Named("psi_star") = psi);
  }

  // compute the second-order properties:
  List target_sop = compute_sop_seas_cpp(nu, phi, kappa, psi, q);

  // now find a completely observed process with the same properties:
  NumericVector nu_star = q*nu; // known for theoretical reasons
  NumericVector phi_plus_kappa_star = phi + kappa;
  NumericVector mu_X_star = target_sop("mu_Y"); // by definition
  NumericVector v_X_star = target_sop("v_Y"); // by definition
  NumericVector cov1_X_star = target_sop("cov1_Y"); // by definition

  // for the remaining one the old values are ususally good starting values:
  NumericVector phi_star = phi;
  NumericVector kappa_star = kappa;
  NumericVector psi_star = psi;

  // given a starting value of psi[L] we can get var_lambda_star[L]
  NumericVector v_lambda_star(L);
  v_lambda_star(L - 1) = (v_X_star(L - 1) - mu_X_star(L - 1) - psi_star(L - 1)*pow(mu_X_star(L - 1), 2))/
    (1 + psi_star(L - 1));

  // now we can iteratively roll through the weeks and correct the values of phi and kappa
  for(int k = 1; k<= 3; k++){
    for(int i = 0; i < L; i++){
      int im1 = i - 1;
      if(i == 0){
        im1 = L - 1;
      }
      // correct phi_star[i]
      phi_star(i) = (cov1_X_star(i) - phi_plus_kappa_star(i)*v_lambda_star(im1) -
        phi_plus_kappa_star(i)*pow(mu_X_star(im1), 2) -
        nu_star(i)*mu_X_star(im1) + mu_X_star(i)*mu_X_star(im1))/
          (v_X_star(im1) - v_lambda_star(im1));
      kappa_star(i) = phi_plus_kappa_star(i) - phi_star(i);
      // update v_lambda_star:
      v_lambda_star(i) = pow(phi_star(i), 2)*v_X_star(im1) +
        (2*phi_star(i)*kappa_star(i) + pow(kappa_star(i), 2))*v_lambda_star(im1);
      // now correct psi_star[i]
      psi_star(i) = (v_X_star(i) - v_lambda_star(i) - mu_X_star(i))/
        (pow(mu_X_star(i), 2) + v_lambda_star(i));
    }
  }
  return Rcpp::List::create(Rcpp::Named("nu") = nu_star,
                            Rcpp::Named("phi") = phi_star,
                            Rcpp::Named("kappa") = kappa_star,
                            Rcpp::Named("psi") = psi_star,
                            Rcpp::Named("q") = 1);
}

// [[Rcpp::export]]
NumericMatrix get_mod_matr_cpp(NumericVector Y, int max_lag){
  int L = Y.size();
  NumericMatrix mod_matr(L, max_lag);
  for(int ro = 0; ro < L; ro++){
    for(int co = 0; co < max_lag; co++){
      if(ro - co >= 1){
        mod_matr(ro, co) = Y(ro - co - 1);
      }else{
        mod_matr(ro, co) = NA_REAL;
      }
    }
  }
  return mod_matr;
}
