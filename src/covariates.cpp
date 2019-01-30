#include <Rcpp.h>
using namespace Rcpp;

NumericVector cond_mean_cov_cpp(double m1, NumericVector nu, NumericVector phi, NumericVector kappa){
  int lgt = nu.size();
  NumericVector mu_X(lgt);
  mu_X(0) = m1;
  for(int i = 1; i < lgt; i++){
    mu_X(i) = nu(i) + (phi(i) + kappa(i))*mu_X(i - 1);
  }
  return mu_X;
}

// [[Rcpp::export]]
List compute_sop_cov_cpp(double m1, double vl1, NumericVector nu,
                             NumericVector phi, NumericVector kappa, NumericVector psi, double q){
  int lgt = nu.size();
  NumericVector mu_X(lgt);
  mu_X(0) = m1;
  NumericVector v_lambda(lgt);
  v_lambda(0) = vl1;
  NumericVector v_X(lgt);
  v_X(0) = (1 + psi(0))*vl1 + mu_X(0) + psi(0)*pow(mu_X(0), 2);
  NumericVector cov1_X(lgt);
  for(int i = 1; i < lgt; i++){
    mu_X(i) = nu(i) + (phi(i) + kappa(i))*mu_X(i - 1);
    v_lambda(i) = pow(phi(i), 2)*v_X(i - 1) + (2*phi(i)*kappa(i) + pow(kappa(i), 2))*v_lambda(i - 1);
    v_X(i) = (1 + psi(i))*v_lambda(i) + mu_X(i) + psi(i)*pow(mu_X(i), 2);
    cov1_X(i) = kappa(i)*v_lambda(i - 1) +
      phi(i)*v_X(i - 1) +
      (phi(i) + kappa(i))*pow(mu_X(i - 1), 2) +
      nu(i)*mu_X(i - 1) -
      mu_X(i - 1)*mu_X(i);
  }

  // thinning:
  NumericVector mu_Y = q*mu_X;
  NumericVector v_Y = pow(q, 2)*v_X + q*(1 - q)*mu_X;
  NumericVector cov1_Y = pow(q, 2)*cov1_X;

  return Rcpp::List::create(Rcpp::Named("mu_X") = mu_X, Rcpp::Named("v_X") = v_X,
                            Rcpp::Named("cov1_X") = cov1_X, Rcpp::Named("decay_cov_X") = phi + kappa,
                            Rcpp::Named("mu_Y") = mu_Y, Rcpp::Named("v_Y") = v_Y,
                            Rcpp::Named("cov1_Y") = cov1_Y, Rcpp::Named("decay_cov_Y") = phi + kappa);
}


// [[Rcpp::export]]
List reparam_cov_cpp(double m1, double vl1, NumericVector nu, NumericVector phi, NumericVector kappa,
                   NumericVector psi, double q){
  int lgt = nu.size();
  // compute target second-order properties:
  List target_sop = compute_sop_cov_cpp(m1, vl1, nu, phi, kappa, psi, q);

  // now find a completely observed process with the same properties:
  NumericVector nu_star = q*nu; // known for theoretical reasons (?)
  NumericVector phi_plus_kappa_star = phi + kappa; // this, too.

  // extract:
  NumericVector mu_X_star = target_sop("mu_Y");
  NumericVector v_X_star = target_sop("v_Y"); // by definition
  NumericVector cov1_X_star = target_sop("cov1_Y"); // by definition

  // initialize:
  NumericVector phi_star(lgt);
  NumericVector kappa_star(lgt);
  NumericVector psi_star(lgt);
  NumericVector v_lambda_star(lgt);

  // first elements are known:
  v_lambda_star(0) = pow(q, 2)*vl1;
  psi_star(0) = psi(0);

  // element [1, 1] is correct by construction (thinning property of NB), the others can be adapted step by step.
  // done in loop:
  for(int i = 1; i < lgt; i++){
    // fill in phi_star[i]
    phi_star(i) = (cov1_X_star(i) - phi_plus_kappa_star(i)*v_lambda_star(i - 1) -
      phi_plus_kappa_star(i)*pow(mu_X_star(i - 1),2) -
      nu_star(i)*mu_X_star(i - 1) + mu_X_star(i)*mu_X_star(i - 1))/
        (v_X_star(i - 1) - v_lambda_star(i - 1));
    kappa_star(i) = phi_plus_kappa_star(i) - phi_star(i);
    // update v_lambda_star:
    v_lambda_star(i) = pow(phi_star(i), 2)*v_X_star(i - 1) +
      (2*phi_star(i)*kappa_star(i) + pow(kappa_star(i), 2))*v_lambda_star(i - 1);
    // now correct psi_star(i)
    psi_star(i) = (v_X_star(i) - v_lambda_star(i) - mu_X_star(i))/
      (pow(mu_X_star(i), 2) + v_lambda_star(i));
  }
  return Rcpp::List::create(Rcpp::Named("m1") = mu_X_star(0), Rcpp::Named("vl1") = v_lambda_star(0),
                            Rcpp::Named("nu") = nu_star, Rcpp::Named("phi") = phi_star, Rcpp::Named("kappa") = kappa_star,
                            Rcpp::Named("psi") = psi_star, Rcpp::Named("q") = 1);
}


// [[Rcpp::export]]
NumericMatrix get_weight_matrix_cov_cpp(NumericVector phi, NumericVector kappa, int max_lag){
  int lgt = phi.size();

  NumericMatrix wgts(lgt, max_lag);
  for(int i = max_lag; i < lgt; i++){
    for(int j = 0; j < max_lag; j++){
      wgts(i, j) = phi(i - j);
    }
    for(int k = 1; k < max_lag; k++){
      for(int l = k; l < max_lag; l++){
        wgts(i, l) = wgts(i, l)*kappa(i - k + 1);
      }
    }
  }

  return wgts;
}

// [[Rcpp::export]]
NumericVector nu_to_nu_tilde_cov_cpp(NumericVector nu, NumericVector kappa){
  return cond_mean_cov_cpp(nu(0), nu, NumericVector(nu.size()), kappa);
}//!!! this does not use max_lag now; takes lags as far as it can get them (approximately); should be good enough.

// a simple copy...
NumericMatrix get_mod_matr_cov_cpp(NumericVector Y, int max_lag){
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

// [[Rcpp::export]]
double nllik_cov_cpp(NumericVector Y, double m1, double vl1, NumericVector nu, NumericVector phi,
                 NumericVector kappa, NumericVector psi, double q, int max_lag = 5){
  int lgt = Y.size();
  // get model matrix:
  NumericMatrix mod_matr = get_mod_matr_cov_cpp(Y, max_lag);
  // get corresponding parameters for unthinned process:
  List pars_star = reparam_cov_cpp(m1, vl1, nu, phi, kappa, psi, q);
  // extract:
  double m1_star = pars_star("m1");
  double vl1_star = pars_star("vl1");
  NumericVector nu_star = pars_star("nu");
  NumericVector phi_star = pars_star("phi");
  NumericVector kappa_star = pars_star("kappa");
  NumericVector psi_star = pars_star("psi");

  // get nu_tilde, i.e. nu for geometric lag formulation
  NumericVector nu_star_tilde  = nu_to_nu_tilde_cov_cpp(nu_star, kappa_star);

  // get weight matrix for lags:
  NumericMatrix weight_matrix = get_weight_matrix_cov_cpp(phi_star, kappa_star, max_lag);

  // get conditional means:
  NumericVector lambda(lgt);
  NumericVector nlliks(lgt);
  for(int i = 0; i < max_lag; i++){
    nlliks(i) = 0;
  }
  for(int i = max_lag; i < lgt; i++){
    lambda(i) = nu_star_tilde(i) + sum(weight_matrix(i,_)*mod_matr(i,_));
    nlliks(i) = -1*dnbinom_mu(Y(i), 1/psi_star(i), lambda(i), true);
   }

  return sum(nlliks);
}

