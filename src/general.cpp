#include <Rcpp.h>
using namespace Rcpp;

NumericVector marg_mean_tv_cpp(double lambda1, NumericVector nu, NumericVector phi, NumericVector kappa){
  int lgt = nu.size();
  NumericVector mu(lgt);
  mu(0) = lambda1;
  for(int i = 1; i < lgt; i++){
    mu(i) = nu(i) + (phi(i) + kappa(i))*mu(i - 1);
  }
  return mu;
}


// [[Rcpp::export]]
List compute_sop_tv_cpp(double lambda1, NumericVector nu,
                        NumericVector phi, NumericVector kappa, NumericVector psi, double q){
  int lgt = nu.size();
  NumericVector mu_X(lgt);
  mu_X(0) = lambda1;
  NumericVector v_lambda(lgt);
  v_lambda(0) = 0;
  NumericVector v_X(lgt);
  v_X(0) = mu_X(0) + psi(0)*pow(mu_X(0), 2);
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
  NumericVector mu_X_tilde = q*mu_X;
  NumericVector v_X_tilde = pow(q, 2)*v_X + q*(1 - q)*mu_X;
  NumericVector cov1_X_tilde = pow(q, 2)*cov1_X;

  return Rcpp::List::create(Rcpp::Named("mu_X") = mu_X, Rcpp::Named("v_X") = v_X,
                            Rcpp::Named("cov1_X") = cov1_X, Rcpp::Named("decay_cov_X") = phi + kappa,
                            Rcpp::Named("mu_X_tilde") = mu_X_tilde, Rcpp::Named("v_X_tilde") = v_X_tilde,
                            Rcpp::Named("cov1_X_tilde") = cov1_X_tilde, Rcpp::Named("decay_cov_X_tilde") = phi + kappa);
}


// [[Rcpp::export]]
List reparam_tv_cpp(double lambda1, NumericVector nu, NumericVector phi, NumericVector kappa,
                    NumericVector psi, double q){
  int lgt = nu.size();
  // compute target second-order properties:
  List target_sop = compute_sop_tv_cpp(lambda1, nu, phi, kappa, psi, q);

  // now find a completely observed process with the same properties:
  NumericVector nu_Y = q*nu; // known for theoretical reasons
  NumericVector phi_plus_kappa_Y = phi + kappa; // this, too.

  // extract:
  NumericVector mu_Y = target_sop("mu_X_tilde");
  NumericVector v_Y = target_sop("v_X_tilde"); // by definition
  NumericVector cov1_Y = target_sop("cov1_X_tilde"); // by definition

  // initialize:
  NumericVector phi_Y(lgt);
  NumericVector kappa_Y(lgt);
  NumericVector psi_Y(lgt);
  NumericVector v_lambda_Y(lgt);

  // first elements are known:
  v_lambda_Y(0) = 0;
  psi_Y(0) = psi(0);

  // element [1, 1] is correct by construction (thinning property of NB), the others can be adapted step by step.
  // done in loop:
  for(int i = 1; i < lgt; i++){
    // fill in phi_Y[i]
    phi_Y(i) = (cov1_Y(i) - phi_plus_kappa_Y(i)*v_lambda_Y(i - 1) -
      phi_plus_kappa_Y(i)*pow(mu_Y(i - 1),2) -
      nu_Y(i)*mu_Y(i - 1) + mu_Y(i)*mu_Y(i - 1))/
        (v_Y(i - 1) - v_lambda_Y(i - 1));
    kappa_Y(i) = phi_plus_kappa_Y(i) - phi_Y(i);
    // update v_lambda_Y:
    v_lambda_Y(i) = pow(phi_Y(i), 2)*v_Y(i - 1) +
      (2*phi_Y(i)*kappa_Y(i) + pow(kappa_Y(i), 2))*v_lambda_Y(i - 1);
    // now correct psi_Y(i)
    psi_Y(i) = (v_Y(i) - v_lambda_Y(i) - mu_Y(i))/
      (pow(mu_Y(i), 2) + v_lambda_Y(i));
  }
  return Rcpp::List::create(Rcpp::Named("lambda1") = mu_Y(0),
                            Rcpp::Named("nu") = nu_Y, Rcpp::Named("phi") = phi_Y, Rcpp::Named("kappa") = kappa_Y,
                                        Rcpp::Named("psi") = psi_Y, Rcpp::Named("q") = 1);
}


// [[Rcpp::export]]
double nllik_tv_cpp(NumericVector observed, double lambda1, NumericVector nu, NumericVector phi,
                    NumericVector kappa, NumericVector psi, double q){
  int lgt = observed.size();
  // get corresponding parameters for unthinned process:
  List pars_Y = reparam_tv_cpp(lambda1, nu, phi, kappa, psi, q);
  // extract:
  double lambda1_Y = pars_Y("lambda1");
  NumericVector nu_Y = pars_Y("nu");
  NumericVector phi_Y = pars_Y("phi");
  NumericVector kappa_Y = pars_Y("kappa");
  NumericVector psi_Y = pars_Y("psi");

  // get conditional means:
  NumericVector lambda_Y(lgt);
  NumericVector nlliks(lgt);

  lambda_Y(0) = lambda1_Y; // indices shifted
  nlliks(0) = -1*dnbinom_mu(observed(0), 1/psi_Y(0), lambda_Y(0), true);

  for(int i = 1; i < lgt; i++){
    lambda_Y(i) = nu_Y(i) + phi_Y(i)*observed(i - 1) + kappa_Y(i)*lambda_Y(i - 1);
    nlliks(i) = -1*dnbinom_mu(observed(i), 1/psi_Y(i), lambda_Y(i), true);
  }

  return sum(nlliks);
}
