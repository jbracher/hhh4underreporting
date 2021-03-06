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
                        NumericVector phi, NumericVector kappa,
                        NumericVector psi, NumericVector q){
  int lgt = nu.size();
  NumericVector mu_X(lgt);
  mu_X(0) = lambda1;
  NumericVector v_lambda(lgt);
  v_lambda(0) = 0;
  NumericVector v_X(lgt);
  v_X(0) = mu_X(0) + psi(0)*pow(mu_X(0), 2);
  NumericVector cov1_X(lgt);
  NumericVector decay_cov_X(lgt);
  for(int i = 1; i < lgt; i++){
    mu_X(i) = nu(i) + (phi(i) + kappa(i))*mu_X(i - 1);
    v_lambda(i) = pow(phi(i), 2)*v_X(i - 1) + (2*phi(i)*kappa(i) + pow(kappa(i), 2))*v_lambda(i - 1);
    v_X(i) = (1 + psi(i))*v_lambda(i) + mu_X(i) + psi(i)*pow(mu_X(i), 2);
    cov1_X(i) = kappa(i)*v_lambda(i - 1) +
      phi(i)*v_X(i - 1);
      // (phi(i) + kappa(i))*pow(mu_X(i - 1), 2) +
      // nu(i)*mu_X(i - 1) -
      // mu_X(i - 1)*mu_X(i); // superfluous, sum up to 0
    decay_cov_X(i) = phi(i) + kappa(i);
  }

  // thinning:
  NumericVector mu_X_tilde = q*mu_X;
  NumericVector v_X_tilde = pow(q, 2)*v_X + q*(1 - q)*mu_X;
  NumericVector cov1_X_tilde(lgt);
  NumericVector decay_cov_X_tilde(lgt);
  for(int i = 1; i < lgt; i++){
    cov1_X_tilde(i) = q(i)*q(i - 1)*cov1_X(i);
    decay_cov_X_tilde(i) = q(i)/q(i - 1)*decay_cov_X(i);
  }

  return Rcpp::List::create(Rcpp::Named("mu_X") = mu_X, Rcpp::Named("v_X") = v_X,
                            Rcpp::Named("cov1_X") = cov1_X, Rcpp::Named("decay_cov_X") = decay_cov_X,
                            Rcpp::Named("mu_X_tilde") = mu_X_tilde, Rcpp::Named("v_X_tilde") = v_X_tilde,
                            Rcpp::Named("cov1_X_tilde") = cov1_X_tilde, Rcpp::Named("decay_cov_X_tilde") = decay_cov_X_tilde);
}

/* old version without decoarsening option
// [[Rcpp::export]]
List reparam_tv_cpp(double lambda1, NumericVector nu, NumericVector phi, NumericVector kappa,
                    NumericVector psi, NumericVector q){
  int lgt = nu.size();
  // compute target second-order properties:
  List target_sop = compute_sop_tv_cpp(lambda1, nu, phi, kappa, psi, q);

  // now find a completely observed process with the same properties:
  NumericVector nu_Y = q*nu; // known for theoretical reasons

  // extract:
  NumericVector mu_Y = target_sop("mu_X_tilde");
  NumericVector v_Y = target_sop("v_X_tilde");
  NumericVector cov1_Y = target_sop("cov1_X_tilde");
  NumericVector decay_cov_Y = target_sop("decay_cov_X_tilde");

  // initialize:
  NumericVector phi_Y(lgt);
  NumericVector kappa_Y(lgt);
  NumericVector psi_Y(lgt);
  NumericVector v_lambda_Y(lgt);
  NumericVector q_Y(10, 1.0);

  // first elements are known:
  v_lambda_Y(0) = 0;
  psi_Y(0) = psi(0);

  // element [1, 1] is correct by construction (thinning property of NB), the others can be adapted step by step.
  // done in loop:
  for(int i = 1; i < lgt; i++){
    // fill in phi_Y[i]
    phi_Y(i) = (cov1_Y(i) - decay_cov_Y(i)*v_lambda_Y(i - 1))/
        (v_Y(i - 1) - v_lambda_Y(i - 1));
    kappa_Y(i) = decay_cov_Y(i) - phi_Y(i);
    // update v_lambda_Y:
    v_lambda_Y(i) = pow(phi_Y(i), 2)*v_Y(i - 1) +
      (2*phi_Y(i)*kappa_Y(i) + pow(kappa_Y(i), 2))*v_lambda_Y(i - 1);
    // now correct psi_Y(i)
    psi_Y(i) = (v_Y(i) - v_lambda_Y(i) - mu_Y(i))/
      (pow(mu_Y(i), 2) + v_lambda_Y(i));
  }
  return Rcpp::List::create(Rcpp::Named("lambda1") = mu_Y(0),
                            Rcpp::Named("nu") = nu_Y, Rcpp::Named("phi") = phi_Y, Rcpp::Named("kappa") = kappa_Y,
                                        Rcpp::Named("psi") = psi_Y, Rcpp::Named("q") = q_Y);
}
*/

// [[Rcpp::export]]
NumericVector coarsen_vector_sum_cpp(NumericVector vect){
  int lgt = vect.size()/2;
  NumericVector vect_coarse(lgt);
  for(int i = 0; i < lgt; i++){
    vect_coarse(i) = vect(2*i) + vect(2*i + 1);
  }
  return(vect_coarse);
}

// [[Rcpp::export]]
NumericVector coarsen_vector_prod_cpp(NumericVector vect){
  int lgt = vect.size()/2;
  NumericVector vect_coarse(lgt);
  for(int i = 0; i < lgt; i++){
    vect_coarse(i) = vect(2*i) * vect(2*i + 1);
  }
  return(vect_coarse);
}

// [[Rcpp::export]]
NumericVector coarsen_vector_only_first_cpp(NumericVector vect){
  int lgt = vect.size()/2;
  NumericVector vect_coarse(lgt);
  for(int i = 0; i < lgt; i++){
    vect_coarse(i) = vect(2*i);
  }
  return(vect_coarse);
}

// [[Rcpp::export]]
NumericVector coarsen_vector_only_sec_cpp(NumericVector vect){
  int lgt = vect.size()/2;
  NumericVector vect_coarse(lgt);
  for(int i = 0; i < lgt; i++){
    vect_coarse(i) = vect(2*i + 1);
  }
  return(vect_coarse);
}


// [[Rcpp::export]]
NumericVector shift_add_zero(NumericVector vect){
  int lgt = vect.size();
  NumericVector shifted_vect(lgt);
  for(int i = 1; i < lgt; i++){
    shifted_vect(i) = vect(i - 1);
  }
  return(shifted_vect);
}


// [[Rcpp::export]]
List coarsen_sop_tv_cpp(List sop){

  // coarsen means:
  NumericVector mu_X = coarsen_vector_sum_cpp(sop["mu_X"]);
  NumericVector mu_X_tilde = coarsen_vector_sum_cpp(sop["mu_X_tilde"]);

  // coarsen variances:
  NumericVector v_X = coarsen_vector_sum_cpp(sop["v_X"]) +
    2 * coarsen_vector_only_sec_cpp(sop["cov1_X"]);
  NumericVector v_X_tilde = coarsen_vector_sum_cpp(sop["v_X_tilde"]) +
    2 * coarsen_vector_only_sec_cpp(sop["cov1_X_tilde"]);

  // coarsen first-order auocovariances:
  NumericVector cov13_X = shift_add_zero(coarsen_vector_only_sec_cpp(sop["cov1_X"])) *
    coarsen_vector_only_first_cpp(sop["decay_cov_X"]); // first with third
  NumericVector cov23_X = coarsen_vector_only_first_cpp(sop["cov1_X"]); // second with third
  NumericVector cov24_X = coarsen_vector_only_first_cpp(sop["cov1_X"])*
    coarsen_vector_only_sec_cpp(sop["decay_cov_X"]); // second with fourth
  NumericVector cov14_X = shift_add_zero(coarsen_vector_only_sec_cpp(sop["cov1_X"]))*
    coarsen_vector_prod_cpp(sop["decay_cov_X"]); // first with fourth
  NumericVector cov1_X = cov13_X + cov14_X + cov23_X + cov24_X;

  NumericVector cov13_X_tilde = shift_add_zero(coarsen_vector_only_sec_cpp(sop["cov1_X_tilde"])) *
    coarsen_vector_only_first_cpp(sop["decay_cov_X_tilde"]); // first with third
  NumericVector cov23_X_tilde = coarsen_vector_only_first_cpp(sop["cov1_X_tilde"]); // second with third
  NumericVector cov24_X_tilde = coarsen_vector_only_first_cpp(sop["cov1_X_tilde"])*
    coarsen_vector_only_sec_cpp(sop["decay_cov_X_tilde"]); // second with fourth
  NumericVector cov14_X_tilde = shift_add_zero(coarsen_vector_only_sec_cpp(sop["cov1_X_tilde"]))*
    coarsen_vector_prod_cpp(sop["decay_cov_X_tilde"]); // first with fourth
  NumericVector cov1_X_tilde = cov13_X_tilde + cov14_X_tilde + cov23_X_tilde + cov24_X_tilde;

  // coarsen decay of autocovariance function:
  NumericVector decay_cov_X = coarsen_vector_prod_cpp(sop["decay_cov_X"]);
  NumericVector decay_cov_X_tilde = coarsen_vector_prod_cpp(sop["decay_cov_X_tilde"]);

  return Rcpp::List::create(Rcpp::Named("mu_X") = mu_X,
                            Rcpp::Named("v_X") = v_X,
                            Rcpp::Named("cov1_X") = cov1_X,
                            Rcpp::Named("decay_cov_X") = decay_cov_X,
                            Rcpp::Named("mu_X_tilde") = mu_X_tilde,
                            Rcpp::Named("v_X_tilde") = v_X_tilde,
                            Rcpp::Named("cov1_X_tilde") = cov1_X_tilde,
                            Rcpp::Named("decay_cov_X_tilde") = decay_cov_X_tilde);
}

// [[Rcpp::export]]
List reparam_tv_cpp(double lambda1, NumericVector nu, NumericVector phi, NumericVector kappa,
                    NumericVector psi, NumericVector q, bool decoarsen){
  int lgt_original = nu.size();
  int lgt_reparam = lgt_original;
  if(decoarsen){
    lgt_reparam = lgt_original/2;
  }

  // compute target second-order properties:
  List target_sop0 = compute_sop_tv_cpp(lambda1, nu, phi, kappa, psi, q);
  List target_sop = target_sop0;
  if(decoarsen){
    target_sop = coarsen_sop_tv_cpp(target_sop0);
  }

  // now find a completely observed process with the same properties:

  // extract:
  NumericVector mu_Y = target_sop("mu_X_tilde");
  NumericVector v_Y = target_sop("v_X_tilde");
  NumericVector cov1_Y = target_sop("cov1_X_tilde");
  NumericVector decay_cov_Y = target_sop("decay_cov_X_tilde");

  // initialize:
  NumericVector nu_Y(lgt_reparam);
  NumericVector phi_Y(lgt_reparam);
  NumericVector kappa_Y(lgt_reparam);
  NumericVector psi_Y(lgt_reparam);
  NumericVector v_lambda_Y(lgt_reparam);
  NumericVector q_Y(lgt_reparam, 1.0);

  // first elements are known:
  v_lambda_Y(0) = 0;
  psi_Y (0) = (v_Y(0) - v_lambda_Y(0) - mu_Y(0))/
    (pow(mu_Y(0), 2) + v_lambda_Y(0));

  // element [1, 1] is correct by construction (thinning property of NB), the others can be adapted step by step.
  // done in loop:
  for(int i = 1; i < lgt_reparam; i++){
    // correct nu_Y[i]
    nu_Y(i) = mu_Y(i) - decay_cov_Y(i)*mu_Y(i - 1);
    // fill in phi_Y[i]
    phi_Y(i) = (cov1_Y(i) - decay_cov_Y(i)*v_lambda_Y(i - 1))/
      (v_Y(i - 1) - v_lambda_Y(i - 1));
    kappa_Y(i) = decay_cov_Y(i) - phi_Y(i);
    // update v_lambda_Y:
    v_lambda_Y(i) = pow(phi_Y(i), 2)*v_Y(i - 1) +
      (2*phi_Y(i)*kappa_Y(i) + pow(kappa_Y(i), 2))*v_lambda_Y(i - 1);
    // now correct psi_Y(i)
    psi_Y(i) = (v_Y(i) - v_lambda_Y(i) - mu_Y(i))/
      (pow(mu_Y(i), 2) + v_lambda_Y(i));
  }
  return Rcpp::List::create(Rcpp::Named("lambda1") = mu_Y(0),
                            Rcpp::Named("nu") = nu_Y, Rcpp::Named("phi") = phi_Y, Rcpp::Named("kappa") = kappa_Y,
                                        Rcpp::Named("psi") = psi_Y, Rcpp::Named("q") = q_Y);
}



// [[Rcpp::export]]
double nllik_tv_cpp(NumericVector observed, double lambda1, NumericVector nu, NumericVector phi,
                    NumericVector kappa, NumericVector psi, NumericVector q, bool decoarsen){
  int lgt = observed.size();
  // get corresponding parameters for unthinned process:
  List pars_Y = reparam_tv_cpp(lambda1, nu, phi, kappa, psi, q, decoarsen);
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

  return(sum(nlliks));
}
