// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// compute_sop_cov_cpp
List compute_sop_cov_cpp(double m1, double vl1, NumericVector nu, NumericVector phi, NumericVector kappa, NumericVector psi, double q);
RcppExport SEXP _hhh4underreporting_compute_sop_cov_cpp(SEXP m1SEXP, SEXP vl1SEXP, SEXP nuSEXP, SEXP phiSEXP, SEXP kappaSEXP, SEXP psiSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type m1(m1SEXP);
    Rcpp::traits::input_parameter< double >::type vl1(vl1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_sop_cov_cpp(m1, vl1, nu, phi, kappa, psi, q));
    return rcpp_result_gen;
END_RCPP
}
// reparam_cov_cpp
List reparam_cov_cpp(double m1, double vl1, NumericVector nu, NumericVector phi, NumericVector kappa, NumericVector psi, double q);
RcppExport SEXP _hhh4underreporting_reparam_cov_cpp(SEXP m1SEXP, SEXP vl1SEXP, SEXP nuSEXP, SEXP phiSEXP, SEXP kappaSEXP, SEXP psiSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type m1(m1SEXP);
    Rcpp::traits::input_parameter< double >::type vl1(vl1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(reparam_cov_cpp(m1, vl1, nu, phi, kappa, psi, q));
    return rcpp_result_gen;
END_RCPP
}
// get_weight_matrix_cov_cpp
NumericMatrix get_weight_matrix_cov_cpp(NumericVector phi, NumericVector kappa, int max_lag);
RcppExport SEXP _hhh4underreporting_get_weight_matrix_cov_cpp(SEXP phiSEXP, SEXP kappaSEXP, SEXP max_lagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< int >::type max_lag(max_lagSEXP);
    rcpp_result_gen = Rcpp::wrap(get_weight_matrix_cov_cpp(phi, kappa, max_lag));
    return rcpp_result_gen;
END_RCPP
}
// nu_to_nu_star_cov_cpp
NumericVector nu_to_nu_star_cov_cpp(NumericVector nu, NumericVector kappa);
RcppExport SEXP _hhh4underreporting_nu_to_nu_star_cov_cpp(SEXP nuSEXP, SEXP kappaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type kappa(kappaSEXP);
    rcpp_result_gen = Rcpp::wrap(nu_to_nu_star_cov_cpp(nu, kappa));
    return rcpp_result_gen;
END_RCPP
}
// nllik_cov_cpp
double nllik_cov_cpp(NumericVector Y, double m1, double vl1, NumericVector nu, NumericVector phi, NumericVector kappa, NumericVector psi, double q, int max_lag);
RcppExport SEXP _hhh4underreporting_nllik_cov_cpp(SEXP YSEXP, SEXP m1SEXP, SEXP vl1SEXP, SEXP nuSEXP, SEXP phiSEXP, SEXP kappaSEXP, SEXP psiSEXP, SEXP qSEXP, SEXP max_lagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< double >::type m1(m1SEXP);
    Rcpp::traits::input_parameter< double >::type vl1(vl1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type max_lag(max_lagSEXP);
    rcpp_result_gen = Rcpp::wrap(nllik_cov_cpp(Y, m1, vl1, nu, phi, kappa, psi, q, max_lag));
    return rcpp_result_gen;
END_RCPP
}
// compute_sop_tv_cpp
List compute_sop_tv_cpp(double lambda1, NumericVector nu, NumericVector phi, NumericVector kappa, NumericVector psi, NumericVector q);
RcppExport SEXP _hhh4underreporting_compute_sop_tv_cpp(SEXP lambda1SEXP, SEXP nuSEXP, SEXP phiSEXP, SEXP kappaSEXP, SEXP psiSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type lambda1(lambda1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_sop_tv_cpp(lambda1, nu, phi, kappa, psi, q));
    return rcpp_result_gen;
END_RCPP
}
// reparam_tv_cpp
List reparam_tv_cpp(double lambda1, NumericVector nu, NumericVector phi, NumericVector kappa, NumericVector psi, NumericVector q);
RcppExport SEXP _hhh4underreporting_reparam_tv_cpp(SEXP lambda1SEXP, SEXP nuSEXP, SEXP phiSEXP, SEXP kappaSEXP, SEXP psiSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type lambda1(lambda1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(reparam_tv_cpp(lambda1, nu, phi, kappa, psi, q));
    return rcpp_result_gen;
END_RCPP
}
// nllik_tv_cpp
double nllik_tv_cpp(NumericVector observed, double lambda1, NumericVector nu, NumericVector phi, NumericVector kappa, NumericVector psi, NumericVector q);
RcppExport SEXP _hhh4underreporting_nllik_tv_cpp(SEXP observedSEXP, SEXP lambda1SEXP, SEXP nuSEXP, SEXP phiSEXP, SEXP kappaSEXP, SEXP psiSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type observed(observedSEXP);
    Rcpp::traits::input_parameter< double >::type lambda1(lambda1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(nllik_tv_cpp(observed, lambda1, nu, phi, kappa, psi, q));
    return rcpp_result_gen;
END_RCPP
}
// coarsen_vector_sum_cpp
NumericVector coarsen_vector_sum_cpp(NumericVector vect);
RcppExport SEXP _hhh4underreporting_coarsen_vector_sum_cpp(SEXP vectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type vect(vectSEXP);
    rcpp_result_gen = Rcpp::wrap(coarsen_vector_sum_cpp(vect));
    return rcpp_result_gen;
END_RCPP
}
// coarsen_vector_prod_cpp
NumericVector coarsen_vector_prod_cpp(NumericVector vect);
RcppExport SEXP _hhh4underreporting_coarsen_vector_prod_cpp(SEXP vectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type vect(vectSEXP);
    rcpp_result_gen = Rcpp::wrap(coarsen_vector_prod_cpp(vect));
    return rcpp_result_gen;
END_RCPP
}
// coarsen_vector_only_first_cpp
NumericVector coarsen_vector_only_first_cpp(NumericVector vect);
RcppExport SEXP _hhh4underreporting_coarsen_vector_only_first_cpp(SEXP vectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type vect(vectSEXP);
    rcpp_result_gen = Rcpp::wrap(coarsen_vector_only_first_cpp(vect));
    return rcpp_result_gen;
END_RCPP
}
// coarsen_vector_only_sec_cpp
NumericVector coarsen_vector_only_sec_cpp(NumericVector vect);
RcppExport SEXP _hhh4underreporting_coarsen_vector_only_sec_cpp(SEXP vectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type vect(vectSEXP);
    rcpp_result_gen = Rcpp::wrap(coarsen_vector_only_sec_cpp(vect));
    return rcpp_result_gen;
END_RCPP
}
// shift_add_zero
NumericVector shift_add_zero(NumericVector vect);
RcppExport SEXP _hhh4underreporting_shift_add_zero(SEXP vectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type vect(vectSEXP);
    rcpp_result_gen = Rcpp::wrap(shift_add_zero(vect));
    return rcpp_result_gen;
END_RCPP
}
// coarsen_sop_tv_cpp
List coarsen_sop_tv_cpp(List sop);
RcppExport SEXP _hhh4underreporting_coarsen_sop_tv_cpp(SEXP sopSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type sop(sopSEXP);
    rcpp_result_gen = Rcpp::wrap(coarsen_sop_tv_cpp(sop));
    return rcpp_result_gen;
END_RCPP
}
// back_seas_cpp
NumericVector back_seas_cpp(NumericVector vect, int to, int lgt);
RcppExport SEXP _hhh4underreporting_back_seas_cpp(SEXP vectSEXP, SEXP toSEXP, SEXP lgtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type vect(vectSEXP);
    Rcpp::traits::input_parameter< int >::type to(toSEXP);
    Rcpp::traits::input_parameter< int >::type lgt(lgtSEXP);
    rcpp_result_gen = Rcpp::wrap(back_seas_cpp(vect, to, lgt));
    return rcpp_result_gen;
END_RCPP
}
// get_weight_matrix_seas_cpp
NumericMatrix get_weight_matrix_seas_cpp(NumericVector phi, NumericVector kappa, int max_lag);
RcppExport SEXP _hhh4underreporting_get_weight_matrix_seas_cpp(SEXP phiSEXP, SEXP kappaSEXP, SEXP max_lagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< int >::type max_lag(max_lagSEXP);
    rcpp_result_gen = Rcpp::wrap(get_weight_matrix_seas_cpp(phi, kappa, max_lag));
    return rcpp_result_gen;
END_RCPP
}
// prod
double prod(NumericVector vect);
RcppExport SEXP _hhh4underreporting_prod(SEXP vectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type vect(vectSEXP);
    rcpp_result_gen = Rcpp::wrap(prod(vect));
    return rcpp_result_gen;
END_RCPP
}
// order_back_in_time
NumericVector order_back_in_time(NumericVector vect, int from);
RcppExport SEXP _hhh4underreporting_order_back_in_time(SEXP vectSEXP, SEXP fromSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type vect(vectSEXP);
    Rcpp::traits::input_parameter< int >::type from(fromSEXP);
    rcpp_result_gen = Rcpp::wrap(order_back_in_time(vect, from));
    return rcpp_result_gen;
END_RCPP
}
// stat_mean_seas_cpp
NumericVector stat_mean_seas_cpp(NumericVector nu, NumericVector phi, NumericVector kappa);
RcppExport SEXP _hhh4underreporting_stat_mean_seas_cpp(SEXP nuSEXP, SEXP phiSEXP, SEXP kappaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type kappa(kappaSEXP);
    rcpp_result_gen = Rcpp::wrap(stat_mean_seas_cpp(nu, phi, kappa));
    return rcpp_result_gen;
END_RCPP
}
// stat_var_seas_cpp
List stat_var_seas_cpp(NumericVector nu, NumericVector phi, NumericVector kappa, NumericVector psi);
RcppExport SEXP _hhh4underreporting_stat_var_seas_cpp(SEXP nuSEXP, SEXP phiSEXP, SEXP kappaSEXP, SEXP psiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type psi(psiSEXP);
    rcpp_result_gen = Rcpp::wrap(stat_var_seas_cpp(nu, phi, kappa, psi));
    return rcpp_result_gen;
END_RCPP
}
// nu_to_nu_star_seas_cpp
NumericVector nu_to_nu_star_seas_cpp(NumericVector nu, NumericVector kappa, int max_lag);
RcppExport SEXP _hhh4underreporting_nu_to_nu_star_seas_cpp(SEXP nuSEXP, SEXP kappaSEXP, SEXP max_lagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< int >::type max_lag(max_lagSEXP);
    rcpp_result_gen = Rcpp::wrap(nu_to_nu_star_seas_cpp(nu, kappa, max_lag));
    return rcpp_result_gen;
END_RCPP
}
// compute_sop_seas_cpp
List compute_sop_seas_cpp(NumericVector nu, NumericVector phi, NumericVector kappa, NumericVector psi, double q);
RcppExport SEXP _hhh4underreporting_compute_sop_seas_cpp(SEXP nuSEXP, SEXP phiSEXP, SEXP kappaSEXP, SEXP psiSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_sop_seas_cpp(nu, phi, kappa, psi, q));
    return rcpp_result_gen;
END_RCPP
}
// reparam_seas_cpp
List reparam_seas_cpp(NumericVector nu, NumericVector phi, NumericVector kappa, NumericVector psi, double q);
RcppExport SEXP _hhh4underreporting_reparam_seas_cpp(SEXP nuSEXP, SEXP phiSEXP, SEXP kappaSEXP, SEXP psiSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(reparam_seas_cpp(nu, phi, kappa, psi, q));
    return rcpp_result_gen;
END_RCPP
}
// get_mod_matr_cpp
NumericMatrix get_mod_matr_cpp(NumericVector observed, int max_lag);
RcppExport SEXP _hhh4underreporting_get_mod_matr_cpp(SEXP observedSEXP, SEXP max_lagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type observed(observedSEXP);
    Rcpp::traits::input_parameter< int >::type max_lag(max_lagSEXP);
    rcpp_result_gen = Rcpp::wrap(get_mod_matr_cpp(observed, max_lag));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_hhh4underreporting_compute_sop_cov_cpp", (DL_FUNC) &_hhh4underreporting_compute_sop_cov_cpp, 7},
    {"_hhh4underreporting_reparam_cov_cpp", (DL_FUNC) &_hhh4underreporting_reparam_cov_cpp, 7},
    {"_hhh4underreporting_get_weight_matrix_cov_cpp", (DL_FUNC) &_hhh4underreporting_get_weight_matrix_cov_cpp, 3},
    {"_hhh4underreporting_nu_to_nu_star_cov_cpp", (DL_FUNC) &_hhh4underreporting_nu_to_nu_star_cov_cpp, 2},
    {"_hhh4underreporting_nllik_cov_cpp", (DL_FUNC) &_hhh4underreporting_nllik_cov_cpp, 9},
    {"_hhh4underreporting_compute_sop_tv_cpp", (DL_FUNC) &_hhh4underreporting_compute_sop_tv_cpp, 6},
    {"_hhh4underreporting_reparam_tv_cpp", (DL_FUNC) &_hhh4underreporting_reparam_tv_cpp, 6},
    {"_hhh4underreporting_nllik_tv_cpp", (DL_FUNC) &_hhh4underreporting_nllik_tv_cpp, 7},
    {"_hhh4underreporting_coarsen_vector_sum_cpp", (DL_FUNC) &_hhh4underreporting_coarsen_vector_sum_cpp, 1},
    {"_hhh4underreporting_coarsen_vector_prod_cpp", (DL_FUNC) &_hhh4underreporting_coarsen_vector_prod_cpp, 1},
    {"_hhh4underreporting_coarsen_vector_only_first_cpp", (DL_FUNC) &_hhh4underreporting_coarsen_vector_only_first_cpp, 1},
    {"_hhh4underreporting_coarsen_vector_only_sec_cpp", (DL_FUNC) &_hhh4underreporting_coarsen_vector_only_sec_cpp, 1},
    {"_hhh4underreporting_shift_add_zero", (DL_FUNC) &_hhh4underreporting_shift_add_zero, 1},
    {"_hhh4underreporting_coarsen_sop_tv_cpp", (DL_FUNC) &_hhh4underreporting_coarsen_sop_tv_cpp, 1},
    {"_hhh4underreporting_back_seas_cpp", (DL_FUNC) &_hhh4underreporting_back_seas_cpp, 3},
    {"_hhh4underreporting_get_weight_matrix_seas_cpp", (DL_FUNC) &_hhh4underreporting_get_weight_matrix_seas_cpp, 3},
    {"_hhh4underreporting_prod", (DL_FUNC) &_hhh4underreporting_prod, 1},
    {"_hhh4underreporting_order_back_in_time", (DL_FUNC) &_hhh4underreporting_order_back_in_time, 2},
    {"_hhh4underreporting_stat_mean_seas_cpp", (DL_FUNC) &_hhh4underreporting_stat_mean_seas_cpp, 3},
    {"_hhh4underreporting_stat_var_seas_cpp", (DL_FUNC) &_hhh4underreporting_stat_var_seas_cpp, 4},
    {"_hhh4underreporting_nu_to_nu_star_seas_cpp", (DL_FUNC) &_hhh4underreporting_nu_to_nu_star_seas_cpp, 3},
    {"_hhh4underreporting_compute_sop_seas_cpp", (DL_FUNC) &_hhh4underreporting_compute_sop_seas_cpp, 5},
    {"_hhh4underreporting_reparam_seas_cpp", (DL_FUNC) &_hhh4underreporting_reparam_seas_cpp, 5},
    {"_hhh4underreporting_get_mod_matr_cpp", (DL_FUNC) &_hhh4underreporting_get_mod_matr_cpp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_hhh4underreporting(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
