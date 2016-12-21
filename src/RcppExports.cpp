// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// cppGeneratePartitions
List cppGeneratePartitions(List adjList, int numBlocks, NumericVector popSizes, int numConstraintLow, int numConstraintHigh, double popConstraintLow, double popConstraintHigh);
RcppExport SEXP redist_cppGeneratePartitions(SEXP adjListSEXP, SEXP numBlocksSEXP, SEXP popSizesSEXP, SEXP numConstraintLowSEXP, SEXP numConstraintHighSEXP, SEXP popConstraintLowSEXP, SEXP popConstraintHighSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type adjList(adjListSEXP);
    Rcpp::traits::input_parameter< int >::type numBlocks(numBlocksSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type popSizes(popSizesSEXP);
    Rcpp::traits::input_parameter< int >::type numConstraintLow(numConstraintLowSEXP);
    Rcpp::traits::input_parameter< int >::type numConstraintHigh(numConstraintHighSEXP);
    Rcpp::traits::input_parameter< double >::type popConstraintLow(popConstraintLowSEXP);
    Rcpp::traits::input_parameter< double >::type popConstraintHigh(popConstraintHighSEXP);
    rcpp_result_gen = Rcpp::wrap(cppGeneratePartitions(adjList, numBlocks, popSizes, numConstraintLow, numConstraintHigh, popConstraintLow, popConstraintHigh));
    return rcpp_result_gen;
END_RCPP
}
// countpartitions
int countpartitions(List aList);
RcppExport SEXP redist_countpartitions(SEXP aListSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type aList(aListSEXP);
    rcpp_result_gen = Rcpp::wrap(countpartitions(aList));
    return rcpp_result_gen;
END_RCPP
}
// calcPWDh
NumericMatrix calcPWDh(NumericMatrix x);
RcppExport SEXP redist_calcPWDh(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(calcPWDh(x));
    return rcpp_result_gen;
END_RCPP
}
// segregationcalc
NumericVector segregationcalc(NumericMatrix distmat, NumericVector grouppop, NumericVector fullpop);
RcppExport SEXP redist_segregationcalc(SEXP distmatSEXP, SEXP grouppopSEXP, SEXP fullpopSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type distmat(distmatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type grouppop(grouppopSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type fullpop(fullpopSEXP);
    rcpp_result_gen = Rcpp::wrap(segregationcalc(distmat, grouppop, fullpop));
    return rcpp_result_gen;
END_RCPP
}
// rsg
List rsg(List adj_list, NumericVector population, int Ndistrict, double target_pop, double thresh, int maxiter);
RcppExport SEXP redist_rsg(SEXP adj_listSEXP, SEXP populationSEXP, SEXP NdistrictSEXP, SEXP target_popSEXP, SEXP threshSEXP, SEXP maxiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type adj_list(adj_listSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type population(populationSEXP);
    Rcpp::traits::input_parameter< int >::type Ndistrict(NdistrictSEXP);
    Rcpp::traits::input_parameter< double >::type target_pop(target_popSEXP);
    Rcpp::traits::input_parameter< double >::type thresh(threshSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    rcpp_result_gen = Rcpp::wrap(rsg(adj_list, population, Ndistrict, target_pop, thresh, maxiter));
    return rcpp_result_gen;
END_RCPP
}
// swMH
List swMH(List aList, NumericVector cdvec, NumericVector cdorigvec, NumericVector popvec, NumericVector grouppopvec, int nsims, double eprob, double pct_dist_parity, NumericVector beta_sequence, NumericVector beta_weights, NumericMatrix ssdmat, int lambda, double beta_population, double beta_compact, double beta_segregation, double beta_similar, int anneal_beta_population, int anneal_beta_compact, int anneal_beta_segregation, int anneal_beta_similar, int adjswap, int exact_mh);
RcppExport SEXP redist_swMH(SEXP aListSEXP, SEXP cdvecSEXP, SEXP cdorigvecSEXP, SEXP popvecSEXP, SEXP grouppopvecSEXP, SEXP nsimsSEXP, SEXP eprobSEXP, SEXP pct_dist_paritySEXP, SEXP beta_sequenceSEXP, SEXP beta_weightsSEXP, SEXP ssdmatSEXP, SEXP lambdaSEXP, SEXP beta_populationSEXP, SEXP beta_compactSEXP, SEXP beta_segregationSEXP, SEXP beta_similarSEXP, SEXP anneal_beta_populationSEXP, SEXP anneal_beta_compactSEXP, SEXP anneal_beta_segregationSEXP, SEXP anneal_beta_similarSEXP, SEXP adjswapSEXP, SEXP exact_mhSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type aList(aListSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cdvec(cdvecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cdorigvec(cdorigvecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type popvec(popvecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type grouppopvec(grouppopvecSEXP);
    Rcpp::traits::input_parameter< int >::type nsims(nsimsSEXP);
    Rcpp::traits::input_parameter< double >::type eprob(eprobSEXP);
    Rcpp::traits::input_parameter< double >::type pct_dist_parity(pct_dist_paritySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta_sequence(beta_sequenceSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta_weights(beta_weightsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ssdmat(ssdmatSEXP);
    Rcpp::traits::input_parameter< int >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type beta_population(beta_populationSEXP);
    Rcpp::traits::input_parameter< double >::type beta_compact(beta_compactSEXP);
    Rcpp::traits::input_parameter< double >::type beta_segregation(beta_segregationSEXP);
    Rcpp::traits::input_parameter< double >::type beta_similar(beta_similarSEXP);
    Rcpp::traits::input_parameter< int >::type anneal_beta_population(anneal_beta_populationSEXP);
    Rcpp::traits::input_parameter< int >::type anneal_beta_compact(anneal_beta_compactSEXP);
    Rcpp::traits::input_parameter< int >::type anneal_beta_segregation(anneal_beta_segregationSEXP);
    Rcpp::traits::input_parameter< int >::type anneal_beta_similar(anneal_beta_similarSEXP);
    Rcpp::traits::input_parameter< int >::type adjswap(adjswapSEXP);
    Rcpp::traits::input_parameter< int >::type exact_mh(exact_mhSEXP);
    rcpp_result_gen = Rcpp::wrap(swMH(aList, cdvec, cdorigvec, popvec, grouppopvec, nsims, eprob, pct_dist_parity, beta_sequence, beta_weights, ssdmat, lambda, beta_population, beta_compact, beta_segregation, beta_similar, anneal_beta_population, anneal_beta_compact, anneal_beta_segregation, anneal_beta_similar, adjswap, exact_mh));
    return rcpp_result_gen;
END_RCPP
}
// genAlConn
List genAlConn(List aList, NumericVector cds);
RcppExport SEXP redist_genAlConn(SEXP aListSEXP, SEXP cdsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type aList(aListSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cds(cdsSEXP);
    rcpp_result_gen = Rcpp::wrap(genAlConn(aList, cds));
    return rcpp_result_gen;
END_RCPP
}
