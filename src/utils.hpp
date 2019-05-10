#pragma once
#include <RcppArmadillo.h>
#if !defined(MYLIB_UTILS)
#define MYLIB_UTILS

const arma::cx_double i(0.0, 1.0);

double sinc( double x );
arma::vec sinc_( arma::vec x );

#endif