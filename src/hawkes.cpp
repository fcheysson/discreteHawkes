#include "hawkes.hpp"
#include "utils.hpp"
#include "expint.h"

//' @export
//[[Rcpp::export]]
double test( double x ) {
  return pkg_expint_E1(x, 0);
};

/////////////////////////////////////////////////////////////// HAWKES ///////////////////////////////////////////////////////////////
// Methods for continuous- and discretized-time spectral densities
double Hawkes::gammaf( double xi ) {
	double term1 = sinc( xi / 2.0 );
	double term2 = 1.0 / std::norm( arma::cx_double(1.0, 0.0) - H( xi / ddata.binsize ) );
	return mean() * ddata.binsize * term1 * term1 * term2; 
};

arma::vec Hawkes::gammaf_( arma::vec xi ) {
	arma::vec term1 = sinc_( xi / 2.0 );
	arma::cx_vec term2 = arma::cx_double(1.0, 0.0) - H_( xi / ddata.binsize );
	return mean() * ddata.binsize * term1 % term1 / arma::conv_to<arma::vec>::from( term2 % arma::conj(term2) );
};

double Hawkes::gammaf1( double xi, int trunc ) {
	arma::vec omega = xi + 2.0 * arma::datum::pi * arma::regspace<arma::vec>(-trunc, trunc);
	return arma::sum( gammaf_(omega) );
};

arma::vec Hawkes::gammaf1_( arma::vec xi, int trunc ) {
	arma::vec omega_ = 2.0 * arma::datum::pi * arma::regspace<arma::vec>(-trunc, trunc);
	arma::vec y(xi.n_elem);
	for (arma::uword k = 0; k < xi.n_elem; k++) {
		y(k) = arma::sum( gammaf_(xi(k) + omega_) );
	}
	return y;
};

double Hawkes::wLik( arma::vec& I, int trunc ) {
	arma::uword n = I.n_elem;
	arma::vec omega = 2.0 * arma::datum::pi * arma::regspace<arma::vec>(0, n-1) / (double)n;
	arma::vec spectrum = gammaf1_( omega, trunc );
	return arma::sum( arma::log(spectrum) + I / spectrum );
};

/////////////////////////////////////////////////////////////// EXPHAWKES ///////////////////////////////////////////////////////////////
double ExpHawkes::mean() {
	return param(0) / ( 1.0 - param(1) );
};

// Virtual methods for time- and frequency-domain excitation functions
double ExpHawkes::h( double x ) {
	return param(1) * param(2) * exp( - param(2) * x );
};

arma::vec ExpHawkes::h_( arma::vec x ) {
	return param(1) * param(2) * arma::exp( - param(2) * x );
};

arma::cx_double ExpHawkes::H( double xi ) {
	double factor = param(1) * param(2) / ( param(2)*param(2) + xi*xi );
	arma::cx_double zeta = arma::cx_double( factor * param(2), - factor * xi );
	return zeta;
}; 

arma::cx_vec ExpHawkes::H_( arma::vec xi ) {
	arma::vec factor = param(1) * param(2) / ( param(2)*param(2) + xi%xi );
	arma::cx_vec zeta = arma::cx_vec( factor * param(2), - factor % xi );
	return zeta;
};


