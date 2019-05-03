#include "hawkes.hpp"
#include "expint.h"

//' @export
//[[Rcpp::export]]
double test( double x ) {
  return pkg_expint_E1(x, 0);
};

double sinc( double x ) {
	if (x == 0.0) return 1.0;
	return sin(x) / x;
};

arma::vec sinc_( arma::vec x ) {
	arma::vec y(x.n_elem);
	for (arma::uword k = 0; k < x.n_elem; k++) {
		if (x(k) == 0.0)
			y(k) = 1.0;
		else
			y(k) = sin(x(k)) / x(k);
		}
	return y;
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

arma::vec Hawkes::gradf( double xi ) {
	//// PLUS TARD
};

double Hawkes::wLik( arma::vec& I, int trunc ) {
	arma::uword n = I.n_elem;
	arma::vec omega = 2.0 * arma::datum::pi * arma::regspace<arma::vec>(0, n-1) / (double)n;
	arma::vec spectrum = gammaf1_( omega, trunc );
	return -arma::sum( arma::log(spectrum) + I / spectrum );
};

/////////////////////////////////////////////////////////////// EXPHAWKES ///////////////////////////////////////////////////////////////
double ExpHawkes::mean() {
	return param(0) / ( 1.0 - param(1) / param(2) );
};

// Virtual methods for time- and frequency-domain excitation functions
double ExpHawkes::h( double x ) {
	return param(1) * exp( - param(2) * x );
};

arma::vec ExpHawkes::h_( arma::vec x ) {
	return param(1) * arma::exp( - param(2) * x );
};

arma::cx_double ExpHawkes::H( double xi ) {
	double factor = param(1) / ( param(2)*param(2) + xi*xi );
	arma::cx_double zeta = arma::cx_double( factor * param(2), - factor * xi );
	return zeta;
}; 

arma::cx_vec ExpHawkes::H_( arma::vec xi ) {
	arma::vec factor = param(1) / ( param(2)*param(2) + xi%xi );
	arma::cx_vec zeta = arma::cx_vec( factor * param(2), - factor % xi );
	return zeta;
};

arma::cx_vec dH( double xi ) {
	
};

// Likelihood methods
double ExpHawkes::loglik() {
	
	// Constants
	const arma::uword n = data.n_elem;
	const double eta = param(0);
	const double alpha = param(1);
	const double beta = param(2);
	
	// Sum log \lambda*
	double A = 0.0;
	double part1 = log(eta);
	for (arma::uword i = 1; i < n; i++) {
		A = exp(- beta * (data(i) - data(i-1))) * (1.0 + A);
		part1 += log(eta + alpha * A);
	}
	
	// Int \lambda* dt
	double part2 = eta * T + (alpha / beta) * ((double)n - exp(-beta * (T - data(n-1))) * (1.0 + A));
	
	return part1 - part2;
};

arma::vec ExpHawkes::gradient() {
	
	// Constants
	const arma::uword n = data.n_elem;
	const double eta = param(0);
	const double alpha = param(1);
	const double beta = param(2);
	const double inv_beta = 1/beta;
	
	// Fill for i = 0
	arma::vec grad = { 1/eta, 0, 0 };
	
	// Iterate on arrival times
	double A = 0.0;
	double C = 0.0;
	double B;
	double denom, expint;
	for (arma::uword i = 1; i < n; i++) {
		expint = exp(- beta * (data(i) - data(i-1)));
		A = expint * (1.0 + A);
		C = expint * (data(i-1) + C);
		B = data(i) * A - C;
		denom = 1.0 / (eta + alpha * A);
		grad(0) += 1.0 * denom;
		grad(1) += A * denom;
		grad(2) -= alpha * B * denom;
	}
	
	expint = exp(- beta * (T - data(n-1)));
	A = expint * (1.0 + A);
	C = expint * (data(n-1) + C);
	B = T * A - C;
	grad(0) -= T;
	grad(1) -= inv_beta * ( (double)n - A );
	grad(2) += alpha * inv_beta * ( inv_beta * ( (double)n - A ) - B );
	
	return grad;
};

arma::mat ExpHawkes::hessian() {
	
	// Constants
	const arma::uword n = data.n_elem;
	const double eta = param(0);
	const double alpha = param(1);
	const double beta = param(2);
	const double inv_beta = 1/beta;
	
	// Fill for i = 0
	arma::mat hess = { {-1/eta, 0, 0},
					   {     0, 0, 0},
					   {     0, 0, 0} };
	
	// Iterate on arrival times
	double A = 0.0;
	double C = 0.0;
	double B;
	double E = 0.0;
	double D;
	double denom, denom2, expint;
	for (arma::uword i = 1; i < n; i++) {
		expint = exp(- beta * (data(i) - data(i-1)));
		A = expint * (1.0 + A);
		C = expint * (data(i-1) + C);
		B = data(i) * A - C;
		E = expint * (data(i-1) * data(i-1) + E);
		D = data(i) * data(i) * A + E - 2.0 * data(i) * C;
		denom = 1.0 / (eta + alpha * A);
		denom2 = denom * denom;
		hess(0,0) -= 1.0 * denom2;
		hess(0,1) -= A * denom2;
		hess(0,2) += alpha * B * denom2;
		hess(1,1) -= A * A * denom2;
		hess(1,2) -= eta * B * denom2;
		hess(2,2) += alpha * denom * (D - alpha * denom * B * B);
	}
	
	expint = exp(- beta * (T - data(n-1)));
	A = expint * (1.0 + A);
	C = expint * (data(n-1) + C);
	B = T * A - C;
	E = expint * (data(n-1) * data(n-1) + E);
	D = T * T * A + E - 2.0 * T * C;
	hess(1,2) += inv_beta * (inv_beta * ((double)n - A) - B);
	hess(2,2) += alpha * inv_beta * (D + 2.0 * inv_beta * (B - inv_beta * ((double)n  - A)));
	
	// Symm
	hess(1,0) = hess(0,1);
	hess(2,0) = hess(0,2);
	hess(2,1) = hess(1,2);
	
	return hess;
};

Rcpp::List ExpHawkes::likngrad() {
	
	// Constants
	const arma::uword n = data.n_elem;
	const double eta = param(0);
	const double alpha = param(1);
	const double beta = param(2);
	const double inv_beta = 1/beta;
	
	// Fill for i = 0
	double lik = log(eta);
	arma::vec grad = { 1/eta, 0, 0 };
	
	// Iterate on arrival times
	double A = 0.0;
	double C = 0.0;
	double B;
	double denom, expint;
	for (arma::uword i = 1; i < n; i++) {
		expint = exp(- beta * (data(i) - data(i-1)));
		A = expint * (1.0 + A);
		C = expint * (data(i-1) + C);
		B = data(i) * A - C;
		denom = 1.0 / (eta + alpha * A);
		lik += log(eta + alpha * A);
		grad(0) += 1.0 * denom;
		grad(1) += A * denom;
		grad(2) -= alpha * B * denom;
	}
	
	// Likelihood of non occurrence
	expint = exp(- beta * (T - data(n-1)));
	A = expint * (1.0 + A);
	C = expint * (data(n-1) + C);
	B = T * A - C;
	lik -= eta * T + (alpha / beta) * ((double)n - A);
	grad(0) -= T;
	grad(1) -= inv_beta * ( (double)n - A );
	grad(2) += alpha * inv_beta * ( inv_beta * ( (double)n - A ) - B );
	
	return Rcpp::List::create(Rcpp::Named("objective") = lik, Rcpp::Named("gradient") = grad);
};

/////////////////////////////////////////////////////////////// MODULE ///////////////////////////////////////////////////////////////
RCPP_MODULE(HawkesModule) {
	using namespace Rcpp;

	class_<Hawkes>("Hawkes")
		.default_constructor() // This exposes the default constructor
		.method("gammaf", &Hawkes::gammaf)
		.method("gammaf_", & Hawkes::gammaf_)
		.method("gammaf1", &Hawkes::gammaf1)
		.method("gammaf1_", &Hawkes::gammaf1_)
		.method("wLik", &Hawkes::wLik)
		.property("param", &Hawkes::getParam, &Hawkes::setParam)
		.property("data", &Hawkes::getData, &Hawkes::setData)
		.property("ddata", &Hawkes::getDData, &Hawkes::setDData)
		.property("binsize", &Hawkes::getBinsize, &Hawkes::setBinsize)	
		.property("T", &Hawkes::getT, &Hawkes::setT)
		.property("vcov", &Hawkes::getVcov, &Hawkes::setVcov)
		.property("opt", &Hawkes::getOpt, &Hawkes::setOpt)
	;
	class_<ExpHawkes>("ExpHawkes")
		.derives<Hawkes>("Hawkes")
		.default_constructor() // This exposes the default constructor
		.method("h", &ExpHawkes::h)
		.method("h_", &ExpHawkes::h_)
		.method("H", &ExpHawkes::H)
		.method("H_", &ExpHawkes::H_)
		.method("loglik", &ExpHawkes::loglik)
		.method("gradient", &ExpHawkes::gradient)
		.method("hessian", &ExpHawkes::hessian)
		.method("likngrad", &ExpHawkes::likngrad)
	;

}