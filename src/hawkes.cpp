#include "hawkes.hpp"

const arma::cx_double i(0.0, 1.0);

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

double Hawkes::whittleLik( arma::vec& I, int trunc ) {
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

// Time-domain covariance functions
// Da Fonseca, J., & Zaatour, R. (2014). Hawkes process: Fast calibration, application to trade clustering, and diffusive limit. Journal of Futures Markets (Vol. 34). http://doi.org/10.1002/fut.21644
double ExpHawkes::var() {
	double kappa = 1.0 / ( 1.0 - param(1) / param(2) );
	double kappa2 = kappa * kappa;
	double gamma = param(2) - param(1);
	return mean() * (ddata.binsize * kappa2 + (1.0 - kappa2) * (1.0 - exp(-ddata.binsize * gamma)) / gamma) ;
};

double ExpHawkes::cov( unsigned int tau ) { // tau = 0 means cov(N[0, binsize], N[binsize, 2*binsize])
	double gamma = param(1) - param(2);
	double gamma2 = gamma * gamma;
	double gamma4 = gamma2 * gamma2;
	double expm1 = exp(gamma * ddata.binsize) - 1.0;
	double expm12 = expm1 * expm1;
	return .5 * param(0) * param(1) * param(2) * (2.0*param(2) - param(1)) * expm12 * exp(gamma * (double)tau) / gamma4;
};

arma::vec ExpHawkes::cov_( arma::uvec tau ) {
	double gamma = param(1) - param(2);
	double gamma2 = gamma * gamma;
	double gamma4 = gamma2 * gamma2;
	double expm1 = exp(gamma * ddata.binsize) - 1.0;
	double expm12 = expm1 * expm1;
	double factor = .5 * param(0) * param(1) * param(2) * (2.0*param(2) - param(1)) * expm12 / gamma4;
	return factor * arma::exp( gamma * arma::conv_to<arma::vec>::from(tau) );
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
		.method("whittleLik", &Hawkes::whittleLik)
		.property("param", &Hawkes::getParam, &Hawkes::setParam)
		.property("data", &Hawkes::getData, &Hawkes::setData)
		.property("ddata", &Hawkes::getDData, &Hawkes::setDData)
		.property("binsize", &Hawkes::getBinsize, &Hawkes::setBinsize)		
	;
	class_<ExpHawkes>("ExpHawkes")
		.derives<Hawkes>("Hawkes")
		.default_constructor() // This exposes the default constructor
		.method("h", &ExpHawkes::h)
		.method("h_", &ExpHawkes::h_)
		.method("H", &ExpHawkes::H)
		.method("H_", &ExpHawkes::H_)
		.method("var", &ExpHawkes::var)
		.method("cov", &ExpHawkes::cov)
		.method("cov_", &ExpHawkes::cov_)
	;

}

// class ExpHawkes: public Hawkes {
	// public:
		
		
		// // Methods for de-biased Whittle
		// // Sykulski, A. M., Olhede, S. C., & Lilly, J. M. (2016). The De-Biased Whittle Likelihood for Second-Order Stationary Stochastic Processes, 1–28. Retrieved from http://arxiv.org/abs/1605.06718
		// double dgammaf( double xi ) {
			// arma::cx_double tmp(0.0, 0.0);
			// arma::uword n = data.n_elem;
			// for (arma::uword tau = 0; tau < n; tau++) {
			// // NEED TO USE FFT HERE TO SPEED UP CALCULATION
				// tmp += (1.0 - (double)tau / (double)n) * cov(tau) * exp(-i * xi * (double)tau);
			// }
			// return 2 * real(tmp) - var();
		// };
		// double wlikCov( arma::vec& I ) {
			// double lik = 0.0;
			// double omega;
			// for (arma::uword k = 1; k < data.n_elem; k++) {
				// omega = 2.0 * (double)k * arma::datum::pi / data.n_elem;
				// lik -= log(dgammaf(omega)) + I(k) / dgammaf(omega);
			// }
			// return lik;
		// };
// };