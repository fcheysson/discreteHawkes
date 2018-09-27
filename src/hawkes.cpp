#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp; 

const arma::cx_double i(0.0, 1.0);

struct DData {
	arma::vec x;
	double binSize;
};

//' @export Hawkes
class Hawkes {
	protected:
		arma::vec data;
		DData ddata;
		arma::vec param;
	public:
		// Constructor
		Hawkes() : {};
		Hawkes( arma::vec data ) : Data(
		Hawkes( arma::vec data, double binSize ) : data(data), binSize(binSize) {};
		virtual ~Hawkes() {};
		
		// Virtual methods for time- and frequency-domain excitation functions
		virtual double h( double x ) { return 0.0; };
		virtual arma::cx_double H( double xi ) { return arma::cx_double(0.0, 0.0); }; 
		
		// Methods for continuous- and discretized-time spectral densities
		double gammaf( double xi ) {
			double term1 = sinc( xi / 2.0 );
			double term2 = 1.0 / std::norm( arma::cx_double(1.0, 0.0) - H( xi / binSize ) );
			return mean() * binSize * term1 * term1 * term2;
		};
		double gammaf1( double xi, int trunc ) {
			double sum = 0.0;
			for (int k = - trunc; k < trunc + 1; k++) {
				sum += gammaf( xi + 2.0 * (double)k * arma::datum::pi );
			}
			return sum;
		};
		
		// Method for Whittle likelihood
		double wlik( arma::vec& I, int trunc ) {
			double lik = 0.0;
			double omega, spec;
			for (arma::uword k = 1; k < data.n_elem; k++) {
				omega = 2.0 * (double)k * arma::datum::pi / data.n_elem;
				spec = gammaf1(omega, trunc);
				lik -= log(spec) + I(k) / spec;
			}
			return lik;
		};
		
		// Get and set methods for param
		void setParam( arma::vec param_ ) {
			param = param_;
		};
		arma::vec getParam() {
			return param;
		};
};

//' @export ExpHawkes
class ExpHawkes: public Hawkes {
	private:
		double mean() {
			return param(0) / ( 1.0 - param(1) / param(2) );
		};
	public:
		// Constructor
		ExpHawkes( double binSize ) : Hawkes(binSize) {};
		ExpHawkes( arma::vec data, double binSize ) : Hawkes(data, binSize) {};
		
		// Time- and frequency-domain excitation functions
		double h( double x ) {
			return param(1) * exp( - param(2) * x );
		};
		arma::cx_double H( double xi ) {
			double factor = param(1) / ( param(2)*param(2) + xi*xi );
			arma::cx_double result = arma::cx_double(factor * param(2), - factor * xi);
			return result;
		};
		
		// Time-domain covariance functions
		// Da Fonseca, J., & Zaatour, R. (2014). Hawkes process: Fast calibration, application to trade clustering, and diffusive limit. Journal of Futures Markets (Vol. 34). http://doi.org/10.1002/fut.21644
		double var() {
			double kappa = 1.0 / ( 1.0 - param(1) / param(2) );
			double kappa2 = kappa * kappa;
			double gamma = param(2) - param(1);
			return mean() * (binSize * kappa2 + (1.0 - kappa2) * (1.0 - exp(-binSize * gamma)) / gamma) ;
		};
		double cov( int tau ) {
			if (tau == 0)
				return var();
			double gamma = param(1) - param(2);
			double gamma2 = gamma * gamma;
			double gamma4 = gamma2 * gamma2;
			double expm1 = exp(gamma * binSize) - 1.0;
			double expm12 = expm1 * expm1;
			return .5 * param(0) * param(1) * param(2) * (2.0*param(2) - param(1)) * expm12 * exp(gamma * (double)(tau - 1)) / gamma4;
		};
		arma::vec covVec( arma::ivec tau ) {
			double gamma = param(1) - param(2);
			double gamma2 = gamma * gamma;
			double gamma4 = gamma2 * gamma2;
			double expm1 = exp(gamma * binSize) - 1.0;
			double expm12 = expm1 * expm1;
			double cmn = .5 * param(0) * param(1) * param(2) * (2.0*param(2) - param(1)) * expm12 / gamma4;
			arma::vec rst(tau.n_elem);
			for (arma::uword k = 0; k < tau.n_elem; k++) {
				if (tau(k) == 0)
					rst(k) = var();
				else
					rst(k) = cmn * exp(gamma * (double)(tau(k) - 1));
			}
			return rst;
		};
		
		// Methods for de-biased Whittle
		// Sykulski, A. M., Olhede, S. C., & Lilly, J. M. (2016). The De-Biased Whittle Likelihood for Second-Order Stationary Stochastic Processes, 1–28. Retrieved from http://arxiv.org/abs/1605.06718
		double dgammaf( double xi ) {
			arma::cx_double tmp(0.0, 0.0);
			arma::uword n = data.n_elem;
			for (arma::uword tau = 0; tau < n; tau++) {
			// NEED TO USE FFT HERE TO SPEED UP CALCULATION
				tmp += (1.0 - (double)tau / (double)n) * cov(tau) * exp(-i * xi * (double)tau);
			}
			return 2 * real(tmp) - var();
		};
		double wlikCov( arma::vec& I ) {
			double lik = 0.0;
			double omega;
			for (arma::uword k = 1; k < data.n_elem; k++) {
				omega = 2.0 * (double)k * arma::datum::pi / data.n_elem;
				lik -= log(dgammaf(omega)) + I(k) / dgammaf(omega);
			}
			return lik;
		};
};

RCPP_MODULE(MyModule) {
	using namespace Rcpp;

	class_<Hawkes>("Hawkes")
		.constructor<double>() // This exposes the default constructor
		.constructor<arma::vec, double>()
		.method("gammaf", &Hawkes::gammaf)
		.method("gammaf1", &Hawkes::gammaf1)
		.method("wlik", &Hawkes::wlik)
		.property("param", &Hawkes::getParam, &Hawkes::setParam)
	;
	class_<ExpHawkes>("ExpHawkes")
		.derives<Hawkes>("Hawkes")
		.constructor<double>() // This exposes the default constructor
		.constructor<arma::vec, double>()
		.method("h", &ExpHawkes::h)
		.method("H", &ExpHawkes::H) // This exposes the H method
		.method("var", &ExpHawkes::var)
		.method("cov", &ExpHawkes::cov)
		.method("covVec", &ExpHawkes::covVec)
		.method("dgammaf", &ExpHawkes::dgammaf)
		.method("wlikCov", &ExpHawkes::wlikCov)
	;

}