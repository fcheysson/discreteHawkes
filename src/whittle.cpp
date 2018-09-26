#include <RcppArmadillo.h>
#include "FftRealPair.hpp"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp; 

//' @export Hawkes
class Hawkes {
	protected:
		// User specified
		arma::vec data;
		double binSize;
		
		// To be estimated
		arma::vec param;
		
		virtual double mean() { return 0.0; };
		
		// Useful function for spectral density
		static double sinc( double x ) {
			if (x == 0) return 1;
			return sin(x) / x;
		};
		
	public:
		// FFT from https://www.nayuki.io/page/free-small-fft-in-multiple-languages
		arma::cx_vec fft() {
			std::vector<double> real = arma::conv_to<std::vector<double>>::from(data);
			std::vector<double> imag = std::vector<double>(data.n_elem, 0.0);
			Fft::transform(real, imag);
			return arma::cx_vec( arma::conv_to<arma::vec>::from(real), arma::conv_to<arma::vec>::from(imag) );
		};
		
		// Constructor
		Hawkes( double binSize ) : binSize(binSize) {};
		Hawkes( arma::vec data, double binSize ) : data(data), binSize(binSize) {};
		virtual ~Hawkes() {};
		
		// Virtual methods for time- and frequency-domain excitation functions
		virtual double h( double x ) { return 0.0; };
		virtual arma::cx_double H( double xi ) { return arma::cx_double(0.0, 0.0); }; 
		
		// Methods for continuous- and discretized-time spectral densities
		double gammaf( double xi ) {
			double term1 = sinc( xi / 2 );
			double term2 = 1 / std::norm( arma::cx_double(1.0, 0.0) - H( xi / binSize ) );
			return mean() * binSize * term1 * term1 * term2;
		};
		double gammaf1( double xi, int trunc ) {
			double sum = 0;
			for (int k = - trunc; k < trunc + 1; k++) {
				sum += gammaf( xi + 2 * k * arma::datum::pi );
			}
			return sum;
		};
		
		// Methods for likelihood
		double wlik( arma::vec& I, int trunc ) {
		  double lik = 0;
		  double omega, spec;
		  for (arma::uword k = 1; k < data.n_elem; k++) {
		    omega = 2 * k * arma::datum::pi / data.n_elem;
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
			return param(0) / ( 1 - param(1) / param(2) );
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
		double var() {
			double kappa = 1.0 / ( 1.0 - param(1) / param(2) );
		  double kappa2 = kappa * kappa;
			double gamma = param(2) - param(1);
			return mean() * (binSize * kappa2 + (1 - kappa2) * (1 - exp(-binSize * gamma)) / gamma) ;
		};
		double cov( int tau ) {
			if (tau == 0)
				return var();
			double gamma = param(1) - param(2);
			double gamma2 = gamma * gamma;
			double gamma4 = gamma2 * gamma2;
			double expm1 = exp(gamma * binSize) - 1;
			double expm12 = expm1 * expm1;
			return .5 * param(0) * param(1) * param(2) * (2*param(2) - param(1)) * expm12 * exp(gamma * (tau - 1)) / gamma4;
		};
};

RCPP_MODULE(MyModule) {
  using namespace Rcpp;

  class_<Hawkes>("Hawkes")
	.constructor<double>() // This exposes the default constructor
	.constructor<arma::vec, double>()
	.method("gammaf", &Hawkes::gammaf)
	.method("gammaf1", &Hawkes::gammaf1)
	.method("fft", &Hawkes::fft)
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
  ;

}