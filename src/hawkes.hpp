#pragma once
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp; 

struct DData {
	arma::vec x;
	double binsize;
};

//' @export Hawkes
class Hawkes {
	protected:
		arma::vec data;
		DData ddata;
		double T;
		arma::vec param;
		arma::mat vcov;
		Rcpp::List opt;
		
	public:
		virtual ~Hawkes() {};
		
		virtual double mean() { return 0.0; };
		
		// Virtual methods for time- and frequency-domain excitation functions
		virtual double h( double x ) { return 0.0; };
		virtual arma::vec h_( arma::vec x ) { return arma::zeros<arma::vec>(x.n_elem); };
		virtual arma::cx_double H( double xi ) { return arma::cx_double(0.0, 0.0); }; 
		virtual arma::cx_vec H_( arma::vec xi ) { return arma::zeros<arma::cx_vec>(xi.n_elem); };
		virtual arma::cx_vec dH( double xi ) { return arma::zeros<arma::cx_vec>(param.n_elem); };
		virtual arma::cx_mat dH_( arma::vec xi ) { return arma::zeros<arma::cx_mat>(param.n_elem, xi.n_elem); };
		
		// Methods for continuous- and discretized-time spectral densities
		double gammaf( double xi );
		arma::vec gammaf_( arma::vec xi );
		double gammaf1( double xi, int trunc );
		arma::vec gammaf1_( arma::vec xi, int trunc );
		arma::vec gradf( double xi );
		
		// Likelihood estimation methods
		virtual double loglik() { return 0.0; };
		virtual arma::vec gradient() { return arma::zeros<arma::vec>(param.n_elem); };
		virtual arma::mat hessian() { return arma::zeros<arma::mat>(param.n_elem, param.n_elem); };
		virtual Rcpp::List likngrad() { return Rcpp::List::create(); };
		
		// Whittle likelihood estimation methods
		double wLik( arma::vec& I, int trunc );
		
		// Get and set methods
		void setParam( arma::vec param_ ) {
			param = param_;
		};
		arma::vec getParam() {
			return param;
		};
		
		void setData( arma::vec data_ ) {
			data = data_;
		};
		arma::vec getData() {
			return data;
		};
		
		void setDData( arma::vec ddata_ ) {
			ddata.x = ddata_;
		};
		arma::vec getDData() {
			return ddata.x;
		};
		
		void setBinsize( double binsize_ ) {
			ddata.binsize = binsize_;
		};
		double getBinsize() {
			return ddata.binsize;
		};
		
		void setT( double T_ ) {
			T = T_;
		};
		double getT() {
			return T;
		};
		
		void setVcov( arma::mat vcov_ ) {
			vcov = vcov_;
		};
		arma::mat getVcov() {
			return vcov;
		};
		
		void setOpt( Rcpp::List opt_ ) {
			opt = opt_;
		};
		Rcpp::List getOpt() {
			return opt;
		};
};

//' @export ExpHawkes
class ExpHawkes: public Hawkes {
	public:
		double mean();
		
		// Methods for time- and frequency-domain excitation functions
		double h( double x );
		arma::vec h_( arma::vec x );
		arma::cx_double H( double xi ); 
		arma::cx_vec H_( arma::vec xi );
		arma::cx_vec dH( double xi );
		arma::cx_mat dH_( arma::vec xi ); 
		
		// Likelihood methods
		double loglik();
		arma::vec gradient();
		arma::mat hessian();
		Rcpp::List likngrad();
};