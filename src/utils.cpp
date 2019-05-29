#include "utils.hpp"

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

arma::uword modulus( arma::sword a, arma::sword b ) {
	return (arma::uword)( (a%b+b)%b );
};

arma::uword modulus( arma::sword a, arma::uword b ) {
	arma::sword c = (arma::sword)b;
	return (arma::uword)( (a%c+c)%c );
};

double HT( double x, double T ) {
	double T4 = std::pow( T, 0.25 );
	if ( x*abs(T4) > arma::datum::pi ) return 0.0;
	else return 0.5 * T4 * inv_pi;
};