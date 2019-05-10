#include "hawkes.hpp"
#include "utils.hpp"

/////////////////////////////////////////////////////////////// HAWKES ///////////////////////////////////////////////////////////////
double Hawkes::G( double xi ) {
	return 1.0 / std::norm( arma::cx_double(1.0, 0.0) - H( xi ) );
};

arma::vec Hawkes::G_( arma::vec xi ) {
	arma::cx_vec temp = arma::cx_double(1.0, 0.0) - H_( xi );
	return 1.0 / arma::conv_to<arma::vec>::from( temp % arma::conj(temp) );
};

arma::vec Hawkes::dG( double xi ) {
	double Gxi = G(xi);
	return 2.0 * Gxi * Gxi * arma::real( ( 1.0 - std::conj(H(xi)) ) * dH(xi) );
};

arma::mat Hawkes::dG_( arma::vec xi ) {
	arma::mat grad( xi.n_elem, param.n_elem );
	arma::cx_mat gradH = dH_(xi);
	arma::cx_vec term1 = 1.0 - arma::conj(H_(xi));
	arma::vec Gxi = G_(xi);
	arma::vec Gxi2 = Gxi % Gxi;
	for (arma::uword k = 0; k < param.n_elem; k++) {
		grad.col(k) = 2.0 * Gxi2 % arma::real( term1 % gradH.col(k) );
	}
	return grad;
};

arma::mat Hawkes::ddG( double xi ) {
	double Gxi = G(xi);
	double Gxi2 = Gxi * Gxi;
	arma::cx_double term0 = 1.0 - std::conj(H(xi));
	arma::mat term1 = arma::real( term0 * ddH(xi) - dH(xi) * dH(xi).t() );
	arma::mat term2 = arma::real( term0 * dH(xi) ) * arma::trans(arma::real( term0 * dH(xi) ));
	return 2 * Gxi2 * (term1 + 4 * Gxi * term2);
};

arma::cube Hawkes::ddG_( arma::vec xi ) {
	arma::cube hess( param.n_elem, param.n_elem, xi.n_elem );
	arma::vec Gxi = G_(xi);
	arma::vec Gxi2 = Gxi % Gxi;
	arma::cx_vec term0 = 1.0 - arma::conj(H_(xi));
	arma::cx_mat gradH = dH_(xi);
	arma::cx_cube hessH = ddH_(xi);
	arma::vec term1(xi.n_elem);
	arma::vec term2(xi.n_elem);
	arma::cx_vec tube(xi.n_elem);
	for (arma::uword i = 0; i < param.n_elem; i++) {
		for (arma::uword j = 0; j < param.n_elem; j++) {
			tube = hessH(arma::span(i),arma::span(j), arma::span::all);
			term1 = arma::real( term0 % tube - arma::conj(gradH.col(i)) % gradH.col(j) );
			term2 = arma::real( term0 % gradH.col(i) ) % arma::real( term0 % gradH.col(j) );
			hess.tube(i, j) = 2 * Gxi2 % (term1 + 4 * Gxi % term2);
		}
	}
	return hess;
};

arma::vec Hawkes::gradf( double xi ) {
	double term0 = sinc( xi / 2.0 );
	return ddata.binsize * term0 * term0 * ( gradmean() * G( xi / ddata.binsize ) + mean() * dG( xi / ddata.binsize ) );
};

arma::mat Hawkes::gradf_( arma::vec xi ) {
	arma::mat grad( xi.n_elem, param.n_elem );
	
	arma::vec term0 = sinc_( xi / 2.0 );
	arma::vec term1 = ddata.binsize * term0 % term0;
	
	double m = mean();
	arma::vec dm = gradmean();
	arma::vec Gxb = G_( xi / ddata.binsize );
	arma::mat gradG = dG_( xi / ddata.binsize );
	
	for (arma::uword k = 0; k < param.n_elem; k++) {
		grad.col(k) = term1 % ( dm(k) * Gxb + m * gradG.col(k) );
	}
	return grad;
};

arma::mat Hawkes::hessf( double xi ) {
	double term0 = sinc( xi / 2.0 );
	return ddata.binsize * term0 * term0 * ( hessmean() * G( xi / ddata.binsize ) + 
											 gradmean() * dG( xi / ddata.binsize ).t() + 
											 dG( xi / ddata.binsize ) * gradmean().t() + 
											 mean() * ddG( xi / ddata.binsize ) );
};

arma::cube Hawkes::hessf_( arma::vec xi ) {
	arma::cube hess( param.n_elem, param.n_elem, xi.n_elem );
	
	arma::vec term0 = sinc_( xi / 2.0 );
	arma::vec term1 = ddata.binsize * term0 % term0;
	
	double m = mean();
	arma::vec dm = gradmean();
	arma::mat ddm = hessmean();
	arma::vec Gxb = G_( xi / ddata.binsize );
	arma::mat gradG = dG_( xi / ddata.binsize );
	arma::cube hessG = ddG_( xi / ddata.binsize );
	
	arma::vec tube(xi.n_elem);
	for (arma::uword i = 0; i < param.n_elem; i++) {
		for (arma::uword j = 0; j < param.n_elem; j++) {
			tube = hessG(arma::span(i),arma::span(j), arma::span::all);
			hess.tube(i, j) = term1 % ( ddm(i, j) * Gxb + dm(i) * gradG.col(j) + dm(j) * gradG.col(i) + m * tube );
		}
	}
	return hess;
};

arma::vec Hawkes::gradf1( double xi, int trunc ) {
	arma::vec omega = xi + 2.0 * arma::datum::pi * arma::regspace<arma::vec>(-trunc, trunc);
	return arma::trans( arma::sum( gradf_(omega), 0 ) );
};

arma::mat Hawkes::gradf1_( arma::vec xi, int trunc ) {
	arma::vec omega_ = 2.0 * arma::datum::pi * arma::regspace<arma::vec>(-trunc, trunc);
	arma::mat y(xi.n_elem, param.n_elem);
	for (arma::uword k = 0; k < xi.n_elem; k++) {
		y.row(k) = arma::sum( gradf_(xi(k) + omega_), 0 );
	}
	return y;
};

arma::mat Hawkes::hessf1( double xi, int trunc ) {
	arma::vec omega = xi + 2.0 * arma::datum::pi * arma::regspace<arma::vec>(-trunc, trunc);
	return arma::sum( hessf_(omega), 2 );
};

arma::cube Hawkes::hessf1_( arma::vec xi, int trunc ) {
	arma::vec omega_ = 2.0 * arma::datum::pi * arma::regspace<arma::vec>(-trunc, trunc);
	arma::cube y(param.n_elem, param.n_elem, xi.n_elem);
	for (arma::uword k = 0; k < xi.n_elem; k++) {
		y.slice(k) = arma::sum( hessf_(xi(k) + omega_), 2 );
	}
	return y;
};

/////////////////////////////////////////////////////////////// EXPHAWKES ///////////////////////////////////////////////////////////////
arma::cx_vec ExpHawkes::dH( double xi ) {
	arma::cx_vec grad = arma::zeros<arma::cx_vec>( param.n_elem );
	arma::cx_double denom = 1.0 / (param(2) - i * xi);
	grad(1) = param(2) * denom;
	grad(2) = - i * param(1) * xi * denom * denom;
	return grad;
};

arma::cx_mat ExpHawkes::dH_( arma::vec xi ) {
	arma::cx_mat grad = arma::zeros<arma::cx_mat>( xi.n_elem, param.n_elem );
	arma::cx_vec denom = 1.0 / (param(2) - i * xi);
	grad.col(1) = param(2) * denom;
	grad.col(2) = - i * param(1) * (xi % denom % denom);
	return grad;
};

arma::cx_mat ExpHawkes::ddH( double xi ) {
	arma::cx_mat hess = arma::zeros<arma::cx_mat>( param.n_elem, param.n_elem );
	arma::cx_double denom = 1.0 / (param(2) - i * xi);
	arma::cx_double grad12 = - i * xi * denom * denom;
	hess(2,1) = grad12;
	hess(1,2) = grad12;
	hess(2,2) = -2.0 * grad12 * denom;
	return hess;
};

arma::cx_cube ExpHawkes::ddH_( arma::vec xi ) {
	arma::cx_cube hess = arma::zeros<arma::cx_cube>( param.n_elem, param.n_elem, xi.n_elem );
	arma::cx_vec denom = 1.0 / (param(2) - i * xi);
	arma::cx_vec grad12 = - i * (xi % denom % denom);
	hess.tube(1,2) = grad12;
	hess.tube(2,1) = grad12;
	hess.tube(2,2) = -2.0 * grad12 % denom;
	return hess;
};