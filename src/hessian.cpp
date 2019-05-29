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
	arma::cx_mat dHxi = dH_(xi);
	arma::cx_vec term1 = 1.0 - arma::conj(H_(xi));
	arma::vec Gxi = G_(xi);
	arma::vec Gxi2 = Gxi % Gxi;
	for (arma::uword k = 0; k < param.n_elem; k++) {
		grad.col(k) = 2.0 * Gxi2 % arma::real( term1 % dHxi.col(k) );
	}
	return grad;
};

arma::mat Hawkes::ddG( double xi ) {
	double Gxi = G(xi);
	double Gxi2 = Gxi * Gxi;
	arma::cx_vec dHxi = dH(xi);
	arma::cx_double term0 = 1.0 - std::conj(H(xi));
	arma::mat term1 = arma::real( term0 * ddH(xi) - dHxi * dHxi.t() );
	arma::mat term2 = arma::real( term0 * dHxi );
	return 2 * Gxi2 * (term1 + 4 * Gxi * term2 * term2.t());
};

arma::cube Hawkes::ddG_( arma::vec xi ) {
	arma::cube hess( param.n_elem, param.n_elem, xi.n_elem );
	arma::vec Gxi = G_(xi);
	arma::vec Gxi2 = Gxi % Gxi;
	arma::cx_vec term0 = 1.0 - arma::conj(H_(xi));
	arma::cx_mat dHxi = dH_(xi);
	arma::cx_cube ddHxi = ddH_(xi);
	arma::vec term1(xi.n_elem);
	arma::vec term2(xi.n_elem);
	arma::cx_vec tube(xi.n_elem);
	for (arma::uword i = 0; i < param.n_elem; i++) {
		for (arma::uword j = 0; j < param.n_elem; j++) {
			tube = ddHxi(arma::span(i),arma::span(j), arma::span::all);
			term1 = arma::real( term0 % tube - dHxi.col(i) % arma::conj(dHxi.col(j)) );
			term2 = arma::real( term0 % dHxi.col(i) ) % arma::real( term0 % dHxi.col(j) );
			hess.tube(i, j) = 2 * Gxi2 % (term1 + 4 * Gxi % term2);
		}
	}
	return hess;
};

arma::vec Hawkes::gradf( double xi ) {
	double term0 = sinc( xi / 2.0 );
	return ddata.binsize * term0 * term0 * ( dmean() * G( xi / ddata.binsize ) + mean() * dG( xi / ddata.binsize ) );
};

arma::mat Hawkes::gradf_( arma::vec xi ) {
	arma::mat grad( xi.n_elem, param.n_elem );
	
	arma::vec term0 = sinc_( xi / 2.0 );
	arma::vec term1 = ddata.binsize * term0 % term0;
	
	double m = mean();
	arma::vec dm = dmean();
	arma::vec Gxb = G_( xi / ddata.binsize );
	arma::mat dGxb = dG_( xi / ddata.binsize );
	
	for (arma::uword k = 0; k < param.n_elem; k++) {
		grad.col(k) = term1 % ( dm(k) * Gxb + m * dGxb.col(k) );
	}
	return grad;
};

arma::mat Hawkes::hessf( double xi ) {
	double term0 = sinc( xi / 2.0 );
	return ddata.binsize * term0 * term0 * ( ddmean() * G( xi / ddata.binsize ) + 
											 dmean() * dG( xi / ddata.binsize ).t() + 
											 dG( xi / ddata.binsize ) * dmean().t() + 
											 mean() * ddG( xi / ddata.binsize ) );
};

arma::cube Hawkes::hessf_( arma::vec xi ) {
	arma::cube hess( param.n_elem, param.n_elem, xi.n_elem );
	
	arma::vec term0 = sinc_( xi / 2.0 );
	arma::vec term1 = ddata.binsize * term0 % term0;
	
	double m = mean();
	arma::vec dm = dmean();
	arma::mat ddm = ddmean();
	arma::vec Gxb = G_( xi / ddata.binsize );
	arma::mat dGxb = dG_( xi / ddata.binsize );
	arma::cube ddGxb = ddG_( xi / ddata.binsize );
	
	arma::vec tube(xi.n_elem);
	for (arma::uword i = 0; i < param.n_elem; i++) {
		for (arma::uword j = 0; j < param.n_elem; j++) {
			tube = ddGxb(arma::span(i),arma::span(j), arma::span::all);
			hess.tube(i, j) = term1 % ( ddm(i, j) * Gxb + dm(i) * dGxb.col(j) + dm(j) * dGxb.col(i) + m * tube );
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

arma::cx_mat Hawkes::wHess( arma::cx_vec& dft, int trunc ) {
	arma::uword n = dft.n_elem;
	arma::sword floor_half_n = floor( .5 * (double)n );
	arma::sword ceil_half_n = ceil( .5 * (double)n );
	double two_pi_over_n = 2 * arma::datum::pi / n;
	arma::vec omega = 2.0 * arma::datum::pi * arma::regspace<arma::vec>(0, n-1) / (double)n;
	
	arma::vec f = gammaf1_( omega, trunc );
	arma::mat df = gradf1_( omega, trunc );
	arma::vec inv_f = 1.0 / f;
	arma::mat dlogf( arma::size(df) );
	for (arma::uword k = 0; k < param.n_elem; k++) {
		dlogf.col(k) = df.col(k) % inv_f;
	}
	
	arma::mat Gamma = arma::zeros<arma::mat>( param.n_elem, param.n_elem );
	for (arma::uword k = 0; k < n; k++) {
		Gamma += dlogf.row(k).t() * dlogf.row(k);
	}
	// Gamma *= .25 * inv_pi / n; // renormalisation non nécessaire ici, car non introduite dans wLik
	arma::cx_mat term1 = arma::zeros<arma::cx_mat>( param.n_elem, param.n_elem ); 
	for (arma::sword j1 = - ceil_half_n + 1; j1 < floor_half_n + 1; j1++) {
		for (arma::sword j2 = - ceil_half_n + 1; j2 < floor_half_n + 1; j2++) {
			for (arma::sword j3 = - ceil_half_n + 1; j3 < floor_half_n + 1; j3++) {
				term1 += HT( two_pi_over_n * (double)(j2+j3), n ) * dft( modulus(j1, n) ) * dft( modulus(j2, n) ) * dft( modulus(j3, n) ) * dft( modulus(-(j1+j2+j3), n) ) * arma::conv_to<arma::cx_mat>::from( dlogf.row( modulus(j1, n) ).t() * dlogf.row( modulus(j2, n) ) );
			}
		}
	}
	term1 /= std::pow( n, 4 );

	arma::mat term2 = arma::zeros<arma::mat>( param.n_elem, param.n_elem );
	double temp;
	for (arma::uword j = 0; j < n; j++) {
		temp = std::norm( dft(j) );
		term2 += dlogf.row(j).t() * dlogf.row(j) * temp * temp;		// dlogf symmetric
	}
	term2 *= .5 * inv_pi / std::pow( n, 3 );
	
	arma::mat term3 = arma::zeros<arma::mat>( param.n_elem, param.n_elem );
	for (arma::uword j1 = 0; j1 < n; j1++) {
		for (arma::uword j2 = 0; j2 < n; j2++) {
			term3 += std::norm( dft(j1) * dft(j2) ) * dlogf.row(j1).t() * dlogf.row(j2);
		}
	}
	term3 *= 0.5 * inv_pi / std::pow( n, 3.75 );
	
	arma::cx_mat W = .5 * (term1 - arma::conv_to<arma::cx_mat>::from(term2 + term3));
	return term1;
};

/////////////////////////////////////////////////////////////// EXPHAWKES ///////////////////////////////////////////////////////////////
arma::vec ExpHawkes::dmean() {
	double denom = 1.0 / ( 1.0 - param(1) );
	arma::vec grad = { denom, param(0) * denom * denom, 0 };
	return grad;
};

arma::mat ExpHawkes::ddmean() {
	double denom = 1.0 / ( 1.0 - param(1) );
	double denom2 = denom * denom;
	arma::mat hess = { {   0.0,                  denom2, 0.0},
					   {denom2, 2*param(0)*denom2*denom, 0.0},
					   {   0.0,                     0.0, 0.0} };
	return hess;
};

arma::cx_vec ExpHawkes::dH( double xi ) {
	arma::cx_vec grad = arma::zeros<arma::cx_vec>( param.n_elem );
	arma::cx_double denom = 1.0 / (param(2) + i * xi);
	grad(1) = param(2) * denom;
	grad(2) = i * param(1) * xi * denom * denom;
	return grad;
};

arma::cx_mat ExpHawkes::dH_( arma::vec xi ) {
	arma::cx_mat grad = arma::zeros<arma::cx_mat>( xi.n_elem, param.n_elem );
	arma::cx_vec denom = 1.0 / (param(2) + i * xi);
	grad.col(1) = param(2) * denom;
	grad.col(2) = i * param(1) * (xi % denom % denom);
	return grad;
};

arma::cx_mat ExpHawkes::ddH( double xi ) {
	arma::cx_mat hess = arma::zeros<arma::cx_mat>( param.n_elem, param.n_elem );
	arma::cx_double denom = 1.0 / (param(2) + i * xi);
	arma::cx_double grad12 = i * xi * denom * denom;
	hess(2,1) = grad12;
	hess(1,2) = grad12;
	hess(2,2) = -2.0 * param(1) * grad12 * denom;
	return hess;
};

arma::cx_cube ExpHawkes::ddH_( arma::vec xi ) {
	arma::cx_cube hess = arma::zeros<arma::cx_cube>( param.n_elem, param.n_elem, xi.n_elem );
	arma::cx_vec denom = 1.0 / (param(2) + i * xi);
	arma::cx_vec grad12 = i * (xi % denom % denom);
	hess.tube(1,2) = grad12;
	hess.tube(2,1) = grad12;
	hess.tube(2,2) = -2.0 * param(1) * grad12 % denom;
	return hess;
};