#include "hawkes.hpp"

RCPP_MODULE(HawkesModule) {
	using namespace Rcpp;

	class_<Hawkes>("Hawkes")
		.default_constructor() // This exposes the default constructor
		.method("gammaf", &Hawkes::gammaf)
		.method("gammaf_", & Hawkes::gammaf_)
		.method("gammaf1", &Hawkes::gammaf1)
		.method("gammaf1_", &Hawkes::gammaf1_)
		.method("G", &Hawkes::G)
		.method("G_", &Hawkes::G_)
		.method("dG", &Hawkes::dG)
		.method("dG_", &Hawkes::dG_)
		.method("ddG", &Hawkes::ddG)
		.method("ddG_", &Hawkes::ddG_)
		.method("gradf", &Hawkes::gradf)
		.method("gradf_", &Hawkes::gradf_)
		.method("hessf", &Hawkes::hessf)
		.method("hessf_", &Hawkes::hessf_)
		.method("gradf1", &Hawkes::gradf1)
		.method("gradf1_", &Hawkes::gradf1_)
		.method("hessf1", &Hawkes::hessf1)
		.method("hessf1_", &Hawkes::hessf1_)
		.method("wLik", &Hawkes::wLik)
		.method("wHess", &Hawkes::wHess)
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
		.method("mean", &ExpHawkes::mean)
		.method("gradmean", &ExpHawkes::gradmean)
		.method("hessmean", &ExpHawkes::hessmean)
		.method("h", &ExpHawkes::h)
		.method("h_", &ExpHawkes::h_)
		.method("H", &ExpHawkes::H)
		.method("H_", &ExpHawkes::H_)
		.method("dH", &ExpHawkes::dH)
		.method("dH_", &ExpHawkes::dH_)
		.method("ddH", &ExpHawkes::ddH)
		.method("ddH_", &ExpHawkes::ddH_)
		.method("loglik", &ExpHawkes::loglik)
		.method("gradient", &ExpHawkes::gradient)
		.method("hessian", &ExpHawkes::hessian)
		.method("likngrad", &ExpHawkes::likngrad)
	;

}