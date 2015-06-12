#include "MEstimator.h"




double
Geman_McClure_rho(double x, double sigma)
{
	return POW2(x) / (sigma + POW2(x));
}


double
Geman_McClure_psi(double x, double sigma)
{
	return 2.0 * x * sigma / POW2(sigma + POW2(x));
}


double
Lorentzian_rho(double x, double sigma)
{
	return log(1.0 + 0.5 * POW2(x / sigma));
}


double
Lorentzian_psi(double x, double sigma)
{
	return 2.0 * x / (2.0 * POW2(sigma) + POW2(x));
}

