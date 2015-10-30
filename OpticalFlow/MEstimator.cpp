#include "MEstimator.h"



#include <cstdio>
double
Geman_McClure_rho(const double &x, const double &sigma)
{
	return POW2(x) / (sigma + POW2(x));
}


double
Geman_McClure_psi(const double &x, const double &sigma)
{
	return 2.0 * x * sigma / POW2(sigma + POW2(x));
}


double
Lorentzian_rho(const double &x, const double &sigma)
{
	return log1p(0.5 * POW2(x / sigma));
}


double
Lorentzian_psi(const double &x, const double &sigma)
{
	return 2.0 * x / (2.0 * POW2(sigma) + POW2(x));
}

