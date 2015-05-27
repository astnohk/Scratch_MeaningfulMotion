#include "MEstimator.h"



double
cal_rhoD(double x, double sigma)
{
	return log(1.0 + 0.5 * POW2(x / sigma));
	//return POW2(x) / (sigma + POW2(x));
}


double
cal_psyD(double x, double sigma)
{
	return 2.0 * x / (2.0 * POW2(sigma) + POW2(x));
	//return 2.0 * x * sigma / POW2(sigma + POW2(x));
}


double
cal_dpsyD(double x, double sigma)
{
	return 4.0 * POW2(sigma) / POW2(2.0 * POW2(sigma) + POW2(x));
	//return 2.0 * sigma * (sigma - 3.0 * POW2(x)) / POW3(sigma + POW2(x));
}


double
cal_rhoS(double x, double sigma)
{
	return log(1.0 + 0.5 * POW2(x / sigma));
	//return POW2(x) / (sigma + POW2(x));
}


double
cal_psyS(double x, double sigma)
{
	return 2.0 * x / (2.0 * POW2(sigma) + POW2(x));
	//return 2.0 * x * sigma / POW2(sigma + POW2(x));
}


double
cal_dpsyS(double x, double sigma)
{
	return 4.0 * POW2(sigma) / POW2(2.0 * POW2(sigma) + POW2(x));
	//return 2.0 * sigma * (sigma - 3.0 * POW2(x)) / POW3(sigma + POW2(x));
}

