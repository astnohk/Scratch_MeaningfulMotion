#define POW2(x) ((x) * (x))
#define POW3(x) ((x) * (x) * (x))
#define SIGN(x) ((x) >= 0 ? (x) > 0 ? 1 : 0 : -1)
#define SIGN_NOZERO(x) ((x) >= 0 ? 1 : -1)
#define SATURATE(x, min, max) (min <= (x) ? (x) <= max ? (x) : max : min)
#define MIN(x, y) ((x) <= (y) ? (x) : (y))
#define MAX(x, y) ((x) >= (y) ? (x) : (y))


#include <math.h>



double cal_rhoD(double x, double sigma);
double cal_psyD(double x, double sigma);
double cal_dpsyD(double x, double sigma);
double cal_rhoS(double x, double sigma);
double cal_psyS(double x, double sigma);
double cal_dpsyS(double x, double sigma);

