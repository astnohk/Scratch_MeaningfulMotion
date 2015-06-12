#define POW2(x) ((x) * (x))
#define POW3(x) ((x) * (x) * (x))
#define SIGN(x) ((x) >= 0 ? (x) > 0 ? 1 : 0 : -1)
#define SIGN_NOZERO(x) ((x) >= 0 ? 1 : -1)
#define SATURATE(x, min, max) (min <= (x) ? (x) <= max ? (x) : max : min)
#define MIN(x, y) ((x) <= (y) ? (x) : (y))
#define MAX(x, y) ((x) >= (y) ? (x) : (y))


#include <math.h>




double Geman_McClure_rho(double x, double sigma);
double Geman_McClure_psi(double x, double sigma);
double Lorentzian_rho(double x, double sigma);
double Lorentzian_psi(double x, double sigma);

