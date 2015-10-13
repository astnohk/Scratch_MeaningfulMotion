/* Robust Estimation of Multiple Motions
 *
 * This motion estimation based on
 * M.J.Black and P.Anandan, "The Robust Estimation of Multiple Motions: Parametric and Piecewise-Smooth Flow Fields," Computer Vision and Image Understanding Vol. 63, No. 1, 1996, pp. 75-104.
 */

#include "MultiResolution.h"
#include "MEstimator.h"




VECTOR_2D* MultipleMotion_OpticalFlow(double *It, double *Itp1, double MaxInt, SIZE size_img, MULTIPLE_MOTION_PARAM MotionParam);

bool IRLS_MultipleMotion_OpticalFlow(VECTOR_2D *u, VECTOR_2D *Img_g, double *Img_t, SIZE size_img, double lambdaD, double lambdaS, double sigmaD, double sigmaS, int IterMax, double ErrorMinThreshold, int level);
void LevelDown(VECTOR_2D *u_l, SIZE size_l, VECTOR_2D *u_lp1, const SIZE &size_lp1);

VECTOR_2D Error_u(int site, VECTOR_2D *u, VECTOR_2D *Img_g, double *Img_t, const SIZE &size_img, const double &lambdaD, const double &lambdaS, const double &sigmaD, const double &sigmaS);
VECTOR_2D sup_Error_uu(VECTOR_2D *Img_g, const SIZE &size, const double &lambdaD, const double &lambdaS, const double &sigmaD, const double &sigmaS);

double Error_MultipleMotion(VECTOR_2D *u, VECTOR_2D *Img_g, double *Img_t, const SIZE &size_img, const double &lambdaD, const double &lambdaS, const double &sigmaD, const double &sigmaS);

bool MultipleMotion_write(VECTOR_2D *u, SIZE size, const char *filename);

