/* Robust Estimation of Multiple Motions
 *
 * This motion estimation based on
 * M.J.Black and P.Anandan, "The Robust Estimation of Multiple Motions: Parametric and Piecewise-Smooth Flow Fields," Computer Vision and Image Understanding Vol. 63, No. 1, 1996, pp. 75-104.
 */

#include "MultiResolution.h"
#include "MEstimator.h"




VECTOR_2D* MultipleMotion_OpticalFlow(double *It, double *Itp1, SIZE size_img, MULTIPLE_MOTION_PARAM MotionParam);
int IRLS_MultipleMotion_OpticalFlow(VECTOR_2D *u, VECTOR_2D *Img_g, double *Img_t, SIZE size_img, double lambdaD, double lambdaS, double sigmaD, double sigmaS);

VECTOR_2D Error_u(int site, VECTOR_2D *u, VECTOR_2D *Img_g, double *Img_t, SIZE size_img, double lambdaD, double lambdaS, double sigmaD, double sigmaS);
VECTOR_2D sup_Error_uu(VECTOR_2D *Img_g, SIZE size, double lambdaD, double lambdaS, double sigmaD, double sigmaS);

double Error_MultipleMotion(VECTOR_2D *u, VECTOR_2D *Img_g, double *Img_t, SIZE size_img, double lambdaD, double lambdaS, double sigmaD, double sigmaS);

int MultipleMotion_write(VECTOR_2D *u, SIZE size, char *filename);

