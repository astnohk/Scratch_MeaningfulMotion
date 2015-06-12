/* Robust Estimation of Multiple Motions
 *
 * This motion estimation based on
 * M.J.Black and P.Anandan, "The Robust Estimation of Multiple Motions: Parametric and Piecewise-Smooth Flow Fields," Computer Vision and Image Understanding Vol. 63, No. 1, 1996, pp. 75-104.
 */

#include "MultiResolution.h"
#include "MEstimator.h"




VECTOR_AFFINE MultipleMotion_Affine(double *It, double *Itp1, double MaxInt, SIZE size_img, MULTIPLE_MOTION_PARAM MotionParam);
int IRLS_MultipleMotion_Affine(VECTOR_AFFINE *u, VECTOR_2D *Img_g, double *Img_t, SIZE size_img, double sigmaD, int IterMax, double ErrorMinThreshold);

VECTOR_AFFINE Error_a(VECTOR_AFFINE *u, VECTOR_2D *Img_g, double *Img_t, SIZE size_img, double sigmaD);
VECTOR_AFFINE sup_Error_aa(VECTOR_2D *Img_g, SIZE size, double sigmaD);

double Error_Affine(VECTOR_AFFINE *u, VECTOR_2D *Img_g, double *Img_t, SIZE size_img, double sigmaD);

int MultipleMotion_Affine_write(VECTOR_AFFINE u, const char *filename);

