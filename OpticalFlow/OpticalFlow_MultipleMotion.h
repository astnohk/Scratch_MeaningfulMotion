/* Robust Estimation of Multiple Motions
 *
 * This motion estimation based on
 * M.J.Black and P.Anandan, "The Robust Estimation of Multiple Motions: Parametric and Piecewise-Smooth Flow Fields," Computer Vision and Image Understanding Vol. 63, No. 1, 1996, pp. 75-104.
 */

#include <string>
#include "MultiResolution.h"
#include "MEstimator.h"
#include "../MotionCompensation/MotionCompensation.h"




ImgVector<VECTOR_2D>* MultipleMotion_OpticalFlow(ImgVector<double> *It, ImgVector<double> *Itp1, double MaxInt, MULTIPLE_MOTION_PARAM MotionParam, int IterMax = 0);

bool IRLS_MultipleMotion_OpticalFlow(ImgVector<VECTOR_2D> *u, ImgVector<VECTOR_2D> *Img_g, ImgVector<double> *Img_t, double lambdaD, double lambdaS, double sigmaD, double sigmaS, int IterMax, double ErrorMinThreshold, int level);
void LevelDown(ImgVector<VECTOR_2D> *u_levels, int level);

VECTOR_2D Error_u(int site, ImgVector<VECTOR_2D> *u, ImgVector<VECTOR_2D> *Img_g, ImgVector<double> *Img_t, const double &lambdaD, const double &lambdaS, const double &sigmaD, const double &sigmaS);
VECTOR_2D sup_Error_uu(ImgVector<VECTOR_2D> *Img_g, const double &lambdaD, const double &lambdaS, const double &sigmaD, const double &sigmaS);

double Error_MultipleMotion(ImgVector<VECTOR_2D> *u, ImgVector<VECTOR_2D> *Img_g, ImgVector<double> *Img_t, const double &lambdaD, const double &lambdaS, const double &sigmaD, const double &sigmaS);

bool MultipleMotion_write(const ImgVector<double> *img_prev, const ImgVector<double> *img_next, const ImgVector<VECTOR_2D> *u, const std::string &filename);

