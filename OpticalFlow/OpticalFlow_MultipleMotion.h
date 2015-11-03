/* Robust Estimation of Multiple Motions
 *
 * This motion estimation based on
 * M.J.Black and P.Anandan, "The Robust Estimation of Multiple Motions: Parametric and Piecewise-Smooth Flow Fields," Computer Vision and Image Understanding Vol. 63, No. 1, 1996, pp. 75-104.
 */

#include <string>
#include <vector>

#include "MEstimator.h"
#include "MultiResolution.h"

#include "../lib/Struct.h"
#include "../ImgClass/Vector.h"
#include "../MotionCompensation/MotionCompensation.h"
#include "../Scratch_MeaningfulMotion.h"




ImgVector<VECTOR_2D<double> >* MultipleMotion_OpticalFlow(ImgVector<double> *It, ImgVector<double> *Itp1, double MaxInt, MULTIPLE_MOTION_PARAM MotionParam, int IterMax = 0);

void IRLS_MultipleMotion_OpticalFlow(ImgVector<VECTOR_2D<double> > *u, ImgVector<VECTOR_2D<double> > *Img_g, ImgVector<double> *Img_t, double lambdaD, double lambdaS, double sigmaD, double sigmaS, int IterMax, double ErrorMinThreshold, int level);
void LevelDown(ImgVector<double> *It_levels, ImgVector<double> *Itp1_levels, ImgVector<double> *I_dt_levels, ImgVector<VECTOR_2D<double> > *u_levels, int level, int MaxLevel);

VECTOR_2D<double> Error_u(int site, ImgVector<VECTOR_2D<double> > *u, ImgVector<VECTOR_2D<double> > *Img_g, ImgVector<double> *Img_t, const double &lambdaD, const double &lambdaS, const double &sigmaD, const double &sigmaS);
VECTOR_2D<double> sup_Error_uu(ImgVector<VECTOR_2D<double> > *Img_g, const double &lambdaD, const double &lambdaS, const double &sigmaD, const double &sigmaS);

double Error_MultipleMotion(ImgVector<VECTOR_2D<double> > *u, ImgVector<VECTOR_2D<double> > *Img_g, ImgVector<double> *Img_t, const double &lambdaD, const double &lambdaS, const double &sigmaD, const double &sigmaS);

void MultipleMotion_write(const ImgVector<double> *img_prev, const ImgVector<double> *img_next, const ImgVector<VECTOR_2D<double> > *u, const std::string &filename);

