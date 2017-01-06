/* Robust Estimation of Multiple Motions
 *
 * This motion estimation based on
 * M.J.Black and P.Anandan, "The Robust Estimation of Multiple Motions: Parametric and Piecewise-Smooth Flow Fields," Computer Vision and Image Understanding Vol. 63, No. 1, 1996, pp. 75-104.
 */
#include <string>

#define OUTPUT_IMG_CLASS

#include "../ImgClass/Color.h"
#include "../ImgClass/Vector.h"
#include "../ImgClass/ImgClass.h"
#include "../ImgClass/MotionCompensation.h"

#include "MEstimator.h"
#include "MultiResolution.h"
#include "../lib/Struct.h"
#include "../Scratch_MeaningfulMotion.h"




ImgVector<VECTOR_2D<double> >* OpticalFlow_Pyramid(ImgVector<double> *It, ImgVector<double> *Itp1, double MaxInt, MULTIPLE_MOTION_PARAM MotionParam, int IterMax = 0);

void LevelDown(ImgVector<double> *I_dt_levels, ImgVector<VECTOR_2D<double> > *u_levels, const ImgVector<double> *It_levels, const ImgVector<double> *Itp1_levels, int level, int MaxLevel);
void Add_VectorOffset(ImgVector<VECTOR_2D<double> > *u_levels, int level, int MaxLevel);

void IRLS_OpticalFlow_Pyramid(ImgVector<VECTOR_2D<double> > *u, const ImgVector<VECTOR_2D<double> > *Img_g, const ImgVector<double> *Img_t, double lambdaD, double lambdaS, double sigmaD, double sigmaS, int IterMax, double ErrorMinThreshold, int level);

VECTOR_2D<double> Error_u(const size_t& site, const ImgVector<VECTOR_2D<double> > *u, const ImgVector<VECTOR_2D<double> > *Img_g, const ImgVector<double> *Img_t, const double &lambdaD, const double &lambdaS, const double &sigmaD, const double &sigmaS);
VECTOR_2D<double> sup_Error_uu(const ImgVector<VECTOR_2D<double> > *Img_g, const double &lambdaD, const double &lambdaS, const double &sigmaD, const double &sigmaS);

double Error_MultipleMotion(const ImgVector<VECTOR_2D<double> > *u, const ImgVector<VECTOR_2D<double> > *Img_g, const ImgVector<double> *Img_t, const double &lambdaD, const double &lambdaS, const double &sigmaD, const double &sigmaS);

void MultipleMotion_write(const ImgVector<double>& img_prev, const ImgVector<double>& img_next, const ImgVector<VECTOR_2D<double> >& u, const std::string& filename);
void MultipleMotion_write(const ImgVector<ImgClass::RGB>& img_prev, const ImgVector<ImgClass::RGB>& img_next, const ImgVector<VECTOR_2D<double> >& u, const std::string& filename);

