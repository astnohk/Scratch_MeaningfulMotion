/* Robust Estimation of Multiple Motions
 *
 * This motion estimation based on
 * M.J.Black and P.Anandan, "The Robust Estimation of Multiple Motions: Parametric and Piecewise-Smooth Flow Fields," Computer Vision and Image Understanding Vol. 63, No. 1, 1996, pp. 75-104.
 */
#include <deque>
#include <string>
#include <vector>

#define OUTPUT_IMG_CLASS

#include "../ImgClass/Color.h"
#include "../ImgClass/Vector.h"
#include "../ImgClass/ImgClass.h"
#include "../ImgClass/Segmentation.h"
#include "../ImgClass/BlockMatching.h"
#include "../ImgClass/MotionCompensation.h"

#include "MEstimator.h"
#include "MultiResolution.h"
#include "Affine_BlockMatching.h"
#include "../lib/Struct.h"
#include "../Scratch_MeaningfulMotion.h"




std::vector<ImgVector<Vector_ST<double> > > OpticalFlow_BlockMatching(const ImgVector<ImgClass::RGB>& It, const ImgVector<ImgClass::RGB>& Itp1, double MaxInt, MULTIPLE_MOTION_PARAM MotionParam, const std::string ofilename, const int Mode, const int IterMax = 2048);


ImgVector<VECTOR_2D<double> > OpticalFlow_GradientMethod(const ImgVector<ImgClass::Lab>* reference, const ImgVector<ImgClass::Lab>* interest, const ImgVector<VECTOR_2D<double> >* MV, const ImgVector<size_t>* region_map, const double& lambdaD, const double& lambdaS, const double& sigmaD, const double& sigmaS, const int IterMax, const double& Error_Min_Threshold);

void IRLS_OpticalFlow_GradientMethod(ImgVector<VECTOR_2D<double> >* u, const ImgVector<size_t>* region_map, const ImgVector<VECTOR_2D<double> >* grad, const ImgVector<double>* dt, const double& lambdaD, const double& lambdaS, const double& sigmaD, const double& sigmaS, const int IterMax, const double& Error_Min_Threshold);

VECTOR_2D<double> Error_u_Block(const size_t& site, const ImgVector<VECTOR_2D<double> >* u, const ImgVector<size_t>* domain_map, const ImgVector<VECTOR_2D<double> >* Img_g, const ImgVector<double>* Img_t, const double& lambdaD, const double& lambdaS, const double& sigmaD, const double& sigmaS);

VECTOR_2D<double> sup_Error_uu_Block(const ImgVector<VECTOR_2D<double> >* Img_g, const double& lambdaD, const double& lambdaS, const double& sigmaD, const double& sigmaS);

double Error_MultipleMotion_Block(const ImgVector<VECTOR_2D<double> >* u, const ImgVector<size_t>* domain_map, const ImgVector<VECTOR_2D<double> >* Img_g, const ImgVector<double>* Img_t, const double& lambdaD, const double& lambdaS, const double& sigmaD, const double& sigmaS);


// Output
void MultipleMotion_write(const ImgVector<double>& img_prev, const ImgVector<double>& img_current, const int MaxInt, const std::vector<ImgVector<Vector_ST<double> > >& u, const std::string& filename);

void MultipleMotion_write(const ImgVector<ImgClass::RGB>& img_prev, const ImgVector<ImgClass::RGB>& img_current, const int MaxInt, const std::vector<ImgVector<Vector_ST<double> > >& u, const std::string& filename);


void MultipleMotion_write(const ImgVector<double>& img_prev, const ImgVector<double>& img_current, const ImgVector<double>& img_next, const int MaxInt, const std::vector<ImgVector<Vector_ST<double> > >& u, const std::string& filename);
void MultipleMotion_write(const ImgVector<ImgClass::RGB>& img_prev, const ImgVector<ImgClass::RGB>& img_current, const ImgVector<ImgClass::RGB>& img_next, const int MaxInt, const std::vector<ImgVector<Vector_ST<double> > >& u, const std::string& filename);

