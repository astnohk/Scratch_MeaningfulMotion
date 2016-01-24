/* Robust Estimation of Multiple Motions
 *
 * This motion estimation based on
 * M.J.Black and P.Anandan, "The Robust Estimation of Multiple Motions: Parametric and Piecewise-Smooth Flow Fields," Computer Vision and Image Understanding Vol. 63, No. 1, 1996, pp. 75-104.
 */
#include <deque>
#include <string>
#include <vector>

#define OUTPUT_IMG_CLASS

#include "../ImgClass/RGB.h"
#include "../ImgClass/Lab.h"
#include "../ImgClass/Vector.h"
#include "../ImgClass/ImgClass.h"
#include "../ImgClass/Segmentation.h"
#include "../ImgClass/BlockMatching.h"
#include "../ImgClass/MotionCompensation.h"

#include "MEstimator.h"
#include "MultiResolution.h"
#include "../lib/Struct.h"
#include "../Scratch_MeaningfulMotion.h"




std::vector<ImgVector<Vector_ST<double> > > OpticalFlow_BlockMatching(const ImgVector<ImgClass::RGB>& It, const ImgVector<ImgClass::RGB>& Itp1, double MaxInt, MULTIPLE_MOTION_PARAM MotionParam, const std::string ofilename, int IterMax = 0);

template<class T = ImgClass::Lab> void LevelDown(ImgVector<double> *I_dt_levels, ImgVector<VECTOR_2D<double> > *u_levels, const ImgVector<double> *It_levels, const ImgVector<double> *Itp1_levels, int level, int MaxLevel, BlockMatching<T>* block_matching = nullptr);
template<class T> void Add_VectorOffset(ImgVector<VECTOR_2D<double> >* u_levels, int level, int MaxLevel, BlockMatching<T>* block_matching);

void IRLS_OpticalFlow_Pyramid_Segment(ImgVector<VECTOR_2D<double> >* u, const ImgVector<size_t>& domain_map, const ImgVector<VECTOR_2D<double> >* Img_g, const ImgVector<double>* Img_t, double lambdaD, double lambdaS, double sigmaD, double sigmaS, int IterMax, double ErrorMinThreshold, int level);

VECTOR_2D<double> Error_u_Block(const size_t& site, const ImgVector<VECTOR_2D<double> >* u, const ImgVector<size_t>& domain_map, const ImgVector<VECTOR_2D<double> >* Img_g, const ImgVector<double>* Img_t, const double& lambdaD, const double& lambdaS, const double& sigmaD, const double& sigmaS);
VECTOR_2D<double> sup_Error_uu_Block(const ImgVector<VECTOR_2D<double> >* Img_g, const double& lambdaD, const double& lambdaS, const double& sigmaD, const double& sigmaS);

double Error_MultipleMotion_Block(const ImgVector<VECTOR_2D<double> >* u, const ImgVector<size_t>& domain_map, const ImgVector<VECTOR_2D<double> >* Img_g, const ImgVector<double>* Img_t, const double& lambdaD, const double& lambdaS, const double& sigmaD, const double& sigmaS);

// Output
void MultipleMotion_write(const ImgVector<double>& img_prev, const ImgVector<double>& img_current, const std::vector<ImgVector<Vector_ST<double> > >& u, const std::string& filename);
void MultipleMotion_write(const ImgVector<ImgClass::RGB>& img_prev, const ImgVector<ImgClass::RGB>& img_current, const std::vector<ImgVector<Vector_ST<double> > >& u, const std::string& filename);

void MultipleMotion_write(const ImgVector<double>& img_prev, const ImgVector<double>& img_current, const ImgVector<double>& img_next, const std::vector<ImgVector<Vector_ST<double> > >& u, const std::string& filename);
void MultipleMotion_write(const ImgVector<ImgClass::RGB>& img_prev, const ImgVector<ImgClass::RGB>& img_current, const ImgVector<ImgClass::RGB>& img_next, const std::vector<ImgVector<Vector_ST<double> > >& u, const std::string& filename);

