/* Robust Estimation of Multiple Motions
 *
 * This motion estimation based on
 * M.J.Black and P.Anandan, "The Robust Estimation of Multiple Motions: Parametric and Piecewise-Smooth Flow Fields," Computer Vision and Image Understanding Vol. 63, No. 1, 1996, pp. 75-104.
 */
#include <string>

#include "../ImgClass/RGB.h"
#include "../ImgClass/Lab.h"
#include "../ImgClass/BlockMatching.h"




ImgVector<VECTOR_2D<double> > *
OpticalFlow_BlockMatching(const ImgVector<ImgClass::RGB<double> >& It, const ImgVector<ImgClass::RGB<double> >& Itp1, double MaxInt, MULTIPLE_MOTION_PARAM MotionParam, const std::string ofilename, int IterMax = 0);

void BM2OpticalFlow(ImgVector<double>* I_dt_levels, ImgVector<VECTOR_2D<double> >* u_levels, const ImgVector<double>* It_levels, const ImgVector<double>* Itp1_levels, const int level, BlockMatching<double>* Motion_Vector);
void Add_VectorOffset(ImgVector<VECTOR_2D<double> >* u_levels, int level, int MaxLevel, BlockMatching<double>* block_matching);

VECTOR_2D<double> Error_u_Block(int site, const ImgVector<VECTOR_2D<double> >* u, const ImgVector<bool>& domain_map, const ImgVector<VECTOR_2D<double> >* Img_g, const ImgVector<double>* Img_t, const double& lambdaD, const double& lambdaS, const double& sigmaD, const double& sigmaS);
VECTOR_2D<double> sup_Error_uu_Block(const ImgVector<VECTOR_2D<double> >* Img_g, const double& lambdaD, const double& lambdaS, const double& sigmaD, const double& sigmaS);

double Error_MultipleMotion_Block(const ImgVector<VECTOR_2D<double> >* u, const ImgVector<bool>& domain_map, const ImgVector<VECTOR_2D<double> >* Img_g, const ImgVector<double>* Img_t, const double& lambdaD, const double& lambdaS, const double& sigmaD, const double& sigmaS);

