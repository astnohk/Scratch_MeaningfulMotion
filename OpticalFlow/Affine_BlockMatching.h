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

#include "MEstimator.h"
#include "../lib/Struct.h"
#include "../Scratch_MeaningfulMotion.h"




ImgVector<VECTOR_2D<double> > AffineParametric(const ImgVector<ImgClass::Lab>& reference, const ImgVector<ImgClass::Lab>& interest, const ImgVector<VECTOR_2D<double> >& MV, const std::vector<std::vector<VECTOR_2D<int> > >& regions, const MULTIPLE_MOTION_PARAM& MotionParam, const int IterMax);


void IRLS_AffineParametric_region(std::vector<double>* u_affine, const std::vector<VECTOR_2D<int> >& region, const ImgVector<VECTOR_2D<double> >& grad, const ImgVector<double>& dt, const double& sigma, const int IterMax, const double& ErrorMinThreshold);

std::vector<double> Error_a_region(const std::vector<double>& u_affine, const std::vector<VECTOR_2D<int> >& region, const ImgVector<VECTOR_2D<double> >& grad, const ImgVector<double>& dt, const double& sigma);

std::vector<double> sup_Error_aa_region(const std::vector<VECTOR_2D<int> >& region, const ImgVector<VECTOR_2D<double> >& grad, const double& sigma);

double Error_Affine_region(const std::vector<double>& u_affine, const std::vector<VECTOR_2D<int> >& region, const ImgVector<VECTOR_2D<double> >& grad, const ImgVector<double>& dt, const double& sigma);

