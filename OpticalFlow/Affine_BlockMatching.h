/* Robust Estimation of Multiple Motions
 *
 * This motion estimation based on
 * M.J.Black and P.Anandan, "The Robust Estimation of Multiple Motions: Parametric and Piecewise-Smooth Flow Fields," Computer Vision and Image Understanding Vol. 63, No. 1, 1996, pp. 75-104.
 */

#include <string>

#include "MEstimator.h"
#include "MultiResolution.h"

#include "../lib/Struct.h"
#include "../Scratch_MeaningfulMotion.h"




ImgVector<VECTOR_2D<double> >* OpticalFlow_Affine_BlockMatching(ImgVector<double> *It, ImgVector<double> *Itp1, double MaxInt, MULTIPLE_MOTION_PARAM MotionParam);

bool IRLS_Affine_Block(std::vector<VECTOR_AFFINE>* u, const std::vector<std::vector<VECTOR_2D<int> > >& connected_domains, const ImgVector<VECTOR_2D<double> > *Img_g, const ImgVector<double> *Img_t, double sigmaD, int IterMax, double ErrorMinThreshold);

VECTOR_AFFINE Error_a_Block(const VECTOR_AFFINE& u, const std::vector<VECTOR_2D<int> >& connected_domain, const ImgVector<VECTOR_2D<double> > *Img_g, const ImgVector<double> *Img_t, double sigmaD);
VECTOR_AFFINE sup_Error_aa_Block(const std::vector<VECTOR_2D<int> >& connected_domain, const ImgVector<VECTOR_2D<double> > *Img_g, double sigmaD);

double Error_Affine_Block(const VECTOR_AFFINE& u, const std::vector<VECTOR_2D<int> >& connected_domain, const ImgVector<VECTOR_2D<double> > *Img_g, const ImgVector<double> *Img_t, double sigmaD);

