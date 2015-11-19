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
bool IRLS_Affine_Block(VECTOR_AFFINE *u, ImgVector<VECTOR_2D<double> > *Img_g, ImgVector<double> *Img_t, double sigmaD, int IterMax, double ErrorMinThreshold);

VECTOR_AFFINE Error_a_Block(VECTOR_AFFINE *u, ImgVector<VECTOR_2D<double> > *Img_g, ImgVector<double> *Img_t, double sigmaD);
VECTOR_AFFINE sup_Error_aa_Block(ImgVector<VECTOR_2D<double> > *Img_g, double sigmaD);

double Error_Affine_Block(const VECTOR_AFFINE *u, ImgVector<VECTOR_2D<double> > *Img_g, ImgVector<double> *Img_t, double sigmaD);

