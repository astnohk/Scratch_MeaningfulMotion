/* Robust Estimation of Multiple Motions
 *
 * This motion estimation based on
 * M.J.Black and P.Anandan, "The Robust Estimation of Multiple Motions: Parametric and Piecewise-Smooth Flow Fields," Computer Vision and Image Understanding Vol. 63, No. 1, 1996, pp. 75-104.
 */

#include "OpticalFlow.h"

#include "../ImgClass/BlockMatching.h"




ImgVector<VECTOR_2D<double> >* OpticalFlow_BlockMatching(const ImgVector<double> *It, const ImgVector<double> *Itp1, double MaxInt, MULTIPLE_MOTION_PARAM MotionParam, int IterMax = 0);

ImgVector<double> OpticalFlow_dt(const ImgVector<double>& It, const ImgVector<double>& Itp1, const ImgVector<VECTOR_2D<double> >& vector);
ImgVector<VECTOR_2D<double> > OpticalFlow_grad(const ImgVector<double>& It);

void IRLS_OpticalFlow(ImgVector<VECTOR_2D<double> >* u, const ImgVector<VECTOR_2D<double> >* Img_g, const ImgVector<double>* Img_t, double lambdaD, double lambdaS, double sigmaD, double sigmaS, int IterMax, double ErrorMinThreshold);

