/* Robust Estimation of Multiple Motions
 *
 * This motion estimation based on
 * M.J.Black and P.Anandan, "The Robust Estimation of Multiple Motions: Parametric and Piecewise-Smooth Flow Fields," Computer Vision and Image Understanding Vol. 63, No. 1, 1996, pp. 75-104.
 */

#include "OpticalFlow.h"

#include "../ImgClass/BlockMatching.h"




ImgVector<VECTOR_2D<double> >* OpticalFlow_BlockMatching(const ImgVector<double> *It, const ImgVector<double> *Itp1, double MaxInt, MULTIPLE_MOTION_PARAM MotionParam, int IterMax = 0);

void BM2OpticalFlow(ImgVector<double>* I_dt_levels, ImgVector<VECTOR_2D<double> >* u_levels, const ImgVector<double>* It_levels, const ImgVector<double>* Itp1_levels, int MaxLevel, const ImgVector<VECTOR_2D<double> >* Motion_Vector);

