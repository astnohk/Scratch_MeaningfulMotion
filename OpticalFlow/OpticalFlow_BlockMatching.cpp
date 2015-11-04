/* Robust Estimation of Multiple Motions
 *
 * This motion estimation based on
 * M.J.Black and P.Anandan, "The Robust Estimation of Multiple Motions: Parametric and Piecewise-Smooth Flow Fields," Computer Vision and Image Understanding, Vol.63, No.1, 1996, pp.75-104.
 */

#include "OpticalFlow_BlockMatching.h"


#define SHOW_IRLS_OPTICALFLOW_E




ImgVector<VECTOR_2D<double> > *
OpticalFlow_BlockMatching(const ImgVector<double>* It, const ImgVector<double>* Itp1, double MaxInt, MULTIPLE_MOTION_PARAM MotionParam, int IterMax)
{
	ImgVector<VECTOR_2D<double> > *u = nullptr; // For RETURN value
	ImgVector<double> It_normalize;
	ImgVector<double> Itp1_normalize;
	ImgVector<double> I_dt;
	ImgVector<VECTOR_2D<double> > grad_It;

	BlockMatching<double> block_matching;
	ImgVector<VECTOR_2D<double> > Motion_Vector;

	// M-estimator parameter
	const double lambdaD = 5.0;
	const double lambdaS = 1.0;
	const double sigmaD_l0 = 0.2 / sqrt(2.0); //4.0 / sqrt(2.0);
	const double sigmaS_l0 = 0.03 / sqrt(2.0);

	if (It == nullptr) {
		throw std::invalid_argument("ImgVector<double>* It");
	} else if (Itp1 == nullptr) {
		throw std::invalid_argument("ImgVector<double>* Itp1");
	}

	// ----- Block Matching -----
	block_matching.reset(It, Itp1, MotionParam.BlockMatching_BlockSize);
	Motion_Vector.copy(block_matching.data());

	// ----- Optical Flow -----
	// Image Normalization
	It_normalize = *It;
	Itp1_normalize = *Itp1;
	for (int i = 0; i < It_normalize.size(); i++) {
		It_normalize[i] /= MaxInt;
		Itp1_normalize[i] /= MaxInt;
	}
	// Multiple Motion Vectors
	try {
		u = new ImgVector<VECTOR_2D<double> >(It->width(), It->height());
	}
	catch (const std::bad_alloc &bad) {
		throw;
	}
	// Derivative about time
	I_dt = OpticalFlow_dt(It_normalize, Itp1_normalize, Motion_Vector);
	grad_It = OpticalFlow_grad(It_normalize);

	u->reset(It->width(), It->height());
	printf("\n(%dx%d)\n  sigmaD = %f\n  sigmaS = %f\n", u->width(), u->height(), sigmaD_l0, sigmaS_l0);
	if (IterMax < 0) {
		IterMax = 10 * MAX(It->width(), It->height());
	}
	printf("IterMax = %d\n", IterMax);
	IRLS_OpticalFlow(
	    u,
	    &grad_It,
	    &I_dt,
	    lambdaD, lambdaS, sigmaD_l0, sigmaS_l0,
	    IterMax,
	    MotionParam.Error_Min_Threshold);
	// Add offset of the motion vector computed by Block Matching
	for (int i = 0; i < u->size(); i++) {
		(*u)[i] += Motion_Vector[i];
	}

	return u;
}


ImgVector<double>
OpticalFlow_dt(const ImgVector<double>& It, const ImgVector<double>& Itp1, const ImgVector<VECTOR_2D<double> >& vector)
{
	ImgVector<double> I_dt;
	int block_width = (int)ceil(It.width() / vector.width());
	int block_height = (int)ceil(It.height() / vector.height());

	I_dt.reset(It.width(), It.height());

	for (int y = 0; y < I_dt.height(); y++) {
		for (int x = 0; x < I_dt.width(); x++) {
			VECTOR_2D<double> v = vector.get(x / block_width, y / block_height);
			I_dt.ref(x, y) =
			    (Itp1.get_zeropad(x + v.x, y + v.y) - It.get_zeropad(x, y)
			    + Itp1.get_zeropad(x + 1 + v.x, y + v.y) - It.get_zeropad(x + 1, y)
			    + Itp1.get_zeropad(x + v.x, y + 1 + v.y) - It.get_zeropad(x, y + 1)
			    + Itp1.get_zeropad(x + 1 + v.x, y + 1 + v.y) - It.get_zeropad(x + 1, y + 1)) / 4.0;
		}
	}
	return I_dt;
}

ImgVector<VECTOR_2D<double> >
OpticalFlow_grad(const ImgVector<double>& It)
{
	ImgVector<VECTOR_2D<double> > grad_It;

	grad_It.reset(It.width(), It.height());

	for (int y = 0; y < grad_It.height(); y++) {
		for (int x = 0; x < grad_It.width(); x++) {
			grad_It.ref(x, y).x =
			    (It.get_zeropad(x + 1, y) - It.get_zeropad(x, y)
			    + It.get_zeropad(x + 1, y + 1) - It.get_zeropad(x, y + 1)) / 2.0;
			grad_It.ref(x, y).y =
			    (It.get_zeropad(x, y + 1) - It.get_zeropad(x, y)
			    + It.get_zeropad(x + 1, y + 1) - It.get_zeropad(x + 1, y)) / 2.0;
		}
	}
	return grad_It;
}


void
IRLS_OpticalFlow(ImgVector<VECTOR_2D<double> >* u, const ImgVector<VECTOR_2D<double> >* Img_g, const ImgVector<double>* Img_t, double lambdaD, double lambdaS, double sigmaD, double sigmaS, int IterMax, double ErrorMinThreshold)
{
	ImgVector<VECTOR_2D<double> > u_np1;
	VECTOR_2D<double> sup;
	VECTOR_2D<double> dE;
	double E = 0.0;
	double E_prev = 0.0;
	int ErrorIncrementCount = 0;
	int site;
	int n;

	if (u == nullptr) {
		throw std::invalid_argument("ImgVector<VECTOR_2D<double> > *u");
	}
	u_np1.copy(u); // Initialize u_np1
	// Reset sup_Error_uu max Img_g
	sup_Error_uu(Img_g, lambdaD, lambdaS, sigmaD, sigmaS);
	sup = sup_Error_uu(nullptr, lambdaD, lambdaS, sigmaD, sigmaS);
	for (n = 0; n < IterMax; n++) {
		// Calc for all sites
#pragma omp parallel for private(dE)
		for (site = 0; site < u->size(); site++) {
			dE = Error_u(site, u, Img_g, Img_t, lambdaD, lambdaS, sigmaD, sigmaS);
			u_np1[site].x = u->get(site).x - dE.x / sup.x;
			u_np1[site].y = u->get(site).y - dE.y / sup.y;
		}
		// Calc for all sites
		for (site = 0; site < u->size(); site++) {
			(*u)[site] = u_np1[site];
		}
		E_prev = E;
		E = Error_MultipleMotion(u, Img_g, Img_t, lambdaD, lambdaS, sigmaD, sigmaS);
		if (E > E_prev) {
			ErrorIncrementCount++;
		} else {
			ErrorIncrementCount = 0;
		}
#ifdef SHOW_IRLS_OPTICALFLOW_E
		if ((n & 0x3F) == 0) {
			printf("E(%4d) = %e\n", n, E);
		}
#endif
		if (E < ErrorMinThreshold || ErrorIncrementCount > 3) {
			break;
		}
	}
}

