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

ImgVector<VECTOR_2D<double> > *
OpticalFlow_BlockMatching(const ImgVector<double>* It, const ImgVector<double>* Itp1, double MaxInt, MULTIPLE_MOTION_PARAM MotionParam, int IterMax)
{
	ImgVector<VECTOR_2D<double> > *u = nullptr; // For RETURN value

	BlockMatching block_matching;
	ImgVector<VECTOR_2D<double> > Motion_Vector;

	ImgVector<double> It_normalize;
	ImgVector<double> Itp1_normalize;
	ImgVector<VECTOR_2D<double> > *u_levels = nullptr;
	ImgVector<double> *I_dt_levels = nullptr;
	ImgVector<double> *It_levels = nullptr;
	ImgVector<double> *Itp1_levels = nullptr;
	ImgVector<VECTOR_2D<double> > *grad_It_levels = nullptr;
	// M-estimator parameter
	const double lambdaD = 5.0;
	const double lambdaS = 1.0;
	double sigmaD;
	const double sigmaD_init = 0.8 / sqrt(2.0); //18.0 / sqrt(2.0);
	const double sigmaD_l0 = 0.2 / sqrt(2.0); //4.0 / sqrt(2.0);
	double sigmaS;
	const double sigmaS_init = 0.3 / sqrt(2.0); //3.0 / sqrt(2.0);
	const double sigmaS_l0 = 0.03 / sqrt(2.0);

	int IterMax_level = 0;
	int MaxLevel = MotionParam.Level;
	int level;
	int i;

	if (It == nullptr) {
		Error.Value("It");
		Error.PointerNull();
		goto ExitError;
	} else if (Itp1 == nullptr) {
		Error.Value("Itp1");
		Error.PointerNull();
		goto ExitError;
	}

	// Image Normalization
	It_normalize = *It;
	Itp1_normalize = *Itp1;
	for (i = 0; i < It_normalize.size(); i++) {
		It_normalize[i] /= MaxInt;
		Itp1_normalize[i] /= MaxInt;
	}

	// Adjust max level to use the Block Matching efficiently
	if (MaxLevel > floor(log((double)MotionParm.BlockMatching_BlockSize) / log(2.0)) {
		MaxLevel = (int)floor(log((double)MotionParm.BlockMatching_BlockSize) / log(2.0);
	}

	// ----- Block Matching -----
	block_matching.reset(It, Itp1, MotionParam.BlockMatching_BlockSize);
	Motion_Vector = block_matching.data();

	// ----- Optical Flow -----
	try {
		u = new ImgVector<VECTOR_2D<double> >(It->width(), It->height());
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("u");
		Error.Malloc();
		goto ExitError;
	}
	try {
		u_levels = new ImgVector<VECTOR_2D<double> >[MaxLevel];
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("u_levels");
		Error.Malloc();
		goto ExitError;
	}
	// Make Pyramid
	if ((It_levels = Pyramider(&It_normalize, MaxLevel)) == nullptr) {
		Error.Function("Pyramider");
		Error.Value("It_levels");
		Error.FunctionFail();
		goto ExitError;
	}
	if ((Itp1_levels = Pyramider(&Itp1_normalize, MaxLevel)) == nullptr) {
		Error.Function("Pyramider");
		Error.Value("Itp1_levels");
		Error.FunctionFail();
		goto ExitError;
	}
	// Derivative about time
	if ((I_dt_levels = dt_Pyramid(It_levels, Itp1_levels, MaxLevel)) == nullptr) {
		Error.Function("dt_Pyramid");
		Error.Value("I_dt_levels");
		Error.FunctionFail();
		goto ExitError;
	}
	// Derivative about space
	if ((grad_It_levels = grad_Pyramid(It_levels, nullptr, MaxLevel)) == nullptr) {
		Error.Function("grad_Pyramid");
		Error.Value("grad_It_levels");
		Error.FunctionFail();
		goto ExitError;
	}

	for (level = MaxLevel - 1; level >= 0; level--) {
		if (MaxLevel > 1) {
			sigmaD = sigmaD_init + (sigmaD_l0 - sigmaD_init) / (MaxLevel - 1.0) * (MaxLevel - 1.0 - level);
			sigmaS = sigmaS_init + (sigmaS_l0 - sigmaS_init) / (MaxLevel - 1.0) * (MaxLevel - 1.0 - level);
		} else {
			sigmaD = sigmaD_l0;
			sigmaS = sigmaS_l0;
		}
		u_levels[level].reset(I_dt_levels[level].width(), I_dt_levels[level].height());
		printf("\nLevel %d : (1 / %d scaled, %dx%d)\n  sigmaD = %f\n  sigmaS = %f\n", level, (int)pow_int(2.0, level), u_levels[level].width(), u_levels[level].height(), sigmaD, sigmaS);
		if (level == MaxLevel - 1) {
			BM2OpticalFlow(I_dt_levels, u_levels, It_levels, Itp1_levels, MaxLevel, &Motion_Vector);
		} else {
			LevelDown(I_dt_levels, u_levels, It_levels, Itp1_levels, level, MaxLevel);
		}
#ifdef DEBUG_STOP_ON_LEVEL_L
		if (level <= 2) {
			continue;
		}
#endif
		IterMax_level = (level + 1) * 10 * MAX(It->width(), It->height());
		if (IterMax < 0 && IterMax_level >= IterMax) {
			IterMax_level = IterMax;
		}
		printf("IterMax = %d\n", IterMax_level);
		IRLS_OpticalFlow_Pyramid(
		    (u_levels + level),
		    (grad_It_levels + level),
		    (I_dt_levels + level),
		    lambdaD, lambdaS, sigmaD, sigmaS,
		    IterMax_level,
		    MotionParam.Error_Min_Threshold,
		    level);
	}
	if (MaxLevel > 1) {
		for (int y = 0; y < u_levels[0].height(); y++) {
			for (int x = 0; x < u_levels[0].width(); x++) {
				u->ref(x, y).x =
				    u_levels[0].get(x, y).x
				    + u_levels[1].get(x / 2, y / 2).x * 2.0;
				u->ref(x, y).y =
				    u_levels[0].get(x, y).y
				    + u_levels[1].get(x / 2, y / 2).y * 2.0;
			}
		}
	} else {
		for (i = 0; i < u->size(); i++) {
			(*u)[i].x = u_levels[0][i].x;
			(*u)[i].y = u_levels[0][i].y;
		}
	}
	delete[] u_levels;
	delete[] grad_It_levels;
	delete[] I_dt_levels;
	delete[] Itp1_levels;
	delete[] It_levels;
	return u;
// Error
ExitError:
	delete[] grad_It_levels;
	delete[] I_dt_levels;
	delete[] Itp1_levels;
	delete[] It_levels;
	delete[] u_levels;
	delete[] u;
	return nullptr;
}

void
BM2OpticalFlow(ImgVector<double>* I_dt_levels, ImgVector<VECTOR_2D<double> >* u_levels, const ImgVector<double>* It_levels, const ImgVector<double>* Itp1_levels, int MaxLevel, const ImgVector<VECTOR_2D<double> >* Motion_Vector)
{
	int level = MaxLevel - 1;
	int W, H;

	W = (int)ceil((double)u_levels[level].width() / Motion_Vector.width());
	H = (int)ceil((double)u_levels[level].height() / Motion_Vector.height());

	for (int y = 0; y < u_levels[level].height(); y++) {
		for (int x = 0; x < u_levels[level].width(); x++) {
			VECTOR_2D<double> u_offset = Motion_Vector.get(x / W, y / W);

			I_dt_levels[level].ref(x, y) =
			    (Itp1_levels[level].get_zeropad(x + (int)floor(2.0 * u_offset.x), y + (int)floor(2.0 * u_offset.y))
			    - It_levels[level].get_zeropad(x, y)
			    + Itp1_levels[level].get_zeropad(x + 1 + (int)floor(2.0 * u_offset.x), y + (int)floor(2.0 * u_offset.y))
			    - It_levels[level].get_zeropad(x + 1, y)
			    + Itp1_levels[level].get_zeropad(x + (int)floor(2.0 * u_offset.x), y + 1 + (int)floor(2.0 * u_offset.y))
			    - It_levels[level].get_zeropad(x, y + 1)
			    + Itp1_levels[level].get_zeropad(x + 1 + (int)floor(2.0 * u_offset.x), y + 1 + (int)floor(2.0 * u_offset.y))
			    - It_levels[level].get_zeropad(x + 1, y + 1)) / 4.0;
			u_levels[level].ref(x, y).x = 0.0;
			u_levels[level].ref(x, y).y = 0.0;
		}
	}
}

