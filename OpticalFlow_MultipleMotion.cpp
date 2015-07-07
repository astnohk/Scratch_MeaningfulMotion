/* Robust Estimation of Multiple Motions
 *
 * This motion estimation based on
 * M.J.Black and P.Anandan, "The Robust Estimation of Multiple Motions: Parametric and Piecewise-Smooth Flow Fields," Computer Vision and Image Understanding, Vol.63, No.1, 1996, pp.75-104.
 */

#include "Scratch_MeaningfulMotion.h"
#include "MultiResolution.h"
#include "OpticalFlow_MultipleMotion.h"


// IRLS parameters
#define IRLS_ITER_MAX 512


// define for debug

#define SHOW_IRLS_MULTIPLEMOTION_OPTICALFLOW_E
//#define DEBUG_STOP_ON_LEVEL_L

// /define for debug




VECTOR_2D *
MultipleMotion_OpticalFlow(double *It, double *Itp1, double MaxInt, SIZE size_img, MULTIPLE_MOTION_PARAM MotionParam)
{
	ERROR Error("MultipleMotion_OpticalFlow");

	// M-estimator parameter
	const double lambdaD = 5.0;
	const double lambdaS = 1.0;
	double sigmaD;
	const double sigmaD_init = 0.8 / sqrt(2.0); //18.0 / sqrt(2.0);
	const double sigmaD_l0 = 0.2 / sqrt(2.0); //4.0 / sqrt(2.0);
	double sigmaS;
	const double sigmaS_init = 0.3 / sqrt(2.0); //3.0 / sqrt(2.0);
	const double sigmaS_l0 = 0.03 / sqrt(2.0);

	double *It_normalize = nullptr;
	double *Itp1_normalize = nullptr;
	VECTOR_2D *u = nullptr; // For RETURN value
	VECTOR_2D **u_levels = nullptr;
	double **I_dt_levels = nullptr;
	double **It_levels = nullptr;
	double **Itp1_levels = nullptr;
	VECTOR_2D **grad_It_levels = nullptr;
	SIZE size_img_l;
	SIZE size_img_lp1;
	int level;
	int IterMax;
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
	try {
		It_normalize = new double[size_img.width * size_img.height];
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("It_normalize");
		Error.Malloc();
		goto ExitError;
	}
	try {
		Itp1_normalize = new double[size_img.width * size_img.height];
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("Itp1_normalize");
		Error.Malloc();
		goto ExitError;
	}
	for (i = 0; i < size_img.width * size_img.height; i++) {
		It_normalize[i] = It[i] / MaxInt;
		Itp1_normalize[i] = Itp1[i] / MaxInt;
	}
	// Multiple Motion Vectors
	try {
		u = new VECTOR_2D[size_img.width * size_img.height];
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("u");
		Error.Malloc();
		goto ExitError;
	}
	try {
		u_levels = new VECTOR_2D*[MotionParam.Level];
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("u_levels");
		Error.Malloc();
		goto ExitError;
	}
	// Make Pyramid
	if ((It_levels = Pyramider(It_normalize, size_img, MotionParam.Level)) == nullptr) {
		Error.Function("Pyramider");
		Error.Value("It_levels");
		Error.FunctionFail();
		goto ExitError;
	}
	if ((Itp1_levels = Pyramider(Itp1_normalize, size_img, MotionParam.Level)) == nullptr) {
		Error.Function("Pyramider");
		Error.Value("Itp1_levels");
		Error.FunctionFail();
		goto ExitError;
	}
	// Derivative about time
	if ((I_dt_levels = dt_Pyramid(It_levels, Itp1_levels, size_img, MotionParam.Level)) == nullptr) {
		Error.Function("dt_Pyramid");
		Error.Value("I_dt_levels");
		Error.FunctionFail();
		goto ExitError;
	}
	// Derivative about space
	if ((grad_It_levels = grad_Pyramid(It_levels, Itp1_levels, size_img, MotionParam.Level)) == nullptr) {
		Error.Function("grad_Pyramid");
		Error.Value("grad_It_levels");
		Error.FunctionFail();
		goto ExitError;
	}

	for (level = MotionParam.Level - 1; level >= 0; level--) {
		if (MotionParam.Level > 1) {
			sigmaD = sigmaD_init + (sigmaD_l0 - sigmaD_init) / (MotionParam.Level - 1.0) * (MotionParam.Level - 1.0 - level);
			sigmaS = sigmaS_init + (sigmaS_l0 - sigmaS_init) / (MotionParam.Level - 1.0) * (MotionParam.Level - 1.0 - level);
		} else {
			sigmaD = sigmaD_l0;
			sigmaS = sigmaS_l0;
		}
		size_img_lp1.width = (int)ceil(size_img.width * pow_int(0.5, level + 1));
		size_img_lp1.height = (int)ceil(size_img.height * pow_int(0.5, level + 1));
		size_img_l.width = (int)ceil(size_img.width * pow_int(0.5, level));
		size_img_l.height = (int)ceil(size_img.height * pow_int(0.5, level));
		printf("\nLevel %d : (1 / %d scaled, %dx%d)\n  sigmaD = %f\n  sigmaS = %f\n", level, (int)pow_int(2.0, level), size_img_l.width, size_img_l.height, sigmaD, sigmaS);
		try {
			u_levels[level] = new VECTOR_2D[size_img_l.width * size_img_l.height];
		}
		catch (const std::bad_alloc &bad) {
			Error.Value("u_levels[level]");
			Error.Malloc();
			goto ExitError;
		}
		if (level < MotionParam.Level - 1) {
			LevelDown(u_levels[level], size_img_l, u_levels[level + 1], size_img_lp1);
		}
#ifdef DEBUG_STOP_ON_LEVEL_L
		if (level <= 2) {
			continue;
		}
#endif
		IterMax = 10 * MAX(size_img_l.width, size_img_l.height);
		IRLS_MultipleMotion_OpticalFlow(u_levels[level], grad_It_levels[level], I_dt_levels[level], size_img_l,
		    lambdaD, lambdaS, sigmaD, sigmaS,
		    IterMax, MotionParam.Error_Min_Threshold,
		    level);
	}
	for (i = 0; i < size_img.width * size_img.height; i++) {
		u[i] = u_levels[0][i];
	}
	for (level = 0; level < MotionParam.Level; level++) {
		delete[] u_levels[level];
		delete[] grad_It_levels[level];
		delete[] I_dt_levels[level];
		delete[] Itp1_levels[level];
		delete[] It_levels[level];
	}
	delete[] u_levels;
	delete[] grad_It_levels;
	delete[] I_dt_levels;
	delete[] Itp1_levels;
	delete[] It_levels;
	delete[] Itp1_normalize;
	delete[] It_normalize;
	return u;
// Error
ExitError:
	for (level = 0; level < MotionParam.Level; level++) {
		delete[] u_levels[level];
		delete[] grad_It_levels[level];
		delete[] I_dt_levels[level];
		delete[] Itp1_levels[level];
		delete[] It_levels[level];
	}
	delete[] grad_It_levels;
	delete[] I_dt_levels;
	delete[] Itp1_levels;
	delete[] It_levels;
	delete[] u_levels;
	delete[] u;
	delete[] Itp1_normalize;
	delete[] It_normalize;
	return nullptr;
}


void
LevelDown(VECTOR_2D *u_l, SIZE size_l, VECTOR_2D *u_lp1, const SIZE &size_lp1)
{
	int x, y;

	for (y = 0; y < size_l.height; y++) {
		for (x = 0; x < size_l.width; x++) {
			u_l[size_l.width * y + x] = u_lp1[size_lp1.width * (y / 2) + (x / 2)];
		}
	}
}


int
IRLS_MultipleMotion_OpticalFlow(VECTOR_2D *u, VECTOR_2D *Img_g, double *Img_t, SIZE size_img, double lambdaD, double lambdaS, double sigmaD, double sigmaS, int IterMax, double ErrorMinThreshold, int level)
{
	ERROR Error("IRLS_MultipleMotion_OpticalFlow");
	VECTOR_2D *u_np1 = nullptr;
	VECTOR_2D sup;
	VECTOR_2D dE;
	double E = 0.0;
	double E_prev = 0.0;
	int ErrorIncrementCount = 0;
	int site;
	int n;

	try {
		u_np1 = new VECTOR_2D[size_img.width * size_img.height];
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("u_np1");
		Error.Malloc();
		goto ExitError;
	}
	// Reset sup_Error_uu max Img_g
	sup_Error_uu(Img_g, size_img, lambdaD, lambdaS, sigmaD, sigmaS);
	sup = sup_Error_uu(nullptr, size_img, lambdaD, lambdaS, sigmaD, sigmaS);
	for (n = 0; n < IterMax; n++) {
#pragma omp parallel for private(dE, sup)
		// Calc for all sites
		for (site = 0; site < size_img.width * size_img.height; site++) {
			dE = Error_u(site, u, Img_g, Img_t, size_img, lambdaD, lambdaS, sigmaD, sigmaS);
			sup = sup_Error_uu(nullptr, size_img, lambdaD, lambdaS, sigmaD, sigmaS);
			u_np1[site].x = u[site].x - dE.x / sup.x;
			u_np1[site].y = u[site].y - dE.y / sup.y;
		}
		// Calc for all sites
		for (site = 0; site < size_img.width * size_img.height; site++) {
			u[site] = u_np1[site];
		}
		if (level == 0) {
			if ((n & 0x3F) == 0) {
				E = Error_MultipleMotion(u, Img_g, Img_t, size_img, lambdaD, lambdaS, sigmaD, sigmaS);
			}
		} else {
			E_prev = E;
			E = Error_MultipleMotion(u, Img_g, Img_t, size_img, lambdaD, lambdaS, sigmaD, sigmaS);
			if (E > E_prev) {
				ErrorIncrementCount++;
			} else {
				ErrorIncrementCount = 0;
			}
		}
#ifdef SHOW_IRLS_MULTIPLEMOTION_OPTICALFLOW_E
		if ((n & 0x3F) == 0) {
			printf("E(%4d) = %e\n", n, E);
		}
#endif
		if (E < ErrorMinThreshold || ErrorIncrementCount > 2) {
			break;
		}
	}
	delete[] u_np1;
	return MEANINGFUL_SUCCESS;
// Error
ExitError:
	delete[] u_np1;
	return MEANINGFUL_FAILURE;
}


VECTOR_2D
Error_u(int site, VECTOR_2D *u, VECTOR_2D *Img_g, double *Img_t, const SIZE &size_img, const double &lambdaD, const double &lambdaS, const double &sigmaD, const double &sigmaS)
{
	double (*psiD)(const double&, const double&) = Geman_McClure_psi;
	double (*psiS)(const double&, const double&) = Geman_McClure_psi;
	VECTOR_2D us;
	double Center;
	VECTOR_2D Neighbor;
	VECTOR_2D E_u;
	int x, y;

	x = site % size_img.width;
	y = (int)floor(site / size_img.width);

	us = u[site];
	Center = (*psiD)(Img_g[site].x * us.x + Img_g[site].y * us.y + Img_t[site], sigmaD);

	Neighbor.x = .0;
	Neighbor.y = .0;
	if (x > 0) {
		Neighbor.x += (*psiS)(us.x - u[size_img.width * y + x - 1].x, sigmaS);
		Neighbor.y += (*psiS)(us.y - u[size_img.width * y + x - 1].y, sigmaS);
	}
	if (x < size_img.width - 1) {
		Neighbor.x += (*psiS)(us.x - u[size_img.width * y + x + 1].x, sigmaS);
		Neighbor.y += (*psiS)(us.y - u[size_img.width * y + x + 1].y, sigmaS);
	}
	if (y > 0) {
		Neighbor.x += (*psiS)(us.x - u[size_img.width * (y - 1) + x].x, sigmaS);
		Neighbor.y += (*psiS)(us.y - u[size_img.width * (y - 1) + x].y, sigmaS);
	}
	if (y < size_img.height - 1) {
		Neighbor.x += (*psiS)(us.x - u[size_img.width * (y + 1) + x].x, sigmaS);
		Neighbor.y += (*psiS)(us.y - u[size_img.width * (y + 1) + x].y, sigmaS);
	}

	E_u.x += lambdaD * Img_g[site].x * Center + lambdaS * Neighbor.x;
	E_u.y += lambdaD * Img_g[site].y * Center + lambdaS * Neighbor.y;
	return E_u;
}


VECTOR_2D
sup_Error_uu(VECTOR_2D *Img_g, const SIZE &size, const double &lambdaD, const double &lambdaS, const double &sigmaD, const double &sigmaS)
{
	static VECTOR_2D Img_g_max;
	VECTOR_2D sup;

	if (Img_g != nullptr) {
		Img_g_max.reset();
		for (int i = 0; i < size.width * size.height; i++) {
			if (Img_g_max.x < POW2(Img_g[i].x)) {
				Img_g_max.x = POW2(Img_g[i].x);
			}
			if (Img_g_max.y < POW2(Img_g[i].y)) {
				Img_g_max.y = POW2(Img_g[i].y);
			}
		}
	}
	sup.x = lambdaD * Img_g_max.x / POW2(sigmaD) + 4.0 * lambdaS / POW2(sigmaS);
	sup.y = lambdaD * Img_g_max.y / POW2(sigmaD) + 4.0 * lambdaS / POW2(sigmaS);
	return sup;
}


double
Error_MultipleMotion(VECTOR_2D *u, VECTOR_2D *Img_g, double *Img_t, const SIZE &size_img, const double &lambdaD, const double &lambdaS, const double &sigmaD, const double &sigmaS)
{
	double (*rhoD)(const double&, const double&) = Geman_McClure_rho;
	double (*rhoS)(const double&, const double&) = Geman_McClure_rho;
	VECTOR_2D us;
	double Center;
	VECTOR_2D Neighbor;
	double E = 0.0;
	int x, y;

#pragma omp parallel for private(x, us, Neighbor, Center) reduction(+:E)
	for (y = 0; y < size_img.height; y++) {
		for (x = 0; x < size_img.width; x++) {
			us = u[size_img.width * y + x];
			Neighbor.x = .0;
			Neighbor.y = .0;
			if (x > 0) {
				Neighbor.x += (*rhoS)(us.x - u[size_img.width * y + x - 1].x, sigmaS);
				Neighbor.y += (*rhoS)(us.y - u[size_img.width * y + x - 1].y, sigmaS);
			}
			if (x < size_img.width - 1) {
				Neighbor.x += (*rhoS)(us.x - u[size_img.width * y + x + 1].x, sigmaS);
				Neighbor.y += (*rhoS)(us.y - u[size_img.width * y + x + 1].y, sigmaS);
			}
			if (y > 0) {
				Neighbor.x += (*rhoS)(us.x - u[size_img.width * (y - 1) + x].x, sigmaS);
				Neighbor.y += (*rhoS)(us.y - u[size_img.width * (y - 1) + x].y, sigmaS);
			}
			if (y < size_img.height - 1) {
				Neighbor.x += (*rhoS)(us.x - u[size_img.width * (y + 1) + x].x, sigmaS);
				Neighbor.y += (*rhoS)(us.y - u[size_img.width * (y + 1) + x].y, sigmaS);
			}
			Center = (*rhoD)(Img_g[size_img.width * y + x].x * us.x
			    + Img_g[size_img.width * y + x].y * us.y
			    + Img_t[size_img.width * y + x],
			    sigmaD);
			E += lambdaD * Center + lambdaS * (Neighbor.x + Neighbor.y);
		}
	}
	return E;
}


int
MultipleMotion_write(VECTOR_2D *u, SIZE size, const char *filename)
{
	ERROR Error("MultipleMotion_write");
	FILE *fp = nullptr;
	int x, y;

	if (u == nullptr) {
		Error.Value("u");
		Error.PointerNull();
		goto ExitError;
	} else if (filename == nullptr) {
		Error.Value("filename");
		Error.PointerNull();
		goto ExitError;
	}

	printf("* Output Optical Flow to '%s'(binary)\n", filename);
	if ((fp = fopen(filename, "wb")) == nullptr) {
		Error.Function("fopen");
		Error.File(filename);
		Error.FileWrite();
		goto ErrorFileOpenFail;
	}
	fprintf(fp, "%d %d\n", size.width, size.height);
	for (y = 0; y < size.height; y++) {
		for (x = 0; x < size.width; x++) {
			if (fwrite(&(u[size.width * y + x].x), sizeof(double), 1, fp) < 1) {
				Error.Function("fwrite");
				Error.Value("u(x, y).x");
				Error.FunctionFail();
				goto ExitError;
			}
			if (fwrite(&(u[size.width * y + x].y), sizeof(double), 1, fp) < 1) {
				Error.Function("fwrite");
				Error.Value("u(x, y).y");
				Error.FunctionFail();
				goto ExitError;
			}
		}
	}
	fclose(fp);
	return MEANINGFUL_SUCCESS;
// Error
ErrorFileOpenFail:
	return MEANINGFUL_FAILURE;
ExitError:
	fclose(fp);
	return MEANINGFUL_FAILURE;
}

