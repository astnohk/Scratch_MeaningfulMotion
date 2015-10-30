/* Robust Estimation of Multiple Motions
 *
 * This motion estimation based on
 * M.J.Black and P.Anandan, "The Robust Estimation of Multiple Motions: Parametric and Piecewise-Smooth Flow Fields," Computer Vision and Image Understanding, Vol.63, No.1, 1996, pp.75-104.
 */

#include "../Scratch_MeaningfulMotion.h"
#include "OpticalFlow_MultipleMotion.h"


// IRLS parameters
#define IRLS_ITER_MAX 512


// define for debug

#define SHOW_IRLS_MULTIPLEMOTION_OPTICALFLOW_E
//#define DEBUG_STOP_ON_LEVEL_L

// /define for debug




ImgVector<VECTOR_2D> *
MultipleMotion_OpticalFlow(ImgVector<double> *It, ImgVector<double> *Itp1, double MaxInt, MULTIPLE_MOTION_PARAM MotionParam, int IterMax)
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

	ImgVector<double> It_normalize;
	ImgVector<double> Itp1_normalize;
	ImgVector<VECTOR_2D> *u = nullptr; // For RETURN value
	ImgVector<VECTOR_2D> *u_levels = nullptr;
	ImgVector<double> *I_dt_levels = nullptr;
	ImgVector<double> *It_levels = nullptr;
	ImgVector<double> *Itp1_levels = nullptr;
	ImgVector<VECTOR_2D> *grad_It_levels = nullptr;
	int IterMax_level = 0;
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
	// Multiple Motion Vectors
	try {
		u = new ImgVector<VECTOR_2D>(It->width(), It->height());
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("u");
		Error.Malloc();
		goto ExitError;
	}
	try {
		u_levels = new ImgVector<VECTOR_2D>[MotionParam.Level];
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("u_levels");
		Error.Malloc();
		goto ExitError;
	}
	// Make Pyramid
	if ((It_levels = Pyramider(&It_normalize, MotionParam.Level)) == nullptr) {
		Error.Function("Pyramider");
		Error.Value("It_levels");
		Error.FunctionFail();
		goto ExitError;
	}
	if ((Itp1_levels = Pyramider(&Itp1_normalize, MotionParam.Level)) == nullptr) {
		Error.Function("Pyramider");
		Error.Value("Itp1_levels");
		Error.FunctionFail();
		goto ExitError;
	}
	// Derivative about time
	if ((I_dt_levels = dt_Pyramid(It_levels, Itp1_levels, MotionParam.Level)) == nullptr) {
		Error.Function("dt_Pyramid");
		Error.Value("I_dt_levels");
		Error.FunctionFail();
		goto ExitError;
	}
	// Derivative about space
	if ((grad_It_levels = grad_Pyramid(It_levels, Itp1_levels, MotionParam.Level)) == nullptr) {
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
		u_levels[level].reset(I_dt_levels[level].width(), I_dt_levels[level].height());
		printf("\nLevel %d : (1 / %d scaled, %dx%d)\n  sigmaD = %f\n  sigmaS = %f\n", level, (int)pow_int(2.0, level), u_levels[level].width(), u_levels[level].height(), sigmaD, sigmaS);
		if (level < MotionParam.Level - 1) {
			LevelDown(u_levels, level);
		}
#ifdef DEBUG_STOP_ON_LEVEL_L
		if (level <= 2) {
			continue;
		}
#endif
//		IterMax_level = MAX(10 * MAX(u_levels[level].width(), u_levels[level].height()), 1000);
		IterMax_level = level * 10 * MAX(It->width(), It->height());
		if (level >= MotionParam.Level - 1) {
		if (IterMax > 0 && IterMax_level >= IterMax) {
			IterMax_level = IterMax;
		}
		} else { IterMax_level = 1; }
		printf("IterMax = %d\n", IterMax_level);
		IRLS_MultipleMotion_OpticalFlow(
		    (u_levels + level),
		    (grad_It_levels + level),
		    (I_dt_levels + level),
		    lambdaD, lambdaS, sigmaD, sigmaS,
		    IterMax_level,
		    MotionParam.Error_Min_Threshold,
		    level);
	}
	for (i = 0; i < u->size(); i++) {
		(*u)[i] = u_levels[0][i];
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
LevelDown(ImgVector<VECTOR_2D> *u_levels, int level)
{
	int x, y;

	for (y = 0; y < u_levels[level].height(); y++) {
		for (x = 0; x < u_levels[level].width(); x++) {
			u_levels[level].ref(x, y).x = 2.0 * u_levels[level + 1].get(x / 2, y / 2).x;
			u_levels[level].ref(x, y).y = 2.0 * u_levels[level + 1].get(x / 2, y / 2).y;
		}
	}
}


bool
IRLS_MultipleMotion_OpticalFlow(ImgVector<VECTOR_2D> *u, ImgVector<VECTOR_2D> *Img_g, ImgVector<double> *Img_t, double lambdaD, double lambdaS, double sigmaD, double sigmaS, int IterMax, double ErrorMinThreshold, int level)
{
	ERROR Error("IRLS_MultipleMotion_OpticalFlow");
	ImgVector<VECTOR_2D> u_np1;
	VECTOR_2D sup;
	VECTOR_2D dE;
	double E = 0.0;
	double E_prev = 0.0;
	int ErrorIncrementCount = 0;
	int site;
	int n;

	if (u == nullptr) {
		Error.Value("u");
		Error.PointerNull();
		goto ExitError;
	} else if (Img_g == nullptr) {
		Error.Value("Img_g");
		Error.PointerNull();
		goto ExitError;
	} else if (Img_t == nullptr) {
		Error.Value("Img_t");
		Error.PointerNull();
		goto ExitError;
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
		if ((n & 0x3F) == 0) {
			printf("u[130, 130]: (%f, %f)\n", u_np1[u->width() * u->width() / 2 + u->width() / 2].x, u_np1[u->width() * u->width() / 2 + u->width() / 2].y);
		}
		if (level == 0) {
			if ((n & 0x3F) == 0) {
				E = Error_MultipleMotion(u, Img_g, Img_t, lambdaD, lambdaS, sigmaD, sigmaS);
			}
		} else {
			E_prev = E;
			E = Error_MultipleMotion(u, Img_g, Img_t, lambdaD, lambdaS, sigmaD, sigmaS);
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
		if (E < ErrorMinThreshold || ErrorIncrementCount > 3) {
			break;
		}
	}
	return MEANINGFUL_SUCCESS;
// Error
ExitError:
	return MEANINGFUL_FAILURE;
}


VECTOR_2D
Error_u(int site, ImgVector<VECTOR_2D> *u, ImgVector<VECTOR_2D> *Img_g, ImgVector<double> *Img_t, const double &lambdaD, const double &lambdaS, const double &sigmaD, const double &sigmaS)
{
	double (*psiD)(const double&, const double&) = Geman_McClure_psi;
	double (*psiS)(const double&, const double&) = Geman_McClure_psi;
	VECTOR_2D us;
	double Center;
	VECTOR_2D Neighbor;
	VECTOR_2D E_u;
	int x, y;

	x = site % u->width();
	y = site / u->width();

	us = u->get(site);
	Center = (*psiD)(Img_g->get(site).x * us.x + Img_g->get(site).y * us.y + Img_t->get(site), sigmaD);

	Neighbor.x = .0;
	Neighbor.y = .0;
	if (x > 0) {
		Neighbor.x += (*psiS)(us.x - u->get(x - 1, y).x, sigmaS);
		Neighbor.y += (*psiS)(us.y - u->get(x - 1, y).y, sigmaS);
	}
	if (x < u->width() - 1) {
		Neighbor.x += (*psiS)(us.x - u->get(x + 1, y).x, sigmaS);
		Neighbor.y += (*psiS)(us.y - u->get(x + 1, y).y, sigmaS);
	}
	if (y > 0) {
		Neighbor.x += (*psiS)(us.x - u->get(x, y - 1).x, sigmaS);
		Neighbor.y += (*psiS)(us.y - u->get(x, y - 1).y, sigmaS);
	}
	if (y < u->height() - 1) {
		Neighbor.x += (*psiS)(us.x - u->get(x, y + 1).x, sigmaS);
		Neighbor.y += (*psiS)(us.y - u->get(x, y + 1).y, sigmaS);
	}

	E_u.x += lambdaD * Img_g->get(site).x * Center + lambdaS * Neighbor.x;
	E_u.y += lambdaD * Img_g->get(site).y * Center + lambdaS * Neighbor.y;
	return E_u;
}


VECTOR_2D
sup_Error_uu(ImgVector<VECTOR_2D> *Img_g, const double &lambdaD, const double &lambdaS, const double &sigmaD, const double &sigmaS)
{
	static VECTOR_2D Img_g_max;
	VECTOR_2D sup;

	if (Img_g != nullptr) {
		Img_g_max.reset();
		for (int i = 0; i < Img_g->size(); i++) {
			if (Img_g_max.x < POW2(Img_g->get(i).x)) {
				Img_g_max.x = POW2(Img_g->get(i).x);
			}
			if (Img_g_max.y < POW2(Img_g->get(i).y)) {
				Img_g_max.y = POW2(Img_g->get(i).y);
			}
		}
	}
	sup.x = lambdaD * Img_g_max.x / POW2(sigmaD) + 4.0 * lambdaS / POW2(sigmaS);
	sup.y = lambdaD * Img_g_max.y / POW2(sigmaD) + 4.0 * lambdaS / POW2(sigmaS);
	return sup;
}


double
Error_MultipleMotion(ImgVector<VECTOR_2D> *u, ImgVector<VECTOR_2D> *Img_g, ImgVector<double> *Img_t, const double &lambdaD, const double &lambdaS, const double &sigmaD, const double &sigmaS)
{
	double (*rhoD)(const double&, const double&) = Geman_McClure_rho;
	double (*rhoS)(const double&, const double&) = Geman_McClure_rho;
	VECTOR_2D us;
	double Center;
	VECTOR_2D Neighbor;
	double E = 0.0;
	int x, y;

#pragma omp parallel for private(x, us, Neighbor, Center) reduction(+:E)
	for (y = 0; y < u->height(); y++) {
		for (x = 0; x < u->width(); x++) {
			us = u->get(x, y);
			Neighbor.x = .0;
			Neighbor.y = .0;
			if (x > 0) {
				Neighbor.x += (*rhoS)(us.x - u->get(x - 1, y).x, sigmaS);
				Neighbor.y += (*rhoS)(us.y - u->get(x - 1, y).y, sigmaS);
			}
			if (x < u->width() - 1) {
				Neighbor.x += (*rhoS)(us.x - u->get(x + 1, y).x, sigmaS);
				Neighbor.y += (*rhoS)(us.y - u->get(x + 1, y).y, sigmaS);
			}
			if (y > 0) {
				Neighbor.x += (*rhoS)(us.x - u->get(x, y - 1).x, sigmaS);
				Neighbor.y += (*rhoS)(us.y - u->get(x, y - 1).y, sigmaS);
			}
			if (y < u->height() - 1) {
				Neighbor.x += (*rhoS)(us.x - u->get(x, y + 1).x, sigmaS);
				Neighbor.y += (*rhoS)(us.y - u->get(x, y + 1).y, sigmaS);
			}
			Center = (*rhoD)(Img_g->get(x, y).x * us.x
			    + Img_g->get(x, y).y * us.y
			    + Img_t->get(x, y),
			    sigmaD);
			E += lambdaD * Center + lambdaS * (Neighbor.x + Neighbor.y);
		}
	}
	return E;
}


bool
MultipleMotion_write(const ImgVector<double> *img_prev, const ImgVector<double> *img_next, const ImgVector<VECTOR_2D> *u, const std::string &filename)
{
	ERROR Error("MultipleMotion_write");
	FILE *fp = nullptr;
	VECTOR_2D v;
	int x, y;
	MotionCompensation compensated(img_prev, img_next, u);
	PNM pnm;
	std::string filename_compensated;

	if (u == nullptr) {
		Error.Value("u");
		Error.PointerNull();
		goto ExitError;
	}

	printf("\n* Output The Optical Flow to '%s'(binary)\n", filename.c_str());
	if ((fp = fopen(filename.c_str(), "wb")) == nullptr) {
		Error.Function("fopen");
		Error.File(filename.c_str());
		Error.FileWrite();
		goto ErrorFileOpenFail;
	}
	fprintf(fp, "%d %d\n", u->width(), u->height());
	for (y = 0; y < u->height(); y++) {
		for (x = 0; x < u->width(); x++) {
			v = u->get(x, y);
			if (fwrite(&v.x, sizeof(double), 1, fp) < 1) {
				Error.Function("fwrite");
				Error.Value("u(x, y).x");
				Error.FunctionFail();
				goto ExitError;
			}
			if (fwrite(&v.y, sizeof(double), 1, fp) < 1) {
				Error.Function("fwrite");
				Error.Value("u(x, y).y");
				Error.FunctionFail();
				goto ExitError;
			}
		}
	}
	fclose(fp);

	compensated.create_image_compensated(); // Make compensated image
	filename_compensated = filename.substr(0, filename.length() - 4) + "compensated" + filename.substr(filename.length() - 4);
	printf("* Output The Compensated Image from Optical Flow to '%s'(binary)\n\n", filename_compensated.c_str());
	pnm.copy(PORTABLE_GRAYMAP_BINARY, compensated.width(), compensated.height(), 255, compensated.ref_image_compensated().data(), 1.0);
	pnm.write(filename_compensated.c_str());
	pnm.free();

	return MEANINGFUL_SUCCESS;
// Error
ErrorFileOpenFail:
	return MEANINGFUL_FAILURE;
ExitError:
	fclose(fp);
	return MEANINGFUL_FAILURE;
}

