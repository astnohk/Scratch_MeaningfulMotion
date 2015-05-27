/* Robust Estimation of Multiple Motions
 *
 * This motion estimation based on
 * M.J.Black and P.Anandan, "The Robust Estimation of Multiple Motions: Parametric and Piecewise-Smooth Flow Fields," Computer Vision and Image Understanding, Vol.63, No.1, 1996, pp.75-104.
 */

#include "Scratch_MeaningfulA.h"
#include "MultiResolution.h"
#include "OpticalFlow_MultipleMotion.h"


/* define for debug */
#define SHOW_IRLS_MULTIPLEMOTION_OPTICALFLOW_ERROR
/* /define for debug */




VECTOR_2D *
MultipleMotion_OpticalFlow(double *It, double *Itp1, SIZE size_img, MULTIPLE_MOTION_PARAM MotionParam)
{
	char *FunctionName = "MultipleMotion_OpticalFlow";
	char *ErrorFunction = "";
	char *ErrorValue = "";

	VECTOR_2D *u = NULL;
	double **I_dt_levels = NULL;
	double **It_levels = NULL;
	double **Itp1_levels = NULL;
	VECTOR_2D **grad_It_levels = NULL;
	SIZE size_img_res;
	int level;
	int i;

	if (It == NULL) {
		ErrorValue = "It";
		goto ErrorPointerNull;
	} else if (Itp1 == NULL) {
		ErrorValue = "Itp1";
		goto ErrorPointerNull;
	}

	if ((u = (VECTOR_2D *)calloc((size_t)size_img.width * size_img.height, sizeof(VECTOR_2D))) == NULL) {
		ErrorValue = "u";
		goto ErrorMalloc;
	}
	/* Make Pyramid */
	if ((It_levels = Pyramider(It, size_img, MotionParam.Level)) == NULL) {
		ErrorFunction = "Pyramider";
		ErrorValue = "It_levels";
		goto ErrorFunctionFail;
	}
	if ((Itp1_levels = Pyramider(Itp1, size_img, MotionParam.Level)) == NULL) {
		ErrorFunction = "Pyramider";
		ErrorValue = "Itp1_levels";
		goto ErrorFunctionFail;
	}
	/* Derivative about time */
	if ((I_dt_levels = dt_Pyramid(It_levels, Itp1_levels, size_img, MotionParam.Level)) == NULL) {
		ErrorFunction = "dt_Pyramid";
		ErrorValue = "I_dt_levels";
		goto ErrorFunctionFail;
	}
	/* Derivative about space */
	if ((grad_It_levels = grad_Pyramid(It_levels, Itp1_levels, size_img, MotionParam.Level)) == NULL) {
		ErrorFunction = "grad_Pyramid";
		ErrorValue = "grad_It_levels";
		goto ErrorFunctionFail;
	}

	for (level = MotionParam.Level - 1; level >= 0; level--) {
		size_img_res.width = floor(size_img.width * pow_int(0.5, level));
		size_img_res.height = floor(size_img.height * pow_int(0.5, level));
		for (i = 0; i < size_img.width * size_img.height; i++) {
			u[i].x *= 2.0;
			u[i].y *= 2.0;
		}
		IRLS_MultipleMotion_OpticalFlow(u, grad_It_levels[level], I_dt_levels[level], size_img_res,
		    MotionParam.lambdaD, MotionParam.lambdaS, MotionParam.sigmaD, MotionParam.sigmaS);
	}
	free(grad_It_levels);
	free(I_dt_levels);
	free(Itp1_levels);
	free(It_levels);
	return u;
/* Error */
ErrorMalloc:
	fprintf(stderr, "*** %s() error - Cannot allocate memory for (*%s) ***\n", FunctionName, ErrorValue);
	goto ErrorReturn;
ErrorPointerNull:
	fprintf(stderr, "*** %s() error - The pointer (*%s) is NULL ***\n", FunctionName, ErrorValue);
	goto ErrorReturn;
ErrorFunctionFail:
	fprintf(stderr, "*** %s() error - %s() failed to compute (%s) ***\n", FunctionName, ErrorFunction, ErrorValue);
ErrorReturn:
	free(grad_It_levels);
	free(I_dt_levels);
	free(Itp1_levels);
	free(It_levels);
	free(u);
	return NULL;
}


int
IRLS_MultipleMotion_OpticalFlow(VECTOR_2D *u, VECTOR_2D *Img_g, double *Img_t, SIZE size_img, double lambdaD, double lambdaS, double sigmaD, double sigmaS)
{
#define IRLS_ITER_MAX 300
	char *FunctionName = "IRLS";
	VECTOR_2D *u_np1 = NULL;
	VECTOR_2D sup;
	VECTOR_2D dE;
	int site;
	double myu_max;
	double omega;
	int n;

	if ((u_np1 = (VECTOR_2D *)calloc(size_img.width * size_img.height, sizeof(VECTOR_2D))) == NULL) {
		fprintf(stderr, "*** %s() error - Cannot allocate memory for (*u_np1) ***\n", FunctionName);
		return MEANINGFUL_FAILURE;
	}
	printf("initialE = %e\n", Error_MultipleMotion(u, Img_g, Img_t, size_img, lambdaD, lambdaS, sigmaD, sigmaS));
	for (n = 0; n < IRLS_ITER_MAX; n++) {
		//myu_max = cos(M_PI / (n + 2.0)); /* Jacobi relaxation determined by the largest eigenvalue of the Jacobi iteration matrix */
		myu_max = 1.0 - pow_int(0.9, n + 1); /* Jacobi relaxation determined by the largest eigenvalue of the Jacobi iteration matrix */
		if (fabs(myu_max) < 1.0E-6) {
			myu_max = 1.0E-6;
		}
		omega = 2.0 * (1.0 - sqrt(1.0 - POW2(myu_max))) / POW2(myu_max);
		for (site = 0; site < size_img.width * size_img.height; site++) { /* Calc for all sites */
			dE = Error_u(site, u, Img_g, Img_t, size_img, lambdaD, lambdaS, sigmaD, sigmaS);
			sup = sup_Error_uu(Img_g, size_img, lambdaD, lambdaS, sigmaD, sigmaS);
			u_np1[site].x = u[site].x - omega / sup.x * dE.x;
			u_np1[site].y = u[site].y - omega / sup.y * dE.y;
		}
		for (site = 0; site < size_img.width * size_img.height; site++) { /* Calc for all sites */
			u[site] = u_np1[site];
		}
#ifdef SHOW_IRLS_MULTIPLEMOTION_OPTICALFLOW_ERROR
		if ((n & 0x0F) == 0) {
			printf("E(%4d) = %e\n", n, Error_MultipleMotion(u, Img_g, Img_t, size_img, lambdaD, lambdaS, sigmaD, sigmaS));
		}
#endif
	}
	free(u_np1);
	return MEANINGFUL_SUCCESS;
}


VECTOR_2D
Error_u(int site, VECTOR_2D *u, VECTOR_2D *Img_g, double *Img_t, SIZE size_img, double lambdaD, double lambdaS, double sigmaD, double sigmaS)
{
	VECTOR_2D us;
	double Center;
	VECTOR_2D Neighbor;
	VECTOR_2D E_u = VECTOR_2D_ZERO;
	int x, y;

	x = site % size_img.width;
	y = (int)floor(site / size_img.width);

	us = u[site];
	Center = cal_psyD(Img_g[site].x * us.x + Img_g[site].y * us.y + Img_t[site], sigmaD);

	Neighbor.x = .0;
	Neighbor.y = .0;
	if (x > 0) {
		Neighbor.x += cal_psyS(us.x - u[size_img.width * y + x - 1].x, sigmaS);
		Neighbor.y += cal_psyS(us.y - u[size_img.width * y + x - 1].y, sigmaS);
	}
	if (x < size_img.width - 1) {
		Neighbor.x += cal_psyS(us.x - u[size_img.width * y + x + 1].x, sigmaS);
		Neighbor.y += cal_psyS(us.y - u[size_img.width * y + x + 1].y, sigmaS);
	}
	if (y > 0) {
		Neighbor.x += cal_psyS(us.x - u[size_img.width * (y - 1) + x].x, sigmaS);
		Neighbor.y += cal_psyS(us.y - u[size_img.width * (y - 1) + x].y, sigmaS);
	}
	if (y < size_img.height - 1) {
		Neighbor.x += cal_psyS(us.x - u[size_img.width * (y + 1) + x].x, sigmaS);
		Neighbor.y += cal_psyS(us.y - u[size_img.width * (y + 1) + x].y, sigmaS);
	}

	E_u.x += lambdaD * Img_g[site].x * Center + lambdaS * Neighbor.x;
	E_u.y += lambdaD * Img_g[site].y * Center + lambdaS * Neighbor.y;
	return E_u;
}


VECTOR_2D
sup_Error_uu(VECTOR_2D *Img_g, SIZE size, double lambdaD, double lambdaS, double sigmaD, double sigmaS)
{
	static int called = 0;
	static VECTOR_2D Img_g_max = VECTOR_2D_ZERO;
	int i;

	if (called == 0) {
		called = 1;
		for (i = 0; i < size.width * size.height; i++) {
			if (Img_g_max.x < POW2(Img_g[i].x)) {
				Img_g_max.x = POW2(Img_g[i].x);
			}
			if (Img_g_max.y < POW2(Img_g[i].y)) {
				Img_g_max.y = POW2(Img_g[i].y);
			}
		}
	}
	return (VECTOR_2D){lambdaD * Img_g_max.x / POW2(sigmaD) + 4.0 * lambdaS / POW2(sigmaS),
	    lambdaD * Img_g_max.y / POW2(sigmaD) + 4.0 * lambdaS / POW2(sigmaS)};
}


double
Error_MultipleMotion(VECTOR_2D *u, VECTOR_2D *Img_g, double *Img_t, SIZE size_img, double lambdaD, double lambdaS, double sigmaD, double sigmaS)
{
	VECTOR_2D us;
	double Center;
	VECTOR_2D Neighbor;
	double E = 0.0;
	int x, y;

	for (y = 0; y < size_img.height; y++) {
		for (x = 0; x < size_img.width; x++) {
			us = u[size_img.width * y + x];
			Neighbor = VECTOR_2D_ZERO;
			if (x > 0) {
				Neighbor.x += cal_rhoS(us.x - u[size_img.width * y + x - 1].x, sigmaS);
				Neighbor.y += cal_rhoS(us.y - u[size_img.width * y + x - 1].y, sigmaS);
			}
			if (x < size_img.width - 1) {
				Neighbor.x += cal_rhoS(us.x - u[size_img.width * y + x + 1].x, sigmaS);
				Neighbor.y += cal_rhoS(us.y - u[size_img.width * y + x + 1].y, sigmaS);
			}
			if (y > 0) {
				Neighbor.x += cal_rhoS(us.x - u[size_img.width * (y - 1) + x].x, sigmaS);
				Neighbor.y += cal_rhoS(us.y - u[size_img.width * (y - 1) + x].y, sigmaS);
			}
			if (y < size_img.height - 1) {
				Neighbor.x += cal_rhoS(us.x - u[size_img.width * (y + 1) + x].x, sigmaS);
				Neighbor.y += cal_rhoS(us.y - u[size_img.width * (y + 1) + x].y, sigmaS);
			}
			Center = cal_rhoD(Img_g[size_img.width * y + x].x * us.x
			    + Img_g[size_img.width * y + x].y * us.y
			    + Img_t[size_img.width * y + x],
			    sigmaD);
			E += lambdaD * Center + lambdaS * (Neighbor.x + Neighbor.y);
		}
	}
	return E;
}


int
MultipleMotion_write(VECTOR_2D *u, SIZE size, char *filename)
{
	char *FunctionName = "MultipleMotion_write";
	char *ErrorFunction = "";
	char *ErrorValue = "";

	FILE *fp;
	int x, y;

	if (u == NULL) {
		ErrorValue = "u";
		goto ErrorPointerNull;
	} else if (filename == NULL) {
		ErrorValue = "filename";
		goto ErrorPointerNull;
	}

	printf("* Output Optical Flow to '%s'\n", filename);
	if ((fp = fopen(filename, "w")) == NULL) {
		ErrorFunction = "fopen";
		ErrorValue = filename;
		goto ErrorFileOpenFail;
	}
	for (y = 0; y < size.width; y++) {
		for (x = 0; x < size.height; x++) {
			if (fprintf(fp, "%0.16e ", u[size.width * y + x].x) < 0) {
				ErrorFunction = "fprintf";
				ErrorValue = "u(x, y).x";
				goto ErrorFunctionFail;
			}
			if (fprintf(fp, "%0.16e ", u[size.width * y + x].y) < 0) {
				ErrorFunction = "fprintf";
				ErrorValue = "u(x, y).y";
				goto ErrorFunctionFail;
			}
		}
		if (fprintf(fp, "\n") < 0) {
			ErrorFunction = "fprintf";
			ErrorValue = "'\n'";
			goto ErrorFunctionFail;
		}
	}
	fclose(fp);
	return MEANINGFUL_SUCCESS;
/* Error */
ErrorPointerNull:
	fprintf(stderr, "*** %s error - The pointer (*%s) is NULL ***\n", FunctionName, ErrorValue);
	goto ErrorReturn;
ErrorFunctionFail:
	fprintf(stderr, "*** %s error - %s() failed to compute (%s) ***\n", FunctionName, ErrorFunction, ErrorValue);
	goto ErrorReturn;
ErrorFileOpenFail:
	fprintf(stderr, "*** %s error - Cannot open the file '%s' ***\n", FunctionName, ErrorValue);
ErrorReturn:
	fclose(fp);
	return MEANINGFUL_FAILURE;
}

