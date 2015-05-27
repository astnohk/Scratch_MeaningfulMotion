/* Robust Estimation of Multiple Motions
 *
 * This motion estimation based on
 * M.J.Black and P.Anandan, "The Robust Estimation of Multiple Motions: Parametric and Piecewise-Smooth Flow Fields," Computer Vision and Image Understanding, Vol.63, No.1, 1996, pp.75-104.
 */

#include "Scratch_MeaningfulA.h"
#include "MultiResolution.h"
#include "Affine_MultipleMotion.h"




VECTOR_AFFINE
MultipleMotion_Affine(double *It, double *Itp1, SIZE size_img, MULTIPLE_MOTION_PARAM MotionParam)
{
	char *FunctionName = "MultipleMotion_Affine";
	char *ErrorFunction = "";
	char *ErrorValue = "";

	VECTOR_AFFINE u;
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

	for (i = 0; i < NUM_AFFINE_PARAMETER; i++) {
		u.a[i] = .0;
	}
	for (level = MotionParam.Level - 1; level >= 0; level--) {
		printf("\nLevel %d :\n", level);
		u.a[0] *= 2;
		u.a[3] *= 2;
		size_img_res.width = floor(size_img.width * pow_int(0.5, level));
		size_img_res.height = floor(size_img.height * pow_int(0.5, level));
		IRLS_MultipleMotion_Affine(&u, grad_It_levels[level], I_dt_levels[level], size_img_res,
		    MotionParam.sigmaD,
		    MotionParam.IRLS_Iter_Max, MotionParam.Error_Min_Threshold);
	}
	free(grad_It_levels);
	free(I_dt_levels);
	free(Itp1_levels);
	free(It_levels);
	return u;
/* Error */
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
	for (i = 0; i < NUM_AFFINE_PARAMETER; i++) {
		u.a[i] = .0;
	}
	return u;
}


int
IRLS_MultipleMotion_Affine(VECTOR_AFFINE *u, VECTOR_2D *Img_g, double *Img_t, SIZE size_img, double sigmaD, int IterMax, double ErrorMinThreshold)
{
	VECTOR_AFFINE u_np1;
	VECTOR_AFFINE sup;
	VECTOR_AFFINE dE;
	double E = 0.0;
	double myu_max;
	double omega;
	int i, n;

	printf("sigmaD = %e\n", sigmaD);
	for (n = 0; n < IterMax; n++) {
		myu_max = cos(M_PI / (n + 1.0)); /* Jacobi relaxation determined by the largest eigenvalue of the Jacobi iteration matrix */
		if (fabs(myu_max) < 1.0E-10) {
			myu_max = 1.0E-10;
		}
		omega = 2.0 * (1.0 - sqrt(1.0 - POW2(myu_max))) / POW2(myu_max) * 1E-4;
		dE = Error_a(u, Img_g, Img_t, size_img, sigmaD);
		sup = sup_Error_aa(Img_g, size_img, sigmaD);
		for (i = 0; i < NUM_AFFINE_PARAMETER; i++) {
			u_np1.a[i] = .0;
		}
		for (i = 0; i < 6; i++) {
			if (fabs(sup.a[i]) < 1.0E-6) {
				u_np1.a[i] = u->a[i] - omega / 1.0E-6 * SIGN_NOZERO(sup.a[i]) * dE.a[i];
			} else {
				u_np1.a[i] = u->a[i] - omega / sup.a[i] * dE.a[i];
			}
		}
		*u = u_np1;
		E = Error_Affine(u, Img_g, Img_t, size_img, sigmaD);
		if ((n & 0x3F) == 0) {
			printf("E(%4d) = %e,  u = [%.4e, %.4e, %.4e, %.4e, %.4e, %.4e]\n", n, E, u->a[0], u->a[1], u->a[2], u->a[3], u->a[4], u->a[5]);
		}
		if (E < ErrorMinThreshold) {
			break;
		}
	}
	return MEANINGFUL_SUCCESS;
}


VECTOR_AFFINE
Error_a(VECTOR_AFFINE *u, VECTOR_2D *Img_g, double *Img_t, SIZE size_img, double sigmaD)
{
	VECTOR_AFFINE E_a;
	VECTOR_2D u_a;
	int site;
	int x, y;
	int i;

	for (i = 0; i < NUM_AFFINE_PARAMETER; i++) {
		E_a.a[i] = .0;
	}
	for (site = 0; site < size_img.width * size_img.height; site++) {
		x = site % size_img.width;
		y = (int)floor(site / size_img.width);
		u_a.x = u->a[0] + u->a[1] * x + u->a[2] * y;
		u_a.y = u->a[3] + u->a[4] * x + u->a[5] * y;
		E_a.a[0] += Img_g[site].x * cal_psyD(Img_g[site].x * u_a.x + Img_g[site].y * u_a.y + Img_t[site], sigmaD);
		E_a.a[1] += Img_g[site].x * x * cal_psyD(Img_g[site].x * u_a.x + Img_g[site].y * u_a.y + Img_t[site], sigmaD);
		E_a.a[2] += Img_g[site].x * y * cal_psyD(Img_g[site].x * u_a.x + Img_g[site].y * u_a.y + Img_t[site], sigmaD);
		E_a.a[3] += Img_g[site].y * cal_psyD(Img_g[site].x * u_a.x + Img_g[site].y * u_a.y + Img_t[site], sigmaD);
		E_a.a[4] += Img_g[site].y * x * cal_psyD(Img_g[site].x * u_a.x + Img_g[site].y * u_a.y + Img_t[site], sigmaD);
		E_a.a[5] += Img_g[site].y * y * cal_psyD(Img_g[site].x * u_a.x + Img_g[site].y * u_a.y + Img_t[site], sigmaD);
	}
	return E_a;
}


VECTOR_AFFINE
sup_Error_aa(VECTOR_2D *Img_g, SIZE size, double sigmaD)
{
	char *FunctionName = "sup_Error_aa";
	char *ErrorValue = "";

	VECTOR_AFFINE sup;
	VECTOR_AFFINE u_aa_max;
	int x, y;
	int i;

	for (i = 0; i < NUM_AFFINE_PARAMETER; i++) {
		u_aa_max.a[i] = .0;
	}
	if (Img_g == NULL) {
		ErrorValue = "Img_g";
		fprintf(stderr, "*** %s() error - The pointer (*%s) is NULL ***\n", FunctionName, ErrorValue);
		return u_aa_max;
	}
	for (i = 0; i < size.width * size.height; i++) {
		x = i % size.width;
		y = (int)floor(i / size.width);
		/* u = a0 + a1 * x + a2 * y */
		if (u_aa_max.a[0] < POW2(Img_g[i].x)) {
			u_aa_max.a[0] = POW2(Img_g[i].x);
		}
		if (u_aa_max.a[1] < POW2(Img_g[i].x * x)) {
			u_aa_max.a[1] = POW2(Img_g[i].x * x);
		}
		if (u_aa_max.a[2] < POW2(Img_g[i].x * y)) {
			u_aa_max.a[2] = POW2(Img_g[i].x * y);
		}
		/* v = a3 + a4 * x + a5 * y */
		if (u_aa_max.a[3] < POW2(Img_g[i].y)) {
			u_aa_max.a[3] = POW2(Img_g[i].y);
		}
		if (u_aa_max.a[4] < POW2(Img_g[i].y * x)) {
			u_aa_max.a[4] = POW2(Img_g[i].y * x);
		}
		if (u_aa_max.a[5] < POW2(Img_g[i].y * y)) {
			u_aa_max.a[5] = POW2(Img_g[i].y * y);
		}
	}
	sup.a[0] = u_aa_max.a[0] * 2.0 / POW2(sigmaD);
	sup.a[1] = u_aa_max.a[1] * 2.0 / POW2(sigmaD);
	sup.a[2] = u_aa_max.a[2] * 2.0 / POW2(sigmaD);
	sup.a[3] = u_aa_max.a[3] * 2.0 / POW2(sigmaD);
	sup.a[4] = u_aa_max.a[4] * 2.0 / POW2(sigmaD);
	sup.a[5] = u_aa_max.a[5] * 2.0 / POW2(sigmaD);
	return sup;
}


double
Error_Affine(VECTOR_AFFINE *u, VECTOR_2D *Img_g, double *Img_t, SIZE size_img, double sigmaD)
{
	double E = 0.0;
	VECTOR_2D u_a;
	int site;
	int x, y;

	for (site = 0; site < size_img.width * size_img.height; site++) {
		x = site % size_img.width;
		y = (int)floor(site / size_img.width);
		u_a.x = u->a[0] + u->a[1] * x + u->a[2] * y;
		u_a.y = u->a[3] + u->a[4] * x + u->a[5] * y;
		E += cal_rhoD(Img_g[site].x * u_a.x + Img_g[site].y * u_a.y + Img_t[site], sigmaD);
	}
	return E;
}


int
MultipleMotion_Affine_write(VECTOR_AFFINE u, char *filename)
{
	char *FunctionName = "MultipleMotion_Affine_write";
	char *ErrorFunction = "";
	char *ErrorValue = "";

	FILE *fp;
	int i;

	if (filename == NULL) {
		ErrorValue = "filename";
		goto ErrorPointerNull;
	}

	printf("* Output Affine Parameter to '%s'\n", filename);
	if ((fp = fopen(filename, "w")) == NULL) {
		ErrorFunction = "fopen";
		ErrorValue = filename;
		goto ErrorFileOpenFail;
	}
	for (i = 0; i < NUM_AFFINE_PARAMETER; i++) {
		if (fprintf(fp, "%0.16e ", u.a[i]) < 0) {
			ErrorFunction = "fprintf";
			ErrorValue = "u(a(i))";
			goto ErrorFunctionFail;
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

