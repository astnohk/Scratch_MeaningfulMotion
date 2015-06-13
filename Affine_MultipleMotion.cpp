/* Robust Estimation of Multiple Motions
 *
 * This motion estimation based on
 * M.J.Black and P.Anandan, "The Robust Estimation of Multiple Motions: Parametric and Piecewise-Smooth Flow Fields," Computer Vision and Image Understanding, Vol.63, No.1, 1996, pp.75-104.
 */

#include "Scratch_MeaningfulMotion.h"
#include "Affine_MultipleMotion.h"
#include "MultiResolution.h"
#include "MEstimator.h"




VECTOR_AFFINE
MultipleMotion_Affine(double *It, double *Itp1, double MaxInt, SIZE size_img, MULTIPLE_MOTION_PARAM MotionParam)
{
	ERROR Error("MultipleMotion_Affine");

	// M-estimator parameter
	const double sigmaD = 0.1 * sqrt(3.0);
	VECTOR_AFFINE u;
	double *It_normalize = nullptr;
	double *Itp1_normalize = nullptr;
	double **I_dt_levels = nullptr;
	double **It_levels = nullptr;
	double **Itp1_levels = nullptr;
	VECTOR_2D **grad_It_levels = nullptr;
	SIZE size_img_l;
	int level;
	int IterMax;

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
	for (int i = 0; i < size_img.width * size_img.height; i++) {
		It_normalize[i] = (double)It[i] / MaxInt;
		Itp1_normalize[i] = (double)Itp1[i] / MaxInt;
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

	for (int i = 0; i < NUM_AFFINE_PARAMETER; i++) {
		u.a[i] = .0;
	}
	for (level = MotionParam.Level - 1; level >= 0; level--) {
		printf("\nLevel %d :\n", level);
		u.a[0] *= 2;
		u.a[3] *= 2;
		size_img_l.width = floor(size_img.width * pow_int(0.5, level));
		size_img_l.height = floor(size_img.height * pow_int(0.5, level));
		IterMax = 2 * MAX(size_img_l.width, size_img_l.height);
		IRLS_MultipleMotion_Affine(&u, grad_It_levels[level], I_dt_levels[level], size_img_l,
		    sigmaD,
		    IterMax, MotionParam.Error_Min_Threshold);
	}
	delete[] grad_It_levels;
	delete[] I_dt_levels;
	delete[] Itp1_levels;
	delete[] It_levels;
	delete[] It_normalize;
	delete[] Itp1_normalize;
	return u;
// Error
ExitError:
	delete[] grad_It_levels;
	delete[] I_dt_levels;
	delete[] Itp1_levels;
	delete[] It_levels;
	delete[] It_normalize;
	delete[] Itp1_normalize;
	for (int i = 0; i < NUM_AFFINE_PARAMETER; i++) {
		u.a[i] = .0;
	}
	return u;
}


int
IRLS_MultipleMotion_Affine(VECTOR_AFFINE *u, VECTOR_2D *Img_g, double *Img_t, SIZE size_img, double sigmaD, int IterMax, double ErrorMinThreshold)
{
	printf("sigmaD = %e\n", sigmaD);
	for (int n = 0; n < IterMax; n++) {
		VECTOR_AFFINE u_np1;
		VECTOR_AFFINE sup;
		VECTOR_AFFINE dE;
		double E = 0.0;
		double omega;
		omega = 1.0E-4;
		dE = Error_a(u, Img_g, Img_t, size_img, sigmaD);
		sup = sup_Error_aa(Img_g, size_img, sigmaD);
		for (int i = 0; i < NUM_AFFINE_PARAMETER; i++) {
			u_np1.a[i] = .0;
		}
		for (int i = 0; i < 6; i++) {
			if (fabs(sup.a[i]) < 1.0E-16) {
				u_np1.a[i] = u->a[i] - omega / 1.0E-16 * SIGN_NOZERO(sup.a[i]) * dE.a[i];
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
	double (*psiD)(double, double) = Geman_McClure_psi;
	VECTOR_AFFINE E_a;
	VECTOR_2D u_a;

	for (int i = 0; i < NUM_AFFINE_PARAMETER; i++) {
		E_a.a[i] = .0;
	}
	for (int site = 0; site < size_img.width * size_img.height; site++) {
		int x, y;
		x = site % size_img.width;
		y = site / size_img.width;
		u_a.x = u->a[0] + u->a[1] * x + u->a[2] * y;
		u_a.y = u->a[3] + u->a[4] * x + u->a[5] * y;
		E_a.a[0] += Img_g[site].x * (*psiD)(Img_g[site].x * u_a.x + Img_g[site].y * u_a.y + Img_t[site], sigmaD);
		E_a.a[1] += Img_g[site].x * x * (*psiD)(Img_g[site].x * u_a.x + Img_g[site].y * u_a.y + Img_t[site], sigmaD);
		E_a.a[2] += Img_g[site].x * y * (*psiD)(Img_g[site].x * u_a.x + Img_g[site].y * u_a.y + Img_t[site], sigmaD);
		E_a.a[3] += Img_g[site].y * (*psiD)(Img_g[site].x * u_a.x + Img_g[site].y * u_a.y + Img_t[site], sigmaD);
		E_a.a[4] += Img_g[site].y * x * (*psiD)(Img_g[site].x * u_a.x + Img_g[site].y * u_a.y + Img_t[site], sigmaD);
		E_a.a[5] += Img_g[site].y * y * (*psiD)(Img_g[site].x * u_a.x + Img_g[site].y * u_a.y + Img_t[site], sigmaD);
	}
	return E_a;
}


VECTOR_AFFINE
sup_Error_aa(VECTOR_2D *Img_g, SIZE size, double sigmaD)
{
	ERROR Error("sup_Error_aa");

	VECTOR_AFFINE sup;
	VECTOR_AFFINE u_aa_max;
	int i;

	for (i = 0; i < NUM_AFFINE_PARAMETER; i++) {
		u_aa_max.a[i] = .0;
	}
	if (Img_g == nullptr) {
		Error.Value("Img_g");
		Error.PointerNull();
		return u_aa_max;
	}
	for (i = 0; i < size.width * size.height; i++) {
		int x, y;
		x = i % size.width;
		y = i / size.width;
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
	double (*rhoD)(double, double) = Geman_McClure_rho;
	double E = 0.0;

	for (int site = 0; site < size_img.width * size_img.height; site++) {
		int x, y;
		VECTOR_2D u_a;
		x = site % size_img.width;
		y = site / size_img.width;
		u_a.x = u->a[0] + u->a[1] * x + u->a[2] * y;
		u_a.y = u->a[3] + u->a[4] * x + u->a[5] * y;
		E += (*rhoD)(Img_g[site].x * u_a.x + Img_g[site].y * u_a.y + Img_t[site], sigmaD);
	}
	return E;
}


int
MultipleMotion_Affine_write(VECTOR_AFFINE u, const char *filename)
{
	ERROR Error("MultipleMotion_Affine_write");

	FILE *fp = nullptr;
	int i;

	if (filename == nullptr) {
		Error.Value("filename");
		Error.PointerNull();
		goto ExitError;
	}

	printf("* Output Affine Parameter to '%s'\n", filename);
	if ((fp = fopen(filename, "w")) == nullptr) {
		Error.Function("fopen");
		Error.Value(filename);
		Error.FileRead();
		goto ExitError;
	}
	for (i = 0; i < NUM_AFFINE_PARAMETER; i++) {
		if (fprintf(fp, "%0.16e ", u.a[i]) < 0) {
			Error.Function("fprintf");
			Error.Value("u(a(i))");
			Error.FunctionFail();
			goto ExitError;
		}
		if (fprintf(fp, "\n") < 0) {
			Error.Function("fprintf");
			Error.Value("'\n'");
			Error.FunctionFail();
			goto ExitError;
		}
	}
	fclose(fp);
	return MEANINGFUL_SUCCESS;
// Error
ExitError:
	fclose(fp);
	return MEANINGFUL_FAILURE;
}

