/* Robust Estimation of Multiple Motions
 *
 * This motion estimation based on
 * M.J.Black and P.Anandan, "The Robust Estimation of Multiple Motions: Parametric and Piecewise-Smooth Flow Fields," Computer Vision and Image Understanding, Vol.63, No.1, 1996, pp.75-104.
 */

#include "Affine_BlockMatching.h"




ImgVector<VECTOR_2D<double> > *
OpticalFlow_Affine_BlockMatching(ImgVector<double> *It, ImgVector<double> *Itp1, double MaxInt, MULTIPLE_MOTION_PARAM MotionParam)
{
	ERROR Error("MultipleMotion_Affine");
	std::bad_alloc bad_alloc;

	// M-estimator parameter
	const double sigmaD = 0.1 * sqrt(3.0);
	ImgVector<VECTOR_2D<double> >* u = nullptr;
	VECTOR_AFFINE* affine = nullptr;
	ImgVector<double> It_normalize;
	ImgVector<double> Itp1_normalize;
	ImgVector<double> *I_dt_levels = nullptr;
	ImgVector<double> *It_levels = nullptr;
	ImgVector<double> *Itp1_levels = nullptr;
	ImgVector<VECTOR_2D<double> > *grad_It_levels = nullptr;
	std::vector<std::vector<VECTOR_2D<int> > > connected_domains;
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
	It_normalize = *It;
	Itp1_normalize = *Itp1;
	for (int i = 0; i < It_normalize.size(); i++) {
		It_normalize[i] /= MaxInt;
		Itp1_normalize[i] /= MaxInt;
	}
	// Make Pyramid
	try {
		It_levels = Pyramider(&It_normalize, MotionParam.Level);
	}
	catch (const std::bad_alloc& bad) {
		Error.Function("Pyramider");
		Error.Value("It_levels");
		Error.FunctionFail();
		bad_alloc = bad;
		goto ExitError;
	}
	if ((Itp1_levels = Pyramider(&Itp1_normalize, MotionParam.Level)) == nullptr) {
		Error.Function("Pyramider");
		Error.Value("Itp1_levels");
		Error.FunctionFail();
		bad_alloc = bad;
		goto ExitError;
	}
	// Derivative about time
	if ((I_dt_levels = dt_Pyramid(It_levels, Itp1_levels, MotionParam.Level)) == nullptr) {
		Error.Function("dt_Pyramid");
		Error.Value("I_dt_levels");
		Error.FunctionFail();
		bad_alloc = bad;
		goto ExitError;
	}
	// Derivative about space
	if ((grad_It_levels = grad_Pyramid(It_levels, Itp1_levels, MotionParam.Level)) == nullptr) {
		Error.Function("grad_Pyramid");
		Error.Value("grad_It_levels");
		Error.FunctionFail();
		bad_alloc = bad;
		goto ExitError;
	}

	for (int i = 0; i < NUM_AFFINE_PARAMETER; i++) {
		u.a[i] = .0;
	}
	for (level = MotionParam.Level - 1; level >= 0; level--) {
		u.a[0] *= 2;
		u.a[3] *= 2;
		IterMax = 2 * MAX(I_dt_levels[level].width(), I_dt_levels[level].height());
		IRLS_Affine_Block(
		    &u,
		    (grad_It_levels + level),
		    (I_dt_levels + level),
		    sigmaD,
		    IterMax,
		    MotionParam.Error_Min_Threshold);
	}
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
	for (int i = 0; i < NUM_AFFINE_PARAMETER; i++) {
		u.a[i] = .0;
	}
	throw bad_alloc;
}


bool
IRLS_Affine_Block(std::vector<VECTOR_AFFINE>* u, const std::vector<std::vector<VECTOR_2D<int> > >& connected_domains, const ImgVector<VECTOR_2D<double> > *Img_g, const ImgVector<double> *Img_t, double sigmaD, int IterMax, double ErrorMinThreshold)
{
	const double omega_initial = 1.0E-4;

	printf("sigmaD = %e\n", sigmaD);
	printf("size (%d, %d)\n", Img_g->width(), Img_g->height());
	printf("E %e\n", Error_Affine(u, Img_g, Img_t, sigmaD));
	for (int R = 0; R < connected_domains.size(); R++) {
		for (int n = 0; n < IterMax; n++) {
			VECTOR_AFFINE u_np1;
			VECTOR_AFFINE dE = Error_a_Block((*u)[R], connected_domains[R], Img_g, Img_t, sigmaD);
			VECTOR_AFFINE sup = sup_Error_aa_Block(connected_domains[R], Img_g, sigmaD);
			double E = 0.0;
			double omega = omega_initial;

			for (int i = 0; i < NUM_AFFINE_PARAMETER; i++) {
				u_np1.a[i] = .0;
			}
			for (int i = 0; i < NUM_AFFINE_PARAMETER; i++) {
				if (fabs(sup.a[i]) < 1.0E-16) {
					u_np1.a[i] = u->at(R).a[i] - omega / 1.0E-16 * SIGN_NOZERO(sup.a[i]) * dE.a[i];
				} else {
					u_np1.a[i] = u->at(R).a[i] - omega / sup.a[i] * dE.a[i];
				}
			}
			u->at(R) = u_np1;
			E = Error_Affine_Block((*u)[R], connected_domains[R], Img_g, Img_t, sigmaD);
			/*if ((n & 0x3F) == 0) {
				printf("E(%4d) = %e,  u = [%.4e, %.4e, %.4e, %.4e, %.4e, %.4e]\n", n, E, u->a[0], u->a[1], u->a[2], u->a[3], u->a[4], u->a[5]);
			}*/
			if (E < ErrorMinThreshold) {
				break;
			}
		}
	}
	return MEANINGFUL_SUCCESS;
}


VECTOR_AFFINE
Error_a_Block(const VECTOR_AFFINE& u, const std::vector<VECTOR_2D<int> >& connected_domain, const ImgVector<VECTOR_2D<double> > *Img_g, const ImgVector<double> *Img_t, double sigmaD)
{
	double (*psiD)(const double&, const double&) = Geman_McClure_psi;
	VECTOR_AFFINE E_a;
	VECTOR_2D<double> u_a;

	for (int i = 0; i < NUM_AFFINE_PARAMETER; i++) {
		E_a.a[i] = .0;
	}
	for (int site = 0; site < connected_domain.size(); site++) {
		int x = connected_domain[site].x;
		int y = connected_domain[site].y;
		u_a.x = u.a[0] + u.a[1] * x + u.a[2] * y;
		u_a.y = u.a[3] + u.a[4] * x + u.a[5] * y;
		E_a.a[0] += Img_g->get(x, y).x * (*psiD)(Img_g->get(x, y).x * u_a.x + Img_g->get(x, y).y * u_a.y + Img_t->get(x, y), sigmaD);
		E_a.a[1] += Img_g->get(x, y).x * x * (*psiD)(Img_g->get(x, y).x * u_a.x + Img_g->get(x, y).y * u_a.y + Img_t->get(x, y), sigmaD);
		E_a.a[2] += Img_g->get(x, y).x * y * (*psiD)(Img_g->get(x, y).x * u_a.x + Img_g->get(x, y).y * u_a.y + Img_t->get(x, y), sigmaD);
		E_a.a[3] += Img_g->get(x, y).y * (*psiD)(Img_g->get(x, y).x * u_a.x + Img_g->get(x, y).y * u_a.y + Img_t->get(x, y), sigmaD);
		E_a.a[4] += Img_g->get(x, y).y * x * (*psiD)(Img_g->get(x, y).x * u_a.x + Img_g->get(x, y).y * u_a.y + Img_t->get(x, y), sigmaD);
		E_a.a[5] += Img_g->get(x, y).y * y * (*psiD)(Img_g->get(x, y).x * u_a.x + Img_g->get(x, y).y * u_a.y + Img_t->get(x, y), sigmaD);
	}
	return E_a;
}


VECTOR_AFFINE
sup_Error_aa_Block(const std::vector<VECTOR_2D<int> >& connected_domain, const ImgVector<VECTOR_2D<double> > *Img_g, double sigmaD)
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
	for (i = 0; i < connected_domain.size(); i++) {
		int x = connected_domain[i].x;
		int y = connected_domain[i].y;
		/* u = a0 + a1 * x + a2 * y */
		if (u_aa_max.a[0] < POW2(Img_g->get(x, y).x)) {
			u_aa_max.a[0] = POW2(Img_g->get(x, y).x);
		}
		if (u_aa_max.a[1] < POW2(Img_g->get(x, y).x * x)) {
			u_aa_max.a[1] = POW2(Img_g->get(x, y).x * x);
		}
		if (u_aa_max.a[2] < POW2(Img_g->get(x, y).x * y)) {
			u_aa_max.a[2] = POW2(Img_g->get(x, y).x * y);
		}
		/* v = a3 + a4 * x + a5 * y */
		if (u_aa_max.a[3] < POW2(Img_g->get(x, y).y)) {
			u_aa_max.a[3] = POW2(Img_g->get(x, y).y);
		}
		if (u_aa_max.a[4] < POW2(Img_g->get(x, y).y * x)) {
			u_aa_max.a[4] = POW2(Img_g->get(x, y).y * x);
		}
		if (u_aa_max.a[5] < POW2(Img_g->get(x, y).y * y)) {
			u_aa_max.a[5] = POW2(Img_g->get(x, y).y * y);
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
Error_Affine_Block(const VECTOR_AFFINE& u, const std::vector<VECTOR_2D<int> >& connected_domain, const ImgVector<VECTOR_2D<double> > *Img_g, const ImgVector<double> *Img_t, double sigmaD)
{
	double (*rhoD)(const double&, const double&) = Geman_McClure_rho;
	double E = 0.0;
	VECTOR_2D<double> u_a;

	for (int site = 0; site < connected_domain.size(); site++) {
		int x = connected_domain[site].x;
		int y = connected_domain[site].y;
		u_a.x = u.a[0] + u.a[1] * x + u.a[2] * y;
		u_a.y = u.a[3] + u.a[4] * x + u.a[5] * y;
		E += (*rhoD)(Img_g->get(x, y).x * u_a.x + Img_g->get(x, y).y * u_a.y + Img_t->get(x, y), sigmaD);
	}
	return E;
}

