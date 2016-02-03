/* Robust Estimation of Multiple Motions
 *
 * This motion estimation based on
 * M.J.Black and P.Anandan, "The Robust Estimation of Multiple Motions: Parametric and Piecewise-Smooth Flow Fields," Computer Vision and Image Understanding, Vol.63, No.1, 1996, pp.75-104.
 */
#include "Affine_BlockMatching.h"




ImgVector<VECTOR_2D<double> >
AffineParametric(const ImgVector<ImgClass::Lab>& reference, const ImgVector<ImgClass::Lab>& interest, const ImgVector<VECTOR_2D<double> >& MV, const std::vector<std::vector<VECTOR_2D<int> > >& regions, const MULTIPLE_MOTION_PARAM& MotionParam, const int IterMax)
{
	// Output MV array
	ImgVector<VECTOR_2D<double> > u(interest.width(), interest.height());
	// M-estimator parameter
	const double sigma = 0.2 / sqrt(2.0); //4.0 / sqrt(2.0);

	// ----- Parametric Motion -----
	std::vector<std::vector<double> > u_affine;
	// Compute gradients
	ImgVector<VECTOR_2D<double> > grad(interest.width(), interest.height());
	for (int y = 0; y < interest.height(); y++) {
		for (int x = 0; x < interest.width(); x++) {
			grad.at(x, y).x =
			    (interest.get_mirror(x + 1, y).L - interest.get(x, y).L
			    + interest.get_mirror(x + 1, y + 1).L - interest.get_mirror(x, y + 1).L)
			    / 2.0;
			grad.at(x, y).y =
			    (interest.get_mirror(x, y + 1).L - interest.get(x, y).L
			    + interest.get_mirror(x + 1, y + 1).L - interest.get_mirror(x + 1, y).L)
			    / 2.0;
		}
	}
	// Initialize affine coefficient vector
	u_affine.resize(regions.size());
	for (size_t n = 0; n < regions.size(); n++) {
		u_affine[n].resize(NUM_AFFINE_PARAMETER);
	}
	// IRLS affine parametric motion estimation
	printf("\n* IRLS\n    sigma = %f, iteration = %d\n", sigma, IterMax);
	ImgVector<double> dt(interest.width(), interest.height());
	for (int y = 0; y < interest.height(); y++) {
		for (int x = 0; x < interest.width(); x++) {
			VECTOR_2D<int> v;
			v.x = int(floor(MV.get(x, y).x));
			v.y = int(floor(MV.get(x, y).y));
			dt.at(x, y) =
			    (reference.get_mirror(x + v.x, y + v.y).L - interest.get(x, y).L
			    + reference.get_mirror(x + v.x + 1, y + v.y).L - interest.get_mirror(x + 1, y).L
			    + reference.get_mirror(x + v.x, y + v.y + 1).L - interest.get_mirror(x, y + 1).L
			    + reference.get_mirror(x + v.x + 1, y + v.y + 1).L - interest.get_mirror(x + 1, y + 1).L)
			    / 4.0;
		}
	}
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for (size_t n = 0; n < regions.size(); n++) {
		IRLS_AffineParametric_region(
		    &u_affine[n],
		    regions[n],
		    grad,
		    dt,
		    sigma,
		    IterMax,
		    MotionParam.Error_Min_Threshold);
		// Compute each pixel's vector
		for (const VECTOR_2D<int>& r : regions[n]) {
			u.at(r.x, r.y).x =
			    u_affine[n][0] + r.x * u_affine[n][1] + r.y * u_affine[n][2];
			u.at(r.x, r.y).y =
			    u_affine[n][3] + r.x * u_affine[n][4] + r.y * u_affine[n][5];
		}
	}
	return u;
}




void
IRLS_AffineParametric_region(std::vector<double>* u_affine, const std::vector<VECTOR_2D<int> >& region, const ImgVector<VECTOR_2D<double> >& grad, const ImgVector<double>& dt, const double& sigma, const int IterMax, const double& ErrorMinThreshold)
{
	const double omega = 1.0E-0;

	// Initialize
	for (size_t i = 0; i < NUM_AFFINE_PARAMETER; i++) {
		u_affine->at(i) = .0;
	}
	// Start IRLS
	for (int n = 0; n < IterMax; n++) {
		std::vector<double> u_np1(NUM_AFFINE_PARAMETER);
		std::vector<double> dE = Error_a_region(*u_affine, region, grad, dt, sigma);
		std::vector<double> sup = sup_Error_aa_region(region, grad, sigma);
		double E = 0.0;

		for (size_t i = 0; i < NUM_AFFINE_PARAMETER; i++) {
			if (fabs(sup[i]) < 1.0E-10) {
				u_np1[i] = u_affine->at(i) - omega * 1.0E+10 * SIGN_NOZERO(sup[i]) * dE[i];
			} else {
				u_np1[i] = u_affine->at(i) - omega / sup[i] * dE[i];
			}
		}
		*u_affine = u_np1;
		E = Error_Affine_region(*u_affine, region, grad, dt, sigma);
		if (E < ErrorMinThreshold) {
			break;
		}
	}
}


std::vector<double>
Error_a_region(const std::vector<double>& u_affine, const std::vector<VECTOR_2D<int> >& region, const ImgVector<VECTOR_2D<double> >& grad, const ImgVector<double>& dt, const double& sigma)
{
	double (*psiD)(const double&, const double&) = Geman_McClure_psi;
	std::vector<double> E_a(NUM_AFFINE_PARAMETER);
	VECTOR_2D<double> u_a;

	for (size_t i = 0; i < NUM_AFFINE_PARAMETER; i++) {
		E_a[i] = .0;
	}
	for (const VECTOR_2D<int>& r : region) {
		u_a.x = u_affine[0] + u_affine[1] * r.x + u_affine[2] * r.y;
		u_a.y = u_affine[3] + u_affine[4] * r.x + u_affine[5] * r.y;
		E_a[0] += grad.get(r.x, r.y).x * (*psiD)(grad.get(r.x, r.y).x * u_a.x + grad.get(r.x, r.y).y * u_a.y + dt.get_zeropad(r.x, r.y), sigma);
		E_a[1] += grad.get(r.x, r.y).x * r.x * (*psiD)(grad.get(r.x, r.y).x * u_a.x + grad.get(r.x, r.y).y * u_a.y + dt.get_zeropad(r.x, r.y), sigma);
		E_a[2] += grad.get(r.x, r.y).x * r.y * (*psiD)(grad.get(r.x, r.y).x * u_a.x + grad.get(r.x, r.y).y * u_a.y + dt.get_zeropad(r.x, r.y), sigma);
		E_a[3] += grad.get(r.x, r.y).y * (*psiD)(grad.get(r.x, r.y).x * u_a.x + grad.get(r.x, r.y).y * u_a.y + dt.get_zeropad(r.x, r.y), sigma);
		E_a[4] += grad.get(r.x, r.y).y * r.x * (*psiD)(grad.get(r.x, r.y).x * u_a.x + grad.get(r.x, r.y).y * u_a.y + dt.get_zeropad(r.x, r.y), sigma);
		E_a[5] += grad.get(r.x, r.y).y * r.y * (*psiD)(grad.get(r.x, r.y).x * u_a.x + grad.get(r.x, r.y).y * u_a.y + dt.get_zeropad(r.x, r.y), sigma);
	}
	return E_a;
}


std::vector<double>
sup_Error_aa_region(const std::vector<VECTOR_2D<int> >& region, const ImgVector<VECTOR_2D<double> >& grad, const double& sigma)
{
	ERROR Error("sup_Error_aa_region");

	std::vector<double> sup(NUM_AFFINE_PARAMETER);
	std::vector<double> u_aa_max(NUM_AFFINE_PARAMETER);

	for (size_t i = 0; i < NUM_AFFINE_PARAMETER; i++) {
		u_aa_max[i] = .0;
	}
	for (const VECTOR_2D<int>& r : region) {
		// u = a0 + a1 * x + a2 * y
		if (u_aa_max[0] < POW2(grad.get(r.x, r.y).x)) {
			u_aa_max[0] = POW2(grad.get(r.x, r.y).x);
		}
		if (u_aa_max[1] < POW2(grad.get(r.x, r.y).x * r.x)) {
			u_aa_max[1] = POW2(grad.get(r.x, r.y).x * r.x);
		}
		if (u_aa_max[2] < POW2(grad.get(r.x, r.y).x * r.y)) {
			u_aa_max[2] = POW2(grad.get(r.x, r.y).x * r.y);
		}
		// v = a3 + a4 * x + a5 * y
		if (u_aa_max[3] < POW2(grad.get(r.x, r.y).y)) {
			u_aa_max[3] = POW2(grad.get(r.x, r.y).y);
		}
		if (u_aa_max[4] < POW2(grad.get(r.x, r.y).y * r.x)) {
			u_aa_max[4] = POW2(grad.get(r.x, r.y).y * r.x);
		}
		if (u_aa_max[5] < POW2(grad.get(r.x, r.y).y * r.y)) {
			u_aa_max[5] = POW2(grad.get(r.x, r.y).y * r.y);
		}
	}
	sup[0] = u_aa_max[0] * 2.0 / POW2(sigma);
	sup[1] = u_aa_max[1] * 2.0 / POW2(sigma);
	sup[2] = u_aa_max[2] * 2.0 / POW2(sigma);
	sup[3] = u_aa_max[3] * 2.0 / POW2(sigma);
	sup[4] = u_aa_max[4] * 2.0 / POW2(sigma);
	sup[5] = u_aa_max[5] * 2.0 / POW2(sigma);
	return sup;
}


double
Error_Affine_region(const std::vector<double>& u_affine, const std::vector<VECTOR_2D<int> >& region, const ImgVector<VECTOR_2D<double> >& grad, const ImgVector<double>& dt, const double& sigma)
{
	double (*rhoD)(const double&, const double&) = Geman_McClure_rho;
	double E = 0.0;
	VECTOR_2D<double> u_a;

	for (const VECTOR_2D<int>& r : region) {
		u_a.x = u_affine[0] + u_affine[1] * r.x + u_affine[2] * r.y;
		u_a.y = u_affine[3] + u_affine[4] * r.x + u_affine[5] * r.y;
		E += (*rhoD)(
		    grad.get(r.x, r.y).x * u_a.x
		    + grad.get(r.x, r.y).y * u_a.y
		    + dt.get(r.x, r.y),
		    sigma);
	}
	return E;
}

