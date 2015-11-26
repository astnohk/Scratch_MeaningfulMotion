/* Robust Estimation of Multiple Motions
 *
 * This motion estimation based on
 * M.J.Black and P.Anandan, "The Robust Estimation of Multiple Motions: Parametric and Piecewise-Smooth Flow Fields," Computer Vision and Image Understanding, Vol.63, No.1, 1996, pp.75-104.
 */

#include "OpticalFlow.h"
#include "OpticalFlow_BlockMatching.h"


#define SHOW_IRLS_OPTICALFLOW_E




// This function will compute INVERSE Optical Flow it points the previous frame which will come to the current (next) frame.
ImgVector<VECTOR_2D<double> > *
OpticalFlow_BlockMatching(const ImgVector<double>* It, const ImgVector<double>* Itp1, double MaxInt, MULTIPLE_MOTION_PARAM MotionParam, int IterMax)
{
	std::bad_alloc except_bad_alloc;

	ImgVector<VECTOR_2D<double> > *u = nullptr; // For RETURN value

	ImgVector<bool> domain_map;
	BlockMatching<double> block_matching;
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

	int BM_Search_Range = 41; // Block Matching search range
	int IterMax_level = 0;
	int MaxLevel = MotionParam.Level;
	int level;
	int i;

	if (It == nullptr) {
		throw std::invalid_argument("OpticalFlow_BlockMatching(const ImgVector<double>*, const ImgVector<double>* double, MULTIPLE_MOTION_PARAM, int) : const ImgVector<double>* It");
	} else if (Itp1 == nullptr) {
		throw std::invalid_argument("OpticalFlow_BlockMatching(const ImgVector<double>*, const ImgVector<double>* double, MULTIPLE_MOTION_PARAM, int) : const ImgVector<double>* Itp1");
	} else if (MaxInt < 0) {
		throw std::invalid_argument("OpticalFlow_BlockMatching(const ImgVector<double>*, const ImgVector<double>* double, MULTIPLE_MOTION_PARAM, int) : double MaxInt");
	}

	// Image Normalization
	It_normalize = *It;
	Itp1_normalize = *Itp1;
	for (i = 0; i < It_normalize.size(); i++) {
		It_normalize[i] /= MaxInt;
		Itp1_normalize[i] /= MaxInt;
	}

	// Adjust max level to use the Block Matching efficiently
	if (MaxLevel > floor(log((double)MotionParam.BlockMatching_BlockSize) / log(2.0))) {
		MaxLevel = (int)floor(log((double)MotionParam.BlockMatching_BlockSize) / log(2.0));
	}

	// ----- Block Matching -----
	int BlockSize = MotionParam.BlockMatching_BlockSize;
	domain_map.reset(It->width(), It->height());
	for (int y = 0; y < It->height(); y++) {
		for (int x = 0; x < It->width(); x++) {
			domain_map.at(x, y) = (int)(BlockSize * floor(y / BlockSize) + floor(x / BlockSize));
		}
	}
	block_matching.reset(It, Itp1, MotionParam.BlockMatching_BlockSize);
	block_matching.block_matching(BM_Search_Range);
	Motion_Vector.copy(block_matching.data());

	// ----- Optical Flow -----
	try {
		u = new ImgVector<VECTOR_2D<double> >(It->width(), It->height());
	}
	catch (const std::bad_alloc &bad) {
		except_bad_alloc = bad;
		goto ExitError;
	}
	try {
		u_levels = new ImgVector<VECTOR_2D<double> >[MaxLevel];
	}
	catch (const std::bad_alloc &bad) {
		except_bad_alloc = bad;
		goto ExitError;
	}
	// Make Pyramid
	try {
		It_levels = Pyramider(&It_normalize, MaxLevel);
	}
	catch (const std::bad_alloc &bad) {
		except_bad_alloc = bad;
		goto ExitError;
	}
	try {
		Itp1_levels = Pyramider(&Itp1_normalize, MaxLevel);
	}
	catch (const std::bad_alloc &bad) {
		except_bad_alloc = bad;
		goto ExitError;
	}
	// Derivative about time
	try {
		// The order reversed along with Block Matching (ordinary It -> Itp1)
		I_dt_levels = dt_Pyramid(Itp1_levels, It_levels, MaxLevel);
	}
	catch (const std::bad_alloc &bad) {
		except_bad_alloc = bad;
		goto ExitError;
	}
	// Derivative about space
	try {
		grad_It_levels = grad_Pyramid(Itp1_levels, nullptr, MaxLevel);
	}
	catch (const std::bad_alloc &bad) {
		except_bad_alloc = bad;
		goto ExitError;
	}

	// Initialize u_levels
	for (level = 0; level < MaxLevel; level++) {
		u_levels[level].reset(It_levels[level].width(), It_levels[level].height());
	}
	// Multi-Resolution IRLS Optical Flow estimation
	for (level = MaxLevel - 1; level >= 0; level--) {
		if (MaxLevel > 1) {
			sigmaD = sigmaD_init + (sigmaD_l0 - sigmaD_init) / (MaxLevel - 1.0) * (MaxLevel - 1.0 - level);
			sigmaS = sigmaS_init + (sigmaS_l0 - sigmaS_init) / (MaxLevel - 1.0) * (MaxLevel - 1.0 - level);
		} else {
			sigmaD = sigmaD_l0;
			sigmaS = sigmaS_l0;
		}
		printf("\nLevel %d : (1 / %d scaled, %dx%d)\n  sigmaD = %f\n  sigmaS = %f\n", level, (int)pow_int(2.0, level), u_levels[level].width(), u_levels[level].height(), sigmaD, sigmaS);
		if (level >= MaxLevel - 1) {
			// The order reversed (ordinary It -> Itp1)
			BM2OpticalFlow(I_dt_levels, u_levels, Itp1_levels, It_levels, level, &block_matching);
		} else {
			// The order reversed (ordinary It -> Itp1)
			LevelDown(I_dt_levels, u_levels, Itp1_levels, It_levels, level, MaxLevel);
		}
		IterMax_level = 4 * MAX(It->width(), It->height());
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
		Add_VectorOffset(u_levels, level, MaxLevel, &block_matching);
	}
	// Copy the lowest vector for output
	for (i = 0; i < u->size(); i++) {
		(*u)[i].x = u_levels[0][i].x;
		(*u)[i].y = u_levels[0][i].y;
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
	throw except_bad_alloc;
}


void
BM2OpticalFlow(ImgVector<double>* I_dt_levels, ImgVector<VECTOR_2D<double> >* u_levels, const ImgVector<double>* It_levels, const ImgVector<double>* Itp1_levels, int level, BlockMatching<double>* block_matching)
{
	double Scale = (double)It_levels[0].width() / It_levels[level].width();

	for (int y = 0; y < u_levels[level].height(); y++) {
		for (int x = 0; x < u_levels[level].width(); x++) {
			VECTOR_2D<double> u_offset = block_matching->get((int)floor(x * Scale / block_matching->block_size()), (int)floor(y * Scale / block_matching->block_size()));
			u_offset.x /= block_matching->vector_width();
			u_offset.y /= block_matching->vector_height();

			I_dt_levels[level].at(x, y) =
			    (Itp1_levels[level].get_zeropad(x + (int)floor(2.0 * u_offset.x), y + (int)floor(2.0 * u_offset.y))
			    - It_levels[level].get_zeropad(x, y)
			    + Itp1_levels[level].get_zeropad(x + 1 + (int)floor(2.0 * u_offset.x), y + (int)floor(2.0 * u_offset.y))
			    - It_levels[level].get_zeropad(x + 1, y)
			    + Itp1_levels[level].get_zeropad(x + (int)floor(2.0 * u_offset.x), y + 1 + (int)floor(2.0 * u_offset.y))
			    - It_levels[level].get_zeropad(x, y + 1)
			    + Itp1_levels[level].get_zeropad(x + 1 + (int)floor(2.0 * u_offset.x), y + 1 + (int)floor(2.0 * u_offset.y))
			    - It_levels[level].get_zeropad(x + 1, y + 1)) / 4.0;
			u_levels[level].at(x, y).x = 0.0;
			u_levels[level].at(x, y).y = 0.0;
		}
	}
}


void
Add_VectorOffset(ImgVector<VECTOR_2D<double> > *u_levels, int level, int MaxLevel, BlockMatching<double>* block_matching)
{
	if (level == MaxLevel - 1) {
		// Add offset calculated by using the motion vector by Block Matching
		double Scale = (double)u_levels[0].width() / u_levels[level].width();

		for (int y = 0; y < u_levels[level].height(); y++) {
			for (int x = 0; x < u_levels[level].width(); x++) {
				u_levels[level].at(x, y).x +=
				    block_matching->get((int)floor(x * Scale / block_matching->block_size()), (int)floor(y * Scale / block_matching->block_size())).x
				    / Scale;
				u_levels[level].at(x, y).y +=
				    block_matching->get((int)floor(x * Scale / block_matching->block_size()), (int)floor(y * Scale / block_matching->block_size())).y
				    / Scale;
			}
		}
	} else {
		// Add offset calculated by using the higher level's motion vector
		for (int y = 0; y < u_levels[level].height(); y++) {
			for (int x = 0; x < u_levels[level].width(); x++) {
				u_levels[level].at(x, y).x += u_levels[level + 1].get(x / 2, y / 2).x * 2.0;
				u_levels[level].at(x, y).y += u_levels[level + 1].get(x / 2, y / 2).y * 2.0;
			}
		}
	}
}


void
IRLS_OpticalFlow_Pyramid_Block(ImgVector<VECTOR_2D<double> > *u, const ImgVector<bool>& domain_map, const ImgVector<VECTOR_2D<double> > *Img_g, const ImgVector<double> *Img_t, double lambdaD, double lambdaS, double sigmaD, double sigmaS, int IterMax, double ErrorMinThreshold, int level)
{
	ERROR Error("IRLS_OpticalFlow_Pyramid_Block");

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
	} else if (Img_g == nullptr) {
		throw std::invalid_argument("ImgVector<VECTOR_2D<double> > *Img_g");
	} else if (Img_t == nullptr) {
		throw std::invalid_argument("ImgVector<double> *Img_t");
	}
	u_np1.copy(u); // Initialize u_np1
	// Reset sup_Error_uu max Img_g
	sup_Error_uu_Block(Img_g, lambdaD, lambdaS, sigmaD, sigmaS);
	sup = sup_Error_uu_Block(nullptr, lambdaD, lambdaS, sigmaD, sigmaS);
	for (n = 0; n < IterMax; n++) {
		// Calc for all sites
#pragma omp parallel for private(dE)
		for (site = 0; site < u->size(); site++) {
			dE = Error_u_Block(site, u, domain_map, Img_g, Img_t, lambdaD, lambdaS, sigmaD, sigmaS);
			u_np1[site].x = u->get(site).x - dE.x / sup.x;
			u_np1[site].y = u->get(site).y - dE.y / sup.y;
		}
		// Calc for all sites
		for (site = 0; site < u->size(); site++) {
			(*u)[site] = u_np1[site];
		}
		if (level == 0) {
			if ((n & 0x3F) == 0) {
				E = Error_MultipleMotion_Block(u, domain_map, Img_g, Img_t, lambdaD, lambdaS, sigmaD, sigmaS);
			}
		} else {
			E_prev = E;
			E = Error_MultipleMotion_Block(u, domain_map, Img_g, Img_t, lambdaD, lambdaS, sigmaD, sigmaS);
			if (E > E_prev) {
				ErrorIncrementCount++;
			} else {
				ErrorIncrementCount = 0;
			}
		}
#ifdef SHOW_IRLS_OPTICALFLOW_PYRAMID_E
		if ((n & 0x3F) == 0) {
			printf("E(%4d) = %e\n", n, E);
		}
#endif
		if (E < ErrorMinThreshold || ErrorIncrementCount > 3) {
			break;
		}
	}
}


VECTOR_2D<double>
Error_u_Block(int site, const ImgVector<VECTOR_2D<double> > *u, const ImgVector<bool>& domain_map, const ImgVector<VECTOR_2D<double> >* Img_g, const ImgVector<double>* Img_t, const double& lambdaD, const double& lambdaS, const double& sigmaD, const double& sigmaS)
{
	double (*psiD)(const double&, const double&) = Geman_McClure_psi;
	double (*psiS)(const double&, const double&) = Geman_McClure_psi;
	VECTOR_2D<double> us;
	double Center;
	VECTOR_2D<double> Neighbor;
	VECTOR_2D<double> E_u;
	int center_domain;
	int x, y;

	x = site % u->width();
	y = site / u->width();
	center_domain = domain_map.get(x, y);

	us = u->get(site);
	Center = (*psiD)(Img_g->get(site).x * us.x + Img_g->get(site).y * us.y + Img_t->get(site), sigmaD);

	Neighbor.x = .0;
	Neighbor.y = .0;
	if (x > 0 && domain_map.get(x - 1, y) == center_domain) {
		Neighbor.x += (*psiS)(us.x - u->get(x - 1, y).x, sigmaS);
		Neighbor.y += (*psiS)(us.y - u->get(x - 1, y).y, sigmaS);
	}
	if (x < u->width() - 1 && domain_map.get(x + 1, y) == center_domain) {
		Neighbor.x += (*psiS)(us.x - u->get(x + 1, y).x, sigmaS);
		Neighbor.y += (*psiS)(us.y - u->get(x + 1, y).y, sigmaS);
	}
	if (y > 0 && domain_map.get(x, y - 1) == center_domain) {
		Neighbor.x += (*psiS)(us.x - u->get(x, y - 1).x, sigmaS);
		Neighbor.y += (*psiS)(us.y - u->get(x, y - 1).y, sigmaS);
	}
	if (y < u->height() - 1 && domain_map.get(x, y + 1) == center_domain) {
		Neighbor.x += (*psiS)(us.x - u->get(x, y + 1).x, sigmaS);
		Neighbor.y += (*psiS)(us.y - u->get(x, y + 1).y, sigmaS);
	}

	E_u.x += lambdaD * Img_g->get(site).x * Center + lambdaS * Neighbor.x;
	E_u.y += lambdaD * Img_g->get(site).y * Center + lambdaS * Neighbor.y;
	return E_u;
}


VECTOR_2D<double>
sup_Error_uu_Block(const ImgVector<VECTOR_2D<double> > *Img_g, const double &lambdaD, const double &lambdaS, const double &sigmaD, const double &sigmaS)
{
	static VECTOR_2D<double> Img_g_max;
	VECTOR_2D<double> sup;

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
Error_MultipleMotion_Block(const ImgVector<VECTOR_2D<double> > *u, const ImgVector<bool>& domain_map, const ImgVector<VECTOR_2D<double> > *Img_g, const ImgVector<double> *Img_t, const double &lambdaD, const double &lambdaS, const double &sigmaD, const double &sigmaS)
{
	double (*rhoD)(const double&, const double&) = Geman_McClure_rho;
	double (*rhoS)(const double&, const double&) = Geman_McClure_rho;
	VECTOR_2D<double> us;
	double Center;
	VECTOR_2D<double> Neighbor;
	double E = 0.0;
	int x, y;

#pragma omp parallel for private(x, us, Neighbor, Center) reduction(+:E)
	for (y = 0; y < u->height(); y++) {
		for (x = 0; x < u->width(); x++) {
			int center_domain = domain_map.get(x, y);
			us = u->get(x, y);
			Neighbor.x = .0;
			Neighbor.y = .0;
			if (x > 0 && domain_map.get(x - 1, y) == center_domain) {
				Neighbor.x += (*rhoS)(us.x - u->get(x - 1, y).x, sigmaS);
				Neighbor.y += (*rhoS)(us.y - u->get(x - 1, y).y, sigmaS);
			}
			if (x < u->width() - 1 && domain_map.get(x + 1, y) == center_domain) {
				Neighbor.x += (*rhoS)(us.x - u->get(x + 1, y).x, sigmaS);
				Neighbor.y += (*rhoS)(us.y - u->get(x + 1, y).y, sigmaS);
			}
			if (y > 0 && domain_map.get(x, y - 1) == center_domain) {
				Neighbor.x += (*rhoS)(us.x - u->get(x, y - 1).x, sigmaS);
				Neighbor.y += (*rhoS)(us.y - u->get(x, y - 1).y, sigmaS);
			}
			if (y < u->height() - 1 && domain_map.get(x, y + 1) == center_domain) {
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

