/* Robust Estimation of Multiple Motions
 *
 * This motion estimation based on
 * M.J.Black and P.Anandan, "The Robust Estimation of Multiple Motions: Parametric and Piecewise-Smooth Flow Fields," Computer Vision and Image Understanding, Vol.63, No.1, 1996, pp.75-104.
 */

#include "OpticalFlow.h"



// define for debug

#define SHOW_IRLS_OPTICALFLOW_PYRAMID_E
//#define DEBUG_STOP_ON_LEVEL_L

// /define for debug




ImgVector<VECTOR_2D<double> > *
OpticalFlow_Pyramid(ImgVector<double> *It, ImgVector<double> *Itp1, double MaxInt, MULTIPLE_MOTION_PARAM MotionParam, int IterMax)
{
	ERROR Error("OpticalFlow_Pyramid");

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
	ImgVector<VECTOR_2D<double> > *u = nullptr; // For RETURN value
	ImgVector<VECTOR_2D<double> > *u_levels = nullptr;
	ImgVector<double> *I_dt_levels = nullptr;
	ImgVector<double> *It_levels = nullptr;
	ImgVector<double> *Itp1_levels = nullptr;
	ImgVector<VECTOR_2D<double> > *grad_It_levels = nullptr;

	int IterMax_level = 0;
	int MaxLevel = MotionParam.Level;
	int level;

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
	for (size_t i = 0; i < It_normalize.size(); i++) {
		It_normalize[i] /= MaxInt;
		Itp1_normalize[i] /= MaxInt;
	}
	// Multiple Motion Vectors
	try {
		u = new ImgVector<VECTOR_2D<double> >(It->width(), It->height());
	}
	catch (const std::bad_alloc &bad) {
		std::cerr << bad.what() << std::endl;
		Error.Value("u");
		Error.Malloc();
		goto ExitError;
	}
	try {
		u_levels = new ImgVector<VECTOR_2D<double> >[MaxLevel + 1];
	}
	catch (const std::bad_alloc &bad) {
		std::cerr << bad.what() << std::endl;
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

	for (level = MaxLevel; level >= 0; level--) {
		if (MaxLevel > 0) {
			sigmaD = sigmaD_init + (sigmaD_l0 - sigmaD_init) / MaxLevel * (MaxLevel - level);
			sigmaS = sigmaS_init + (sigmaS_l0 - sigmaS_init) / MaxLevel * (MaxLevel - level);
		} else {
			sigmaD = sigmaD_l0;
			sigmaS = sigmaS_l0;
		}
		u_levels[level].reset(I_dt_levels[level].width(), I_dt_levels[level].height());
		printf("\nLevel %d : (1 / %d scaled, %dx%d)\n  sigmaD = %f\n  sigmaS = %f\n", level, int(pow_int(2.0, level)), u_levels[level].width(), u_levels[level].height(), sigmaD, sigmaS);
		if (level < MaxLevel) {
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
		Add_VectorOffset(u_levels, level, MaxLevel);
	}
	// Copy the lowest vector for output
	for (size_t i = 0; i < u->size(); i++) {
		u->at(i).x = u_levels[0][i].x;
		u->at(i).y = u_levels[0][i].y;
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
	delete u;
	return nullptr;
}


void
LevelDown(ImgVector<double> *I_dt_levels, ImgVector<VECTOR_2D<double> > *u_levels, const ImgVector<double> *It_levels, const ImgVector<double> *Itp1_levels, int level, int MaxLevel)
{
	if (level == MaxLevel) {
		// Do NOT need projection from higher level
		return;
	}
	for (int y = 0; y < u_levels[level].height(); y++) {
		for (int x = 0; x < u_levels[level].width(); x++) {
			VECTOR_2D<double> u_offset = u_levels[level + 1].get(x / 2, y / 2);

			I_dt_levels[level].at(x, y) =
			    (Itp1_levels[level].get_zeropad(x + int(floor(2.0 * u_offset.x)), y + int(floor(2.0 * u_offset.y)))
			    - It_levels[level].get_zeropad(x, y)
			    + Itp1_levels[level].get_zeropad(x + 1 + int(floor(2.0 * u_offset.x)), y + int(floor(2.0 * u_offset.y)))
			    - It_levels[level].get_zeropad(x + 1, y)
			    + Itp1_levels[level].get_zeropad(x + int(floor(2.0 * u_offset.x)), y + 1 + int(floor(2.0 * u_offset.y)))
			    - It_levels[level].get_zeropad(x, y + 1)
			    + Itp1_levels[level].get_zeropad(x + 1 + int(floor(2.0 * u_offset.x)), y + 1 + int(floor(2.0 * u_offset.y)))
			    - It_levels[level].get_zeropad(x + 1, y + 1)) / 4.0;
			u_levels[level].at(x, y).x = 0.0;
			u_levels[level].at(x, y).y = 0.0;
		}
	}
}


void
Add_VectorOffset(ImgVector<VECTOR_2D<double> > *u_levels, int level, int MaxLevel)
{
	if (level == MaxLevel) {
		// Do NOT need projection from higher level
		return;
	}
	// Add offset calculated by using the higher level's motion vector
	for (int y = 0; y < u_levels[level].height(); y++) {
		for (int x = 0; x < u_levels[level].width(); x++) {
			u_levels[level].at(x, y).x += u_levels[level + 1].get(x / 2, y / 2).x * 2.0;
			u_levels[level].at(x, y).y += u_levels[level + 1].get(x / 2, y / 2).y * 2.0;
		}
	}
}


void
IRLS_OpticalFlow_Pyramid(ImgVector<VECTOR_2D<double> > *u, const ImgVector<VECTOR_2D<double> > *Img_g, const ImgVector<double> *Img_t, double lambdaD, double lambdaS, double sigmaD, double sigmaS, int IterMax, double ErrorMinThreshold, int level)
{
	ERROR Error("IRLS_OpticalFlow_Pyramid");

	if (u == nullptr) {
		throw std::invalid_argument("ImgVector<VECTOR_2D<double> > *u");
	} else if (Img_g == nullptr) {
		throw std::invalid_argument("ImgVector<VECTOR_2D<double> > *Img_g");
	} else if (Img_t == nullptr) {
		throw std::invalid_argument("ImgVector<double> *Img_t");
	}
	ImgVector<VECTOR_2D<double> > u_np1(*u);
	// Reset sup_Error_uu max Img_g
	sup_Error_uu(Img_g, lambdaD, lambdaS, sigmaD, sigmaS);
	VECTOR_2D<double> sup;
	sup = sup_Error_uu(nullptr, lambdaD, lambdaS, sigmaD, sigmaS);
	double E = 0.0;
	int ErrorIncrementCount = 0;
	for (int n = 0; n < IterMax; n++) {
		// Calc for all sites
		size_t site;
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (site = 0; site < u->size(); site++) {
			VECTOR_2D<double> dE;
			dE = Error_u(site, u, Img_g, Img_t, lambdaD, lambdaS, sigmaD, sigmaS);
			u_np1[site].x = u->get(site).x - dE.x / sup.x;
			u_np1[site].y = u->get(site).y - dE.y / sup.y;
		}
		// Calc for all sites
		for (site = 0; site < u->size(); site++) {
			u->at(site) = u_np1[site];
		}
		if (level == 0) {
			if ((n & 0x3F) == 0) {
				E = Error_MultipleMotion(u, Img_g, Img_t, lambdaD, lambdaS, sigmaD, sigmaS);
			}
		} else {
			double E_prev = E;
			E = Error_MultipleMotion(u, Img_g, Img_t, lambdaD, lambdaS, sigmaD, sigmaS);
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
Error_u(const size_t& site, const ImgVector<VECTOR_2D<double> > *u, const ImgVector<VECTOR_2D<double> > *Img_g, const ImgVector<double> *Img_t, const double &lambdaD, const double &lambdaS, const double &sigmaD, const double &sigmaS)
{
	double (*psiD)(const double&, const double&) = Geman_McClure_psi;
	double (*psiS)(const double&, const double&) = Geman_McClure_psi;
	VECTOR_2D<double> us;
	double Center;
	VECTOR_2D<double> E_u;

	int x = static_cast<int>(site % size_t(u->width()));
	int y = static_cast<int>(site / size_t(u->width()));

	us = u->get(site);
	Center = (*psiD)(Img_g->get(site).x * us.x + Img_g->get(site).y * us.y + Img_t->get(site), sigmaD);

	VECTOR_2D<double> Neighbor(0.0, 0.0);
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


VECTOR_2D<double>
sup_Error_uu(const ImgVector<VECTOR_2D<double> > *Img_g, const double &lambdaD, const double &lambdaS, const double &sigmaD, const double &sigmaS)
{
	static VECTOR_2D<double> Img_g_max;
	VECTOR_2D<double> sup;

	if (Img_g != nullptr) {
		Img_g_max.reset();
		for (size_t i = 0; i < Img_g->size(); i++) {
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
Error_MultipleMotion(const ImgVector<VECTOR_2D<double> > *u, const ImgVector<VECTOR_2D<double> > *Img_g, const ImgVector<double> *Img_t, const double &lambdaD, const double &lambdaS, const double &sigmaD, const double &sigmaS)
{
	double (*rhoD)(const double&, const double&) = Geman_McClure_rho;
	double (*rhoS)(const double&, const double&) = Geman_McClure_rho;
	VECTOR_2D<double> us;
	double Center;
	VECTOR_2D<double> Neighbor;
	double E = 0.0;
	int x, y;

#ifdef _OPENMP
#pragma omp parallel for private(x, us, Neighbor, Center) reduction(+:E)
#endif
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


void
MultipleMotion_write(const ImgVector<double>& img_prev, const ImgVector<double>& img_next, const ImgVector<VECTOR_2D<double> >& u, const std::string& filename)
{
	ERROR Error("MultipleMotion_write");

	FILE *fp = nullptr;
	VECTOR_2D<double> v;
	int x, y;
	MotionCompensation<double> compensated(img_prev, img_next, u);
	PNM pnm;
	std::string filename_compensated;

	printf("\n* Output The Optical Flow to '%s'(binary)\n", filename.c_str());
	if ((fp = fopen(filename.c_str(), "wb")) == nullptr) {
		Error.Function("fopen");
		Error.File(filename.c_str());
		Error.FileWrite();
		throw std::logic_error("fopen");
	}
	fprintf(fp, "%d %d\n", u.width(), u.height());
	for (y = 0; y < u.height(); y++) {
		for (x = 0; x < u.width(); x++) {
			v = u.get(x, y);
			if (fwrite(&v.x, sizeof(double), 1, fp) < 1) {
				Error.Function("fwrite");
				Error.Value("u(x, y).x");
				Error.FunctionFail();
				throw std::logic_error("fwrite");
			}
			if (fwrite(&v.y, sizeof(double), 1, fp) < 1) {
				Error.Function("fwrite");
				Error.Value("u(x, y).y");
				Error.FunctionFail();
				throw std::logic_error("fwrite");
			}
		}
	}
	fclose(fp);

	compensated.create_image_compensated(); // Make compensated image
	std::string::size_type found = filename.find_last_of("/\\");
	filename_compensated = filename.substr(0, found + 1) + "compensated_" + filename.substr(found + 1);
	printf("* Output The Compensated Image from Optical Flow to '%s'(binary)\n\n", filename_compensated.c_str());
	pnm.copy(PORTABLE_GRAYMAP_BINARY, compensated.width(), compensated.height(), 255, compensated.ref_image_compensated().data(), 1.0);
	pnm.write(filename_compensated.c_str());
	pnm.free();
}

void
MultipleMotion_write(const ImgVector<ImgClass::RGB>& img_prev, const ImgVector<ImgClass::RGB>& img_next, const ImgVector<VECTOR_2D<double> >& u, const std::string &filename)
{
	ERROR Error("MultipleMotion_write");

	FILE *fp = nullptr;
	VECTOR_2D<double> v;
	int x, y;
	MotionCompensation<ImgClass::RGB> compensated(img_prev, img_next, u);
	PNM pnm;
	std::string filename_compensated;

	printf("\n* Output The Optical Flow to '%s'(binary)\n", filename.c_str());
	if ((fp = fopen(filename.c_str(), "wb")) == nullptr) {
		Error.Function("fopen");
		Error.File(filename.c_str());
		Error.FileWrite();
		throw std::logic_error("fopen");
	}
	fprintf(fp, "%d %d\n", u.width(), u.height());
	for (y = 0; y < u.height(); y++) {
		for (x = 0; x < u.width(); x++) {
			v = u.get(x, y);
			if (fwrite(&v.x, sizeof(double), 1, fp) < 1) {
				Error.Function("fwrite");
				Error.Value("u(x, y).x");
				Error.FunctionFail();
				throw std::logic_error("fwrite");
			}
			if (fwrite(&v.y, sizeof(double), 1, fp) < 1) {
				Error.Function("fwrite");
				Error.Value("u(x, y).y");
				Error.FunctionFail();
				throw std::logic_error("fwrite");
			}
		}
	}
	fclose(fp);

	compensated.create_image_compensated(); // Make compensated image
	int* compensated_image = nullptr;
	size_t size = compensated.ref_image_compensated().size();
	std::string::size_type found = filename.find_last_of("/\\");
	filename_compensated = filename.substr(0, found + 1) + "compensated_" + filename.substr(found + 1);
	printf("* Output The Compensated Image from Optical Flow to '%s'(binary)\n\n", filename_compensated.c_str());
	try {
		compensated_image = new int[size * 3];
	}
	catch (const std::bad_alloc& bad) {
		std::cerr << bad.what() << std::endl;
		fprintf(stderr, "void MultipleMotion_write(const ImgVector<ImgClass::RGB>*, const ImgVector<ImgClass::RGB>*, const ImgVector<VECTOR_2D<double> >*, const std::string&) : Cannot allocate memory\n");
		throw;
	}
	for (size_t n = 0; n < size; n++) {
		compensated_image[n] = int(compensated.ref_image_compensated().get(n).R);
		compensated_image[n + size] = int(compensated.ref_image_compensated().get(n).G);
		compensated_image[n + 2 * size] = int(compensated.ref_image_compensated().get(n).B);
	}
	pnm.copy(PORTABLE_PIXMAP_BINARY, compensated.width(), compensated.height(), 255, compensated_image);
	pnm.write(filename_compensated.c_str());
	pnm.free();
}

