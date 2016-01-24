/* Robust Estimation of Multiple Motions
 *
 * This motion estimation based on
 * M.J.Black and P.Anandan, "The Robust Estimation of Multiple Motions: Parametric and Piecewise-Smooth Flow Fields," Computer Vision and Image Understanding, Vol.63, No.1, 1996, pp.75-104.
 */
#include "OpticalFlow_BlockMatching.h"




// This function will compute INVERSE Optical Flow it points the previous frame which will come to the current (next) frame.
std::vector<ImgVector<Vector_ST<double> > >
OpticalFlow_BlockMatching(const ImgVector<ImgClass::RGB>& It_color, const ImgVector<ImgClass::RGB>& Itp1_color, double MaxInt, MULTIPLE_MOTION_PARAM MotionParam, const std::string ofilename, int IterMax)
{
	const bool Bidirectional_with_Time = true; // on almost all cases it is true
	const size_t History_Max = 4;
	static std::deque<ImgVector<ImgClass::RGB> > sequence_sRGB;
	static std::deque<ImgVector<double> > sequence_Grayscale;
	static std::deque<ImgVector<ImgClass::Lab> > sequence_Lab;
	static std::deque<Segmentation<ImgClass::Lab> > segmentations;

	std::bad_alloc except_bad_alloc;

	std::vector<ImgVector<Vector_ST<double> > > u; // For RETURN value

	ImgVector<size_t> domain_map;
	const double coeff_MAD = 1.0;
	const double coeff_ZNCC = 0.0;
	BlockMatching<ImgClass::Lab> block_matching;
	int BM_Search_Range = 61; // Block Matching search range
	int Subpixel_Scale = 2;

	ImgVector<double> It;
	ImgVector<double> Itp1;
	ImgVector<ImgClass::RGB> It_sRGB_normalize;
	ImgVector<ImgClass::RGB> Itp1_sRGB_normalize;
	ImgVector<ImgClass::Lab> It_Lab_normalize;
	ImgVector<ImgClass::Lab> Itp1_Lab_normalize;
	ImgVector<double> It_normalize;
	ImgVector<double> Itp1_normalize;
	ImgVector<VECTOR_2D<double> >* u_levels = nullptr;
	ImgVector<double>* I_dt_levels = nullptr;
	ImgVector<double>* It_levels = nullptr;
	ImgVector<double>* Itp1_levels = nullptr;
	ImgVector<VECTOR_2D<double> >* grad_It_levels = nullptr;
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

	if (It_color.isNULL()) {
		throw std::invalid_argument("OpticalFlow_BlockMatching(const ImgVector<double>*, const ImgVector<double>* double, MULTIPLE_MOTION_PARAM, int) : const ImgVector<double>* It");
	} else if (Itp1_color.isNULL()) {
		throw std::invalid_argument("OpticalFlow_BlockMatching(const ImgVector<double>*, const ImgVector<double>* double, MULTIPLE_MOTION_PARAM, int) : const ImgVector<double>* Itp1");
	} else if (MaxInt < 0) {
		throw std::invalid_argument("OpticalFlow_BlockMatching(const ImgVector<double>*, const ImgVector<double>* double, MULTIPLE_MOTION_PARAM, int) : double MaxInt");
	}

	// sRGB image
	It_sRGB_normalize.copy(It_color);
	Itp1_sRGB_normalize.copy(Itp1_color);
	// Grayscale
	It.cast_copy(It_color);
	Itp1.cast_copy(Itp1_color);
	// Image Normalization
	It_normalize = It;
	Itp1_normalize = Itp1;
	for (size_t i = 0; i < It_normalize.size(); i++) {
		// sRGB
		It_sRGB_normalize[i] /= MaxInt;
		Itp1_sRGB_normalize[i] /= MaxInt;
		// Grayscale
		It_normalize[i] /= MaxInt;
		Itp1_normalize[i] /= MaxInt;
	}
	// Convert sRGB to CIE Lab
	It_Lab_normalize.reset(It_sRGB_normalize.width(), It_sRGB_normalize.height());
	Itp1_Lab_normalize.reset(It_sRGB_normalize.width(), It_sRGB_normalize.height());
	for (size_t i = 0; i < It_sRGB_normalize.size(); i++) {
		It_Lab_normalize[i].set(It_sRGB_normalize[i]);
		Itp1_Lab_normalize[i].set(Itp1_sRGB_normalize[i]);
	}

	// Shift image sequence by 1 and assign current image
	if (sequence_sRGB.empty()) {
		sequence_sRGB.push_front(It_sRGB_normalize);
		sequence_Grayscale.push_front(It_normalize);
		sequence_Lab.push_front(It_Lab_normalize);
	}
	sequence_sRGB.push_front(Itp1_sRGB_normalize);
	sequence_Grayscale.push_front(Itp1_normalize);
	sequence_Lab.push_front(Itp1_Lab_normalize);
	if (sequence_sRGB.size() >= History_Max) {
		sequence_sRGB.pop_back();
		sequence_Grayscale.pop_back();
		sequence_Lab.pop_back();
	}

	// Adjust max level to use the Block Matching efficiently
	if (MaxLevel > floor(log(double(MotionParam.BlockMatching_BlockSize)) / log(2.0))) {
		MaxLevel = int(floor(log(double(MotionParam.BlockMatching_BlockSize)) / log(2.0)));
	}

	// ----- Block Matching -----
#if 0
	// Normal Block Matching
	printf("* * Compute Block Matching\n");
	int BlockSize = MotionParam.BlockMatching_BlockSize;
	domain_map.reset(It.width(), It.height());
	for (int y = 0; y < It.height(); y++) {
		for (int x = 0; x < It.width(); x++) {
			domain_map.at(x, y) = size_t(BlockSize * floor(y / BlockSize) + floor(x / BlockSize));
		}
	}
	int BlockMatching_BlockSize = 8;
	if (sequence_Lab.size() <= 2) {
		block_matching.reset(It_Lab_normalize, Itp1_Lab_normalize, BlockMatching_BlockSize, Subpixel_Scale);
	} else {
		block_matching.reset(sequence_Lab[2], sequence_Lab[1], sequence_Lab[0], BlockMatching_BlockSize, Subpixel_Scale);
	}
	block_matching.block_matching(BM_Search_Range, coeff_MAD, coeff_ZNCC);
#else
	{
		// Segmentation
		printf("* * Compute Segmentation by Mean Shift\n");

#ifdef MEANSHIFT_KERNEL_SPATIAL
		double kernel_spatial = MEANSHIFT_KERNEL_SPATIAL, kernel_intensity = 8.0 / 255.0; // for images under about HD resolution
#else
		//double kernel_spatial = 64.0, kernel_intensity = 12.0 / 255.0; // for 4K Film kernel(spatial = 64.0, intensity = 12.0 / 255.0)
		double kernel_spatial = 8.0, kernel_intensity = 8.0 / 255.0; // for images under about HD resolution
#endif

		if (segmentations.empty()) {
			segmentations.push_front(Segmentation<ImgClass::Lab>(It_Lab_normalize, kernel_spatial, kernel_intensity));
		}
		segmentations.push_front(Segmentation<ImgClass::Lab>(Itp1_Lab_normalize, kernel_spatial, kernel_intensity));
		if (segmentations.size() >= History_Max) {
			segmentations.pop_back();
		}

		PNM pnm;
		std::string::size_type found = 1 + ofilename.find_last_not_of("0123456789", ofilename.find_last_of("0123456789"));
		if (found == std::string::npos) {
			found = ofilename.find_last_of(".");
		}
		std::string ofilename_segmentation = ofilename.substr(0, found) + "segmentation" + ofilename.substr(found);
		printf("* Output The Segmentation result to '%s'(binary)\n\n", ofilename_segmentation.c_str());
		{
			ImgVector<int> tmp_vector(segmentations[0].width(), segmentations[0].height());
			for (size_t i = 0; i < segmentations[0].size(); i++) {
				tmp_vector[i] = static_cast<int>(segmentations[0][i]);
			}
			pnm.copy(PORTABLE_GRAYMAP_BINARY, segmentations[0].width(), segmentations[0].height(), int(tmp_vector.max()), tmp_vector.data());
			pnm.write(ofilename_segmentation.c_str());
			pnm.free();
		}

		ImgVector<int> quantized(segmentations[0].width(), segmentations[0].height());
		for (size_t i = 0; i < segmentations[0].ref_color_quantized_image().size(); i++) {
			quantized.at(i) = int(round(double(segmentations[0].ref_color_quantized_image().get(i)) / 100.0));
			if (quantized.get(i) < 0) {
				quantized.at(i) = 0;
			}
		}
		std::string ofilename_quantized = ofilename.substr(0, found) + "color-quantized_" + ofilename.substr(found);
		printf("* Output The color quantized image '%s'(binary)\n\n", ofilename_quantized.c_str());
		pnm.copy(PORTABLE_GRAYMAP_BINARY, segmentations[0].width(), segmentations[0].height(), quantized.max(), quantized.data());
		pnm.write(ofilename_quantized.c_str());
		pnm.free();
		// Output vectors
		std::string ofilename_vector = ofilename.substr(0, found) + "shift-vector" + ofilename.substr(found);
		FILE *fp;
		fp = fopen(ofilename_vector.c_str(), "w");
		fprintf(fp, "%d %d\n", segmentations[0].width(), segmentations[0].height());
		for (int y = 0; y < segmentations[0].height(); y++) {
			for (int x = 0; x < segmentations[0].width(); x++) {
				VECTOR_2D<double> v;
				v.x = segmentations[0].ref_shift_vector().get(x, y).x - x;
				v.y = segmentations[0].ref_shift_vector().get(x, y).y - y;
				fwrite(&v.x, sizeof(double), 1, fp);
				fwrite(&v.y, sizeof(double), 1, fp);
			}
		}
		fclose(fp);
		// Arbitrary shaped Block Matching
		printf("* * Compute Block Matching\n");
		//block_matching.reset(segmentations.begin()->ref_segmentation_map(), It, Itp1);
		if (sequence_Lab.size() <= 2) {
			block_matching.reset(
			    It_Lab_normalize, segmentations[1].ref_segmentation_map(),
			    Itp1_Lab_normalize, segmentations[0].ref_segmentation_map(),
			    Subpixel_Scale);
		} else {
			block_matching.reset(
			    sequence_Lab[2], segmentations[2].ref_segmentation_map(),
			    sequence_Lab[1], segmentations[1].ref_segmentation_map(),
			    sequence_Lab[0], segmentations[0].ref_segmentation_map(),
			    Subpixel_Scale);
		}
		block_matching.block_matching(BM_Search_Range, coeff_MAD, coeff_ZNCC);
		if (MaxLevel > 0) {
			MaxLevel = 0;
		}
	}
#endif

	// ----- Optical Flow -----
	if (MaxLevel >= 0) {
		try {
			u_levels = new ImgVector<VECTOR_2D<double> >[MaxLevel + 1];
		}
		catch (const std::bad_alloc &bad) {
			std::cerr << "error : ImgVector<VECTOR_2D<double> >* OpticalFlow_BlockMatching(const ImgVector<ImgClass::RGB>&, const ImgVector<ImgClass::RGB>&, double, const MULTIPLE_MOTION_PARAM&, const std::string, int)" << std::endl
			    << bad.what() << std::endl;
			throw;
		}
		// Make Pyramid
		try {
			if (sequence_Grayscale.size() >= 3) {
				It_levels = Pyramider(&(sequence_Grayscale[2]), MaxLevel);
				Itp1_levels = Pyramider(&(sequence_Grayscale[1]), MaxLevel);
			} else {
				It_levels = Pyramider(&It_normalize, MaxLevel);
				Itp1_levels = Pyramider(&Itp1_normalize, MaxLevel);
			}
			// The order reversed along with Block Matching (ordinary (It, Itp1))
			// Derivative about time
			I_dt_levels = dt_Pyramid(Itp1_levels, It_levels, MaxLevel);
			// Derivative about space
			grad_It_levels = grad_Pyramid(Itp1_levels, nullptr, MaxLevel);
		}
		catch (const std::bad_alloc& bad) {
			std::cerr << "error : ImgVector<VECTOR_2D<double> >* OpticalFlow_BlockMatching(const ImgVector<ImgClass::RGB>&, const ImgVector<ImgClass::RGB>&, double, const MULTIPLE_MOTION_PARAM&, const std::string, int)" << std::endl
			    << bad.what() << std::endl;
			delete[] grad_It_levels;
			delete[] I_dt_levels;
			delete[] Itp1_levels;
			delete[] It_levels;
			delete[] u_levels;
			throw;
		}
		// Initialize u_levels
		for (int level = 0; level <= MaxLevel; level++) {
			u_levels[level].reset(It_levels[level].width(), It_levels[level].height());
		}
		// Multi-Resolution IRLS Optical Flow estimation
		for (int level = MaxLevel; level >= 0; level--) {
			if (MaxLevel > 0) {
				sigmaD = sigmaD_init + (sigmaD_l0 - sigmaD_init) / MaxLevel * (MaxLevel - level);
				sigmaS = sigmaS_init + (sigmaS_l0 - sigmaS_init) / MaxLevel * (MaxLevel - level);
			} else {
				sigmaD = sigmaD_l0;
				sigmaS = sigmaS_l0;
			}
			printf("\nLevel %d : (1 / %d scaled, %dx%d)\n  sigmaD = %f\n  sigmaS = %f\n", level, int(pow_int(2.0, level)), u_levels[level].width(), u_levels[level].height(), sigmaD, sigmaS);
			if (level >= MaxLevel) {
				// The order of It_levels and Itp1_levels are reversed (ordinary It -> Itp1)
				LevelDown(I_dt_levels, u_levels, Itp1_levels, It_levels, level, MaxLevel, &block_matching);
			} else {
				// The order of It_levels and Itp1_levels are reversed (ordinary It -> Itp1)
				LevelDown(I_dt_levels, u_levels, Itp1_levels, It_levels, level, MaxLevel);
			}
			IterMax_level = 4 * MAX(It.width(), It.height());
			if (IterMax < 0 && IterMax_level >= IterMax) {
				IterMax_level = IterMax;
			}
			printf("IterMax = %d\n", IterMax_level);
			IRLS_OpticalFlow_Pyramid_Segment(
			    (u_levels + level),
			    segmentations.begin()->ref_segmentation_map(),
			    (grad_It_levels + level),
			    (I_dt_levels + level),
			    lambdaD, lambdaS, sigmaD, sigmaS,
			    IterMax_level,
			    MotionParam.Error_Min_Threshold,
			    level);
			Add_VectorOffset(u_levels, level, MaxLevel, &block_matching);
		}
		delete[] grad_It_levels;
		delete[] I_dt_levels;
		delete[] Itp1_levels;
		delete[] It_levels;
	}
	// Copy the lowest vector for output
	if (u_levels != nullptr) {
		u.resize(1);
		u[0].reset(It.width(), It.height());
		for (size_t i = 0; i < u[0].size(); i++) {
			u[0][i].x = u_levels[0][i].x;
			u[0][i].y = u_levels[0][i].y;
		}
	} else if (Bidirectional_with_Time) {
		u.resize(1);
		u[0].reset(It.width(), It.height());
		for (int y = 0; y < block_matching.height(); y++) {
			for (int x = 0; x < block_matching.width(); x++) {
				u[0].at(x, y) = block_matching.get(x, y);
			}
		}
	} else {
		u.resize(2);
		u[0].reset(It.width(), It.height());
		u[1].reset(It.width(), It.height());
		for (int y = 0; y < block_matching.height(); y++) {
			for (int x = 0; x < block_matching.width(); x++) {
				u[0].at(x, y) = block_matching.get_prev(x, y);
				u[0].at(x, y).t = -1;
				if (u.size() > 1) { // bi-directional
					u[1].at(x, y) = block_matching.get_next(x, y);
					u[1].at(x, y).t = 1;
				}
			}
		}
	}
	delete[] u_levels;
	return u;
}


template <class T>
void
Add_VectorOffset(ImgVector<VECTOR_2D<double> >* u_levels, int level, int MaxLevel, BlockMatching<T>* block_matching)
{
	if (level == MaxLevel) {
		// Add offset calculated by using the motion vector by Block Matching
		double Scale = double(u_levels[0].width()) / u_levels[level].width();

		for (int y = 0; y < u_levels[level].height(); y++) {
			for (int x = 0; x < u_levels[level].width(); x++) {
				u_levels[level].at(x, y).x +=
				    block_matching->get_prev(int(round(x * Scale)), int(round(y * Scale))).x
				    / Scale;
				u_levels[level].at(x, y).y +=
				    block_matching->get_prev(int(round(x * Scale)), int(round(y * Scale))).y
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

template <class T>
void
LevelDown(ImgVector<double> *I_dt_levels, ImgVector<VECTOR_2D<double> > *u_levels, const ImgVector<double> *It_levels, const ImgVector<double> *Itp1_levels, int level, int MaxLevel, BlockMatching<T>* block_matching)
{
	if (level == MaxLevel) {
		// Do NOT need projection from higher level
		return;
	}
	double Scale;
	if (block_matching == nullptr) {
		Scale = double(It_levels[0].width()) / double(It_levels[level].width());
	} else {
		Scale = 0.5;
	}

	for (int y = 0; y < u_levels[level].height(); y++) {
		for (int x = 0; x < u_levels[level].width(); x++) {
			VECTOR_2D<double> u_offset;
			if (block_matching == nullptr) {
				u_offset = u_levels[level + 1].get(int(round(x * Scale)), int(round(y * Scale)));
			} else {
				u_offset = block_matching->get_prev(int(round(x * Scale)), int(round(y * Scale)));
			}
			u_offset /= Scale;

			I_dt_levels[level].at(x, y) =
			    (Itp1_levels[level].get_zeropad(x + int(floor(u_offset.x)), y + int(floor(u_offset.y)))
			    - It_levels[level].get_zeropad(x, y)
			    + Itp1_levels[level].get_zeropad(x + int(floor(u_offset.x)) + 1, y + int(floor(u_offset.y)))
			    - It_levels[level].get_zeropad(x + 1, y)
			    + Itp1_levels[level].get_zeropad(x + int(floor(u_offset.x)), y + int(floor(u_offset.y)) + 1)
			    - It_levels[level].get_zeropad(x, y + 1)
			    + Itp1_levels[level].get_zeropad(x + int(floor(u_offset.x)) + 1, y + int(floor(u_offset.y)) + 1)
			    - It_levels[level].get_zeropad(x + 1, y + 1)) / 4.0;
			u_levels[level].at(x, y).x = 0.0;
			u_levels[level].at(x, y).y = 0.0;
		}
	}
}



void
IRLS_OpticalFlow_Pyramid_Segment(ImgVector<VECTOR_2D<double> >* u, const ImgVector<size_t>& domain_map, const ImgVector<VECTOR_2D<double> >* Img_g, const ImgVector<double>* Img_t, double lambdaD, double lambdaS, double sigmaD, double sigmaS, int IterMax, double ErrorMinThreshold, int level)
{
	ERROR Error("IRLS_OpticalFlow_Pyramid_Block");

	ImgVector<VECTOR_2D<double> > u_np1;
	VECTOR_2D<double> sup;
	double E = 0.0;
	double E_prev = 0.0;
	int ErrorIncrementCount = 0;

	if (u == nullptr) {
		throw std::invalid_argument("ImgVector<VECTOR_2D<double> > *u");
	} else if (Img_g == nullptr) {
		throw std::invalid_argument("ImgVector<VECTOR_2D<double> > *Img_g");
	} else if (Img_t == nullptr) {
		throw std::invalid_argument("ImgVector<double> *Img_t");
	}
	u_np1.copy(*u); // Initialize u_np1
	// Reset sup_Error_uu max Img_g
	sup_Error_uu_Block(Img_g, lambdaD, lambdaS, sigmaD, sigmaS);
	sup = sup_Error_uu_Block(nullptr, lambdaD, lambdaS, sigmaD, sigmaS);
	for (int n = 0; n < IterMax; n++) {
		// Calc for all sites
		size_t site;
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (site = 0; site < u->size(); site++) {
			VECTOR_2D<double> dE;
			dE = Error_u_Block(site, u, domain_map, Img_g, Img_t, lambdaD, lambdaS, sigmaD, sigmaS);
			u_np1[site].x = u->get(site).x - dE.x / sup.x;
			u_np1[site].y = u->get(site).y - dE.y / sup.y;
		}
		// Calc for all sites
		for (site = 0; site < u->size(); site++) {
			u->at(site) = u_np1.get(site);
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
Error_u_Block(const size_t& site, const ImgVector<VECTOR_2D<double> >* u, const ImgVector<size_t>& domain_map, const ImgVector<VECTOR_2D<double> >* Img_g, const ImgVector<double>* Img_t, const double& lambdaD, const double& lambdaS, const double& sigmaD, const double& sigmaS)
{
	double (*psiD)(const double&, const double&) = Geman_McClure_psi;
	double (*psiS)(const double&, const double&) = Geman_McClure_psi;
	VECTOR_2D<double> us;
	double Center;
	VECTOR_2D<double> Neighbor;
	VECTOR_2D<double> E_u;

	int x = static_cast<int>(site % size_t(u->width()));
	int y = static_cast<int>(site / size_t(u->width()));
	size_t center_domain = domain_map.get(x, y);

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
sup_Error_uu_Block(const ImgVector<VECTOR_2D<double> >* Img_g, const double& lambdaD, const double& lambdaS, const double& sigmaD, const double& sigmaS)
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
Error_MultipleMotion_Block(const ImgVector<VECTOR_2D<double> >* u, const ImgVector<size_t>& domain_map, const ImgVector<VECTOR_2D<double> >* Img_g, const ImgVector<double>* Img_t, const double& lambdaD, const double& lambdaS, const double& sigmaD, const double& sigmaS)
{
	double (*rhoD)(const double&, const double&) = Geman_McClure_rho;
	double (*rhoS)(const double&, const double&) = Geman_McClure_rho;
	double E = 0.0;
	int y;

#ifdef _OPENMP
#pragma omp parallel for reduction(+:E)
#endif
	for (y = 0; y < u->height(); y++) {
		for (int x = 0; x < u->width(); x++) {
			size_t center_domain = domain_map.get(x, y);
			VECTOR_2D<double> us(u->get(x, y));
			VECTOR_2D<double> Neighbor(0.0, 0.0);
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
			double Center = (*rhoD)(Img_g->get(x, y).x * us.x
			    + Img_g->get(x, y).y * us.y
			    + Img_t->get(x, y),
			    sigmaD);
			E += lambdaD * Center + lambdaS * (Neighbor.x + Neighbor.y);
		}
	}
	return E;
}




void
MultipleMotion_write(const ImgVector<double>& img_prev, const ImgVector<double>& img_current, const std::vector<ImgVector<Vector_ST<double> > >& u, const std::string& filename)
{
	ERROR Error("MultipleMotion_write");
	FILE *fp = nullptr;
	MotionCompensation<double> compensated(img_prev, img_current, u[0]);
	PNM pnm;
	std::string filename_compensated;

	printf("\n* Output The Optical Flow to '%s'(binary)\n", filename.c_str());
	if ((fp = fopen(filename.c_str(), "wb")) == nullptr) {
		Error.Function("fopen");
		Error.File(filename.c_str());
		Error.FileWrite();
		throw std::logic_error("fopen");
	}
	fprintf(fp, "%d %d\n", u[0].width(), u[0].height());
	for (int y = 0; y < u[0].height(); y++) {
		for (int x = 0; x < u[0].width(); x++) {
			Vector_ST<double> v = u[0].get(x, y);
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
MultipleMotion_write(const ImgVector<ImgClass::RGB>& img_prev, const ImgVector<ImgClass::RGB>& img_current, const std::vector<ImgVector<Vector_ST<double> > >& u, const std::string &filename)
{
	ERROR Error("MultipleMotion_write");
	FILE *fp = nullptr;
	MotionCompensation<ImgClass::RGB> compensated(img_prev, img_current, u[0]);
	PNM pnm;
	std::string filename_compensated;

	printf("\n* Output The Optical Flow to '%s'(binary)\n", filename.c_str());
	if ((fp = fopen(filename.c_str(), "wb")) == nullptr) {
		Error.Function("fopen");
		Error.File(filename.c_str());
		Error.FileWrite();
		throw std::logic_error("fopen");
	}
	fprintf(fp, "%d %d\n", u[0].width(), u[0].height());
	for (int y = 0; y < u[0].height(); y++) {
		for (int x = 0; x < u[0].width(); x++) {
			Vector_ST<double> v = u[0].get(x, y);
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

void
MultipleMotion_write(const ImgVector<double>& img_prev, const ImgVector<double>& img_current, const ImgVector<double>& img_next, const std::vector<ImgVector<Vector_ST<double> > >& u, const std::string& filename)
{
	ERROR Error("MultipleMotion_write");

	FILE *fp = nullptr;
	int x, y;
	MotionCompensation<double> compensated;
	PNM pnm;
	std::string filename_compensated;

	if (u.size() == 1) {
		compensated.set(img_prev, img_current, img_next, u[0]);
	} else {
		compensated.set(img_prev, img_current, img_next, u);
	}
	printf("\n* Output The Optical Flow to '%s'(binary)\n", filename.c_str());
	if ((fp = fopen(filename.c_str(), "wb")) == nullptr) {
		Error.Function("fopen");
		Error.File(filename.c_str());
		Error.FileWrite();
		throw std::logic_error("fopen");
	}
	fprintf(fp, "%d %d\n", u[0].width(), u[0].height());
	for (y = 0; y < u[0].height(); y++) {
		for (x = 0; x < u[0].width(); x++) {
			Vector_ST<double> v = u[0].get(x, y);
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
MultipleMotion_write(const ImgVector<ImgClass::RGB>& img_prev, const ImgVector<ImgClass::RGB>& img_current, const ImgVector<ImgClass::RGB>& img_next, const std::vector<ImgVector<Vector_ST<double> > >& u, const std::string &filename)
{
	ERROR Error("MultipleMotion_write");

	FILE *fp = nullptr;
	VECTOR_2D<double> v;
	int x, y;
	MotionCompensation<ImgClass::RGB> compensated;
	PNM pnm;
	std::string filename_compensated;

	if (u.size() == 1) {
		compensated.set(img_prev, img_current, img_next, u[0]);
	} else {
		compensated.set(img_prev, img_current, img_next, u);
	}
	printf("\n* Output The Optical Flow to '%s'(binary)\n", filename.c_str());
	if ((fp = fopen(filename.c_str(), "wb")) == nullptr) {
		Error.Function("fopen");
		Error.File(filename.c_str());
		Error.FileWrite();
		throw std::logic_error("fopen");
	}
	fprintf(fp, "%d %d\n", u[0].width(), u[0].height());
	for (y = 0; y < u[0].height(); y++) {
		for (x = 0; x < u[0].width(); x++) {
			v = u[0].get(x, y);
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

