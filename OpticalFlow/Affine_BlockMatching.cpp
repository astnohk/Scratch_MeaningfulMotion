/* Robust Estimation of Multiple Motions
 *
 * This motion estimation based on
 * M.J.Black and P.Anandan, "The Robust Estimation of Multiple Motions: Parametric and Piecewise-Smooth Flow Fields," Computer Vision and Image Understanding, Vol.63, No.1, 1996, pp.75-104.
 */
#include "Affine_BlockMatching.h"




std::vector<ImgVector<Vector_ST<double> > >
AffineParametric(const ImgVector<ImgClass::RGB>& It_color, const ImgVector<ImgClass::RGB>& Itp1_color, double MaxInt, MULTIPLE_MOTION_PARAM MotionParam, const std::string ofilename, int IterMax)
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
	int Subpixel_Scale = 1;

	ImgVector<double> It;
	ImgVector<double> Itp1;
	ImgVector<ImgClass::RGB> It_sRGB_normalize;
	ImgVector<ImgClass::RGB> Itp1_sRGB_normalize;
	ImgVector<ImgClass::Lab> It_Lab_normalize;
	ImgVector<ImgClass::Lab> Itp1_Lab_normalize;
	ImgVector<double> It_normalize;
	ImgVector<double> Itp1_normalize;

	// M-estimator parameter
	const double sigma = 0.2 / sqrt(2.0); //4.0 / sqrt(2.0);

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
		double kernel_spatial = MEANSHIFT_KERNEL_SPATIAL, kernel_intensity = 9.0 / 255.0; // for images under about HD resolution
#else
		//double kernel_spatial = 64.0, kernel_intensity = 12.0 / 255.0; // for 4K Film kernel(spatial = 64.0, intensity = 12.0 / 255.0)
		double kernel_spatial = 10.0, kernel_intensity = 9.0 / 255.0; // for images under about HD resolution
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
		std::string ofilename_segmentation = ofilename.substr(0, found) + "segmentation_" + ofilename.substr(found);
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
		std::string ofilename_vector = ofilename.substr(0, found) + "shift-vector_" + ofilename.substr(found);
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
	}
#endif

	// ----- Optical Flow -----
	std::vector<std::vector<double> > u_affine;
	if (MotionParam.Level > 0) {
		const std::vector<std::vector<VECTOR_2D<int> > >* regions;
		ImgVector<ImgClass::Lab>* interest = nullptr;
		std::vector<ImgVector<ImgClass::Lab>*> references;

		if (sequence_Grayscale.size() <= 2) {
			regions = &(segmentations[0].ref_regions());
			interest = &Itp1_Lab_normalize;
			references.push_back(&It_Lab_normalize);
		} else {
			regions = &(segmentations[1].ref_regions());
			interest = &sequence_Lab[1];
			references.push_back(&sequence_Lab[2]);
			references.push_back(&sequence_Lab[0]);
		}
		// Gradient
		ImgVector<VECTOR_2D<double> > grad(It_color.width(), It_color.height());
		for (int y = 0; y < It_color.height(); y++) {
			for (int x = 0; x < It_color.width(); x++) {
				grad.at(x, y).x =
				    (interest->get_mirror(x + 1, y).L - interest->get(x, y).L
				    + interest->get_mirror(x + 1, y + 1).L - interest->get_mirror(x, y + 1).L)
				    / 2.0;
				grad.at(x, y).y =
				    (interest->get_mirror(x, y + 1).L - interest->get(x, y).L
				    + interest->get_mirror(x + 1, y + 1).L - interest->get_mirror(x + 1, y).L)
				    / 2.0;
			}
		}
		// Initialize affine coefficient vector
		u_affine.resize(regions->size());
		for (size_t n = 0; n < regions->size(); n++) {
			u_affine[n].resize(NUM_AFFINE_PARAMETER);
		}
		// Gradient-based estimation
		for (size_t ref = 0; ref < references.size(); ref++) {
			ImgVector<double> dt_I(interest->width(), interest->height());
			for (int y = 0; y < interest->height(); y++) {
				for (int x = 0; x < interest->width(); x++) {
					dt_I.at(x, y) =
					    (references[ref]->get(x, y).L - interest->get(x, y).L
					    + references[ref]->get_mirror(x + 1, y).L - interest->get_mirror(x + 1, y).L
					    + references[ref]->get_mirror(x, y + 1).L - interest->get_mirror(x, y + 1).L
					    + references[ref]->get_mirror(x + 1, y + 1).L - interest->get_mirror(x + 1, y + 1).L)
					    / 4.0;
				}
			}
			// IRLS Optical Flow estimation
			printf("\n* IRLS\n    sigma = %f, iteration = %d\n", sigma, IterMax);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
			for (size_t n = 0; n < regions->size(); n++) {
				IRLS_AffineParametric_region(
				    &u_affine[n],
				    regions->at(n),
				    &grad,
				    &dt_I,
				    sigma,
				    IterMax,
				    MotionParam.Error_Min_Threshold);
			}
		}
	}

	// Copy the lowest vector for output
	if (Bidirectional_with_Time) {
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
	return u;
}




void
IRLS_AffineParametric_region(std::vector<double>* u_affine, const std::vector<VECTOR_2D<int> >& region, const ImgVector<VECTOR_2D<double> >* grad, const ImgVector<double>* dt, const double& sigma, const int IterMax, const double& ErrorMinThreshold)
{
	const double omega = 1.0E-3;

	// Initialize
	for (size_t i = 0; i < NUM_AFFINE_PARAMETER; i++) {
		u_affine->at(i) = .0;
	}
	// Start IRLS
	for (int n = 0; n < IterMax; n++) {
		std::vector<double> u_np1(NUM_AFFINE_PARAMETER);
		std::vector<double> dE = Error_a_region(u_affine, region, grad, dt, sigma);
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
		E = Error_Affine_region(u_affine, region, grad, dt, sigma);
		if (E < ErrorMinThreshold) {
			break;
		}
	}
}


std::vector<double>
Error_a_region(const std::vector<double>* u_affine, const std::vector<VECTOR_2D<int> >& region, const ImgVector<VECTOR_2D<double> >* grad, const ImgVector<double>* dt, double sigma)
{
	double (*psiD)(const double&, const double&) = Geman_McClure_psi;
	std::vector<double> E_a(NUM_AFFINE_PARAMETER);
	VECTOR_2D<double> u_a;

	for (size_t i = 0; i < NUM_AFFINE_PARAMETER; i++) {
		E_a[i] = .0;
	}
	for (const VECTOR_2D<int>& element : region) {
		int x = element.x;
		int y = element.y;
		u_a.x = u_affine->at(0) + u_affine->at(1) * x + u_affine->at(2) * y;
		u_a.y = u_affine->at(3) + u_affine->at(4) * x + u_affine->at(5) * y;
		E_a[0] += grad->get(x, y).x * (*psiD)(grad->get(x, y).x * u_a.x + grad->get(x, y).y * u_a.y + dt->get_zeropad(x, y), sigma);
		E_a[1] += grad->get(x, y).x * x * (*psiD)(grad->get(x, y).x * u_a.x + grad->get(x, y).y * u_a.y + dt->get_zeropad(x, y), sigma);
		E_a[2] += grad->get(x, y).x * y * (*psiD)(grad->get(x, y).x * u_a.x + grad->get(x, y).y * u_a.y + dt->get_zeropad(x, y), sigma);
		E_a[3] += grad->get(x, y).y * (*psiD)(grad->get(x, y).x * u_a.x + grad->get(x, y).y * u_a.y + dt->get_zeropad(x, y), sigma);
		E_a[4] += grad->get(x, y).y * x * (*psiD)(grad->get(x, y).x * u_a.x + grad->get(x, y).y * u_a.y + dt->get_zeropad(x, y), sigma);
		E_a[5] += grad->get(x, y).y * y * (*psiD)(grad->get(x, y).x * u_a.x + grad->get(x, y).y * u_a.y + dt->get_zeropad(x, y), sigma);
	}
	return E_a;
}


std::vector<double>
sup_Error_aa_region(const std::vector<VECTOR_2D<int> >& region, const ImgVector<VECTOR_2D<double> >* grad, double sigma)
{
	ERROR Error("sup_Error_aa_region");

	std::vector<double> sup(NUM_AFFINE_PARAMETER);
	std::vector<double> u_aa_max(NUM_AFFINE_PARAMETER);

	for (size_t i = 0; i < NUM_AFFINE_PARAMETER; i++) {
		u_aa_max[i] = .0;
	}
	if (grad == nullptr) {
		Error.Value("grad");
		Error.PointerNull();
		return u_aa_max;
	}
	for (const VECTOR_2D<int>& element : region) {
		int x = element.x;
		int y = element.y;
		// u = a0 + a1 * x + a2 * y
		if (u_aa_max[0] < POW2(grad->get(x, y).x)) {
			u_aa_max[0] = POW2(grad->get(x, y).x);
		}
		if (u_aa_max[1] < POW2(grad->get(x, y).x * x)) {
			u_aa_max[1] = POW2(grad->get(x, y).x * x);
		}
		if (u_aa_max[2] < POW2(grad->get(x, y).x * y)) {
			u_aa_max[2] = POW2(grad->get(x, y).x * y);
		}
		// v = a3 + a4 * x + a5 * y
		if (u_aa_max[3] < POW2(grad->get(x, y).y)) {
			u_aa_max[3] = POW2(grad->get(x, y).y);
		}
		if (u_aa_max[4] < POW2(grad->get(x, y).y * x)) {
			u_aa_max[4] = POW2(grad->get(x, y).y * x);
		}
		if (u_aa_max[5] < POW2(grad->get(x, y).y * y)) {
			u_aa_max[5] = POW2(grad->get(x, y).y * y);
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
Error_Affine_region(const std::vector<double>* u_affine, const std::vector<VECTOR_2D<int> >& region, const ImgVector<VECTOR_2D<double> >* grad, const ImgVector<double>* dt, double sigma)
{
	double (*rhoD)(const double&, const double&) = Geman_McClure_rho;
	double E = 0.0;
	VECTOR_2D<double> u_a;

	for (const VECTOR_2D<int>& element : region) {
		int x = element.x;
		int y = element.y;
		u_a.x = u_affine->at(0) + u_affine->at(1) * x + u_affine->at(2) * y;
		u_a.y = u_affine->at(3) + u_affine->at(4) * x + u_affine->at(5) * y;
		E += (*rhoD)(grad->get(x, y).x * u_a.x + grad->get(x, y).y * u_a.y + dt->get(x, y), sigma);
	}
	return E;
}

