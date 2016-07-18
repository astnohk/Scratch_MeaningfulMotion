/* Robust Estimation of Multiple Motions
 *
 * This motion estimation based on
 * M.J.Black and P.Anandan, "The Robust Estimation of Multiple Motions: Parametric and Piecewise-Smooth Flow Fields," Computer Vision and Image Understanding, Vol.63, No.1, 1996, pp.75-104.
 */
#include "OpticalFlow_BlockMatching.h"

using namespace ImgClass;



// This function will compute INVERSE Optical Flow it points the previous frame which will come to the current (next) frame.
std::vector<ImgVector<Vector_ST<double> > >
OpticalFlow_BlockMatching(const ImgVector<RGB>& It_color, const ImgVector<RGB>& Itp1_color, double MaxInt, MULTIPLE_MOTION_PARAM MotionParam, const std::string newest_filename, const int Mode, const int IterMax)
{
	const bool Bidirectional = true;
	const bool Bidirectional_with_Time = true; // on almost all cases it is true
	const size_t History_Max = 4;

	static std::deque<ImgVector<RGB> > sequence_sRGB;
	static std::deque<ImgVector<Lab> > sequence_Lab;
	static std::deque<Segmentation<Lab> > segmentations;

	std::bad_alloc except_bad_alloc;

	std::vector<ImgVector<Vector_ST<double> > > u; // For RETURN value

	ImgVector<size_t> domain_map;
	const double coeff_MAD = 1.0;
	const double coeff_ZNCC = 0.5;
	BlockMatching<Lab> block_matching;
	int BM_Search_Range = 61; // Block Matching search range
	int Subpixel_Scale = 2;

	ImgVector<double> It;
	ImgVector<double> Itp1;
	ImgVector<RGB> It_sRGB_normalize;
	ImgVector<RGB> Itp1_sRGB_normalize;
	ImgVector<Lab> It_Lab_normalize;
	ImgVector<Lab> Itp1_Lab_normalize;
	ImgVector<double> It_normalize;
	ImgVector<double> Itp1_normalize;
	ImgVector<RGB> It_sRGB_quantized;
	// M-estimator parameter
	const double lambdaD = 5.0;
	const double lambdaS = 1.0;
	const double sigmaD = 0.2 / sqrt(2.0); //4.0 / sqrt(2.0);
	const double sigmaS = 0.03 / sqrt(2.0);

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
		sequence_Lab.push_front(It_Lab_normalize);
	}
	sequence_sRGB.push_front(Itp1_sRGB_normalize);
	sequence_Lab.push_front(Itp1_Lab_normalize);
	if (sequence_sRGB.size() >= History_Max) {
		sequence_sRGB.pop_back();
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
	if (sequence_Lab.size() >= 3) {
		if (Bidirectional) {
			block_matching.reset(sequence_Lab[2], sequence_Lab[1], sequence_Lab[0], BlockMatching_BlockSize, Subpixel_Scale);
		} else {
			block_matching.reset(sequence_Lab[2], sequence_Lab[1], BlockMatching_BlockSize, Subpixel_Scale);
		}
	} else {
		block_matching.reset(It_Lab_normalize, Itp1_Lab_normalize, BlockMatching_BlockSize, Subpixel_Scale);
	}
	block_matching.block_matching(BM_Search_Range, coeff_MAD, coeff_ZNCC);
#else
	{
		// Segmentation
		printf("* * Compute Segmentation by Mean Shift\n");

#ifdef MEANSHIFT_KERNEL_SPATIAL
		double kernel_spatial = MEANSHIFT_KERNEL_SPATIAL, kernel_intensity = 16.0 / 255.0; // for images under about HD resolution
#else
		//double kernel_spatial = 64.0, kernel_intensity = 16.0 / 255.0; // for 4K Film kernel(spatial = 64.0, intensity = 12.0 / 255.0)
		double kernel_spatial = 20.0, kernel_intensity = 16.0 / 255.0; // for images under about HD resolution
#endif

		if (segmentations.empty()) {
			segmentations.push_front(Segmentation<Lab>(It_Lab_normalize, kernel_spatial, kernel_intensity));
		}
		segmentations.push_front(Segmentation<Lab>(Itp1_Lab_normalize, kernel_spatial, kernel_intensity));
		if (segmentations.size() >= History_Max) {
			segmentations.pop_back();
		}

		PNM pnm;
		std::string::size_type found = 1 + newest_filename.find_last_not_of("0123456789", newest_filename.find_last_of("0123456789"));
		if (found == std::string::npos) {
			found = newest_filename.find_last_of(".");
		}
		std::string newest_filename_segmentation = newest_filename.substr(0, found) + "segmentation_" + newest_filename.substr(found);
		printf("* Output The Segmentation result to '%s'(binary)\n\n", newest_filename_segmentation.c_str());
		{
			ImgVector<int> tmp_vector(segmentations[0].width(), segmentations[0].height());
			for (size_t i = 0; i < segmentations[0].size(); i++) {
				tmp_vector[i] = static_cast<int>(segmentations[0][i]);
			}
			pnm.copy(PORTABLE_GRAYMAP_BINARY, segmentations[0].width(), segmentations[0].height(), int(tmp_vector.max()), tmp_vector.data());
			pnm.write(newest_filename_segmentation.c_str());
			pnm.free();
		}

		{
			int width = sequence_sRGB[0].width();
			int height = sequence_sRGB[0].height();

			int *quantized = new int[3 * segmentations[0].width() * segmentations[0].height()];
			for (const std::vector<VECTOR_2D<int> >& region : segmentations[0].ref_regions()) {
				RGB sum_sRGB(.0, .0, .0);
				for (const VECTOR_2D<int>& r : region) {
					sum_sRGB += sequence_sRGB[0].get(r.x, r.y);
				}
				sum_sRGB *= 255.0 / region.size();
				sum_sRGB.R = sum_sRGB.R > 255.0 ? 255 : sum_sRGB.R;
				sum_sRGB.G = sum_sRGB.G > 255.0 ? 255 : sum_sRGB.G;
				sum_sRGB.B = sum_sRGB.B > 255.0 ? 255 : sum_sRGB.B;
				for (const VECTOR_2D<int>& r : region) {
					quantized[width * r.y + r.x] = int(sum_sRGB.R);
					quantized[width * height + width * r.y + r.x] = int(sum_sRGB.G);
					quantized[2 * width * height + width * r.y + r.x] = int(sum_sRGB.B);
				}
			}
			std::string newest_filename_quantized = newest_filename.substr(0, found) + "color-quantized_" + newest_filename.substr(found);
			printf("* Output The color quantized image '%s'(binary)\n\n", newest_filename_quantized.c_str());
			pnm.copy(PORTABLE_PIXMAP_BINARY, segmentations[0].width(), segmentations[0].height(), 255, quantized);
			delete[] quantized;
			quantized = nullptr;
			pnm.write(newest_filename_quantized.c_str());
			pnm.free();
		}
		// Output vectors
		std::string newest_filename_vector = newest_filename.substr(0, found) + "shift-vector_" + newest_filename.substr(found);
		FILE *fp;
		fp = fopen(newest_filename_vector.c_str(), "w");
		fprintf(fp, "%d %d\n", segmentations[0].width(), segmentations[0].height());
		for (int y = 0; y < segmentations[0].height(); y++) {
			for (int x = 0; x < segmentations[0].width(); x++) {
				VECTOR_2D<double> v;
				v.x = segmentations[0].ref_shift_vector_spatial().get(x, y).x - x;
				v.y = segmentations[0].ref_shift_vector_spatial().get(x, y).y - y;
				fwrite(&v.x, sizeof(double), 1, fp);
				fwrite(&v.y, sizeof(double), 1, fp);
			}
		}
		fclose(fp);
		// Arbitrary shaped Block Matching
		printf("* * Compute Block Matching\n");
		//block_matching.reset(segmentations.begin()->ref_segmentation_map(), It, Itp1);
		if (sequence_Lab.size() >= 3) {
			if (Bidirectional) {
				block_matching.reset(
				    sequence_Lab[2], segmentations[2].ref_segmentation_map(),
				    sequence_Lab[1], segmentations[1].ref_segmentation_map(),
				    sequence_Lab[0], segmentations[0].ref_segmentation_map(),
				    Subpixel_Scale);
			} else {
				block_matching.reset(
				    sequence_Lab[2], segmentations[2].ref_segmentation_map(),
				    sequence_Lab[1], segmentations[1].ref_segmentation_map(),
				    Subpixel_Scale);
			}
		} else {
			block_matching.reset(
			    It_Lab_normalize, segmentations[1].ref_segmentation_map(),
			    Itp1_Lab_normalize, segmentations[0].ref_segmentation_map(),
			    Subpixel_Scale);
		}
		block_matching.block_matching(BM_Search_Range, coeff_MAD, coeff_ZNCC);
	}
#endif

	// ----- Optical Flow -----
	std::vector<ImgVector<VECTOR_2D<double> > > u_optical;
	if (Bidirectional && sequence_Lab.size() >= 3) {
		u_optical.resize(2);
		u_optical[0].resize(It_color.width(), It_color.height());
		u_optical[1].resize(It_color.width(), It_color.height());
	} else {
		u_optical.resize(1);
		u_optical[0].resize(It_color.width(), It_color.height());
	}
	if (MotionParam.Level >= 0) {
		const ImgVector<Lab>* interest = nullptr;
		const std::vector<std::vector<VECTOR_2D<int> > >* regions = nullptr;
		const ImgVector<size_t>* region_map = nullptr;
		std::vector<const ImgVector<Lab>*> references;
		std::vector<ImgVector<VECTOR_2D<double> > > MVs;
		if (Bidirectional && sequence_Lab.size() >= 3) {
			interest = &sequence_Lab[1];
			regions = &(segmentations[1].ref_regions());
			region_map = &(segmentations[1].ref_segmentation_map());
			references.push_back(&sequence_Lab[2]);
			references.push_back(&sequence_Lab[0]);
			if (block_matching.vector_field_width() == It_color.width()
			    && block_matching.vector_field_height() == It_color.height()) {
				MVs.push_back(block_matching.ref_motion_vector_prev());
				MVs.push_back(block_matching.ref_motion_vector_next());
			} else {
				MVs.resize(2);
				MVs[0].reset(It_color.width(), It_color.height());
				MVs[1].reset(It_color.width(), It_color.height());
				for (int y = 0; y < It_color.height(); y++) {
					for (int x = 0; x < It_color.width(); x++) {
						MVs[0].at(x, y) = block_matching.get_prev(x, y);
						MVs[1].at(x, y) = block_matching.get_next(x, y);
					}
				}
			}
		} else {
			interest = &Itp1_Lab_normalize;
			regions = &(segmentations[0].ref_regions());
			region_map = &(segmentations[0].ref_segmentation_map());
			references.push_back(&It_Lab_normalize);
			if (block_matching.vector_field_width() == It_color.width()
			    && block_matching.vector_field_height() == It_color.height()) {
				MVs.push_back(block_matching.ref_motion_vector_prev());
			} else {
				MVs.resize(1);
				MVs[0].reset(It_color.width(), It_color.height());
				for (int y = 0; y < It_color.height(); y++) {
					for (int x = 0; x < It_color.width(); x++) {
						MVs[0].at(x, y) = block_matching.get_prev(x, y);
					}
				}
			}
		}
		switch (Mode) {
			case MODE_OUTPUT_AFFINE_BLOCKMATCHING: // Affine Parametric motion
				for (size_t ref = 0; ref < references.size(); ref++) {
					u_optical[ref] = AffineParametric(
					    *references[ref],
					    *interest,
					    MVs[ref],
					    *regions,
					    MotionParam,
					    IterMax);
				}
				break;
			default: // A Gradient-based method
				for (size_t ref = 0; ref < references.size(); ref++) {
					for (unsigned int i = 0; i < MVs[ref].size(); i++) { // for DEBUG
						MVs[ref][i] = 0.0;
					}
					u_optical[ref] = OpticalFlow_GradientMethod(
					    references[ref],
					    interest,
					    &MVs[ref],
					    region_map,
					    lambdaD, lambdaS, sigmaD, sigmaS,
					    IterMax,
					    MotionParam.Error_Min_Threshold);
				}
		}
	}
	// Copy the lowest vector for output
	if (Bidirectional) {
		if (Bidirectional_with_Time) {
			u.resize(1);
			u[0].reset(It.width(), It.height());
			for (int y = 0; y < It.height(); y++) {
				for (int x = 0; x < It.width(); x++) {
					u[0].at(x, y) = block_matching.get(x, y);
					if (u_optical.size() > 0) {
						if (u[0].get(x, y).t < 0) {
							u[0].at(x, y).x += u_optical[0].get(x, y).x;
							u[0].at(x, y).y += u_optical[0].get(x, y).y;
						} else {
							u[0].at(x, y).x += u_optical[1].get(x, y).x;
							u[0].at(x, y).y += u_optical[1].get(x, y).y;
						}
					}
				}
			}
		} else {
			u.resize(2);
			u[0].reset(It.width(), It.height());
			u[1].reset(It.width(), It.height());
			for (int y = 0; y < It.height(); y++) {
				for (int x = 0; x < It.width(); x++) {
					u[0].at(x, y) = block_matching.get_prev(x, y);
					u[0].at(x, y).t = -1;
					if (u_optical.size() > 0) {
						u[0].at(x, y).x += u_optical[0].get(x, y).x;
						u[0].at(x, y).y += u_optical[0].get(x, y).y;
					}
					if (u.size() > 1) { // bi-directional
						u[1].at(x, y) = block_matching.get_next(x, y);
						u[1].at(x, y).t = 1;
						if (u_optical.size() > 1) {
							u[1].at(x, y).x += u_optical[1].get(x, y).x;
							u[1].at(x, y).y += u_optical[1].get(x, y).y;
						}
					}
				}
			}
		}
	} else {
		u.resize(1);
		u[0].reset(It.width(), It.height());
		for (int y = 0; y < It.height(); y++) {
			for (int x = 0; x < It.width(); x++) {
				u[0].at(x, y) = block_matching.get(x, y);
				if (u_optical.size() > 0) {
					u[0].at(x, y).x += u_optical[0].get(x, y).x;
					u[0].at(x, y).y += u_optical[0].get(x, y).y;
				}
			}
		}
	}
	return u;
}




ImgVector<VECTOR_2D<double> >
OpticalFlow_GradientMethod(const ImgVector<Lab>* reference, const ImgVector<Lab>* interest, const ImgVector<VECTOR_2D<double> >* MV, const ImgVector<size_t>* region_map, const double& lambdaD, const double& lambdaS, const double& sigmaD, const double& sigmaS, const int IterMax, const double& Error_Min_Threshold)
{
	ImgVector<VECTOR_2D<double> > u(reference->width(), reference->height());
	// Compute gradient and derivation
	ImgVector<VECTOR_2D<double> > grad(reference->width(), reference->height());
	for (int y = 0; y < reference->height(); y++) {
		for (int x = 0; x < reference->width(); x++) {
			grad.at(x, y).x =
			    (interest->get_mirror(x + 1, y).L - interest->get_mirror(x, y).L
			    + interest->get_mirror(x + 1, y + 1).L - interest->get_mirror(x, y + 1).L)
			    / 2.0;
			grad.at(x, y).y =
			    (interest->get_mirror(x, y + 1).L - interest->get_mirror(x, y).L
			    + interest->get_mirror(x + 1, y + 1).L - interest->get_mirror(x + 1, y).L)
			    / 2.0;
		}
	}
	ImgVector<double> dt(reference->width(), reference->height());
	for (int y = 0; y < reference->height(); y++) {
		for (int x = 0; x < reference->width(); x++) {
			int x_t = x + int(floor(MV->get(x, y).x));
			int y_t = y + int(floor(MV->get(x, y).y));
			dt.at(x, y) =
			    (reference->get_mirror(x_t, y_t).L - interest->get_mirror(x, y).L
			    + reference->get_mirror(x_t + 1, y_t).L - interest->get_mirror(x + 1, y).L
			    + reference->get_mirror(x_t, y_t + 1).L - interest->get_mirror(x, y + 1).L
			    + reference->get_mirror(x_t + 1, y_t + 1).L - interest->get_mirror(x + 1, y + 1).L)
			    / 4.0;
		}
	}
	// Optical Flow estimation with IRLS
	printf("\n    lambdaD = %f, lambdaS = %f\n    sigmaD = %f\n    sigmaS = %f\n    iteration = %d\n", lambdaD, lambdaS, sigmaD, sigmaS, IterMax);
	IRLS_OpticalFlow_GradientMethod(
	    &u,
	    region_map,
	    &grad,
	    &dt,
	    lambdaD, lambdaS, sigmaD, sigmaS,
	    IterMax,
	    Error_Min_Threshold);
	return u;
}


void
IRLS_OpticalFlow_GradientMethod(ImgVector<VECTOR_2D<double> >* u, const ImgVector<size_t>* region_map, const ImgVector<VECTOR_2D<double> >* grad, const ImgVector<double>* dt, const double& lambdaD, const double& lambdaS, const double& sigmaD, const double& sigmaS, const int IterMax, const double& Error_Min_Threshold)
{
	ERROR Error("IRLS_OpticalFlow_Pyramid_Block");

	ImgVector<VECTOR_2D<double> > u_np1;
	VECTOR_2D<double> sup;
	double E = 0.0;
	double E_prev = 0.0;
	int Error_IncrementCount = 0;

	if (u == nullptr) {
		throw std::invalid_argument("ImgVector<VECTOR_2D<double> >* u");
	}
	u_np1.copy(*u); // Initialize u_np1
	// Reset sup_Error_uu max Img_g
	sup_Error_uu_Block(grad, lambdaD, lambdaS, sigmaD, sigmaS);
	sup = sup_Error_uu_Block(nullptr, lambdaD, lambdaS, sigmaD, sigmaS);
	for (int n = 0; n < IterMax; n++) {
		// Calc for all sites
		size_t site;
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (site = 0; site < u->size(); site++) {
			VECTOR_2D<double> dE;
			dE = Error_u_Block(site, u, region_map, grad, dt, lambdaD, lambdaS, sigmaD, sigmaS);
			u_np1[site].x = u->get(site).x - dE.x / sup.x;
			u_np1[site].y = u->get(site).y - dE.y / sup.y;
		}
		// Calc for all sites
		for (site = 0; site < u->size(); site++) {
			u->at(site) = u_np1.get(site);
		}
		if ((n & 0x3F) == 0) {
			E_prev = E;
			E = Error_MultipleMotion_Block(u, region_map, grad, dt, lambdaD, lambdaS, sigmaD, sigmaS);
#ifdef SHOW_IRLS_OPTICALFLOW_PYRAMID_E
			printf("E(%4d) = %e\n", n, E);
#endif
			if (E > E_prev) {
				Error_IncrementCount++;
			} else {
				Error_IncrementCount = 0;
			}
			if (E < Error_Min_Threshold || Error_IncrementCount > 3) {
				break;
			}
		}
	}
}


VECTOR_2D<double>
Error_u_Block(const size_t& site, const ImgVector<VECTOR_2D<double> >* u, const ImgVector<size_t>* domain_map, const ImgVector<VECTOR_2D<double> >* Img_g, const ImgVector<double>* Img_t, const double& lambdaD, const double& lambdaS, const double& sigmaD, const double& sigmaS)
{
	double (*psiD)(const double&, const double&) = Geman_McClure_psi;
	double (*psiS)(const double&, const double&) = Geman_McClure_psi;
	VECTOR_2D<double> us;
	double Center;
	VECTOR_2D<double> Neighbor;
	VECTOR_2D<double> E_u;
	double inner_prod;
	double coeff;

	int x = static_cast<int>(site % size_t(u->width()));
	int y = static_cast<int>(site / size_t(u->width()));
	size_t center_domain = domain_map->get(x, y);

	us = u->get(site);
	Center = (*psiD)(Img_g->get(site).x * us.x + Img_g->get(site).y * us.y + Img_t->get(site), sigmaD);

	Neighbor.x = .0;
	Neighbor.y = .0;
	if (x > 0 && domain_map->get(x - 1, y) == center_domain) {
		inner_prod = us * u->get(x - 1, y);
		coeff = 0.5 * (1.0 + inner_prod / (norm(us) * norm(u->get(x - 1, y))));
		Neighbor.x += coeff * (*psiS)(us.x - u->get(x - 1, y).x, sigmaS);
		Neighbor.y += coeff * (*psiS)(us.y - u->get(x - 1, y).y, sigmaS);
	}
	if (x < u->width() - 1 && domain_map->get(x + 1, y) == center_domain) {
		inner_prod = us * u->get(x + 1, y);
		coeff = 0.5 * (1.0 + inner_prod / (norm(us) * norm(u->get(x + 1, y))));
		Neighbor.x += coeff * (*psiS)(us.x - u->get(x + 1, y).x, sigmaS);
		Neighbor.y += coeff * (*psiS)(us.y - u->get(x + 1, y).y, sigmaS);
	}
	if (y > 0 && domain_map->get(x, y - 1) == center_domain) {
		inner_prod = us * u->get(x, y - 1);
		coeff = 0.5 * (1.0 + inner_prod / (norm(us) * norm(u->get(x, y - 1))));
		Neighbor.x += coeff * (*psiS)(us.x - u->get(x, y - 1).x, sigmaS);
		Neighbor.y += coeff * (*psiS)(us.y - u->get(x, y - 1).y, sigmaS);
	}
	if (y < u->height() - 1 && domain_map->get(x, y + 1) == center_domain) {
		inner_prod = us * u->get(x, y + 1);
		coeff = 0.5 * (1.0 + inner_prod / (norm(us) * norm(u->get(x, y + 1))));
		Neighbor.x += coeff * (*psiS)(us.x - u->get(x, y + 1).x, sigmaS);
		Neighbor.y += coeff * (*psiS)(us.y - u->get(x, y + 1).y, sigmaS);
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
Error_MultipleMotion_Block(const ImgVector<VECTOR_2D<double> >* u, const ImgVector<size_t>* domain_map, const ImgVector<VECTOR_2D<double> >* Img_g, const ImgVector<double>* Img_t, const double& lambdaD, const double& lambdaS, const double& sigmaD, const double& sigmaS)
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
			size_t center_domain = domain_map->get(x, y);
			VECTOR_2D<double> us(u->get(x, y));
			VECTOR_2D<double> Neighbor(0.0, 0.0);
			double inner_prod;
			double coeff;
			if (x > 0 && domain_map->get(x - 1, y) == center_domain) {
				inner_prod = us * u->get(x - 1, y);
				coeff = 0.5 * (1.0 + inner_prod / (norm(us) * norm(u->get(x - 1, y))));
				Neighbor.x += coeff * (*rhoS)(us.x - u->get(x - 1, y).x, sigmaS);
				Neighbor.y += coeff * (*rhoS)(us.y - u->get(x - 1, y).y, sigmaS);
			}
			if (x < u->width() - 1 && domain_map->get(x + 1, y) == center_domain) {
				inner_prod = us * u->get(x + 1, y);
				coeff = 0.5 * (1.0 + inner_prod / (norm(us) * norm(u->get(x + 1, y))));
				Neighbor.x += coeff * (*rhoS)(us.x - u->get(x + 1, y).x, sigmaS);
				Neighbor.y += coeff * (*rhoS)(us.y - u->get(x + 1, y).y, sigmaS);
			}
			if (y > 0 && domain_map->get(x, y - 1) == center_domain) {
				inner_prod = us * u->get(x, y - 1);
				coeff = 0.5 * (1.0 + inner_prod / (norm(us) * norm(u->get(x, y - 1))));
				Neighbor.x += coeff * (*rhoS)(us.x - u->get(x, y - 1).x, sigmaS);
				Neighbor.y += coeff * (*rhoS)(us.y - u->get(x, y - 1).y, sigmaS);
			}
			if (y < u->height() - 1 && domain_map->get(x, y + 1) == center_domain) {
				inner_prod = us * u->get(x, y + 1);
				coeff = 0.5 * (1.0 + inner_prod / (norm(us) * norm(u->get(x, y + 1))));
				Neighbor.x += coeff * (*rhoS)(us.x - u->get(x, y + 1).x, sigmaS);
				Neighbor.y += coeff * (*rhoS)(us.y - u->get(x, y + 1).y, sigmaS);
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
MultipleMotion_write(const ImgVector<double>& img_prev, const ImgVector<double>& img_current, const int MaxInt, const std::vector<ImgVector<Vector_ST<double> > >& u, const std::string& filename)
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
	pnm.copy(PORTABLE_GRAYMAP_BINARY, compensated.width(), compensated.height(), MaxInt, compensated.ref_image_compensated().data(), 1.0);
	pnm.write(filename_compensated.c_str());
	pnm.free();
}

void
MultipleMotion_write(const ImgVector<RGB>& img_prev, const ImgVector<RGB>& img_current, const int MaxInt, const std::vector<ImgVector<Vector_ST<double> > >& u, const std::string &filename)
{
	ERROR Error("MultipleMotion_write");
	FILE *fp = nullptr;
	MotionCompensation<RGB> compensated(img_prev, img_current, u[0]);
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
		fprintf(stderr, "void MultipleMotion_write(const ImgVector<RGB>*, const ImgVector<RGB>*, const ImgVector<VECTOR_2D<double> >*, const std::string&) : Cannot allocate memory\n");
		throw;
	}
	for (size_t n = 0; n < size; n++) {
		compensated_image[n] = int(compensated.ref_image_compensated().get(n).R);
		compensated_image[n + size] = int(compensated.ref_image_compensated().get(n).G);
		compensated_image[n + 2 * size] = int(compensated.ref_image_compensated().get(n).B);
	}
	pnm.copy(PORTABLE_PIXMAP_BINARY, compensated.width(), compensated.height(), MaxInt, compensated_image);
	delete[] compensated_image;
	pnm.write(filename_compensated.c_str());
	pnm.free();
}

void
MultipleMotion_write(const ImgVector<double>& img_prev, const ImgVector<double>& img_current, const ImgVector<double>& img_next, const int MaxInt, const std::vector<ImgVector<Vector_ST<double> > >& u, const std::string& filename)
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
	pnm.copy(PORTABLE_GRAYMAP_BINARY, compensated.width(), compensated.height(), MaxInt, compensated.ref_image_compensated().data(), 1.0);
	pnm.write(filename_compensated.c_str());
	pnm.free();
}

void
MultipleMotion_write(const ImgVector<RGB>& img_prev, const ImgVector<RGB>& img_current, const ImgVector<RGB>& img_next, const int MaxInt, const std::vector<ImgVector<Vector_ST<double> > >& u, const std::string &filename)
{
	ERROR Error("MultipleMotion_write");

	FILE *fp = nullptr;
	VECTOR_2D<double> v;
	int x, y;
	MotionCompensation<RGB> compensated;
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
	// Saturation
	ImgVector<RGB> compensated_image(compensated.ref_image_compensated());
	compensated_image.saturate(
	    static_cast<RGB (*)(const RGB&, const double&, const double&)>(saturate),
	    0.0, double(MaxInt));
	size_t size = compensated_image.size();
	std::string::size_type found = filename.find_last_of("/\\");
	filename_compensated = filename.substr(0, found + 1) + "compensated_" + filename.substr(found + 1);
	printf("* Output The Compensated Image from Optical Flow to '%s'(binary)\n\n", filename_compensated.c_str());
	{
		int* image_tmp = nullptr;
		try {
			image_tmp = new int[size * 3];
		}
		catch (const std::bad_alloc& bad) {
			std::cerr << bad.what() << std::endl;
			fprintf(stderr, "void MultipleMotion_write(const ImgVector<RGB>*, const ImgVector<RGB>*, const ImgVector<VECTOR_2D<double> >*, const std::string&) : Cannot allocate memory\n");
			throw;
		}
		// Copy for output
		for (size_t n = 0; n < size; n++) {
			image_tmp[n] = int(compensated_image.get(n).R);
			image_tmp[n + size] = int(compensated_image.get(n).G);
			image_tmp[n + 2 * size] = int(compensated_image.get(n).B);
		}
		pnm.copy(PORTABLE_PIXMAP_BINARY, compensated_image.width(), compensated_image.height(), MaxInt, image_tmp);
		delete[] image_tmp;
		image_tmp = nullptr;
	}
	pnm.write(filename_compensated.c_str());
	pnm.free();
}

