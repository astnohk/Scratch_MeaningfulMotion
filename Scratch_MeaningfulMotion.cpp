#include "Scratch_MeaningfulMotion.h"
#include "OpticalFlow/Affine_MultipleMotion.h"
#include "OpticalFlow/OpticalFlow.h"
#include "OpticalFlow/OpticalFlow_BlockMatching.h"
#include "HOG/HOG.h"




int
Scratch_MeaningfulMotion(char *OutputName, char *InputName, unsigned int OutputNameLength, unsigned int InputNameLength, int Start, int End, OPTIONS Options, FILTER_PARAM FilterParam)
{
	ERROR Error("Scratch_MeaningfulMotion");
	const char *Bars = "------------------------------------------------";

	const std::string FilterNames[] = {"undefined", "Epsilon", "Gaussian"};
	std::string InputNameNums;
	std::string OutputNameNums;
	char *char_tmp = nullptr;
	int CurrentFileNum;

	PNM pnm_in;
	PNM pnm_out;
	PNM pnm_orig;
	PNM pnm_res;
	PNM_DOUBLE pnmd_in;
	PNM_DOUBLE pnmd_out;
	ImgVector<double> imgd_prev;
	ImgVector<double> imgd_in;
	ImgVector<int> img_orig; // For X11 plotting

	// Optical Flow
	VECTOR_AFFINE MultipleMotion_AffineCoeff;
	ImgVector<VECTOR_2D<double> > *MultipleMotion_u = nullptr;
	// HOG
	HOG hog_raw;
	HOG hog;
	HOG hog_raw_prv;
	HOG hog_prv;
	VECTOR_2D_W_SCORE *hog_vector = nullptr;

	int *k_list = nullptr;
	ImgVector<double> *Pr_table = nullptr;

	ImgVector<double> *filtered = nullptr;
	ImgVector<double> *scratches = nullptr;
	ImgVector<double> *angles = nullptr;

	SEGMENT *MaximalSegments = nullptr;
	SEGMENT *EPSegments = nullptr;
	int *segments = nullptr;
	int Num_Segments = 0;

	int Initialize = 0;

	SIZE size;
	SIZE size_prev;
	SIZE size_orig;
	SIZE size_res;
	SIZE size_out;
	int maxMN = 0;
	int l_min = 1;
	int i, k, L;
	double progress;
	int count;

	for (CurrentFileNum = Start; CurrentFileNum <= End; CurrentFileNum++) {
		// Read PNM Files
		if (strchr(InputName, '%') == nullptr) {
			InputNameNums = InputName;
		} else {
			try {
				char_tmp = new char[InputNameLength + 1];
			}
			catch (const std::bad_alloc &bad) {
				Error.Value("char_tmp");
				Error.Malloc();
				goto ExitError;
			}
			sprintf(char_tmp, InputName, CurrentFileNum);
			InputNameNums = char_tmp;
			delete[] char_tmp;
			char_tmp = nullptr;
		}
		if (pnm_orig.read(InputNameNums.c_str()) == PNM_FUNCTION_ERROR) {
			Error.Function("pnm_orig.read");
			Error.File(InputNameNums.c_str());
			Error.FileRead();
			goto ExitError;
		}
		// END Read
		if (size.height == 0 || size.width == 0) { // Initialize size.height (width) and size.width (height)
			size.height = pnm_orig.Height();
			size.width = pnm_orig.Width();
			size_orig.height = pnm_orig.Height();
			size_orig.width = pnm_orig.Width();
			size_prev = size_orig;
			maxMN = (size.height > size.width ? size.height : size.width);
			size_res = Options.ResampleSize;
			if (size_res.width == 0) {
				size_res.width = pnm_orig.Width();
			}
			if (size_res.height == 0) {
				size_res.height = pnm_orig.Height();
			}
			if ((Options.PlotOptions & PLOT_AS_RESAMPLE) != 0) {
				size_out = size_res;
			} else {
				size_out = size;
			}
		}
		if (size_prev.height != pnm_orig.Height() || size_prev.width != pnm_orig.Width()) {
			Error.Others("Image size are not match with previous one");
			goto ExitError;
		}
		printf("- The input image size is %dx%d\n- and bit depth is %d\n", size.width, size.height, (int)round(log2((double)pnm_orig.MaxInt())));
		if (size_res.width > 0 || size_res.height > 0) { // Resample
			size = size_res;
			if (pnmd_in.copy(pnm_orig, 1.0) != PNM_FUNCTION_SUCCESS) {
				Error.Function("pnmd_in.copy");
				Error.Value("pnmd_in <- pnmd_orig");
				Error.FunctionFail();
				goto ExitError;
			}
			if (pnm_resize(&pnmd_out, pnmd_in, size_res.width, size_res.height, Options.ResampleMethod) != PNM_FUNCTION_SUCCESS) {
				Error.Function("pnm_resize");
				Error.Value("(pnmd_in -> pnmd_out)");
				Error.FunctionFail();
				goto ExitError;
			}
			if (pnm_res.copy(pnmd_out, 1.0, "round") != PNM_FUNCTION_SUCCESS) {
				Error.Function("pnm_double2int");
				Error.Value("(pnmd_out -> pnm_res)");
				Error.FunctionFail();
				goto ExitError;
			}
			if (pnm_in.copy(pnmd_out, 1.0, "round") != PNM_FUNCTION_SUCCESS) {
				Error.Function("pnm_double2int");
				Error.Value("(pnmd_out -> pnm_in)");
				Error.FunctionFail();
				goto ExitError;
			}
			pnmd_in.free();
			pnmd_out.free();
		} else {
			if (pnm_in.copy(pnm_orig) != PNM_FUNCTION_SUCCESS) {
				Error.Function("pnm_in.copy");
				Error.Value("pnm_in <- pnm_orig");
				Error.FunctionFail();
				goto ExitError;
			}
		}
		if ((Options.PlotOptions & PLOT_RESAMPLED_IMG_ONLY) != 0) {
			// Just output only the resampled image
			if (pnm_out.copy(pnm_in) != PNM_FUNCTION_SUCCESS) {
				Error.Function("pnm_out.copy");
				Error.Value("pnm_out <- pnm_in");
				Error.FunctionFail();
				goto ExitError;
			}
			goto Write;
		}
		if ((pnm_in.Desc() % 3) == 0) { // Convert to Grayscale
			printf("- The input image is color data\n");
			printf("* Convert the image to grayscale before applying Meaningful Alignments.\n");
			printf("Convert...   ");
			if (pnmd_in.copy(pnm_in, 1.0) != PNM_FUNCTION_SUCCESS) {
				Error.Function("pnmd_in.copy");
				Error.Value("pnmd_in <- pnm_in");
				Error.FunctionFail();
				goto ExitError;
			}
			if (pnmd_out.RGB2Gray(pnmd_in) != PNM_FUNCTION_SUCCESS) {
				Error.Function("pnmd_out.RGB2Gray");
				Error.Value("pnmd_out <- pnmd_in");
				Error.FunctionFail();
				goto ExitError;
			}
			pnmd_in.free();
			pnm_in.free();
			if (pnm_in.copy(pnmd_out, 1.0, "round") != PNM_FUNCTION_SUCCESS) {
				Error.Function("pnm_in.copy");
				Error.Value("pnm_in <- pnmd_out");
				Error.FunctionFail();
				goto ExitError;
			}
			pnmd_out.free();
			printf("Finished\n\n");
		}
		imgd_in.reset(pnm_in.Width(), pnm_in.Height());
		for (i = 0; i < imgd_in.size(); i++) {
			imgd_in[i] = (double)pnm_in[i];
		}

		// Show Parameters
		printf("\n      --- Parameters ---\n  %s\n", Bars);
		if (size_res.width > 0 || size_res.height > 0) {
			printf("  | Resample (%d), 0:z-hold, 1:bicubic\n", Options.ResampleMethod);
			printf("  |   %dx%d -> %dx%d\n", pnm_orig.Width(), pnm_orig.Width(), size_res.width, size_res.height);
		}
		printf("  | filter type = %s\n", FilterNames[FilterParam.type < NUM_FILTER_TYPE ? FilterParam.type : 0].c_str());
		printf("  | filter size = %dx%d\n", FilterParam.size.width, FilterParam.size.height);
		switch (FilterParam.type) {
			case FILTER_ID_EPSILON:
				printf("  | filter epsilon = %.0f\n", FilterParam.epsilon);
				break;
			case FILTER_ID_GAUSSIAN:
				printf("  | Gaussian filter's standard deviation = %f\n", FilterParam.std_deviation);
		}
		printf("  | Scratch Detection :\n");
		printf("  |   s_med = %d\n", Options.s_med);
		printf("  |   s_avg = %d\n", Options.s_avg);
		printf("  |   p = %f\n", Options.p);
		printf("  | Meaningful Alignments :\n");
		printf("  |   precision = %f [deg]\n", 360.0 * Options.p);
		printf("  |   epsilon = %e\n", Options.ep);
		printf("  |   searching angle = [%f, %f] [deg] (from Vertical line)\n", -0.5 * 180.0 / DIV_ANGLE_VERTICAL, 0.5 * 180.0 / DIV_ANGLE_VERTICAL);
		printf("  |   exclusive radius = %f [px]\n", Options.Exclusive_Max_Radius);
		if (Options.Max_Length > 0) {
			printf("  |   Max length = %d [px]\n", Options.Max_Length);
		} else {
			printf("  |   NO search length limit\n");
		}
		if (Options.Max_Output_Length > 0) {
			printf("  |   Max output length = %d [px]\n", Options.Max_Output_Length);
		} else {
			printf("  |   NO output length limit\n");
		}
		if (Options.ExclusivePrinciple != 0) {
			printf("  |   Turn ON Exclusive Principle\n");
		}
		printf("  %s\n\n", Bars);


		if ((Options.mode & MODE_OUTPUT_FILTERED_IMAGE) != 0) {
			// Output filtered image
			printf("* Filtering\n");
			filtered = DetectScratch(pnm_in, Options.s_med, Options.s_avg, FilterParam, DO_NOT_DETECTION);
			if (filtered == nullptr) {
				Error.Function("DetectScratch");
				Error.Value("filtered");
				Error.FunctionFail();
				goto ExitError;
			}
			printf("* Output Filtered Image\n");
			if (pnm_out.copy(PORTABLE_GRAYMAP_BINARY, pnm_in.Width(), pnm_in.Height(), pnm_in.MaxInt(), filtered->data(), 1.0) != PNM_FUNCTION_SUCCESS) {
				Error.Function("pnm_out.copy");
				Error.Value("pnm_out");
				Error.FunctionFail();
				goto ExitError;
			}
			delete filtered;
			filtered = nullptr;
		} else if ((Options.mode & MODE_OUTPUT_MULTIPLE_MOTIONS_AFFINE) != 0) {
			// Computte and output Multiple Motion Affine Parameters by method of M.J.Black
			if (imgd_prev.isNULL() != false) {
				printf("* Skip Calculate Multiple Motions by Affine while there is NOT any previous frame\n");
			} else {
				printf("* Compute Multiple Motions Affine Parameters by method of M.J.Black\n");
				MultipleMotion_AffineCoeff = MultipleMotion_Affine(&imgd_prev, &imgd_in, pnm_in.MaxInt(), Options.MultipleMotion_Param);
			}
		} else if ((Options.mode & MODE_OUTPUT_OPTICALFLOW) != 0) {
			// Computte and output Multiple Motion Optical Flow by method of M.J.Black
			if (imgd_prev.isNULL() != false) {
				printf("* Skip Calculate Optical Flow while there is NOT any previous frame\n");
			} else {
				printf("* Compute Optical Flow by method of M.J.Black\n");
				//MultipleMotion_u = OpticalFlow_Pyramid(&imgd_prev, &imgd_in, pnm_in.MaxInt(), Options.MultipleMotion_Param);
				MultipleMotion_u = OpticalFlow_BlockMatching(&imgd_prev, &imgd_in, pnm_in.MaxInt(), Options.MultipleMotion_Param);
			}
		} else if ((Options.mode & MODE_OUTPUT_HISTOGRAMS_OF_ORIENTED_GRADIENTS) != 0
		    || (Options.mode & MODE_OUTPUT_HISTOGRAMS_OF_ORIENTED_GRADIENTS_RAW_HOG) != 0
		    || (Options.mode & MODE_OUTPUT_HISTOGRAMS_OF_ORIENTED_GRADIENTS_MATCHING_VECTOR) != 0) {
			printf("* Compute HOG\n");
			pnmd_in.copy(pnm_in, 1.0 / pnm_in.MaxInt());
			hog_raw_prv.copy(hog_raw);
			hog_prv.copy(hog);
			hog_raw.free();
			hog.free();
			HistogramsOfOrientedGradients(&hog_raw, &hog, pnmd_in, Options.HOG_Param);
		} else {
			// Scratch Detection
			printf("* Detect Scratch like vertical lines\n");
			scratches = DetectScratch(pnm_in, Options.s_med, Options.s_avg, FilterParam, DO_DETECTION);
			if (scratches == nullptr) {
				Error.Function("DetectScratch");
				Error.Value("scratches");
				Error.FunctionFail();
				goto ExitError;
			}
			pnm_in.copy(pnm_in.Desc(), pnm_in.Width(), pnm_in.Height(), pnm_in.MaxInt(), scratches->data(), 1.0);
			if ((Options.mode & MODE_OUTPUT_BINARY_IMAGE) != 0) {
				// Output Scratch Map without Meaningful Alignments
				printf("* Output Scratch binary image\n");
				if (pnm_out.copy(pnm_in) != PNM_FUNCTION_SUCCESS) {
					Error.Function("pnm_out.copy");
					Error.Value("(pnm_out <- pnm_in)");
					Error.FunctionFail();
					goto ExitError;
				}
				for (int i = 0; i < pnm_out.Width() * pnm_out.Height(); i++) {
					pnm_out[i] = scratches->get(i);
				}
			} else {
				// A Contrario Method : Meaningful Alignments
				if (Initialize == 0) {
					Initialize = 1;
					try {
						Pr_table = new ImgVector<double>(maxMN + 1, maxMN + 1);
					}
					catch (const std::bad_alloc &bad) {
						Error.Function("new");
						Error.Value("Pr_table");
						Error.Malloc();
						goto ExitError;
					}
					printf("* Calculate Pr(k, L) table :\n[L =     0]   0.0%% |%s\x1b[1A\n", Progress_End.c_str());
					progress = 0;
					count = 0;
#pragma omp parallel for schedule(dynamic) private(k)
					for (L = 1; L <= maxMN; L++) {
						for (k = 0; k <= L; k++) {
							Pr_table->ref(k, L) = Pr(k, L, Options.p);
						}
#pragma omp critical
						{
							count++;
							if (round((double)count / maxMN * 1000.0) > progress) {
								progress = round((double)count / maxMN * 1000.0);
								printf("\r[L = %5d] %5.1f%% |%s#\x1b[1A\n", count, progress * 0.1, Progress[NUM_PROGRESS * count / (1 + maxMN)].c_str());
							}
						}
					}
					printf("\nComplete!\n");
					l_min = (int)ceil((log(Options.ep) - (log(DIV_ANGLE) + log((double)size.height) + 2.0 * log((double)size.width))) / log(Options.p));
					if (l_min < 1) {
						l_min = 1;
					}
					printf("* Compute k_list\n");
					k_list = Calc_k_l(size, Options.p, Options.ep);
					if (k_list == nullptr) {
						Error.Function("Calc_k_l");
						Error.Value("k_list");
						Error.FunctionFail();
						goto ExitError;
					}
				}

				printf("* Compute Direction Field\n");
				angles = DerivativeAngler(scratches);
				if (angles == nullptr) {
					Error.Function("Derivation");
					Error.Value("angles");
					Error.FunctionFail();
					goto ExitError;
				}
				printf("* Compute Segments and Maximal Meaningfulness\n");
				if (MaximalSegments == nullptr) {
					delete[] MaximalSegments;
				}
				MaximalSegments = AlignedSegment_vertical(angles, k_list, l_min, Pr_table, &Num_Segments, Options.Max_Length, Options.Max_Output_Length);
				if (MaximalSegments == nullptr) {
					Error.Function("AlignedSegment_vertical");
					Error.Value("MaximalSegments");
					Error.FunctionFail();
					goto ExitError;
				}
				printf("- Found (%d) Maximal Meaningful Segments\n", Num_Segments);
				if (Options.ExclusivePrinciple != 0) {
					printf("* Delete Redundant Segments by Exclusive Principle\n");
					EPSegments = ExclusivePrinciple(angles, k_list, Pr_table, MaximalSegments, &Num_Segments, Options.Exclusive_Max_Radius);
					if (EPSegments == nullptr) {
						Error.Function("ExclusivePrinciple");
						Error.Value("EPSegments");
						Error.FunctionFail();
						goto ExitError;
					}
					printf("- Reduced to (%d) EP-Maximal Meaningful Segments\n", Num_Segments);
					delete[] MaximalSegments;
					MaximalSegments = EPSegments;
					EPSegments = nullptr;
				}
				printf("* Plot The Segments");
				if (Options.Max_Output_Length > 0) {
					printf(" that satisfy (length < %d)", Options.Max_Output_Length);
				}
				printf("\n");
				segments = PlotSegment(MaximalSegments, Num_Segments, size, size_out, Options.PlotOptions & PLOT_NEGATE);
				if (segments == nullptr) {
					Error.Function("PlotSegment");
					Error.Value("segments");
					Error.FunctionFail();
					goto ExitError;
				}
				if (Options.Superimpose != 0) {
					printf("* Superimpose plot image on original image\n");
					if ((Options.PlotOptions & PLOT_AS_RESAMPLE) != 0) {
						if (Superimposer(&pnm_out, pnm_res, segments, size_out, Options.Superimpose, Options.PlotOptions & PLOT_NEGATE) != MEANINGFUL_SUCCESS) {
							Error.Function("Superimposer");
							Error.Value("pnm_out");
							Error.FunctionFail();
							goto ExitError;
						}
					} else {
						if (Superimposer(&pnm_out, pnm_orig, segments, size_out, Options.Superimpose, Options.PlotOptions & PLOT_NEGATE) != MEANINGFUL_SUCCESS) {
							Error.Function("Superimposer");
							Error.Value("pnm_out");
							Error.FunctionFail();
							goto ExitError;
						}
					}
				} else {
					pnm_out.copy(PORTABLE_GRAYMAP_BINARY, size_out.width, size_out.height, pnm_orig.MaxInt(), segments);
				}
				delete[] segments;
				segments = nullptr;
				delete angles;
				angles = nullptr;
			}
			delete scratches;
			scratches = nullptr;
		}
		// X11 Plotting
		if (Options.x11_plot) {
			img_orig.reset(pnm_orig.Width(), pnm_orig.Height());
			for (i = 0; i < img_orig.size(); i++) {
				img_orig[i] = (int)pnm_orig[i];
			}
			ShowSegments_X11(&img_orig, size, pnm_orig.MaxInt(), MaximalSegments, Num_Segments);
		}
		delete[] MaximalSegments;
		MaximalSegments = nullptr;

Write:
		if (strchr(OutputName, '%') == nullptr) {
			OutputNameNums = OutputName;
		} else {
			try {
				char_tmp = new char[OutputNameLength + 1];
			}
			catch (const std::bad_alloc &bad) {
				Error.Value("char_tmp");
				Error.Malloc();
				goto ExitError;
			}
			sprintf(char_tmp, OutputName, CurrentFileNum);
			OutputNameNums = char_tmp;
			delete[] char_tmp;
			char_tmp = nullptr;
		}
		if ((Options.mode & MODE_OUTPUT_MULTIPLE_MOTIONS_AFFINE) != 0) {
			if (imgd_prev.isNULL() == false) {
				MultipleMotion_Affine_write(MultipleMotion_AffineCoeff, OutputNameNums);
			}
		} else if ((Options.mode & MODE_OUTPUT_OPTICALFLOW) != 0) {
			if (imgd_prev.isNULL() == false) {
				MultipleMotion_write(&imgd_prev, &imgd_in, MultipleMotion_u, OutputNameNums);
			}
		} else if ((Options.mode & MODE_OUTPUT_HISTOGRAMS_OF_ORIENTED_GRADIENTS_RAW_HOG) != 0) {
			if (HOG_write(hog_raw, OutputNameNums) == false) {
				Error.Function("HOG_write");
				Error.Value("hog");
				Error.FunctionFail();
				goto ExitError;
			}
		} else if ((Options.mode & MODE_OUTPUT_HISTOGRAMS_OF_ORIENTED_GRADIENTS) != 0) {
			if (HOG_write(hog, OutputNameNums) == false) {
				Error.Function("HOG_write");
				Error.Value("hog");
				Error.FunctionFail();
				goto ExitError;
			}
		} else if ((Options.mode & MODE_OUTPUT_HISTOGRAMS_OF_ORIENTED_GRADIENTS_MATCHING_VECTOR) != 0) {
			if (imgd_prev.isNULL() == false
			    && hog_prv.Bins() == hog.Bins()) {
				printf("* Compute matching each images HOG feature\n");
				hog_vector = HOG_Matching(&hog_prv, &hog);
				if (HOG_vector_write(hog_vector, hog.Width(), hog.Height(), OutputNameNums) == false) {
					Error.Function("HOG_vector_write");
					Error.Value("hog_vector");
					Error.FunctionFail();
					goto ExitError;
				}
				HOG_vector_compensated_write(&imgd_prev, &imgd_in, hog_vector, hog.Width(), hog.Height(), OutputNameNums);
				delete[] hog_vector;
				hog_vector = nullptr;
			} else {
				printf("There are NO previous image or HOG data\n");
			}
		} else {
			if (pnm_out.write(OutputNameNums.c_str()) == PNM_FUNCTION_ERROR) {
				Error.Function("pnm_out.write");
				Error.File(OutputNameNums.c_str());
				Error.FileWrite();
				goto ExitError;
			}
		}
		delete MultipleMotion_u;
		MultipleMotion_u = nullptr;
		if (imgd_in.isNULL() == false) {
			imgd_prev = imgd_in;
		}
		pnmd_in.free();
		pnm_out.free();
		pnm_in.free();
		pnm_res.free();
		pnm_orig.free();
	}
	delete[] hog_vector;
	delete[] k_list;
	delete Pr_table;
	return MEANINGFUL_SUCCESS;
// Exit Error
ExitError:
	delete[] hog_vector;
	delete[] segments;
	delete[] EPSegments;
	delete[] MaximalSegments;
	delete[] k_list;
	delete MultipleMotion_u;
	delete angles;
	delete scratches;
	delete filtered;
	delete Pr_table;
	return MEANINGFUL_FAILURE;
}

