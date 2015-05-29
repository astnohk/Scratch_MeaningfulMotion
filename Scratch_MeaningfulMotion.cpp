#include "Scratch_MeaningfulMotion.h"
#include "Affine_MultipleMotion.h"
#include "OpticalFlow_AffineParamet.h"
#include "OpticalFlow_MultipleMotion.h"




int
Scratch_MeaningfulMotion(char *OutputName, char *InputName, unsigned int OutputNameLength, unsigned int InputNameLength, int Start, int End, OPTIONS Options, FILTER_PARAM FilterParam)
{
	ERROR Error("Scratch_MeaningfulMotion()");
	char Bars[] = "------------------------------------------------";

	const std::string FilterNames[] = {"undefined", "Epsilon", "Gaussian"};
	char *InputNameNums = NULL;
	char *OutputNameNums = NULL;
	int CurrentFileNum;
	PNM pnm_in = PNM_NULL;
	PNM pnm_out = PNM_NULL;
	PNM pnm_orig = PNM_NULL;
	PNM pnm_res = PNM_NULL;
	PNM_DOUBLE pnmd_in = PNM_DOUBLE_NULL;
	PNM_DOUBLE pnmd_out = PNM_DOUBLE_NULL;
	PNM_DOUBLE pnmd_prev = PNM_DOUBLE_NULL;
	VECTOR_AFFINE MultipleMotion_AffineCoeff;
	VECTOR_2D *MultipleMotion_u;
	TUPLE_VEC_SCALAR *OpticalFlow_Affine = NULL;

	int Initialize = 0;

	int *k_list = NULL;
	double *Pr_table = NULL;
	double *filtered = NULL;
	double *scratches = NULL;
	double *binary = NULL;
	double *angles = NULL;
	SEGMENT *MaximalSegments = NULL;
	SEGMENT *EPSegments = NULL;
	int *segments = NULL;
	int Num_Segments = 0;
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

	try {
		InputNameNums = new char[InputNameLength + 1u];
	}
	catch (std::bad_alloc bad) {
		Error.Function("new");
		Error.Value("InputNameNums");
		Error.Malloc();
		goto ExitError;
	}
	try {
		OutputNameNums = new char[OutputNameLength + 1u];
	}
	catch (std::bad_alloc bad) {
		Error.Function("new");
		Error.Value("OutputNameNums");
		Error.Malloc();
		goto ExitError;
	}
	for (CurrentFileNum = Start; CurrentFileNum <= End; CurrentFileNum++) {
		// Read PNM Files
		if (strchr(InputName, '%') == NULL) {
			sprintf(InputNameNums, InputName);
		} else {
			sprintf(InputNameNums, InputName, CurrentFileNum);
		}
		if (pnmread(&pnm_orig, InputNameNums) == PNM_FUNCTION_ERROR) {
			fprintf(stderr, "*** Scratch_MeaningfulA error - Failed to read the PNM file \"%s\" ***\n", InputNameNums);
			Error.Function("pnmread");
			Error.File(InputNameNums);
			Error.FileRead();
			goto ExitError;
		}
		// END Read

		if (((int)pnm_orig.height < 0) || (int)pnm_orig.width < 0) {
			fprintf(stderr, "The size of the image is too large to processing (image size OVERFLOWED)");
			goto ExitError;
		}
		if (size.height == 0 || size.width == 0) { // Initialize size.height (width) and size.width (height)
			size.height = (int)pnm_orig.height;
			size.width = (int)pnm_orig.width;
			size_orig.height = (int)pnm_orig.height;
			size_orig.width = (int)pnm_orig.width;
			size_prev = size_orig;
			maxMN = (size.height > size.width ? size.height : size.width);
			size_res = Options.ResampleSize;
			if (size_res.width == 0) {
				size_res.width = (int)pnm_orig.width;
			}
			if (size_res.height == 0) {
				size_res.height = (int)pnm_orig.height;
			}
			if ((Options.PlotOptions & PLOT_AS_RESAMPLE) != 0) {
				size_out = size_res;
			} else {
				size_out = size;
			}
			printf("* Initialize atan2_div_pi_table()\n");
			atan2_div_pi_table(0, 0, &size);
		}
		if (size_prev.height != (int)pnm_orig.height || size_prev.width != (int)pnm_orig.width) {
			Error.Others("Image size are not match with previous one");
			goto ExitError;
		}
		printf("- The input image size is %dx%d\n- and bit depth is %d\n", size.width, size.height, (int)round(log2((double)pnm_orig.maxint)));
		if (size_res.width > 0 || size_res.height > 0) { /* Resample */
			size = size_res;
			if (pnm_int2double(&pnmd_in, &pnm_orig, 1.0, NULL) != PNM_FUNCTION_SUCCESS) {
				Error.Function("pnm_int2double");
				Error.Value("(pnmd_orig -> pnmd_in)");
				Error.FunctionFail();
				goto ExitError;
			}
			if (pnm_resize(&pnmd_out, &pnmd_in, size_res.width, size_res.height, Options.ResampleMethod.c_str()) != PNM_FUNCTION_SUCCESS) {
				Error.Function("pnm_resize");
				Error.Value("(pnmd_in -> pnmd_out)");
				Error.FunctionFail();
				goto ExitError;
			}
			pnmdouble_free(&pnmd_in);
			if (pnm_double2int(&pnm_res, &pnmd_out, 1.0, "round", NULL) != PNM_FUNCTION_SUCCESS) {
				Error.Function("pnm_double2int");
				Error.Value("(pnmd_out -> pnm_res)");
				Error.FunctionFail();
				goto ExitError;
			}
			if (pnm_double2int(&pnm_in, &pnmd_out, 1.0, "round", NULL) != PNM_FUNCTION_SUCCESS) {
				Error.Function("pnm_double2int");
				Error.Value("(pnmd_out -> pnm_in)");
				Error.FunctionFail();
				goto ExitError;
			}
			pnmdouble_free(&pnmd_out);
		} else {
			if (pnmcp(&pnm_in, &pnm_orig) != PNM_FUNCTION_SUCCESS) {
				Error.Function("pnmcp");
				Error.Value("(pnm_orig -> pnm_in)");
				Error.FunctionFail();
				goto ExitError;
			}
		}
		if ((Options.PlotOptions & PLOT_RESAMPLED_IMG_ONLY) != 0) {
			/* Just output only the resampled image */
			if (pnmcp(&pnm_out, &pnm_in) != PNM_FUNCTION_SUCCESS) {
				Error.Function("pnmcp");
				Error.Value("(pnm_in -> pnm_out)");
				Error.FunctionFail();
				goto ExitError;
			}
			goto Write;
		}
		if ((pnm_in.desc % 3) == 0) { /* Convert to Grayscale */
			printf("- The input image is color data\n");
			printf("* Convert the image to grayscale before applying Meaningful Alignments.\n");
			printf("Convert...   ");
			if (pnm_int2double(&pnmd_in, &pnm_in, 1.0, NULL) != PNM_FUNCTION_SUCCESS) {
				Error.Function("pnm_int2double");
				Error.Value("(pnm_in -> pnmd_in)");
				Error.FunctionFail();
				goto ExitError;
			}
			if (pnm_RGB2Gray(&pnmd_out, &pnmd_in) != PNM_FUNCTION_SUCCESS) {
				Error.Function("pnm_RGB2Gray");
				Error.Value("(pnmd_in -> pnmd_out)");
				Error.FunctionFail();
				goto ExitError;
			}
			pnmdouble_free(&pnmd_in);
			pnmfree(&pnm_in);
			if (pnm_double2int(&pnm_in, &pnmd_out, 1.0, "round", NULL) != PNM_FUNCTION_SUCCESS) {
				Error.Function("pnm_double2int");
				Error.Value("(pnmd_out -> pnm_in)");
				Error.FunctionFail();
				goto ExitError;
			}
			pnmdouble_free(&pnmd_out);
			printf("Finished\n\n");
		}

		/* Show Parameters */
		printf("\n      --- Parameters ---\n  %s\n", Bars);
		if (size_res.width > 0 || size_res.height > 0) {
			printf("  | Resample (%s)\n", Options.ResampleMethod.c_str());
			printf("  |   %ux%u -> %dx%d\n", pnm_orig.width, pnm_orig.width, size_res.width, size_res.height);
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
			/* Output filtered image */
			printf("* Filtering\n");
			filtered = DetectScratch(&pnm_in, Options.s_med, Options.s_avg, FilterParam, DO_NOT_DETECTION);
			if (filtered == NULL) {
				Error.Function("DetectScratch");
				Error.Value("filtered");
				Error.FunctionFail();
				goto ExitError;
			}
			printf("* Output Filtered Image\n");
			if (pnmnew(&pnm_out, PORTABLE_GRAYMAP_BINARY, size_res.width, size_res.height, pnm_in.maxint) != PNM_FUNCTION_SUCCESS) {
				Error.Function("pnmnew");
				Error.Value("pnm_out");
				Error.Malloc();
				goto ExitError;
			}
			for (i = 0; i < size_res.width * size_res.height; i++) {
				pnm_out.img[i] = (int)round(filtered[i]);
			}
			delete[] filtered;
			filtered = NULL;
		} else if ((Options.mode & MODE_OUTPUT_MULTIPLE_MOTIONS_AFFINE) != 0) {
			/* Computte and output Multiple Motion Affine Parameters by method of M.J.Black */
			pnm_int2double(&pnmd_in, &pnm_orig, 1.0, NULL);
			if (pnmdouble_isNULL(&pnmd_prev) != PNM_FALSE) {
				printf("* Skip Calculate Multiple Motions by Affine while there is NOT any previous frame\n");
			} else {
				printf("* Compute Multiple Motions Affine Parameters by method of M.J.Black\n");
				MultipleMotion_AffineCoeff = MultipleMotion_Affine(pnmd_prev.imgd, pnmd_in.imgd, size_orig, Options.MultipleMotion_Param);
			}
		} else if ((Options.mode & MODE_OUTPUT_MULTIPLE_MOTIONS_OPTICALFLOW) != 0) {
			/* Computte and output Multiple Motion Optical Flow by method of M.J.Black */
			pnm_int2double(&pnmd_in, &pnm_orig, 1.0, NULL);
			if (pnmdouble_isNULL(&pnmd_prev) != PNM_FALSE) {
				printf("* Skip Calculate Multiple Motions while there is NOT any previous frame\n");
			} else {
				printf("* Compute Multiple Motions Optical Flow by method of M.J.Black\n");
				MultipleMotion_u = MultipleMotion_OpticalFlow(pnmd_prev.imgd, pnmd_in.imgd, size_orig, Options.MultipleMotion_Param);
			}
		} else if ((Options.mode & MODE_OUTPUT_OPTICALFLOW_AFFINE_PARAMETER) != 0) {
			/* Computte and output Affine vectors which represent Optical Flow by method of J.M.Odobez */
			pnm_int2double(&pnmd_in, &pnm_orig, 1.0, NULL);
			if (pnmdouble_isNULL(&pnmd_prev) != PNM_FALSE) {
				printf("* Skip Calculate Optical Flow while there is NOT any previous frame\n");
			} else {
				printf("* Compute Optical Flow Affine Parameters by method of J.M.Odobez\n");
				OpticalFlow_Affine = OpticalFlow_AffineParamet(pnmd_prev.imgd, pnmd_in.imgd, size_orig, Options.OpticalFlow_Param);
				if (OpticalFlow_Affine == NULL) {
					Error.Function("OpticalFlow_RMR");
					Error.Value("OpticalFlow_Affine");
					Error.FunctionFail();
					goto ExitError;
				}
			}
		} else {
			/* Scratch Detection */
			printf("* Detect Scratch like vertical lines\n");
			scratches = DetectScratch(&pnm_in, Options.s_med, Options.s_avg, FilterParam, DO_DETECTION);
			if (scratches == NULL) {
				Error.Function("DetectScratch");
				Error.Value("scratches");
				Error.FunctionFail();
				goto ExitError;
			}
			for (i = 0; i < size.height * size.width; i++) {
				pnm_in.img[i] = scratches[i];
			}

			if ((Options.mode & MODE_OUTPUT_BINARY_IMAGE) != 0) {
				// Output Scratch Map without Meaningful Alignments
				if (pnmcp(&pnm_out, &pnm_in) != PNM_FUNCTION_SUCCESS) {
					Error.Function("pnmcp");
					Error.Value("(pnm_in -> pnm_out)");
					Error.FunctionFail();
					goto ExitError;
				}
			} else {
				/* A Contrario Method : Meaningful Alignments */
				if (Initialize == 0) {
					Initialize = 1;
					try {
						Pr_table = new double[(maxMN + 1) * (maxMN + 1)];
					}
					catch (std::bad_alloc bad) {
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
							Pr_table[maxMN * k + L] = Pr(k, L, Options.p);
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
					if (k_list == NULL) {
						Error.Function("Calc_k_l");
						Error.Value("k_list");
						Error.FunctionFail();
						goto ExitError;
					}
				}

				printf("* Compute Direction Field\n");
				angles = DerivativeAngler(scratches, size);
				if (angles ==NULL) {
					Error.Function("Derivation");
					Error.Value("angles");
					Error.FunctionFail();
					goto ExitError;
				}
				printf("* Compute Segments and Maximal Meaningfulness\n");
				MaximalSegments = AlignedSegment_vertical(angles, size, k_list, l_min, Pr_table, &Num_Segments, Options.Max_Length, Options.Max_Output_Length);
				if (MaximalSegments == NULL) {
					Error.Function("AlignedSegment_vertical");
					Error.Value("MaximalSegments");
					Error.FunctionFail();
					goto ExitError;
				}
				printf("- Found (%d) Maximal Meaningful Segments\n", Num_Segments);
				if (Options.ExclusivePrinciple != 0) {
					printf("* Delete Redundant Segments by Exclusive Principle\n");
					EPSegments = ExclusivePrinciple(angles, size, k_list, Pr_table, MaximalSegments, &Num_Segments, Options.Exclusive_Max_Radius);
					if (EPSegments == NULL) {
						Error.Function("ExclusivePrinciple");
						Error.Value("EPSegments");
						Error.FunctionFail();
						goto ExitError;
					}
					printf("- Reduced to (%d) EP-Maximal Meaningful Segments\n", Num_Segments);
					delete[] MaximalSegments;
					MaximalSegments = EPSegments;
					EPSegments = NULL;
				}
				printf("* Plot The Segments");
				if (Options.Max_Output_Length > 0) {
					printf(" that satisfy (length < %d)", Options.Max_Output_Length);
				}
				printf("\n");
				segments = PlotSegment(MaximalSegments, Num_Segments, size, size_out, Options.PlotOptions & PLOT_NEGATE);
				if (segments == NULL) {
					Error.Function("PlotSegment");
					Error.Value("segments");
					Error.FunctionFail();
					goto ExitError;
				}
				if (Options.Superimpose != 0) {
					printf("* Superimpose plot image on original image\n");
					if ((Options.PlotOptions & PLOT_AS_RESAMPLE) != 0) {
						if (Superimposer(&pnm_out, &pnm_res, segments, size_out, Options.Superimpose, Options.PlotOptions & PLOT_NEGATE) != MEANINGFUL_SUCCESS) {
							Error.Function("Superimposer");
							Error.Value("pnm_out");
							Error.FunctionFail();
							goto ExitError;
						}
					} else {
						if (Superimposer(&pnm_out, &pnm_orig, segments, size_out, Options.Superimpose, Options.PlotOptions & PLOT_NEGATE) != MEANINGFUL_SUCCESS) {
							Error.Function("Superimposer");
							Error.Value("pnm_out");
							Error.FunctionFail();
							goto ExitError;
						}
					}
				} else {
					if (pnmnew(&pnm_out, PORTABLE_GRAYMAP_BINARY, size_out.width, size_out.height, pnm_orig.maxint) != PNM_FUNCTION_SUCCESS) {
						Error.Function("pnmnew");
						Error.Value("pnm_out");
						Error.Malloc();
						goto ExitError;
					}
					for (i = 0; i < size_out.height * size_out.width; i++) {
						pnm_out.img[i] = segments[i];
					}
				}

				/* X11 Plotting */
				ShowSegments_X11(pnm_orig.img, size_orig, size, pnm_orig.maxint, MaximalSegments, Num_Segments);
				/* /X11 Plotting */

				delete[] segments;
				segments = NULL;
				delete[] angles;
				angles = NULL;
				delete[] MaximalSegments;
				MaximalSegments = NULL;
			}
			delete[] scratches;
			scratches = NULL;
		}
		// Write
Write:
		if (strchr(OutputName, '%') == NULL) {
			sprintf(OutputNameNums, OutputName);
		} else {
			sprintf(OutputNameNums, OutputName, CurrentFileNum);
		}
		if ((Options.mode & MODE_OUTPUT_MULTIPLE_MOTIONS_AFFINE) != 0) {
			if (pnmdouble_isNULL(&pnmd_prev) == PNM_FALSE
			    && MultipleMotion_Affine_write(MultipleMotion_AffineCoeff, OutputNameNums) == MEANINGFUL_FAILURE) {
				Error.Function("MultipleMotion_Affine_write");
				Error.Value("MultipleMotions_Affine");
				Error.FunctionFail();
				goto ExitError;
			}
		} else if ((Options.mode & MODE_OUTPUT_MULTIPLE_MOTIONS_OPTICALFLOW) != 0) {
			if (pnmdouble_isNULL(&pnmd_prev) == PNM_FALSE
			    && MultipleMotion_write(MultipleMotion_u, size_orig, OutputNameNums) == MEANINGFUL_FAILURE) {
				Error.Function("MultipleMotion_write");
				Error.Value("MultipleMotions_u");
				Error.FunctionFail();
				goto ExitError;
			}
		} else if ((Options.mode & MODE_OUTPUT_OPTICALFLOW_AFFINE_PARAMETER) != 0) {
			if (pnmdouble_isNULL(&pnmd_prev) == PNM_FALSE
			    && OpticalFlow_write(OpticalFlow_Affine, size_orig, OutputNameNums) == MEANINGFUL_FAILURE) {
				Error.Function("OpticalFlow_write");
				Error.Value("OpticalFlow_Affine");
				Error.FunctionFail();
				goto ExitError;
			}
		} else {
			if (pnmwrite(&pnm_out, OutputNameNums) == PNM_FUNCTION_ERROR) {
				Error.Function("pnmwrite");
				Error.File(OutputNameNums);
				Error.FileWrite();
				goto ExitError;
			}
		}
		delete[] MultipleMotion_u;
		MultipleMotion_u = NULL;
		delete[] OpticalFlow_Affine;
		OpticalFlow_Affine = NULL;
		pnmdouble_free(&pnmd_prev);
		if (pnmdouble_isNULL(&pnmd_in) == 0) {
			pnmdouble_cp(&pnmd_prev, &pnmd_in);
		}
		pnmdouble_free(&pnmd_in);
		pnmfree(&pnm_out);
		pnmfree(&pnm_in);
		pnmfree(&pnm_res);
		pnmfree(&pnm_orig);
	}
	pnmdouble_free(&pnmd_prev);
	pnmdouble_free(&pnmd_in);
	pnmfree(&pnm_out);
	pnmfree(&pnm_in);
	pnmfree(&pnm_res);
	pnmfree(&pnm_orig);
	delete[] k_list;
	delete[] Pr_table;
	delete[] OutputNameNums;
	delete[] InputNameNums;
	atan2_div_pi_table(0, 0, &size);
	return MEANINGFUL_SUCCESS;
// Exit Error
ExitError:
	delete[] segments;
	delete[] EPSegments;
	delete[] MaximalSegments;
	delete[] angles;
	delete[] binary;
	delete[] scratches;
	delete[] k_list;
	delete[] Pr_table;
	delete[] OutputNameNums;
	delete[] InputNameNums;
	delete[] filtered;
	atan2_div_pi_table(0, 0, &size);
	pnmfree(&pnm_out);
	pnmfree(&pnm_in);
	pnmfree(&pnm_res);
	pnmfree(&pnm_orig);
	pnmdouble_free(&pnmd_in);
	pnmdouble_free(&pnmd_out);
	return MEANINGFUL_FAILURE;
}

