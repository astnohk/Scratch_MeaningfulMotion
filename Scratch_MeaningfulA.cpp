#include "Scratch_MeaningfulA.h"
#include "Affine_MultipleMotion.h"
#include "OpticalFlow_MultipleMotion.h"
#include "OpticalFlow_AffineParamet.h"




int
Scratch_MeaningfulA(char *OutputName, char *InputName, unsigned int OutputNameLength, unsigned int InputNameLength, int Start, int End, OPTIONS Options, FILTER_PARAM FilterParam)
{
	char *ErrorFunctionName = "";
	char *ErrorValueName = "";
	char *ErrorFileName = "";
	char *ErrorDescription = "";
	char Bars[] = "------------------------------------------------";

	char *FilterNames[] = {"undefined", "Epsilon", "Gaussian"};
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
	VECTOR_AFFINE MultipleMotions_Affine;
	VECTOR_2D *MultipleMotions_u = NULL;
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
	SIZE size = SIZE_ZERO;
	SIZE size_prev = SIZE_ZERO;
	SIZE size_orig = SIZE_ZERO;
	SIZE size_res = SIZE_ZERO;
	SIZE size_out = SIZE_ZERO;
	int maxMN = 0;
	int l_min = 1;
	int i, k, L;
	double progress;
	int count;

	if ((InputNameNums = (char *)calloc((size_t)InputNameLength + 1u, sizeof(char))) == NULL) {
		ErrorFunctionName = "calloc";
		ErrorValueName = "InputNameNums";
		goto ErrorMalloc;
	}
	if ((OutputNameNums = (char *)calloc((size_t)OutputNameLength + 1u, sizeof(char))) == NULL) {
		ErrorFunctionName = "calloc";
		ErrorValueName = "OutputNameNums";
		goto ErrorMalloc;
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
			ErrorFunctionName = "pnmread";
			ErrorFileName = InputNameNums;
			goto ErrorFileRead;
		}
		// END Read

		if (((int)pnm_orig.height < 0) || (int)pnm_orig.width < 0) {
			ErrorDescription = "The size of the image is too large to processing (image size OVERFLOWED)";
			goto ErrorOthers;
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
			goto ImageSizeNotMatch;
		}
		printf("- The input image size is %dx%d\n- and bit depth is %d\n", size.width, size.height, (int)round(log2((double)pnm_orig.maxint)));
		if (size_res.width > 0 || size_res.height > 0) { /* Resample */
			size = size_res;
			if (pnm_int2double(&pnmd_in, &pnm_orig, 1.0, NULL) != PNM_FUNCTION_SUCCESS) {
				ErrorFunctionName = "pnm_int2double";
				goto ErrorConvert;
			}
			if (pnm_resize(&pnmd_out, &pnmd_in, size_res.width, size_res.height, Options.ResampleMethod) != PNM_FUNCTION_SUCCESS) {
				ErrorFunctionName = "pnm_resize";
				ErrorValueName = "(pnmd_in -> pnmd_out)";
				goto ErrorFunctionFailed;
			}
			pnmdouble_free(&pnmd_in);
			if (pnm_double2int(&pnm_res, &pnmd_out, 1.0, "round", NULL) != PNM_FUNCTION_SUCCESS) {
				ErrorFunctionName = "pnm_double2int";
				ErrorValueName = "(pnmd_out -> pnm_res)";
				goto ErrorFunctionFailed;
			}
			if (pnm_double2int(&pnm_in, &pnmd_out, 1.0, "round", NULL) != PNM_FUNCTION_SUCCESS) {
				ErrorFunctionName = "pnm_double2int";
				ErrorValueName = "(pnmd_out -> pnm_in)";
				goto ErrorConvert;
			}
			pnmdouble_free(&pnmd_out);
		} else {
			if (pnmcp(&pnm_in, &pnm_orig) != PNM_FUNCTION_SUCCESS) {
				ErrorFunctionName = "pnmcp";
				ErrorValueName = "(pnm_orig -> pnm_in)";
				goto ErrorConvert;
			}
		}
		if ((Options.PlotOptions & PLOT_RESAMPLED_IMG_ONLY) != 0) {
			/* Just output only the resampled image */
			if (pnmcp(&pnm_out, &pnm_in) != PNM_FUNCTION_SUCCESS) {
				ErrorFunctionName = "pnmcp";
				ErrorValueName = "(pnm_in -> pnm_out)";
				goto ErrorFunctionFailed;
			}
			goto Write;
		}
		if ((pnm_in.desc % 3) == 0) { /* Convert to Grayscale */
			printf("- The input image is color data\n");
			printf("* Convert the image to grayscale before applying Meaningful Alignments.\n");
			printf("Convert...   ");
			if (pnm_int2double(&pnmd_in, &pnm_in, 1.0, NULL) != PNM_FUNCTION_SUCCESS) {
				ErrorFunctionName = "pnm_int2double";
				goto ErrorConvert;
			}
			if (pnm_RGB2Gray(&pnmd_out, &pnmd_in) != PNM_FUNCTION_SUCCESS) {
				ErrorFunctionName = "pnm_RGB2Gray";
				goto ErrorConvert;
			}
			pnmdouble_free(&pnmd_in);
			pnmfree(&pnm_in);
			if (pnm_double2int(&pnm_in, &pnmd_out, 1.0, "round", NULL) != PNM_FUNCTION_SUCCESS) {
				ErrorFunctionName = "pnm_double2int";
				goto ErrorConvert;
			}
			pnmdouble_free(&pnmd_out);
			printf("Finished\n\n");
		}

		/* Show Parameters */
		printf("\n      --- Parameters ---\n  %s\n", Bars);
		if (size_res.width > 0 || size_res.height > 0) {
			printf("  | Resample (%s)\n", Options.ResampleMethod);
			printf("  |   %ux%u -> %dx%d\n", pnm_orig.width, pnm_orig.width, size_res.width, size_res.height);
		}
		printf("  | filter type = %s\n", FilterNames[FilterParam.type < NUM_FILTER_TYPE ? FilterParam.type : 0]);
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
				ErrorFunctionName = "DetectScratch";
				ErrorValueName = "filtered";
				goto ErrorFunctionFailed;
			}
			printf("* Output Filtered Image\n");
			if (pnmnew(&pnm_out, PORTABLE_GRAYMAP_BINARY, size_res.width, size_res.height, pnm_in.maxint) != PNM_FUNCTION_SUCCESS) {
				ErrorFunctionName = "pnmnew";
				ErrorValueName = "pnm_out";
				goto ErrorMalloc;
			}
			for (i = 0; i < size_res.width * size_res.height; i++) {
				pnm_out.img[i] = (int)round(filtered[i]);
			}
			free(filtered);
			filtered = NULL;
		} else if ((Options.mode & MODE_OUTPUT_MULTIPLE_MOTIONS_AFFINE) != 0) {
			/* Computte and output Multiple Motion Affine Parameters by method of M.J.Black */
			pnm_int2double(&pnmd_in, &pnm_orig, 1.0, NULL);
			if (pnmdouble_isNULL(&pnmd_prev) != PNM_FALSE) {
				printf("* Skip Calculate Multiple Motions by Affine while there is NOT any previous frame\n");
			} else {
				printf("* Compute Multiple Motions Affine Parameters by method of M.J.Black\n");
				MultipleMotions_Affine = MultipleMotion_Affine(pnmd_prev.imgd, pnmd_in.imgd, size_orig, Options.MultipleMotion_Param);
			}
		} else if ((Options.mode & MODE_OUTPUT_MULTIPLE_MOTIONS_OPTICALFLOW) != 0) {
			/* Computte and output Multiple Motion Optical Flow by method of M.J.Black */
			pnm_int2double(&pnmd_in, &pnm_orig, 1.0, NULL);
			if (pnmdouble_isNULL(&pnmd_prev) != PNM_FALSE) {
				printf("* Skip Calculate Multiple Motions while there is NOT any previous frame\n");
			} else {
				printf("* Compute Multiple Motions Optical Flow by method of M.J.Black\n");
				MultipleMotions_u = MultipleMotion_OpticalFlow(pnmd_prev.imgd, pnmd_in.imgd, size_orig, Options.MultipleMotion_Param);
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
					ErrorFunctionName = "OpticalFlow_AffineParamet";
					ErrorValueName = "OpticalFlow_Affine";
					goto ErrorFunctionFailed;
				}
			}
		} else {
			/* Scratch Detection */
			printf("* Detect Scratch like vertical lines\n");
			scratches = DetectScratch(&pnm_in, Options.s_med, Options.s_avg, FilterParam, DO_DETECTION);
			if (scratches == NULL) {
				ErrorFunctionName = "DetectScratch";
				goto ErrorFunctionFailed;
			}
			for (i = 0; i < size.height * size.width; i++) {
				pnm_in.img[i] = scratches[i];
			}

			if ((Options.mode & MODE_OUTPUT_BINARY_IMAGE) != 0) {
				// Output Scratch Map without Meaningful Alignments
				if (pnmcp(&pnm_out, &pnm_in) != PNM_FUNCTION_SUCCESS) {
					ErrorFunctionName = "pnmcp";
					ErrorValueName = "(pnm_in -> pnm_out)";
					goto ErrorFunctionFailed;
				}
			} else {
				/* A Contrario Method : Meaningful Alignments */
				if (Initialize == 0) {
					Initialize = 1;
					if ((Pr_table = (double *)calloc((size_t)((maxMN + 1) * (maxMN + 1)), sizeof(double))) == NULL) {
						ErrorFunctionName = "calloc";
						ErrorValueName = "Pr_table";
						goto ErrorMalloc;
					}
					printf("* Calculate Pr(k, L) table :\n[L =     0]   0.0%% |%s\x1b[1A\n", Progress_End);
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
								printf("\r[L = %5d] %5.1f%% |%s#\x1b[1A\n", count, progress * 0.1, Progress[NUM_PROGRESS * count / (1 + maxMN)]);
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
						ErrorFunctionName = "Calc_k_l";
						goto ErrorFunctionFailed;
					}
				}

				printf("* Compute Direction Field\n");
				angles = DerivativeAngler(scratches, size);
				if (angles ==NULL) {
					ErrorFunctionName = "Derivation";
					goto ErrorFunctionFailed;
				}
				printf("* Compute Segments and Maximal Meaningfulness\n");
				MaximalSegments = AlignedSegment_vertical(angles, size, k_list, l_min, Pr_table, &Num_Segments, Options.Max_Length, Options.Max_Output_Length);
				if (MaximalSegments == NULL) {
					ErrorFunctionName = "AlignedSegment_vertical";
					goto ErrorFunctionFailed;
				}
				printf("- Found (%d) Maximal Meaningful Segments\n", Num_Segments);
				if (Options.ExclusivePrinciple != 0) {
					printf("* Delete Redundant Segments by Exclusive Principle\n");
					EPSegments = ExclusivePrinciple(angles, size, k_list, Pr_table, MaximalSegments, &Num_Segments, Options.Exclusive_Max_Radius);
					if (EPSegments == NULL) {
						ErrorFunctionName = "ExclusivePrinciple";
						goto ErrorFunctionFailed;
					}
					printf("- Reduced to (%d) EP-Maximal Meaningful Segments\n", Num_Segments);
					free(MaximalSegments);
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
					ErrorFunctionName = "PlotSegment";
					goto ErrorFunctionFailed;
				}
				if (Options.Superimpose != 0) {
					printf("* Superimpose plot image on original image\n");
					if ((Options.PlotOptions & PLOT_AS_RESAMPLE) != 0) {
						if (Superimposer(&pnm_out, &pnm_res, segments, size_out, Options.Superimpose, Options.PlotOptions & PLOT_NEGATE) != MEANINGFUL_SUCCESS) {
							ErrorFunctionName = "Superimposer";
							ErrorValueName = "pnm_out";
							goto ErrorFunctionFailed;
						}
					} else {
						if (Superimposer(&pnm_out, &pnm_orig, segments, size_out, Options.Superimpose, Options.PlotOptions & PLOT_NEGATE) != MEANINGFUL_SUCCESS) {
							ErrorFunctionName = "Superimposer";
							ErrorValueName = "pnm_out";
							goto ErrorFunctionFailed;
						}
					}
				} else {
					if (pnmnew(&pnm_out, PORTABLE_GRAYMAP_BINARY, size_out.width, size_out.height, pnm_orig.maxint) != PNM_FUNCTION_SUCCESS) {
						ErrorFunctionName = "pnmnew";
						ErrorValueName = "pnm_out";
						goto ErrorMalloc;
					}
					for (i = 0; i < size_out.height * size_out.width; i++) {
						pnm_out.img[i] = segments[i];
					}
				}

				/* X11 Plotting */
				ShowSegments_X11(pnm_orig.img, size_orig, size, pnm_orig.maxint, MaximalSegments, Num_Segments);
				/* /X11 Plotting */

				free(segments);
				segments = NULL;
				free(angles);
				angles = NULL;
				free(MaximalSegments);
				MaximalSegments = NULL;
			}
			free(scratches);
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
			    && MultipleMotion_Affine_write(MultipleMotions_Affine, OutputNameNums) == MEANINGFUL_FAILURE) {
				ErrorFunctionName = "MultipleMotion_Affine_write";
				ErrorValueName = "MultipleMotions_Affine";
				goto ErrorFunctionFailed;
			}
		} else if ((Options.mode & MODE_OUTPUT_MULTIPLE_MOTIONS_OPTICALFLOW) != 0) {
			if (pnmdouble_isNULL(&pnmd_prev) == PNM_FALSE
			    && MultipleMotion_write(MultipleMotions_u, size_orig, OutputNameNums) == MEANINGFUL_FAILURE) {
				ErrorFunctionName = "MultipleMotion_write";
				ErrorValueName = "MultipleMotions_u";
				goto ErrorFunctionFailed;
			}
		} else if ((Options.mode & MODE_OUTPUT_OPTICALFLOW_AFFINE_PARAMETER) != 0) {
			if (pnmdouble_isNULL(&pnmd_prev) == PNM_FALSE
			    && OpticalFlow_write(OpticalFlow_Affine, size_orig, OutputNameNums) == MEANINGFUL_FAILURE) {
				ErrorFunctionName = "OpticalFlow_write";
				ErrorValueName = "OpticalFlow_Affine";
				goto ErrorFunctionFailed;
			}
		} else {
			if (pnmwrite(&pnm_out, OutputNameNums) == PNM_FUNCTION_ERROR) {
				ErrorFunctionName = "pnmwrite";
				ErrorFileName = OutputNameNums;
				goto ErrorFileWrite;
			}
		}
		free(MultipleMotions_u);
		MultipleMotions_u = NULL;
		free(OpticalFlow_Affine);
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
	free(OpticalFlow_Affine);
	pnmdouble_free(&pnmd_prev);
	pnmdouble_free(&pnmd_in);
	pnmfree(&pnm_out);
	pnmfree(&pnm_in);
	pnmfree(&pnm_res);
	pnmfree(&pnm_orig);
	free(k_list);
	free(Pr_table);
	free(OutputNameNums);
	free(InputNameNums);
	size = SIZE_ZERO;
	atan2_div_pi_table(0, 0, &size);
	return MEANINGFUL_SUCCESS;
// Errors
ErrorMalloc:
	fprintf(stderr, "*** Scratch_MeaningfulA() error - Cannot allocate memory for (*%s) by %s() ***\n", ErrorValueName, ErrorFunctionName);
	goto ErrorReturn;
ImageSizeNotMatch:
	fprintf(stderr, "*** Scratch_MeaningfulA() error - The size of image (%u x %u) is varied from First Frame (%d x %d) ***\n", pnm_orig.width, pnm_orig.height, size.width, size.height);
	goto ErrorReturn;
ErrorConvert:
	fprintf(stderr, "*** Scratch_MeaningfulA() error - Failed to convert by %s() ***\n", ErrorFunctionName);
	goto ErrorReturn;
ErrorFileRead:
	fprintf(stderr, "*** Scratch_MeaningfulA() error - Failed to read the image from file \"%s\" by %s() ***\n", ErrorFileName, ErrorFunctionName);
	goto ErrorReturn;
ErrorFileWrite:
	fprintf(stderr, "*** Scratch_MeaningfulA() error - Failed to write the image to file \"%s\" by %s() ***\n", ErrorFileName, ErrorFunctionName);
	goto ErrorReturn;
ErrorFunctionFailed:
	fprintf(stderr, "*** Scratch_MeaningfulA() error - %s() exited with FAILURE signal (value : %s) ***\n", ErrorFunctionName, ErrorValueName);
	goto ErrorReturn;
ErrorOthers:
	fprintf(stderr, "*** Scratch_MeaningfulA() error - %s ***\n", ErrorDescription);
ErrorReturn:
	free(MultipleMotions_u);
	free(OpticalFlow_Affine);
	free(segments);
	free(EPSegments);
	free(MaximalSegments);
	free(angles);
	free(binary);
	free(scratches);
	free(k_list);
	free(Pr_table);
	free(OutputNameNums);
	free(InputNameNums);
	free(filtered);
	size = SIZE_ZERO;
	atan2_div_pi_table(0, 0, &size);
	pnmfree(&pnm_out);
	pnmfree(&pnm_in);
	pnmfree(&pnm_res);
	pnmfree(&pnm_orig);
	pnmdouble_free(&pnmd_in);
	pnmdouble_free(&pnmd_out);
	return MEANINGFUL_FAILURE;
}

