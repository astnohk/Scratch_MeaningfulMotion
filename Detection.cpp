#include "Scratch_MeaningfulA.h"

#define DEBUG_FILTER



double*
DetectScratch(PNM *pnm, double s_med, double s_avg, FILTER_PARAM FilterParam, int Do_Detection)
{
	char *ErrorFunctionName = "";
	char *ErrorValueName = "";

	double *img = NULL;
	double *scratches = NULL;
	double *Ig = NULL;
	SIZE size = SIZE_ZERO;
	int x, y, m, n;
	double Im, Il, Ir;
#if defined(DEBUG_FILTER)
	PNM pnm_out = PNM_NULL;
#endif

	size.width = (int)pnm->width;
	size.height = (int)pnm->height;
	if ((img = (double *)calloc((size_t)size.height * size.width, sizeof(double))) == NULL) {
		ErrorFunctionName = "calloc";
		ErrorValueName = "scratches";
		goto ErrorMalloc;
	}
	for (m = 0; m < size.height * size.width; m++) {
		img[m] = (double)pnm->img[m];
	}

	switch (FilterParam.type) {
		case FILTER_ID_EPSILON: /* Epsilon Filter */
			printf("* Epsilon filtering on input\n");
			if ((Ig = EpsilonFilter(img, size, FilterParam)) == NULL) {
				ErrorFunctionName = "EpsilonFilter";
				ErrorValueName = "*Ig";
				goto ErrorFunctionFailed;
			}
			break;
		case FILTER_ID_GAUSSIAN: /* Gaussian Filter */
			printf("* Gaussian filtering on input\n");
			if ((Ig = Gaussian(img, size, FilterParam)) == NULL) {
				ErrorFunctionName = "Gaussian";
				ErrorValueName = "*Ig";
				goto ErrorFunctionFailed;
			}
			break;
		default:
			printf("* Do NOT any filtering\n");
	}
#if defined(DEBUG_FILTER)
	if (pnmnew(&pnm_out, PORTABLE_GRAYMAP_BINARY, size.width, size.height, pnm->maxint) != PNM_FUNCTION_SUCCESS) {
		ErrorFunctionName = "pnmnew";
		ErrorValueName = "pnm_out";
		goto ErrorMalloc;
	}
	for (m = 0; m < size.height * size.width; m++)
		pnm_out.img[m] = (int)Ig[m];
	pnmwrite(&pnm_out, "filtered.pgm");
	pnmfree(&pnm_out);
	pnm_out = PNM_NULL;
#endif

	if (Do_Detection == 0) {
		scratches = Ig;
		Ig = NULL;
	} else {
		if ((scratches = (double *)calloc((size_t)size.height * size.width, sizeof(double))) == NULL) {
			ErrorFunctionName = "calloc";
			ErrorValueName = "scratches";
			goto ErrorMalloc;
		}
#pragma omp parallel for private(x, m, n, Im, Il, Ir)
		for (y = 0; y < size.height; y++) {
			for (x = 0; x < size.width; x++) {
				Im = HorizontalMedian(Ig, size.width, x, y, MEAN_WIDTH);
				if (fabs(Ig[size.width * y + x] - Im) >= s_med) {
					Il = m = 0;
					for (n = ((x - AVE_FAR) < 0) ? 0 : (x - AVE_FAR); (n < size.width) && (n < (x - SCRATCH_WIDTH / 2)); n++) {
						Il += Ig[size.width * y + n];
						m++;
					}
					Il /= (double)m;
					Ir = m = 0;
					for (n = x + SCRATCH_WIDTH / 2 + 1; n < size.width && n <= x + AVE_FAR; n++) {
						Ir += Ig[size.width * y + n];
						m++;
					}
					Ir /= (double)m;
					if (fabs(Il - Ir) <= s_avg) {
						scratches[size.width * y + x] = PLOT_INTENSITY_MAX;
					}
				}
			}
		}
	}
	free(Ig);
	Ig = NULL;
	free(img);
	img = NULL;
	return scratches;
// Errors
ErrorMalloc:
	fprintf(stderr, "*** DetectScratch error - Cannot allocate memory for (*%s) by %s() ***\n", ErrorValueName, ErrorFunctionName);
	goto ErrorReturn;
ErrorFunctionFailed:
	fprintf(stderr, "*** DetectScratch() error - %s() exited with FAILURE signal ***\n", ErrorFunctionName);
ErrorReturn:
	free(Ig);
	Ig = NULL;
	free(scratches);
	scratches = NULL;
	free(img);
	img = NULL;
	return NULL;
}


SEGMENT*
AlignedSegment_vertical(double *angles, SIZE size, int *k_list, int l_min, double *Pr_table, int *Num_Segments, int Max_Length, int Max_Output_Length)
{
	char *ErrorFunctionName = "";
	char *ErrorValueName = "";
	char *ErrorDescription = "";
	double rad_offset = M_PI * (0.5 - 0.5 / (double)DIV_ANGLE_VERTICAL);
	double tan_list[DIV_ANGLE];
	int r, m, n, x, y;
	int L;
	double dx, dy;
	SEGMENT *coord_array = NULL;
	SEGMENTS_LIST coord_list = SEGMENTS_LIST_NULL;
	SEGMENTS_LIST *plist = NULL;
	SEGMENTS_LIST *plist_del = NULL;
	double progress;
	int current_count;

	printf("- Search only the segments that satisfy\n\tlength > %d", l_min);
	if (Max_Length > 0) {
		printf(" AND length <= %d", Max_Length);
	}
	printf("\n");

	for (r = 0; r < DIV_ANGLE; r++) {
		if (r == DIV_ANGLE / 2) {
			tan_list[r] = 2.0 * (size.height > size.width ? size.height : size.width);
		} else {
			tan_list[r] = tan((M_PI / (double)DIV_ANGLE_VERTICAL) * r / (double)DIV_ANGLE + rad_offset);
		}
	}
	progress = .0;
	current_count = 0;
	printf("* Search segments starts from Upper or Bottom edge :\n   0.0%% |%s\x1b[1A\n", Progress_End);
#pragma omp parallel for schedule(dynamic) private(r, x, y, dx, dy)
	for (n = 0; n < size.width; n++) {
		for (r = 0; r < DIV_ANGLE; r++) {
			// Upper side to Other 3 sides (n, 0) -> (x, y)
			dx = n + round((size.height - 1) / tan_list[r]);
			x = (dx >= 0.0) ? (dx < (double)size.width) ? (int)dx : size.width - 1 : 0;
			if (tan_list[r] >= 0.0) {
				dy = round((size.width - 1 - n) * tan_list[r]);
			} else {
				dy = round(-n * tan_list[r]);
			}
			y = (dy >= 0.0) ? (dy < (double)size.height) ? (int)dy : size.height - 1 : 0;
			if (AlignedCheck(angles, size, k_list, Pr_table, &coord_list, l_min, 0, n, x, y, Max_Length, Max_Output_Length)
			    != MEANINGFUL_SUCCESS) {
				ErrorDescription = "Occured at Upper side to Other 3 sides";
				n = size.width;
				break;
			}
			// Bottom side to Other 3 sides (n, size.height) -> (x, y)
			dx = n + round(-(size.height - 1) / tan_list[r]);
			x = (dx >= 0.0) ? (dx < (double)size.width) ? (int)dx : size.width - 1 : 0;
			if (tan_list[r] >= 0.0) {
				dy = size.height - 1 + round(-n * tan_list[r]);
			} else {
				dy = size.height - 1 + round((size.width - 1 - n) * tan_list[r]);
			}
			y = (dy >= 0.0) ? (dy < (double)size.height) ? (int)dy : size.height - 1 : 0;
			if (AlignedCheck(angles, size, k_list, Pr_table, &coord_list, l_min, size.height - 1, n, x, y, Max_Length, Max_Output_Length)
			    != MEANINGFUL_SUCCESS) {
				ErrorDescription = "Occured at Bottom side to Other 3 sides";
				n = size.width;
				break;
			}
#pragma omp critical
			{
				current_count++;
				if (round((double)current_count / (DIV_ANGLE * size.width) * 1000.0) > progress) {
					progress = round((double)current_count / (DIV_ANGLE * size.width) * 1000.0); // Take account of Overflow
					printf("\r %5.1f%% |%s#\x1b[1A\n", progress * 0.1, Progress[NUM_PROGRESS * current_count / (DIV_ANGLE * size.width + 1)]);
				}
			}
		}
	}
	if (ErrorDescription[0] !='\0') {
		goto ErrorOthers;
	}
	printf("\nComplete!\n");
	plist = coord_list.next;
	L = 0;
	while (plist) {
		L++;
		plist = plist->next;
	}
	if ((coord_array = (SEGMENT *)calloc((size_t)L, sizeof(SEGMENT))) == NULL) {
		ErrorFunctionName = "calloc";
		ErrorValueName = "coord_array";
		goto ErrorMalloc;
	}
	plist = coord_list.next;
	for (m = 0; m < L; m++) {
		coord_array[m] = (SEGMENT){plist->n, plist->m, plist->x, plist->y, plist->Pr};
		plist_del = plist;
		plist = plist->next;
		free(plist_del);
	}
	plist_del = NULL;
	*Num_Segments = L;
	return coord_array;
// Errors
ErrorMalloc:
	fprintf(stderr, "*** AlignedSegment_vertical() error - Cannot allocate memory for (*%s) by %s() ***\n", ErrorValueName, ErrorFunctionName);
	goto ErrorReturn;
ErrorOthers:
	fprintf(stderr, "*** AlignedSegment_vertical() error - %s ***\n", ErrorDescription);
ErrorReturn:
	free(coord_array);
	coord_array = NULL;
	list_free(coord_list.next);
	coord_list.next = NULL;
	return NULL;
}


int
AlignedCheck(double *angles, SIZE size, int *k_list, double *Pr_table, SEGMENTS_LIST *list_start, int l_min, int m, int n, int x, int y, int Max_Length, int Max_Output_Length)
{
	char *ErrorFunctionName = "";
	char *ErrorValueName = "";

	FRAGMENT List_Segments = {0, 0, 0.0, NULL};
	FRAGMENT *psegments = NULL;
	FRAGMENT *seg_ref_prev = NULL;
	FRAGMENT *seg_ref = NULL;
	FRAGMENT *seg_current_prev = NULL;
	FRAGMENT *seg_current = NULL;
	FRAGMENT *pseg_del = NULL;
	SEGMENTS_LIST *plist = NULL;
	double Pr_max, Pr_new;
	double aligned_angle;
	int L;
	double dx, dy;
	int t_start, t_end, t_end_max;
	int k, t;
	int tmpx, tmpy;

	if (abs(x - n) > abs(y - m)) {
		L = abs(x - n) + 1;
	} else {
		L = abs(y - m) + 1;
	}
	aligned_angle = atan2_div_pi_table(y - m, x - n, NULL);
	if (aligned_angle < 0.0) {
		aligned_angle += ANGLE_MAX;
	}
	dx = (x - n) / (double)(L - 1.0);
	dy = (y - m) / (double)(L - 1.0);
	psegments = NULL;
	for (t_start = 0; t_start <= L - l_min; t_start++) {
		tmpx = (int)round(dx * t_start + n);
		tmpx = tmpx >= 0 ? (tmpx < size.width ? tmpx : size.width - 1) : 0;
		tmpy = (int)round(dy * t_start + m);
		tmpy = tmpy >= 0 ? (tmpy < size.height ? tmpy : size.height - 1) : 0;
		if (fabs(angles[size.width * tmpy + tmpx] - aligned_angle) <= DIR_PROBABILITY
		    || fabs(angles[size.width * tmpy + tmpx] - ANGLE_MAX - aligned_angle) <= DIR_PROBABILITY
		    || fabs(angles[size.width * tmpy + tmpx] + ANGLE_MAX - aligned_angle) <= DIR_PROBABILITY) {
			Pr_max = 1;
			t_end_max = 0;
			for (t_end = (l_min > 1) ? t_start + l_min - 1 : t_start + 1; t_end < L; t_end++) {
				if (Max_Length > 0) {
					if (t_end_max > 0 && t_end_max - t_start + 1 <= Max_Length && t_end - t_start + 1 > Max_Length) {
						// The length of the segment EXCEEDS the length limit
						if ((psegments = (FRAGMENT *)malloc(sizeof(FRAGMENT))) == NULL) {
							ErrorFunctionName = "malloc";
							ErrorValueName = "psegments";
							goto ErrorMalloc;
						}
						*psegments = (FRAGMENT){t_start, t_end_max, Pr_max, List_Segments.next};
						List_Segments.next = psegments;
						// Reset "t_end_max" but keep "Pr_max" to search the Longer Segments as if the process not limited by Length.
						t_end_max = 0;
					}
				}
				tmpx = (int)round(dx * t_end + n);
				tmpx = tmpx >= 0 ? (tmpx < size.width ? tmpx : size.width - 1) : 0;
				tmpy = (int)round(dy * t_end + m);
				tmpy = tmpy >= 0 ? (tmpy < size.height ? tmpy : size.height - 1) : 0;
				if (fabs(angles[size.width * tmpy + tmpx] - aligned_angle) <= DIR_PROBABILITY
				    || fabs(angles[size.width * tmpy + tmpx] - ANGLE_MAX - aligned_angle) <= DIR_PROBABILITY
				    || fabs(angles[size.width * tmpy + tmpx] + ANGLE_MAX - aligned_angle) <= DIR_PROBABILITY) {
					k = 2;
					// Check if the segment which has aligned point on both ends is epsilon-Meaningful
					for (t = t_start + 1; t <= t_end - 1; t++) {
						tmpx = (int)round(dx * t + n);
						tmpx = tmpx >= 0 ? (tmpx < size.width ? tmpx : size.width - 1) : 0;
						tmpy = (int)round(dy * t + m);
						tmpy = tmpy >= 0 ? (tmpy < size.height ? tmpy : size.height - 1) : 0;
						if (fabs(angles[size.width * tmpy + tmpx] - aligned_angle) <= DIR_PROBABILITY
						    || fabs(angles[size.width * tmpy + tmpx] - ANGLE_MAX - aligned_angle) <= DIR_PROBABILITY
						    || fabs(angles[size.width * tmpy + tmpx] + ANGLE_MAX - aligned_angle) <= DIR_PROBABILITY) {
							k++;
						}
					}
					if (k >= k_list[t_end - t_start + 1]) {
						Pr_new = Pr_table[(size.height > size.width ? size.height : size.width) * k + t_end - t_start + 1];
						if (Pr_new <= Pr_max) {
							Pr_max = Pr_new;
							t_end_max = t_end;
						}
					}
				}
			}
			if (t_end_max > 0) {
				if ((psegments = (FRAGMENT *)malloc(sizeof(FRAGMENT))) == NULL) {
					ErrorFunctionName = "malloc";
					ErrorValueName = "psegments";
					goto ErrorMalloc;
				}
				*psegments = (FRAGMENT){t_start, t_end_max, Pr_max, List_Segments.next};
				List_Segments.next = psegments;
			}
		}
	}
	// Maximal Meaningfulness
	pseg_del = NULL;
	seg_ref_prev = &List_Segments;
	seg_ref = List_Segments.next;
	while (seg_ref != NULL) {
		seg_current_prev = &List_Segments;
		seg_current = List_Segments.next;
		while (seg_ref != NULL && seg_current != NULL) {
			if (seg_ref == seg_current) {
				seg_current_prev = seg_current;
				seg_current = seg_current->next;
				continue;
			}
			if (seg_ref->start <= seg_current->start && seg_current->end <= seg_ref->end) {
				if (seg_ref->Pr <= seg_current->Pr) {
					pseg_del = seg_current;
					seg_current_prev->next = seg_current = pseg_del->next;
					if (pseg_del == seg_ref_prev) {
						seg_ref_prev = seg_current_prev;
					}
				} else {
					pseg_del = seg_ref;
					seg_ref_prev->next = seg_ref = pseg_del->next;
					if (pseg_del == seg_current_prev) {
						seg_current_prev = seg_ref_prev;
					}
				}
				free(pseg_del);
				pseg_del = NULL;
			} else if (seg_current->start <= seg_ref->start && seg_ref->end <= seg_current->end) {
				if (seg_current->Pr <= seg_ref->Pr) {
					pseg_del = seg_ref;
					seg_ref_prev->next = seg_ref = pseg_del->next;
					if (pseg_del == seg_current_prev) {
						seg_current_prev = seg_ref_prev;
					}
				} else {
					pseg_del = seg_current;
					seg_current_prev->next = seg_current = pseg_del->next;
					if (pseg_del == seg_ref_prev) {
						seg_ref_prev = seg_current_prev;
					}
				}
				free(pseg_del);
				pseg_del = NULL;
			} else {
				seg_current_prev = seg_current;
				seg_current = seg_current->next;
			}
		}
		seg_current_prev = NULL;
		seg_current = NULL;
		if (seg_ref != NULL) {
			seg_ref_prev = seg_ref;
			seg_ref = seg_ref->next;
		}
	}
	seg_ref_prev = NULL;
	seg_ref = NULL;
	// Write out to coord_list
	psegments = List_Segments.next;
	while (psegments != NULL) {
		// If Max_Output_Length <= 0 then Do NOT Limit The Length of Segments
		if (Max_Output_Length <= 0
		    || ((Max_Output_Length > 0) && (psegments->end - psegments->start + 1) <= Max_Output_Length)) {
			if ((plist = (SEGMENTS_LIST *)malloc(sizeof(SEGMENTS_LIST))) == NULL) {
				ErrorFunctionName = "malloc";
				ErrorValueName = "plist";
				goto ErrorMalloc;
			}
#pragma omp critical
			{
				*plist = (SEGMENTS_LIST){(int)round(n + dx * psegments->start), (int)round(m + dy * psegments->start), (int)round(n + dx * psegments->end), (int)round(m + dy * psegments->end), psegments->Pr, list_start->next};
				list_start->next = plist;
			}
		}
		pseg_del = psegments;
		psegments = psegments->next;
		List_Segments.next = psegments;
		free(pseg_del);
		pseg_del = NULL;
	}
	return MEANINGFUL_SUCCESS;
// Errors
ErrorMalloc:
	fprintf(stderr, "*** AlignedCheck() error - Cannot allocate memory for (*%s) by %s() ***\n", ErrorValueName, ErrorFunctionName);
	segments_free(List_Segments.next);
	List_Segments.next = NULL;
	return MEANINGFUL_FAILURE;
}

