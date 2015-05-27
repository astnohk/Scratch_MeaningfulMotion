#include "Scratch_MeaningfulA.h"



SEGMENT*
ExclusivePrinciple(double *angles, SIZE size, int *k_list, double *Pr_table, SEGMENT *MaximalSegments, int *Num_Segments, double Exclusive_max_radius)
{
	char *ErrorFunctionName = "";
	char *ErrorPointerName = "";
	int *IndexMap = NULL;
	LINEPOLE *Lines = NULL;
	double *EPSegments_Pr = NULL;
	int Num_EPSegments = 0;
	SEGMENT *MaxEPSegments = NULL;
	int maxMN = (size.height > size.width ? size.height : size.width);
	int line_index;
	double Pr_min;
	double aligned_angle;
	double d, d_triangle, d_max;
	int m, n;
	int x, y;
	int t;
	int x_t, y_t;
	double dx, dy;
	int k, L;
	int n_seg;
	PNM pnm = PNM_NULL;
	int progress_count;
	int present_count;

	if ((IndexMap = (int *)calloc((size_t)(size.height * size.width), sizeof(int))) == NULL) {
		ErrorFunctionName = "calloc";
		ErrorPointerName = "IndexMap";
		goto ErrorMalloc;
	}
	if ((Lines = (LINEPOLE *)calloc((size_t)(*Num_Segments), sizeof(LINEPOLE))) == NULL) {
		ErrorFunctionName = "calloc";
		ErrorPointerName = "Lines";
		goto ErrorMalloc;
	}

	/* Convert the segments Cartesian coordinates to Polar coordinates Expression */
	for (n_seg = 0; n_seg < (*Num_Segments); n_seg++) {
		Lines[n_seg].theta = M_PI * atan2_div_pi_table(MaximalSegments[n_seg].n - MaximalSegments[n_seg].x, MaximalSegments[n_seg].y - MaximalSegments[n_seg].m, NULL);
		if (Lines[n_seg].theta >= M_PI) {
			Lines[n_seg].theta -= M_PI;
		} else if (Lines[n_seg].theta < 0.0) {
			Lines[n_seg].theta += M_PI;
		}
		Lines[n_seg].cos = cos(Lines[n_seg].theta);
		Lines[n_seg].sin = sin(Lines[n_seg].theta);
		Lines[n_seg].r = MaximalSegments[n_seg].x * Lines[n_seg].cos + MaximalSegments[n_seg].y * Lines[n_seg].sin;
	}
	printf("* Add all pixels to the Segments Exclusively :\n  0%% |%s\x1b[1A\n", Progress_End);
	/* Select the segments each Pixel exclusively belongs to */
	progress_count = 0;
	present_count = 0;
#pragma omp parallel for schedule(dynamic) private(y, L, d, d_triangle, d_max, Pr_min, n_seg, line_index)
	for (x = 0; x < size.width; x++) {
		for (y = 0; y < size.height; y++) {
			Pr_min = 1.0; // Reset
			line_index = -1;
			for (n_seg = 0; n_seg < (*Num_Segments); n_seg++) {
				d = fabs(Lines[n_seg].r - (x * Lines[n_seg].cos + y * Lines[n_seg].sin)); // Calc distance
				d_triangle =
				    sqrt((double)POW2(x - MaximalSegments[n_seg].x) + (double)POW2(y - MaximalSegments[n_seg].y))
				    + sqrt((double)POW2(x - MaximalSegments[n_seg].n) + (double)POW2(y - MaximalSegments[n_seg].m));
				d_max = d + sqrt(
				    (double)POW2(MaximalSegments[n_seg].x - MaximalSegments[n_seg].n)
				    + (double)POW2(MaximalSegments[n_seg].y - MaximalSegments[n_seg].m)
				    + (double)d * d);
				if ((d < Exclusive_max_radius)
				    && (d_triangle <= d_max)
				    && (MaximalSegments[n_seg].Pr < Pr_min)) {
					line_index = n_seg;
					Pr_min = MaximalSegments[n_seg].Pr;
				}
			}
			IndexMap[size.width * y + x] = line_index;
		}
#pragma omp critical
		{
			present_count++;
			if (NUM_PROGRESS * present_count / size.width > progress_count) {
				progress_count = NUM_PROGRESS * (present_count - 1) / size.width; // Take account of Overflow
				printf("\r%3d%% |%s#\x1b[1A\n", 100 * present_count / size.width, Progress[progress_count]);
			}
		}
	}
	printf("\nComplete!\n");
	free(Lines);
	Lines = NULL;
	if ((EPSegments_Pr = (double *)calloc((size_t)(*Num_Segments), sizeof(double))) == NULL) {
		ErrorFunctionName = "calloc";
		ErrorPointerName = "EPSegments_Pr";
		goto ErrorMalloc;
	}
	progress_count = 0;
	present_count = 0;
	printf("* Delete Redundant Segments by Exclusive Principle :\n  0%% |%s\x1b[1A\n", Progress_End);
#pragma omp parallel for schedule(dynamic) private(m, n, x, y, k, L, dx, dy, d, d_triangle, aligned_angle, t, x_t, y_t) reduction(+:Num_EPSegments)
	for (n_seg = 0; n_seg < (*Num_Segments); n_seg++) {
		/* Re-meaningful segments */
		n = MaximalSegments[n_seg].n;
		m = MaximalSegments[n_seg].m;
		x = MaximalSegments[n_seg].x;
		y = MaximalSegments[n_seg].y;
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
		k = 0;
		for (t = 0; t < L; t++) {
			x_t = (int)round(dx * t + n);
			y_t = (int)round(dy * t + m);
			if (x_t < 0 || size.width <= x_t || y_t < 0 || size.height <= y_t) {
				break;
			}
			if (IndexMap[size.width * y_t + x_t] == n_seg) {
				if (fabs(angles[size.width * y_t + x_t] - aligned_angle) <= DIR_PROBABILITY
				    || fabs(angles[size.width * y_t + x_t] - ANGLE_MAX - aligned_angle) <= DIR_PROBABILITY
				    || fabs(angles[size.width * y_t + x_t] + ANGLE_MAX - aligned_angle) <= DIR_PROBABILITY) {
					k++;
				}
			}
		}
#pragma omp critical
		if (k >= k_list[L]) {
			Num_EPSegments++;
			EPSegments_Pr[n_seg] = Pr_table[maxMN * k + L];
		}
#pragma omp critical
		{
			present_count++;
			if (NUM_PROGRESS * present_count / (*Num_Segments) > progress_count) {
				progress_count = NUM_PROGRESS * (present_count - 1) / (*Num_Segments); // Take account of Overflow
				printf("\r%3d%% |%s#\x1b[1A\n", 100 * present_count / (*Num_Segments), Progress[progress_count]);
			}
		}
	}
	printf("\nComplete!\n");
	if (pnmnew(&pnm, PORTABLE_GRAYMAP_ASCII, (unsigned int)size.width, (unsigned int)size.height, (unsigned int)(*Num_Segments - 1)) != PNM_FUNCTION_SUCCESS) {
		ErrorFunctionName = "pnmnew";
		ErrorPointerName = "pnm";
		goto ErrorMalloc;
	}
	for (x = 0; x < size.height * size.width; x++) {
		pnm.img[x] = IndexMap[x];
	}
	if (pnmwrite(&pnm, "IndexMap.pgm") != PNM_FUNCTION_SUCCESS) {
		fprintf(stderr, "*** ExclusivePrinciple error - CanNOT write out the IndexMap to \"IndexMap.pgm\" ***\n");
		// Do NOT EXIT because it is NOT FATAL ERROR
	}
	pnmfree(&pnm);
	free(IndexMap);
	IndexMap = NULL;
	if ((MaxEPSegments = (SEGMENT *)calloc((size_t)(Num_EPSegments), sizeof(SEGMENT))) == NULL) {
		ErrorFunctionName = "calloc";
		ErrorPointerName = "MaxEPSegments";
		goto ErrorMalloc;
	}
	for (n_seg = k = 0; n_seg < (*Num_Segments); n_seg++) {
		if (EPSegments_Pr[n_seg] > 0.0) {
			MaxEPSegments[k].n = MaximalSegments[n_seg].n;
			MaxEPSegments[k].m = MaximalSegments[n_seg].m;
			MaxEPSegments[k].x = MaximalSegments[n_seg].x;
			MaxEPSegments[k].y = MaximalSegments[n_seg].y;
			MaxEPSegments[k].Pr = EPSegments_Pr[n_seg];
			k++;
			if (k >= Num_EPSegments) {
				break;
			}
		}
	}
	*Num_Segments = Num_EPSegments;
	free(EPSegments_Pr);
	EPSegments_Pr = NULL;
	return MaxEPSegments;
// Errors
ErrorMalloc:
	fprintf(stderr, "*** ExclusivePrinciple error - Cannot allocate memory for (*%s) by %s() ***\n", ErrorPointerName, ErrorFunctionName);
	free(EPSegments_Pr);
	EPSegments_Pr = NULL;
	free(IndexMap);
	IndexMap = NULL;
	free(Lines);
	Lines = NULL;
	free(MaxEPSegments);
	MaxEPSegments = NULL;
	return NULL;
}

