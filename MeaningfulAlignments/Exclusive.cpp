#include "../Scratch_MeaningfulMotion.h"



SEGMENT *
ExclusivePrinciple(ImgVector<double> *angles, int *k_list, ImgVector<double> *Pr_table, SEGMENT *MaximalSegments, int *Num_Segments, double Exclusive_max_radius)
{
	ERROR Error("ExclusivePrinciple");

	ImgVector<int> *IndexMap = nullptr;
	SEGMENT *MaxEPSegments = nullptr;
	SIZE size;

	PNM pnm;

	// Compute Exclusive Index Map
	size.width = angles->width();
	size.height = angles->height();
	IndexMap = ExclusiveIndexMap(size, MaximalSegments, Num_Segments, Exclusive_max_radius);
	if (IndexMap == nullptr) {
		Error.Function("ExclusiveIndexMap");
		Error.Value("IndexMap");
		Error.FunctionFail();
		goto ExitError;
	}
	printf("\nComplete!\n");
	// DEBUG : Output IndexMap
	pnm.copy(PORTABLE_GRAYMAP_ASCII, angles->width(), angles->height(), (*Num_Segments) - 1, IndexMap->data());
	if (pnm.write("IndexMap.pgm") == PNM_FUNCTION_ERROR) {
		fprintf(stderr, "*** ExclusivePrinciple error - CanNOT write out the IndexMap to \"IndexMap.pgm\" ***\n");
	}
	pnm.free();
	// Select the segments which satisfy The Maximal Exclusive Principle
	MaxEPSegments = ExclusiveSegments(IndexMap, angles, MaximalSegments, Num_Segments, k_list, Pr_table);
	if (MaxEPSegments == nullptr) {
		Error.Function("ExclusiveSegments");
		Error.Value("MaxEPSegments");
		Error.FunctionFail();
		goto ExitError;
	}
	printf("\nComplete!\n");
	delete IndexMap;
	IndexMap = nullptr;
	return MaxEPSegments;
// Error
ExitError:
	delete[] MaxEPSegments;
	delete IndexMap;
	return nullptr;
}


ImgVector<int> *
ExclusiveIndexMap(SIZE size, SEGMENT *MaximalSegments, int *Num_Segments, double Exclusive_max_radius)
{
	ERROR Error("ExclusiveIndexMap");

	ImgVector<int> *IndexMap = nullptr;

	LINEPOLE *Lines = nullptr;

	ATAN2_DIV_PI atan2_div_pi(size.width, size.height);

	int line_index;
	double Pr_min;
	double d, d_triangle, d_max;
	int n_seg;
	int x, y;

	int progress_count;
	int present_count;

	try {
		IndexMap = new ImgVector<int>(size.width, size.height);
	}
	catch (const std::bad_alloc &bad) {
		Error.Function("new");
		Error.Value("IndexMap");
		Error.Malloc();
		goto ExitError;
	}
	try {
		Lines = new LINEPOLE[(*Num_Segments)];
	}
	catch (const std::bad_alloc &bad) {
		Error.Function("new");
		Error.Value("Lines");
		Error.Malloc();
		goto ExitError;
	}

	// Convert the segments Cartesian coordinates to Polar coordinates Expression
	for (n_seg = 0; n_seg < (*Num_Segments); n_seg++) {
		Lines[n_seg].theta = M_PI * atan2_div_pi.val(
		    MaximalSegments[n_seg].n - MaximalSegments[n_seg].x,
		    MaximalSegments[n_seg].y - MaximalSegments[n_seg].m);
		if (Lines[n_seg].theta >= M_PI) {
			Lines[n_seg].theta -= M_PI;
		} else if (Lines[n_seg].theta < 0.0) {
			Lines[n_seg].theta += M_PI;
		}
		Lines[n_seg].cos = cos(Lines[n_seg].theta);
		Lines[n_seg].sin = sin(Lines[n_seg].theta);
		Lines[n_seg].r =
		    MaximalSegments[n_seg].x * Lines[n_seg].cos
		    + MaximalSegments[n_seg].y * Lines[n_seg].sin;
	}
	printf("* Add all pixels to the Segments Exclusively :\n  0%% |%s\x1b[1A\n", Progress_End.c_str());
	// Select the segments each Pixel exclusively belongs to
	progress_count = 0;
	present_count = 0;
#pragma omp parallel for schedule(dynamic) private(y, d, d_triangle, d_max, Pr_min, n_seg, line_index)
	for (x = 0; x < size.width; x++) {
		for (y = 0; y < size.height; y++) {
			Pr_min = 1.0; // Reset
			line_index = -1;
			for (n_seg = 0; n_seg < (*Num_Segments); n_seg++) {
				d = fabs(Lines[n_seg].r - (x * Lines[n_seg].cos + y * Lines[n_seg].sin)); // Calculate distance
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
			IndexMap->set(x, y, line_index);
		}
#pragma omp critical
		{
			present_count++;
			if (NUM_PROGRESS * present_count / size.width > progress_count) {
				progress_count = NUM_PROGRESS * (present_count - 1) / size.width; // Take account of Overflow
				printf("\r%3d%% |%s#\x1b[1A\n", 100 * present_count / size.width, Progress[progress_count].c_str());
			}
		}
	}
	delete[] Lines;
	Lines = nullptr;
	return IndexMap;
// Error
ExitError:
	delete IndexMap;
	delete[] Lines;
	return nullptr;
}


SEGMENT *
ExclusiveSegments(ImgVector<int> *IndexMap, ImgVector<double> *angles, SEGMENT *MaximalSegments, int *Num_Segments, int *k_list, ImgVector<double> *Pr_table)
{
	ERROR Error("ExclusiveSegments");

	SEGMENT *MaxEPSegments = nullptr;

	double *EPSegments_Pr = nullptr;
	int Num_EPSegments = 0;

	ATAN2_DIV_PI atan2_div_pi(angles->width(), angles->height());
	double aligned_angle;

	int n_seg;
	int m, n;
	int x, y, t;
	int x_t, y_t;
	double dx, dy;
	int k, L;
	
	int progress_count;
	int present_count;

	try {
		EPSegments_Pr = new double[(*Num_Segments)];
	}
	catch (const std::bad_alloc &bad) {
		Error.Function("new");
		Error.Value("EPSegments_Pr");
		Error.Malloc();
		goto ExitError;
	}
	progress_count = 0;
	present_count = 0;
	printf("* Delete Redundant Segments by Exclusive Principle :\n  0%% |%s\x1b[1A\n", Progress_End.c_str());
#pragma omp parallel for schedule(dynamic) private(m, n, x, y, t, k, L, dx, dy, x_t, y_t, aligned_angle) reduction(+:Num_EPSegments)
	for (n_seg = 0; n_seg < (*Num_Segments); n_seg++) {
		// Re-meaningful segments
		n = MaximalSegments[n_seg].n;
		m = MaximalSegments[n_seg].m;
		x = MaximalSegments[n_seg].x;
		y = MaximalSegments[n_seg].y;
		if (abs(x - n) > abs(y - m)) {
			L = abs(x - n) + 1;
		} else {
			L = abs(y - m) + 1;
		}
		aligned_angle = atan2_div_pi.val(y - m, x - n);
		if (aligned_angle < 0.0) {
			aligned_angle += ANGLE_MAX;
		}
		dx = (x - n) / (double)(L - 1.0);
		dy = (y - m) / (double)(L - 1.0);
		k = 0;
		for (t = 0; t < L; t++) {
			x_t = (int)round(dx * t + n);
			y_t = (int)round(dy * t + m);
			if (x_t < 0 || angles->width() <= x_t || y_t < 0 || angles->height() <= y_t) {
				break;
			}
			if (IndexMap->get(x_t, y_t) == n_seg) {
				if (fabs(angles->get(x_t, y_t) - aligned_angle) <= DIR_PROBABILITY
				    || fabs(angles->get(x_t, y_t) - ANGLE_MAX - aligned_angle) <= DIR_PROBABILITY
				    || fabs(angles->get(x_t, y_t) + ANGLE_MAX - aligned_angle) <= DIR_PROBABILITY) {
					k++;
				}
			}
		}
#pragma omp critical
		if (k >= k_list[L]) {
			Num_EPSegments++;
			EPSegments_Pr[n_seg] = Pr_table->get(k, L);
		}
#pragma omp critical
		{
			present_count++;
			if (NUM_PROGRESS * present_count / (*Num_Segments) > progress_count) {
				progress_count = NUM_PROGRESS * (present_count - 1) / (*Num_Segments); // Take account of Overflow
				printf("\r%3d%% |%s#\x1b[1A\n", 100 * present_count / (*Num_Segments), Progress[progress_count].c_str());
			}
		}
	}
	// Make the maximal exclusive principle segments list
	try {
		MaxEPSegments = new SEGMENT[Num_EPSegments];
	}
	catch (const std::bad_alloc &bad) {
		Error.Function("new");
		Error.Value("MaxEPSegments");
		Error.Malloc();
		goto ExitError;
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
	delete[] EPSegments_Pr;
	EPSegments_Pr = nullptr;
	return MaxEPSegments;
// Error
ExitError:
	delete[] EPSegments_Pr;
	delete[] MaxEPSegments;
	return nullptr;
}

