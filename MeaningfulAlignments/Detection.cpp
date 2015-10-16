#include "../Scratch_MeaningfulMotion.h"

#define DEBUG_FILTER



ImgVector<double> *
DetectScratch(const PNM &pnm, double s_med, double s_avg, FILTER_PARAM FilterParam, int Do_Detection)
{
	ERROR Error("DetectScratch");

	ImgVector<double> *scratches = nullptr;
	ImgVector<double> *img = nullptr;
	ImgVector<double> *img_filtered = nullptr;
	SIZE size;
	int x, y, m, n;
	double Im, Il, Ir;

#if defined(DEBUG_FILTER)
	PNM pnm_out;
#endif

	// Initialize
	size.width = pnm.Width();
	size.height = pnm.Height();
	try {
		img = new ImgVector<double>(pnm.Width(), pnm.Height(), pnm.get_double());
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("img");
		Error.Malloc();
		goto ExitError;
	}
	switch (FilterParam.type) {
		case FILTER_ID_EPSILON: // Epsilon Filter
			printf("* Epsilon filtering on input\n");
			if ((img_filtered = EpsilonFilter(img, FilterParam)) == nullptr) {
				Error.Function("EpsilonFilter");
				Error.Value("*img_filtered");
				Error.FunctionFail();
				goto ExitError;
			}
			break;
		case FILTER_ID_GAUSSIAN: // Gaussian Filter
			printf("* Gaussian filtering on input\n");
			if ((img_filtered = Gaussian(img, FilterParam)) == nullptr) {
				Error.Function("Gaussian");
				Error.Value("*img_filtered");
				Error.FunctionFail();
				goto ExitError;
			}
			break;
		default:
			printf("* Do NOT any filtering\n");
			try {
				img_filtered = new ImgVector<double>(*img);
			}
			catch (const std::bad_alloc &bad) {
				Error.Value("img_filtered");
				Error.Malloc();
				goto ExitError;
			}
	}
#if defined(DEBUG_FILTER)
	if (pnm_out.copy(PORTABLE_GRAYMAP_BINARY, size.width, size.height, pnm.MaxInt(), img_filtered->data(), 1.0) == PNM_FUNCTION_ERROR) {
		Error.Function("pnm_out.copy");
		Error.FunctionFail();
		goto ExitError;
	}
	if (pnm_out.write("filtered.pgm") == PNM_FUNCTION_ERROR) {
		Error.Function("pnm_out.write");
		Error.FunctionFail();
		goto ExitError;
	}
	pnm_out.free();
#endif

	if (Do_Detection == 0) {
		scratches = img_filtered;
		img_filtered = nullptr;
	} else {
		try {
			scratches = new ImgVector<double>(size.width, size.height);
		}
		catch (const std::bad_alloc &bad) {
			Error.Function("new");
			Error.Value("scratches");
			Error.Malloc();
			goto ExitError;
		}
#pragma omp parallel for private(x, m, n, Im, Il, Ir)
		for (y = 0; y < size.height; y++) {
			for (x = 0; x < size.width; x++) {
				Im = HorizontalMedian(img_filtered, x, y, MEAN_WIDTH);
				if (fabs(img_filtered->get(x, y) - Im) >= s_med) {
					Il = m = 0;
					for (n = ((x - AVE_FAR) < 0) ? 0 : (x - AVE_FAR); (n < size.width) && (n < (x - SCRATCH_WIDTH / 2)); n++) {
						Il += img_filtered->get(n, y);
						m++;
					}
					Il /= (double)m;
					Ir = m = 0;
					for (n = x + SCRATCH_WIDTH / 2 + 1; n < size.width && n <= x + AVE_FAR; n++) {
						Ir += img_filtered->get(n, y);
						m++;
					}
					Ir /= (double)m;
					if (fabs(Il - Ir) <= s_avg) {
						scratches->ref(x, y) = PLOT_INTENSITY_MAX;
					}
				}
			}
		}
	}
	delete img_filtered;
	delete img;
	return scratches;
// Errors
ExitError:
	delete img_filtered;
	delete img;
	delete scratches;
	return nullptr;
}


SEGMENT *
AlignedSegment_vertical(ImgVector<double> *angles, int *k_list, int l_min, ImgVector<double> *Pr_table, int *Num_Segments, int Max_Length, int Max_Output_Length)
{
	ERROR Error("AlignedSegment_vertical");
	std::string ErrorDescription;

	double rad_offset = M_PI * (0.5 - 0.5 / (double)DIV_ANGLE_VERTICAL);
	double tan_list[DIV_ANGLE];
	int r, m, n, x, y;
	int L;
	double dx, dy;
	SEGMENT *array_segment = nullptr;
	std::list<FRAGMENT>* list_fragment = nullptr;
	std::list<SEGMENT> list_segment;
	std::list<SEGMENT>::iterator itr_segment;
	double progress;
	int current_count;

	if (angles == nullptr) {
		Error.Value("angles");
		Error.PointerNull();
		goto ExitError;
	} else if (k_list == nullptr) {
		Error.Value("k_list");
		Error.PointerNull();
		goto ExitError;
	} else if (Pr_table == nullptr) {
		Error.Value("Pr_table");
		Error.PointerNull();
		goto ExitError;
	} else if (Num_Segments == nullptr) {
		Error.Value("Num_Segments");
		Error.PointerNull();
		goto ExitError;
	}
	printf("- Search only the segments that satisfy\n\tlength > %d", l_min);
	if (Max_Length > 0) {
		printf(" AND length <= %d", Max_Length);
	}
	printf("\n");

	for (r = 0; r < DIV_ANGLE; r++) {
		if (r == DIV_ANGLE / 2) {
			// Set the value out of range instead of infinity
			tan_list[r] = 2.0 * (angles->height() > angles->width() ? angles->height() : angles->width());
		} else {
			tan_list[r] = tan((M_PI / (double)DIV_ANGLE_VERTICAL) * r / (double)DIV_ANGLE + rad_offset);
		}
	}
	progress = .0;
	current_count = 0;
	printf("* Search segments starts from Upper or Bottom edge :\n   0.0%% |%s\x1b[1A\n", Progress_End.c_str());
#pragma omp parallel for schedule(dynamic) private(list_fragment, r, x, y, dx, dy)
	for (n = 0; n < angles->width(); n++) {
		for (r = 0; r < DIV_ANGLE; r++) {
			// Upper side to Other 3 sides (n, 0) -> (x, y)
			dx = n + round((angles->height() - 1) / tan_list[r]);
			x = (dx >= 0.0) ? (dx < (double)angles->width()) ? (int)dx : angles->width() - 1 : 0;
			if (tan_list[r] >= 0.0) {
				dy = round((angles->width() - 1 - n) * tan_list[r]);
			} else {
				dy = round(-n * tan_list[r]);
			}
			y = (dy >= 0.0) ? (dy < (double)angles->height()) ? (int)dy : angles->height() - 1 : 0;
			list_fragment = AlignedCheck(angles, k_list, Pr_table, l_min, 0, n, x, y, Max_Length);
			if (list_fragment == nullptr) {
				printf("error\n");
				ErrorDescription = "Occured at Upper side to Other 3 sides";
				break;
			}
			if (MaximalMeaningfulness(&list_segment, list_fragment, 0, n, x, y, Max_Output_Length) != MEANINGFUL_SUCCESS) {
				printf("error\n");
				ErrorDescription = "Occured at Bottom side to Other 3 sides";
				break;
			}
			if (list_fragment != nullptr) {
				list_fragment->clear();
				delete list_fragment;
				list_fragment = nullptr;
			}
			// Bottom side to Other 3 sides (n, size.height) -> (x, y)
			dx = n + round(-(angles->height() - 1) / tan_list[r]);
			x = (dx >= 0.0) ? (dx < (double)angles->width()) ? (int)dx : angles->width() - 1 : 0;
			if (tan_list[r] >= 0.0) {
				dy = angles->height() - 1 + round(-n * tan_list[r]);
			} else {
				dy = angles->height() - 1 + round((angles->width() - 1 - n) * tan_list[r]);
			}
			y = (dy >= 0.0) ? (dy < (double)angles->height()) ? (int)dy : angles->height() - 1 : 0;
			list_fragment = AlignedCheck(angles, k_list, Pr_table, l_min, angles->height() - 1, n, x, y, Max_Length);
			if (list_fragment == nullptr) {
				printf("error\n");
				ErrorDescription = "Occured at Bottom side to Other 3 sides";
				break;
			}
			if (MaximalMeaningfulness(&list_segment, list_fragment, angles->height() - 1, n, x, y, Max_Output_Length) != MEANINGFUL_SUCCESS) {
				printf("error\n");
				ErrorDescription = "Occured at Bottom side to Other 3 sides";
				break;
			}
			if (list_fragment != nullptr) {
				list_fragment->clear();
				delete list_fragment;
				list_fragment = nullptr;
			}
#pragma omp critical
			{
				current_count++;
				if (round((double)current_count / (DIV_ANGLE * angles->width()) * 1000.0) > progress) {
					progress = round((double)current_count / (DIV_ANGLE * angles->width()) * 1000.0); // Take account of Overflow
					printf("\r %5.1f%% |%s#\x1b[1A\n", progress * 0.1, Progress[NUM_PROGRESS * current_count / (DIV_ANGLE * angles->width() + 1)].c_str());
				}
			}
		}
	}
	if (ErrorDescription.empty() == false) {
		Error.Others(ErrorDescription.c_str());
		goto ExitError;
	}
	printf("\nComplete!\n");
	// Count the number of segments
	L = 0;
	for (itr_segment = list_segment.begin();
	    itr_segment != list_segment.end();
	    ++itr_segment) {
		L++;
	}
	try {
		array_segment = new SEGMENT[L];
	}
	catch (const std::bad_alloc &bad) {
		Error.Function("new");
		Error.Value("array_coord");
		Error.Malloc();
		goto ExitError;
	}
	itr_segment = list_segment.begin();
	for (m = 0; m < L; m++) {
		array_segment[m] = *itr_segment;
		++itr_segment;
	}
	*Num_Segments = L;
	list_segment.clear();
	return array_segment;
// Error
ExitError:
	delete[] array_segment;
	list_segment.clear();
	return nullptr;
}


std::list<FRAGMENT> *
AlignedCheck(ImgVector<double> *angles, int *k_list, ImgVector<double> *Pr_table, int l_min, int m, int n, int x, int y, int Max_Length)
{
	ERROR Error("AlignedCheck");
	ATAN2_DIV_PI atan2_div_pi(angles->width(), angles->height());
	std::list<FRAGMENT>* list_fragment = nullptr;
	FRAGMENT fragment_data;
	SEGMENT segment_data;
	double Pr_max, Pr_new;
	double aligned_angle;
	int L;
	double dx, dy;
	int t_start, t_end, t_end_max;
	int k, t;

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
	try {
		list_fragment = new std::list<FRAGMENT>;
	}
	catch (const std::bad_alloc &bad) {
		Error.Function("new");
		Error.Value("list_fragment");
		Error.Malloc();
		goto ExitError;
	}
	for (t_start = 0; t_start <= L - l_min; t_start++) {
		int tmpx, tmpy;
		tmpx = (int)round(dx * t_start + n);
		tmpx = tmpx >= 0 ? (tmpx < angles->width() ? tmpx : angles->width() - 1) : 0;
		tmpy = (int)round(dy * t_start + m);
		tmpy = tmpy >= 0 ? (tmpy < angles->height() ? tmpy : angles->height() - 1) : 0;
		if (fabs(angles->get(tmpx, tmpy) - aligned_angle) <= DIR_PROBABILITY
		    || fabs(angles->get(tmpx, tmpy) - ANGLE_MAX - aligned_angle) <= DIR_PROBABILITY
		    || fabs(angles->get(tmpx, tmpy) + ANGLE_MAX - aligned_angle) <= DIR_PROBABILITY) {
			Pr_max = 1;
			t_end_max = 0;
			for (t_end = (l_min > 1) ? t_start + l_min - 1 : t_start + 1; t_end < L; t_end++) {
				if (Max_Length > 0) {
					if (t_end_max > 0 && t_end_max - t_start + 1 <= Max_Length && t_end - t_start + 1 > Max_Length) {
						// The length of the segment EXCEEDS the length limit
						fragment_data = (FRAGMENT){t_start, t_end_max, Pr_max};
						list_fragment->push_front(fragment_data);
						// Reset "t_end_max" but keep "Pr_max" to search the Longer Segments as if the process not limited by Length.
						t_end_max = 0;
					}
				}
				tmpx = (int)round(dx * t_end + n);
				tmpx = tmpx >= 0 ? (tmpx < angles->width() ? tmpx : angles->width() - 1) : 0;
				tmpy = (int)round(dy * t_end + m);
				tmpy = tmpy >= 0 ? (tmpy < angles->height() ? tmpy : angles->height() - 1) : 0;
				if (fabs(angles->get(tmpx, tmpy) - aligned_angle) <= DIR_PROBABILITY
				    || fabs(angles->get(tmpx, tmpy) - ANGLE_MAX - aligned_angle) <= DIR_PROBABILITY
				    || fabs(angles->get(tmpx, tmpy) + ANGLE_MAX - aligned_angle) <= DIR_PROBABILITY) {
					k = 2;
					// Check if the segment which has aligned point on both ends is epsilon-Meaningful
					for (t = t_start + 1; t <= t_end - 1; t++) {
						tmpx = (int)round(dx * t + n);
						tmpx = tmpx >= 0 ? (tmpx < angles->width() ? tmpx : angles->width() - 1) : 0;
						tmpy = (int)round(dy * t + m);
						tmpy = tmpy >= 0 ? (tmpy < angles->height() ? tmpy : angles->height() - 1) : 0;
						if (fabs(angles->get(tmpx, tmpy) - aligned_angle) <= DIR_PROBABILITY
						    || fabs(angles->get(tmpx, tmpy) - ANGLE_MAX - aligned_angle) <= DIR_PROBABILITY
						    || fabs(angles->get(tmpx, tmpy) + ANGLE_MAX - aligned_angle) <= DIR_PROBABILITY) {
							k++;
						}
					}
					if (k >= k_list[t_end - t_start + 1]) {
						Pr_new = Pr_table->get(k, t_end - t_start + 1);
						if (Pr_new <= Pr_max) {
							Pr_max = Pr_new;
							t_end_max = t_end;
						}
					}
				}
			}
			if (t_end_max > 0) {
				fragment_data = (FRAGMENT){t_start, t_end_max, Pr_max};
				list_fragment->push_front(fragment_data);
			}
		}
	}
	return list_fragment;
// Error
ExitError:
	if (list_fragment != nullptr) {
		list_fragment->clear();
		delete list_fragment;
	}
	return nullptr;
}


bool
MaximalMeaningfulness(std::list<SEGMENT>* list_segment, std::list<FRAGMENT>* list_fragment, int m, int n, int x, int y, int Max_Output_Length)
{
	std::list<FRAGMENT>::iterator itr_fragment;
	std::list<FRAGMENT>::iterator fragment_ref;
	std::list<FRAGMENT>::iterator fragment_cur;
	int L;
	double dx, dy;
	SEGMENT segment_data;

	if (list_fragment == nullptr) {
		return MEANINGFUL_SUCCESS;
	}
	if (abs(x - n) > abs(y - m)) {
		L = abs(x - n) + 1;
	} else {
		L = abs(y - m) + 1;
	}
	dx = (x - n) / (double)(L - 1.0);
	dy = (y - m) / (double)(L - 1.0);
	// Maximal Meaningfulness
	for (fragment_ref = list_fragment->begin();
	    fragment_ref != list_fragment->end();
	    ++fragment_ref) {
		fragment_cur = list_fragment->begin();
		while (fragment_ref != list_fragment->end()
		    && fragment_cur != list_fragment->end()) {
			if (fragment_cur == fragment_ref) {
				++fragment_cur;
			} else if (fragment_ref->start <= fragment_cur->start && fragment_cur->end <= fragment_ref->end) {
				// Delete the fragment which has higher probability
				if (fragment_ref->Pr <= fragment_cur->Pr) {
					fragment_cur = list_fragment->erase(fragment_cur);
				} else {
					fragment_ref = list_fragment->erase(fragment_ref);
				}
			} else if (fragment_cur->start <= fragment_ref->start && fragment_ref->end <= fragment_cur->end) {
				// Delete the fragment which has higher probability
				if (fragment_cur->Pr <= fragment_ref->Pr) {
					fragment_ref = list_fragment->erase(fragment_ref);
				} else {
					fragment_cur = list_fragment->erase(fragment_cur);
				}
			} else {
				++fragment_cur;
			}
		}
	}
	// Write out to coord_list
	for (itr_fragment = list_fragment->begin();
	    itr_fragment != list_fragment->end();
	    ++itr_fragment) {
		// If Max_Output_Length <= 0 then Do NOT Limit The Length of Segments
		if (Max_Output_Length <= 0
		    || (itr_fragment->end - itr_fragment->start + 1) <= Max_Output_Length) {
#pragma omp critical
			{
				segment_data = (SEGMENT){
				    (int)round(n + dx * itr_fragment->start),
				    (int)round(m + dy * itr_fragment->start),
				    (int)round(n + dx * itr_fragment->end),
				    (int)round(m + dy * itr_fragment->end),
				    itr_fragment->Pr};
				list_segment->push_front(segment_data);
			}
		}
	}
	return MEANINGFUL_SUCCESS;
}

