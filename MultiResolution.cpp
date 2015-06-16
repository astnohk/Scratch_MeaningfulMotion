#include "Scratch_MeaningfulMotion.h"




double **
Pyramider(double *img, SIZE size, int Level)
{
#define WEIGHTED_FILTER_SIZE 5
	ERROR Error("Pyramider");

	double **Pyramid = nullptr;
	double w[WEIGHTED_FILTER_SIZE];
	double a = 0.4; // Typically a = {z | 0.3 <= z <= 0.5}
	SIZE size_l;
	SIZE size_lm1;
	int x, y;
	int m, n;
	int xn, ym;
	int l;
	double sum;

	PNM pnm;
	char filename[128];

	try {
		Pyramid = new double*[Level];
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("Pyramid");
		Error.Malloc();
		goto ExitError;
	}
	// Lowpass Filter
	w[0] = a / 2.0;
	w[1] = 0.5;
	w[2] = a;
	w[3] = 0.5;
	w[4] = a / 2.0;
	// Scaling for sum equal to zero
	sum = .0;
	for (x = 0; x < WEIGHTED_FILTER_SIZE; x++) {
		sum += w[x];
	}
	for (x = 0; x < WEIGHTED_FILTER_SIZE; x++) {
		w[x] /= sum;
	}

	// Level 0 (Just same as original)
	size_l = size;
	try {
		Pyramid[0] = new double[size.width * size.height];
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("Pyramid[0]");
		Error.Malloc();
		goto ExitError;
	}
	for (x = 0; x < size.width * size.height; x++) {
		Pyramid[0][x] = img[x];
	}

	/* Test Output */
	pnm.copy(PORTABLE_GRAYMAP_BINARY, size.width, size.height, 255, Pyramid[0], 256.0);
	sprintf(filename, "Pyramid_0000.pgm");
	pnm.write(filename);
	pnm.free();

	/* Make image pyramid (l > 0) */
	for (l = 1; l < Level; l++) {
		size_lm1 = size_l;
		size_l.width = (int)floor(size.width * pow_int(0.5, l));
		size_l.height = (int)floor(size.height * pow_int(0.5, l));
		if (size_l.width <= 0 || size_l.height <= 0) {
			Error.OthersWarning("The image reaches minimum size in Pyramid");
			break;
		}
		try {
			Pyramid[l] = new double[size_l.width * size_l.height];
		}
		catch (const std::bad_alloc &bad) {
			Error.Value("Pyramid[l]");
			Error.Malloc();
			goto ExitError;
		}
		for (x = 0; x < size_l.width; x++) {
			for (y = 0; y < size_l.height; y++) {
				Pyramid[l][size_l.width * y + x] = .0;
				for (m = 0; m < WEIGHTED_FILTER_SIZE; m++) {
					ym = 2 * y + m - (int)floor(WEIGHTED_FILTER_SIZE / 2.0);
					if (ym < 0) {
						ym = abs(ym) - 1;
					} else if (ym >= size_lm1.height) {
						ym = 2 * size_lm1.height - ym - 1;
					}
					for (n = 0; n < WEIGHTED_FILTER_SIZE; n++) {
						xn = 2 * x + n - (int)floor(WEIGHTED_FILTER_SIZE / 2.0);
						if (xn < 0) {
							xn = abs(xn) - 1;
						} else if (xn >= size_lm1.width) {
							xn = 2 * size_lm1.width - xn - 1;
						}
						Pyramid[l][size_l.width * y + x] += w[m] * w[n] * Pyramid[l - 1][size_lm1.width * ym + xn];
					}
				}
			}
		}
		pnm.copy(PORTABLE_GRAYMAP_BINARY, size_l.width, size_l.height, 255, Pyramid[l], 256.0);
		sprintf(filename, "Pyramid_%04d.pgm", l);
		pnm.write(filename);
		pnm.free();
	}
	return Pyramid;
/* Error */
ExitError:
	for (l = 0; Pyramid != nullptr && l < Level; l++) {
		delete[] Pyramid[l];
	}
	delete[] Pyramid;
	return nullptr;
}


VECTOR_2D **
grad_Pyramid(double **img_t, double **img_tp1, SIZE size, int Level)
{
	ERROR Error("grad_Pyramid");

	/* Note that convolution invert *Filter */
	VECTOR_2D **grad_levels = nullptr;
	SIZE size_l;
	int x, y;
	int m, n;
	int l;

	/* img_tp1 is allowed to be NULL */
	if (img_t == nullptr) {
		Error.Value("img_t");
		Error.PointerNull();
		goto ExitError;
	}
	try {
		grad_levels = new VECTOR_2D*[Level];
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("grad_levels");
		Error.Malloc();
		goto ExitError;
	}
	for (l = 0; l < Level; l++) {
		size_l.width = (int)floor((double)size.width * pow_int(0.5, l));
		size_l.height = (int)floor((double)size.height * pow_int(0.5, l));
		try {
			grad_levels[l] = new VECTOR_2D[size_l.width * size_l.height];
		}
		catch (const std::bad_alloc &bad) {
			Error.Value("grad_levels[l]");
			Error.Malloc();
			goto ExitError;
		}
		for (m = 0; m < size_l.height; m++) {
			y = SATURATE(m, 0, size_l.height - 2);
			for (n = 0; n < size_l.width; n++) {
				x = SATURATE(n, 0, size_l.width - 2);
				/* dx */
				grad_levels[l][size_l.width * m + n].x =
				    (img_t[l][size_l.width * y + x + 1] - img_t[l][size_l.width * y + x]
				    + img_t[l][size_l.width * (y + 1) + x + 1] - img_t[l][size_l.width * (y + 1) + x])
				    / 2.0;
				/* dy */
				grad_levels[l][size_l.width * m + n].y =
				    (img_t[l][size_l.width * (y + 1) + x] - img_t[l][size_l.width * y + x]
				    + img_t[l][size_l.width * (y + 1) + x + 1] - img_t[l][size_l.width * y + x + 1])
				    / 2.0;
				if (img_tp1 != NULL) {
					/* dx */
					grad_levels[l][size_l.width * m + n].x +=
					    (img_tp1[l][size_l.width * y + x + 1] - img_tp1[l][size_l.width * y + x]
					    + img_tp1[l][size_l.width * (y + 1) + x + 1] - img_tp1[l][size_l.width * (y + 1) + x])
					    / 2.0;
					/* dy */
					grad_levels[l][size_l.width * m + n].y +=
					    (img_tp1[l][size_l.width * (y + 1) + x] - img_tp1[l][size_l.width * y + x]
					    + img_tp1[l][size_l.width * (y + 1) + x + 1] - img_tp1[l][size_l.width * y + x + 1])
					    / 2.0;
				}
			}
		}
	}
	return grad_levels;
// Error
ExitError:
	if (grad_levels != NULL) {
		for (l = 0; l < Level; l++) {
			delete[] grad_levels[l];
		}
	}
	delete[] grad_levels;
	return nullptr;
}


double **
dt_Pyramid(double **img_t, double **img_tp1, SIZE size, int Level)
{
	ERROR Error("grad_Pyramid");

	// Note that convolution invert *Filter
	double **dt_levels = nullptr;
	SIZE size_l;
	int x, y;
	int m, n;
	int l;

	if (img_t == nullptr) {
		Error.Value("img_t");
		Error.PointerNull();
		goto ExitError;
	} else if (img_tp1 == nullptr) {
		Error.Value("img_tp1");
		Error.PointerNull();
		goto ExitError;
	}
	try {
		dt_levels = new double*[Level];
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("dt_levels");
		Error.Malloc();
		goto ExitError;
	}
	for (l = 0; l < Level; l++) {
		size_l.width = (int)floor((double)size.width * pow_int(0.5, l));
		size_l.height = (int)floor((double)size.height * pow_int(0.5, l));
		try {
			dt_levels[l] = new double[size_l.width * size_l.height];
		}
		catch (const std::bad_alloc &bad) {
			Error.Value("dt_levels[l]");
			Error.Malloc();
			goto ExitError;
		}
		for (m = 0; m < size_l.height; m++) {
			y = SATURATE(m, 0, size_l.height - 2);
			for (n = 0; n < size_l.width; n++) {
				x = SATURATE(n, 0, size_l.width - 2);
				// dt
				dt_levels[l][size_l.width * m + n] =
				    (img_tp1[l][size_l.width * y + x] - img_t[l][size_l.width * y + x]
				    + img_tp1[l][size_l.width * y + x + 1] - img_t[l][size_l.width * y + x + 1]
				    + img_tp1[l][size_l.width * (y + 1) + x] - img_t[l][size_l.width * (y + 1) + x]
				    + img_tp1[l][size_l.width * (y + 1) + x + 1] - img_t[l][size_l.width * (y + 1) + x + 1])
				    / 4.0;
			}
		}
	}
	return dt_levels;
// Error
ExitError:
	if (dt_levels != nullptr) {
		for (l = 0; l < Level; l++) {
			delete[] dt_levels[l];
		}
	}
	delete[] dt_levels;
	return nullptr;
}

