#include "Scratch_MeaningfulMotion.h"




double **
Pyramider(double *img, SIZE size, int Level)
{
#define WEIGHTED_FILTER_SIZE 5
	char *FunctionName = "Pyramider";
	char *ErrorFunctionName = "";
	char *ErrorValueName = "";

	double **Pyramid = NULL;
	double w[WEIGHTED_FILTER_SIZE];
	double a = 0.4; /* Typically a = {z | 0.3 <= z <= 0.5} */
	SIZE size_l = SIZE_ZERO;
	SIZE size_lm1 = SIZE_ZERO;
	int x, y;
	int m, n;
	int xn, ym;
	int l;
	double sum;

	PNM pnm = PNM_NULL;
	char filename[128];

	if ((Pyramid = (double **)calloc(Level, sizeof(double *))) == NULL) {
		ErrorFunctionName = "calloc";
		ErrorValueName = "Pyramid";
		goto ErrorMalloc;
	}
	/* Lowpass Filter */
	w[0] = a / 2.0;
	w[1] = 0.5;
	w[2] = a;
	w[3] = 0.5;
	w[4] = a / 2.0;
	/* Scaling for sum equal to zero */
	sum = .0;
	for (x = 0; x < WEIGHTED_FILTER_SIZE; x++) {
		sum += w[x];
	}
	for (x = 0; x < WEIGHTED_FILTER_SIZE; x++) {
		w[x] /= sum;
	}

	/* Level 0 (Just same as original) */
	size_l = size;
	if ((Pyramid[0] = (double *)calloc((size_t)size.width * size.height, sizeof(double))) == NULL) {
		ErrorFunctionName = "calloc";
		ErrorValueName = "Pyramid[0]";
		goto ErrorMalloc;
	}
	for (x = 0; x < size.width * size.height; x++) {
		Pyramid[0][x] = img[x];
	}

	/* Test Output */
	pnmnew(&pnm, PORTABLE_GRAYMAP_BINARY, size.width, size.height, 255);
	for (x = 0; x < size.width * size.height; x++) {
		pnm.img[x] = SATURATE(Pyramid[0][x], 0, 255);
	}
	sprintf(filename, "Pyramid_0000.pgm");
	pnmwrite(&pnm, filename);
	pnmfree(&pnm);

	/* Make image pyramid (l > 0) */
	for (l = 1; l < Level; l++) {
		size_lm1 = size_l;
		size_l.width = (int)floor(size_lm1.width / 2.0);
		size_l.height = (int)floor(size_lm1.height / 2.0);
		if (size_l.width <= 0 || size_l.height <= 0) {
			fprintf(stderr, "*** %s() warning - The image reaches minimum size in Pyramid ***\n", FunctionName);
			break;
		}
		if ((Pyramid[l] = (double *)calloc((size_t)size_l.width * size_l.height, sizeof(double))) == NULL) {
			ErrorFunctionName = "calloc";
			ErrorValueName = "Pyramid[l]";
			goto ErrorMalloc;
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
		pnmnew(&pnm, PORTABLE_GRAYMAP_BINARY, size_l.width, size_l.height, 255);
		for (x = 0; x < size_l.width * size_l.height; x++) {
			pnm.img[x] = SATURATE(Pyramid[l][x], 0, 255);
		}
		sprintf(filename, "Pyramid_%04d.pgm", l);
		pnmwrite(&pnm, filename);
		pnmfree(&pnm);
	}
	return Pyramid;
/* Error */
ErrorMalloc:
	fprintf(stderr, "*** %s() error - Cannot allocate memory for (*%s) by %s() ***\n", FunctionName, ErrorValueName, ErrorFunctionName);
	for (l = 0; Pyramid != NULL && l < Level; l++) {
		free(Pyramid[l]);
	}
	free(Pyramid[l]);
	return NULL;
}


VECTOR_2D **
grad_Pyramid(double **img_t, double **img_tp1, SIZE size, int Level)
{
	char *FunctionName = "grad_Pyramid";
	char *ErrorFunction = "";
	char *ErrorValue = "";

	/* Note that convolution invert *Filter */
	VECTOR_2D **grad_levels = NULL;
	SIZE size_l;
	int x, y;
	int m, n;
	int l;

	/* img_tp1 is allowed to be NULL */
	if (img_t == NULL) {
		ErrorValue = "img_t";
		goto ErrorPointerNull;
	}
	if ((grad_levels = (VECTOR_2D **)calloc(Level, sizeof(VECTOR_2D *))) == NULL) {
		ErrorFunction = "calloc";
		ErrorValue = "grad_levels";
		goto ErrorMalloc;
	}
	for (l = 0; l < Level; l++) {
		size_l.width = (int)floor((double)size.width / pow_int(2.0, l));
		size_l.height = (int)floor((double)size.height / pow_int(2.0, l));
		if ((grad_levels[l] = (VECTOR_2D *)calloc(size_l.width * size_l.height, sizeof(VECTOR_2D))) == NULL) {
			ErrorFunction = "calloc";
			ErrorValue = "grad_levels[l]";
			goto ErrorMalloc;
		}
		for (m = 0; m < size_l.height; m++) {
			y = SATURATE(m, 0, size_l.height - 2);
			for (n = 0; n < size_l.width; n++) {
				x = SATURATE(n, 0, size_l.width - 2);
				/* dx */
				grad_levels[l][size_l.width * y + x].x =
				    (img_t[l][size_l.width * y + x + 1] - img_t[l][size_l.width * y + x]
				    + img_t[l][size_l.width * (y + 1) + x + 1] - img_t[l][size_l.width * (y + 1) + x])
				    / 2.0;
				/* dy */
				grad_levels[l][size_l.width * y + x].y =
				    (img_t[l][size_l.width * (y + 1) + x] - img_t[l][size_l.width * y + x]
				    + img_t[l][size_l.width * (y + 1) + x + 1] - img_t[l][size_l.width * y + x + 1])
				    / 2.0;
				if (img_tp1 != NULL) {
					/* dx */
					grad_levels[l][size_l.width * y + x].x +=
					    (img_tp1[l][size_l.width * y + x + 1] - img_tp1[l][size_l.width * y + x]
					    + img_tp1[l][size_l.width * (y + 1) + x + 1] - img_tp1[l][size_l.width * (y + 1) + x])
					    / 2.0;
					/* dy */
					grad_levels[l][size_l.width * y + x].y +=
					    (img_tp1[l][size_l.width * (y + 1) + x] - img_tp1[l][size_l.width * y + x]
					    + img_tp1[l][size_l.width * (y + 1) + x + 1] - img_tp1[l][size_l.width * y + x + 1])
					    / 2.0;
				}
			}
		}
	}
	return grad_levels;
/* Error */
ErrorPointerNull:
	fprintf(stderr, "*** %s() error - The pointer (*%s) is NULL ***\n", FunctionName, ErrorValue);
	return NULL;
ErrorMalloc:
	fprintf(stderr, "*** %s() error - Cannot allocate memory for (*%s) by %s() ***\n", FunctionName, ErrorValue, ErrorFunction);
	if (grad_levels != NULL) {
		for (l = 0; l < Level; l++) {
			free(grad_levels[l]);
		}
	}
	free(grad_levels);
	return NULL;
}


double **
dt_Pyramid(double **img_t, double **img_tp1, SIZE size, int Level)
{
	char *FunctionName = "grad_Pyramid";
	char *ErrorFunction = "";
	char *ErrorValue = "";

	/* Note that convolution invert *Filter */
	double **dt_levels = NULL;
	SIZE size_l;
	int x, y;
	int m, n;
	int l;

	if (img_t == NULL) {
		ErrorValue = "img_t";
		goto ErrorPointerNull;
	} else if (img_tp1 == NULL) {
		ErrorValue = "img_tp1";
		goto ErrorPointerNull;
	}
	if ((dt_levels = (double **)calloc(Level, sizeof(double *))) == NULL) {
		ErrorFunction = "calloc";
		ErrorValue = "dt_levels";
		goto ErrorMalloc;
	}
	for (l = 0; l < Level; l++) {
		size_l.width = (int)floor((double)size.width / pow_int(2.0, l));
		size_l.height = (int)floor((double)size.height / pow_int(2.0, l));
		if ((dt_levels[l] = (double *)calloc(size_l.width * size_l.height, sizeof(double))) == NULL) {
			ErrorFunction = "calloc";
			ErrorValue = "dt_levels[l]";
			goto ErrorMalloc;
		}
		for (m = 0; m < size_l.height - 1; m++) {
			y = SATURATE(m, 0, size_l.height - 2);
			for (n = 0; n < size_l.width - 1; n++) {
				x = SATURATE(n, 0, size_l.width - 2);
				/* dt */
				dt_levels[l][size_l.width * y + x] =
				    (img_tp1[l][size_l.width * y + x] - img_t[l][size_l.width * y + x]
				    + img_tp1[l][size_l.width * y + x + 1] - img_t[l][size_l.width * y + x + 1]
				    + img_tp1[l][size_l.width * (y + 1) + x] - img_t[l][size_l.width * (y + 1) + x]
				    + img_tp1[l][size_l.width * (y + 1) + x + 1] - img_t[l][size_l.width * (y + 1) + x + 1])
				    / 4.0;
			}
		}
	}
	return dt_levels;
/* Error */
ErrorPointerNull:
	fprintf(stderr, "*** %s() error - The pointer (*%s) is NULL ***\n", FunctionName, ErrorValue);
	return NULL;
ErrorMalloc:
	fprintf(stderr, "*** %s() error - Cannot allocate memory for (*%s) by %s() ***\n", FunctionName, ErrorValue, ErrorFunction);
	if (dt_levels != NULL) {
		for (l = 0; l < Level; l++) {
			free(dt_levels[l]);
		}
	}
	free(dt_levels);
	return NULL;
}

