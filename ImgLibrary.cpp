#include "Scratch_MeaningfulA.h"

//#define SHOW_GAUSSIAN_FILTER


double
HorizontalMedian(double *img, int img_width, int x, int y, int width)
{
	double *array = NULL;
	int L = 0;
	int m_s, m_e;
	int m, n;
	double tmp;

	if (x < (width - 1) / 2) {
		m_s = 0;
		m_e = width / 2;
	} else if (x + width / 2 >= img_width) {
		m_s = -(width - 1) / 2;
		m_e = img_width - 1 - x;
	} else {
		m_s = -(width - 1) / 2;
		m_e = width / 2;
	}
	if ((array = (double *)calloc((size_t)(m_e - m_s + 1), sizeof(double))) == NULL) {
		fprintf(stderr, "HorizontalMedian() error - Cannot allocate memory for (*array) by calloc() ***\n");
		return 0;
	}
	for (m = m_s; m < m_e; m++) {
		array[m - m_s] = img[img_width * y + x + m];
	}
	L = m_e - m_s + 1;
	for (m = L - 1; m > 0; m--) { // Sort
		for (n = 0; n < m; n++) {
			if (array[n] < array[n + 1]) {
				tmp = array[n + 1];
				array[n + 1] = array[n];
				array[n] = tmp;
			}
		}
	}
	if (L & 1) {
		tmp = array[L / 2];
	} else {
		tmp = (array[L / 2 - 1] + array[L / 2]) / 2.0;
	}
	free(array);
	array = NULL;
	return tmp;
}


double*
EpsilonFilter(double *img, SIZE size, FILTER_PARAM Param)
{
	char *ErrorFunctionName = "";
	char *ErrorValueName = "";
	int ErrorValue = 0;
	double ErrorValueDouble = .0;

	double *Epsilon = NULL;
	int filter_width_2 = (int)floor(Param.size.width / 2.0);
	int filter_height_2 = (int)floor(Param.size.height / 2.0);
	double Div_coeff = 1.0 / (Param.size.width * Param.size.height);
	int center;
	int x, y;
	int fx, fy;
	int m, n;
	int tmp;

	if (Param.epsilon < .0) {
		ErrorValueName = "epsilon";
		ErrorValueDouble = Param.epsilon;
		goto ErrorIncorrectValueDouble;
	} else if (Param.size.width <= 0 || Param.size.width % 2 == 0) {
		ErrorValueName = "Param.size.width";
		ErrorValue = Param.size.width;
		goto ErrorIncorrectValue;
	} else if (Param.size.height <= 0 || Param.size.height % 2 == 0) {
		ErrorValueName = "Param.size.height";
		ErrorValue = Param.size.height;
		goto ErrorIncorrectValue;
	} else if (img == NULL) {
		ErrorValueName = "img";
		goto ErrorPointerNull;
	}

	if ((Epsilon = (double *)calloc((size_t)(size.width * size.height), sizeof(double))) == NULL) {
		ErrorFunctionName = "calloc";
		ErrorValueName = "Epsilon";
		goto ErrorMalloc;
	}
#pragma omp parallel for private(x, fx, fy, n, m, center, tmp)
	for (y = 0; y < size.height; y++) {
		for (x = 0; x < size.width; x++) {
			center = img[size.width * y + x];
			Epsilon[size.width * y + x] = .0;
			for (fy = -filter_height_2; fy <= filter_height_2; fy++) {
				if (y + fy < 0) {
					tmp = abs(y + fy) / size.height;
					if ((tmp & 1) == 1) {
						m = size.height - (abs(y + fy) - size.height * tmp);
					} else {
						m = abs(y + fy) - size.height * tmp;
					}
				} else if (y + fy > size.height - 1) {
					tmp = (y + fy) / size.height;
					if ((tmp & 1) == 1) {
						m = size.height - 1 - (y + fy - size.height * tmp);
					} else {
						m = y + fy - size.height * tmp;
					}
				} else {
					m = y + fy;
				}
				for (fx = -filter_width_2; fx <= filter_width_2; fx++) {
					if (x + fx < 0) {
						tmp = abs(x + fx) / size.width;
						if ((tmp & 1) == 1) {
							n = size.width - (abs(x + fx) - size.width * tmp);
						} else {
							n = abs(x + fx) - size.width * tmp;
						}
					} else if (x + fx > size.width - 1) {
						tmp = abs(x + fx) / size.width;
						if ((tmp & 1) == 1) {
							n = size.width - 1 - (x + fx - size.width * tmp);
						} else {
							n = x + fx - size.width * tmp;
						}
					} else {
						n = x + fx;
					}
					if (fabs(center - img[size.width * m + n]) <= Param.epsilon) {
						Epsilon[size.width * y + x] += img[size.width * m + n];
					} else {
						Epsilon[size.width * y + x] += center;
					}
				}
			}
			Epsilon[size.width * y + x] *= Div_coeff;
		}
	}
	return Epsilon;
/* Error */
ErrorMalloc:
	fprintf(stderr, "*** EpsilonFilter() error - Cannot allocate memory for (*%s) by %s() ***\n", ErrorValueName, ErrorFunctionName);
	goto ErrorReturn;
ErrorPointerNull:
	fprintf(stderr, "*** EpsilonFilter() error - The pointer (*%s) is NULL ***\n", ErrorValueName);
	goto ErrorReturn;
ErrorIncorrectValue:
	fprintf(stderr, "*** EpsilonFilter() error - The variable (%s) has incorrect value (%d) ***\n", ErrorValueName, ErrorValue);
	goto ErrorReturn;
ErrorIncorrectValueDouble:
	fprintf(stderr, "*** EpsilonFilter() error - The variable (%s) has incorrect value (%f) ***\n", ErrorValueName, ErrorValueDouble);
ErrorReturn:
	free(Epsilon);
	return NULL;
}


double*
Gaussian(double *img, SIZE size, FILTER_PARAM Param)
{
	char *ErrorFunctionName = "";
	char *ErrorValueName = "";
	double *blurred = NULL;
	double *Gauss = NULL;
	int filter_width_2 = 0;
	int filter_height_2 = 0;
	int m, n, y, x, i, j;
	int tmp;
	double sum = 0;
	int norm = 0;
	
	if ((Param.size.width & 1) == 0) {
		norm = 1;
		Param.size.width++;
	}
	if ((Param.size.height & 1) == 0) {
		norm = 1;
		Param.size.height++;
	}
	filter_width_2 = Param.size.width / 2;
	filter_height_2 = Param.size.height / 2;
	if ((Gauss = (double *)calloc((size_t)(Param.size.width * Param.size.height), sizeof(double))) == NULL) {
		ErrorFunctionName = "calloc";
		ErrorValueName = "Gauss";
		goto ErrorMalloc;
	}
	sum = 0;
#if defined(SHOW_GAUSSIAN_FILTER)
	printf("\n  Gaussian Filter :\n");
#endif
	if (norm == 1) {
		// Diamond
		for (m = 0; m < Param.size.height; m++) {
#if defined(SHOW_GAUSSIAN_FILTER)
			printf("    |");
#endif
			for (n = 0; n < Param.size.width; n++) {
				if (filter_width_2 * abs(m - filter_height_2) + filter_height_2 * abs(n - filter_width_2) <= filter_width_2 * filter_height_2) {
					Gauss[Param.size.width * m + n] = exp(-(
					    POW2(m - filter_height_2) + POW2(n - filter_width_2)
					    ) / (2.0 * POW2(Param.std_deviation)));
					sum += Gauss[Param.size.width * m + n];
				}
#if defined(SHOW_GAUSSIAN_FILTER)
				printf(" %.5f", Gauss[Param.size.width * m + n]);
#endif
			}
#if defined(SHOW_GAUSSIAN_FILTER)
			printf(" |\n");
#endif
		}
	} else {
		// Square
		for (m = 0; m < Param.size.height; m++) {
#if defined(SHOW_GAUSSIAN_FILTER)
			printf("    |");
#endif
			for (n = 0; n < Param.size.width; n++) {
				Gauss[Param.size.width * m + n] = exp(-(
				    POW2(m - filter_height_2) + POW2(n - filter_width_2)
				    ) / (2.0 * POW2(Param.std_deviation)));
				sum += Gauss[Param.size.width * m + n];
#if defined(SHOW_GAUSSIAN_FILTER)
				printf(" %.5f", Gauss[Param.size.width * m + n]);
#endif
			}
#if defined(SHOW_GAUSSIAN_FILTER)
			printf(" |\n");
#endif
		}
	}
#if defined(SHOW_GAUSSIAN_FILTER)
	printf("\n");
#endif
	for (m = 0; m < Param.size.width * Param.size.height; m++) {
		Gauss[m] /= sum;
	}
	if ((blurred = (double *)calloc((size_t)(size.height * size.width), sizeof(double))) == NULL) {
		ErrorFunctionName = "calloc";
		ErrorValueName = "blurred";
		goto ErrorMalloc;
	}
#pragma omp parallel for private(n, x, y, i, j, tmp)
	for (m = 0; m < size.height; m++) {
		for (n = 0; n < size.width; n++) {
			for (y = -filter_height_2; y <= filter_height_2; y++) {
				if (m + y < 0) {
					tmp = abs(m + y) / size.height;
					if ((tmp & 1) == 1) {
						i = size.height - (abs(m + y) - size.height * tmp);
					} else {
						i = abs(m + y) - size.height * tmp;
					}
				} else if (m + y > size.height - 1) {
					tmp = (m + y) / size.height;
					if ((tmp & 1) == 1) {
						i = size.height - 1 - (m + y - size.height * tmp);
					} else {
						i = m + y - size.height * tmp;
					}
				} else {
					i = m + y;
				}
				for (x = -filter_width_2; x <= filter_width_2; x++) {
					if (n + x < 0) {
						tmp = abs(n + x) / size.width;
						if ((tmp & 1) == 1) {
							j = size.width - (abs(n + x) - size.width * tmp);
						} else {
							j = abs(n + x) - size.width * tmp;
						}
					} else if (n + x > size.width - 1) {
						tmp = abs(n + x) / size.width;
						if ((tmp & 1) == 1) {
							j = size.width - 1 - (n + x - size.width * tmp);
						} else {
							j = n + x - size.width * tmp;
						}
					} else {
						j = n + x;
					}
					blurred[size.width * m + n] += img[size.width * i + j] * Gauss[Param.size.width * (y + filter_height_2) + x + filter_width_2];
				}
			}
		}
	}
	free(Gauss);
	Gauss = NULL;
	return blurred;
// Errors
ErrorMalloc:
	fprintf(stderr, "*** Gaussian() error - Cannot allocate memory for (*%s) by %s() ***\n", ErrorValueName, ErrorFunctionName);
	free(blurred);
	blurred = NULL;
	free(Gauss);
	Gauss = NULL;
	return NULL;
}


double *
DerivativeAngler(double *img, SIZE size)
{
	char *ErrorFunctionName = "";
	char *ErrorValueName = "";

	VECTOR_2D *Derivative = NULL;
	double *angles = NULL;
	int n;
	double dx, dy;
 
	if (img == NULL) {
		ErrorValueName = "img";
		goto ErrorPointerNull;
	} else if (size.width < 0 || size.height < 0) {
		ErrorValueName = "size";
		goto ErrorIncorrectValue;
	}

	if ((Derivative = Derivator(img, size, "Sobel")) == NULL) {
		ErrorFunctionName = "Derivator";
		ErrorValueName = "Derivative";
		goto ErrorFunctionFailed;
	}
	if ((angles = (double *)calloc((size_t)(size.height * size.width), sizeof(double))) == NULL) {
		ErrorFunctionName = "calloc";
		ErrorValueName = "angles";
		goto ErrorMalloc;
	}
	for (n = 0; n < size.width * size.height; n++) {
		dx = Derivative[n].x;
		dy = Derivative[n].y;
		if (fabs(dx) <= DERIVATIVE_MINIMUM && fabs(dy) <= DERIVATIVE_MINIMUM) {
			angles[n] = -2.0 * ANGLE_MAX;
		} else {
			angles[n] = atan2(dy, dx) / M_PI + 0.5; // Add 0.5 for pi/2 rotation
			if (angles[n] > ANGLE_MAX) {
				angles[n] -= ANGLE_MAX;
			} else if (angles[n] < 0.0) {
				angles[n] += ANGLE_MAX;
			}
		}
	}
	free(Derivative);
	return angles;
// Errors
ErrorMalloc:
	fprintf(stderr, "*** DeriveAngler() error - Cannot allocate memory for (*%s) by %s() ***\n", ErrorValueName, ErrorFunctionName);
	goto ErrorReturn;
ErrorPointerNull:
	fprintf(stderr, "*** DeriveAngler() error - The pointer (*%s) is NULL ***\n", ErrorValueName);
	goto ErrorReturn;
ErrorIncorrectValue:
	fprintf(stderr, "*** DeriveAngler() error - The variable (%s) has incorrect value ***\n", ErrorValueName);
	goto ErrorReturn;
ErrorFunctionFailed:
	fprintf(stderr, "*** DeriveAngler() error - %s() failed to compute (%s) ***\n", ErrorFunctionName, ErrorValueName);
ErrorReturn:
	free(angles);
	angles = NULL;
	return NULL;
}


VECTOR_2D *
Derivator(double *Image, SIZE size, char *Type)
{
	char *FunctionName = "Derivation";
	char *ErrorFunctionName = "";
	char *ErrorValueName = "";

	/* Note that convolution invert *Filter */
	double DiffFilter_x[4] = {0.5, -0.5, 0.5, -0.5};
	double DiffFilter_y[4] = {0.5, 0.5, -0.5, -0.5};
	double SobelFilter_x[9] = {0.25, 0, -0.25, 0.5, 0, -0.5, 0.25, 0, -0.25};
	double SobelFilter_y[9] = {0.25, 0.5, 0.25, 0, 0, 0, -0.25, -0.5, -0.25};
	double *Dx = NULL;
	double *Dy = NULL;
	SIZE size_f = SIZE_ZERO;
	double *Image_dx = NULL;
	double *Image_dy = NULL;
	VECTOR_2D *Derivative = NULL;
	int x;
 
	if (Image == NULL) {
		ErrorValueName = "Image";
		goto ErrorPointerNull;
	} else if (size.width < 0 || size.height < 0) {
		ErrorValueName = "size";
		goto ErrorIncorrectValue;
	}

	if (Type == NULL || strcmp(Type, "Normal") == 0) {
		/* use four pixels differential filter */
		size_f.width = 2;
		size_f.height = 2;
		Dx = DiffFilter_x;
		Dy = DiffFilter_y;
	} else if (strcmp(Type, "Sobel") == 0) {
		/* use Sobel filter */
		size_f.width = 3;
		size_f.height = 3;
		Dx = SobelFilter_x;
		Dy = SobelFilter_y;
	} else {
		fprintf(stderr, "*** %s() error - Filter options incorrect '%s' ***\n", FunctionName, Type);
		return NULL;
	}

	Image_dx = Filterer(Image, size, Dx, size_f, MEANINGFUL_FALSE);
	Image_dy = Filterer(Image, size, Dy, size_f, MEANINGFUL_FALSE);

	if ((Derivative = (VECTOR_2D *)calloc((size_t)size.width * size.height, sizeof(VECTOR_2D))) == NULL) {
		ErrorFunctionName = "calloc";
		ErrorValueName = "Derivative";
		goto ErrorMalloc;
	}
	for (x = 0; x < size.width * size.height; x++) {
		Derivative[x].x = Image_dx[x];
		Derivative[x].y = Image_dy[x];
	}
	free(Image_dx);
	free(Image_dy);
	return Derivative;
/* Errors */
ErrorMalloc:
	fprintf(stderr, "*** %s() error - Cannot allocate memory for (*%s) by %s() ***\n", FunctionName, ErrorValueName, ErrorFunctionName);
	goto ErrorReturn;
ErrorPointerNull:
	fprintf(stderr, "*** %s() error - The pointer (*%s) is NULL ***\n", FunctionName, ErrorValueName);
	goto ErrorReturn;
ErrorIncorrectValue:
	fprintf(stderr, "*** %s() error - The variable (%s) has incorrect value ***\n", FunctionName, ErrorValueName);
ErrorReturn:
	free(Derivative);
	free(Image_dy);
	free(Image_dx);
	return NULL;
}


double *
Derivation_abs(VECTOR_2D *Derivative_2D, SIZE size)
{
	char *FunctionName = "Derivation_abs";
	char *ErrorFunctionName = "";
	char *ErrorValueName = "";

	double *Derivative = NULL;
	int x;

	if (Derivative_2D == NULL) {
		ErrorValueName = "Derivative_2D";
		goto ErrorPointerNull;
	}
	if ((Derivative = (double *)calloc((size_t)size.width * size.height, sizeof(double))) == NULL) {
		ErrorFunctionName = "calloc";
		ErrorValueName = "Derivative";
		goto ErrorMalloc;
	}

	for (x = 0; x < size.width * size.height; x++) {
		Derivative[x] = sqrt(POW2(Derivative_2D[x].x) + POW2(Derivative_2D[x].y));
	}
	return Derivative;
/* Error */
ErrorMalloc:
	fprintf(stderr, "*** %s error - Cannot allocate memory for (*%s) by %s() ***\n", FunctionName, ErrorValueName, ErrorFunctionName);
	return NULL;
ErrorPointerNull:
	fprintf(stderr, "*** %s error - The pointer (*%s) is NULL ***\n", FunctionName, ErrorValueName);
	return NULL;
}


double *
Filterer(double *Image, SIZE size, double *Filter, SIZE size_f, int Mirroring)
{
	char *FunctionName = "Filterer";
	char *ErrorValueName = "";

	double *Filtered = NULL;
	SIZE center;
	int x, y;
	int m, n;
	int xx, yy;
	double sum;

	if (Image == NULL) {
		ErrorValueName = "Image";
		goto ErrorPointerNull;
	} else if (Filter == NULL) {
		ErrorValueName = "Filter";
		goto ErrorPointerNull;
	} else if (size.width < 0 || size.height < 0) {
		ErrorValueName = "size";
		goto ErrorIncorrectValue;
	} else if (size_f.width < 0 || size_f.height < 0) {
		ErrorValueName = "size_f";
		goto ErrorIncorrectValue;
	}

	if ((Filtered = (double *)calloc(size.width * size.height, sizeof(double))) == NULL) {
		fprintf(stderr, "*** %s() error - Cannot allocate memory for (*Filtered) by calloc() ***\n", FunctionName);
		return NULL;
	}
	center.width = floor(size_f.width / 2.0);
	center.height = floor(size_f.height / 2.0);
#pragma omp parallel for private(x, m, n, xx, yy, sum)
	for (y = 0; y < size.height; y++) {
		for (x = 0; x < size.width; x++) {
			sum = .0;
			for (m = 0; m < size_f.height; m++) {
				if (Mirroring == MEANINGFUL_TRUE) {
					yy = IndexOfMirroring(y + center.height - m, size.height);
				} else {
					yy = y + center.height - m;
				}
				for (n = 0; n < size_f.width; n++) {
					if (Mirroring == MEANINGFUL_TRUE) {
						xx = IndexOfMirroring(x + center.width - n, size.width);
					} else {
						xx = x + center.width - n;
					}
					if (0 <= xx && xx < size.width && 0 <= yy && yy < size.height) {
						sum += Image[size.width * yy + xx] * Filter[size_f.width * m + n];
					}
				}
			}
			Filtered[size.width * y + x] = sum;
		}
	}
	return Filtered;
/* Error */
ErrorPointerNull:
	fprintf(stderr, "*** %s() error - The pointer (*%s) is NULL ***\n", FunctionName, ErrorValueName);
	goto ErrorReturn;
ErrorIncorrectValue:
	fprintf(stderr, "*** %s() error - The variable (%s) has incorrect value ***\n", FunctionName, ErrorValueName);
ErrorReturn:
	free(Filtered);
	return NULL;
}


int
IndexOfMirroring(int x, int size)
{
	int index = 0;
	int x_mod = 0;

	if (x < 0) {
		x_mod = -x % (2 * size);
		if (x_mod < size) {
			index = x_mod;
		} else {
			index = 2 * size - 1 - x_mod;
		}
	} else if (x >= size) {
		x_mod = x % (2 * size);
		if (x_mod < size) {
			index = size - 1 - x_mod;
		} else {
			index = x_mod - size;
		}
	}
	return index;
}

