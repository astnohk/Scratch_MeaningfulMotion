#include "Scratch_MeaningfulMotion.h"

//#define SHOW_GAUSSIAN_FILTER


double
HorizontalMedian(double *img, int img_width, int x, int y, int width)
{
	ERROR Error("HorizontalMedian");

	double *array = nullptr;
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
	try {
		array = new double[m_e - m_s + 1];
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("array");
		Error.Malloc();
		return 0.0;
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
	delete[] array;
	array = nullptr;
	return tmp;
}


double*
EpsilonFilter(double *img, SIZE size, FILTER_PARAM Param)
{
	ERROR Error("EpsilonFilter");

	double *Epsilon = nullptr;
	int filter_width_2 = (int)floor(Param.size.width / 2.0);
	int filter_height_2 = (int)floor(Param.size.height / 2.0);
	double Div_coeff = 1.0 / (Param.size.width * Param.size.height);
	int center;
	int x, y;
	int fx, fy;
	int m, n;
	int tmp;

	if (Param.epsilon < .0) {
		Error.Value("epsilon");
		Error.ValueIncorrect();
		goto ExitError;
	} else if (Param.size.width <= 0 || Param.size.width % 2 == 0) {
		Error.Value("Param.size.width");
		Error.ValueIncorrect();
		goto ExitError;
	} else if (Param.size.height <= 0 || Param.size.height % 2 == 0) {
		Error.Value("Param.size.height");
		Error.ValueIncorrect();
		goto ExitError;
	} else if (img == nullptr) {
		Error.Value("img");
		Error.PointerNull();
		goto ExitError;
	}

	try {
		Epsilon = new double[size.width * size.height];
	}
	catch (const std::bad_alloc &bad) {
		Error.Function("new");
		Error.Value("Epsilon");
		Error.Malloc();
		goto ExitError;
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
// Error
ExitError:
	delete[] Epsilon;
	return nullptr;
}


double*
Gaussian(double *img, SIZE size, FILTER_PARAM Param)
{
	ERROR Error("Gaussian");

	double *blurred = nullptr;
	double *Gauss = nullptr;
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
	try {
		Gauss = new double[Param.size.width * Param.size.height];
	}
	catch (const std::bad_alloc &bad) {
		Error.Function("new");
		Error.Value("Gauss");
		Error.Malloc();
		goto ExitError;
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
	try {
		blurred = new double[size.height * size.width];
	}
	catch (const std::bad_alloc &bad) {
		Error.Function("new");
		Error.Value("blurred");
		Error.Malloc();
		goto ExitError;
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
	delete[] Gauss;
	Gauss = nullptr;
	return blurred;
// Errors
ExitError:
	delete[] blurred;
	delete[] Gauss;
	return nullptr;
}


double *
DerivativeAngler(double *img, SIZE size)
{
	ERROR Error("DerivativeAngler");

	VECTOR_2D *Derivative = nullptr;
	double *angles = nullptr;
	int n;
	double dx, dy;
 
	if (img == nullptr) {
		Error.Value("img");
		Error.PointerNull();
		goto ExitError;
	} else if (size.width < 0 || size.height < 0) {
		Error.Value("size");
		Error.ValueIncorrect();
		goto ExitError;
	}

	if ((Derivative = Derivator(img, size, "Sobel")) == nullptr) {
		Error.Function("Derivator");
		Error.Value("Derivative");
		Error.FunctionFail();
		goto ExitError;
	}
	try {
		angles = new double[size.height * size.width];
	}
	catch (const std::bad_alloc &bad) {
		Error.Function("new");
		Error.Value("angles");
		Error.Malloc();
		goto ExitError;
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
	delete[] Derivative;
	return angles;
// Error
ExitError:
	delete[] angles;
	return nullptr;
}


VECTOR_2D *
Derivator(double *Image, SIZE size, const char *Type)
{
	ERROR Error("Derivator");

	/* Note that convolution invert *Filter */
	double DiffFilter_x[4] = {0.5, -0.5, 0.5, -0.5};
	double DiffFilter_y[4] = {0.5, 0.5, -0.5, -0.5};
	double SobelFilter_x[9] = {0.25, 0, -0.25, 0.5, 0, -0.5, 0.25, 0, -0.25};
	double SobelFilter_y[9] = {0.25, 0.5, 0.25, 0, 0, 0, -0.25, -0.5, -0.25};
	double *Dx = nullptr;
	double *Dy = nullptr;
	SIZE size_f;
	double *Image_dx = nullptr;
	double *Image_dy = nullptr;
	VECTOR_2D *Derivative = nullptr;
	int x;
 
	if (Image == nullptr) {
		Error.Value("Image");
		Error.PointerNull();
		goto ExitError;
	} else if (size.width < 0 || size.height < 0) {
		Error.Value("size");
		Error.ValueIncorrect();
		goto ExitError;
	}

	if (Type == nullptr || strcmp(Type, "Normal") == 0) {
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
		Error.Others("Filter options incorrect");
		goto ExitError;
	}

	Image_dx = Filterer(Image, size, Dx, size_f, MEANINGFUL_FALSE);
	Image_dy = Filterer(Image, size, Dy, size_f, MEANINGFUL_FALSE);

	try {
		Derivative = new VECTOR_2D[size.width * size.height];
	}
	catch (const std::bad_alloc &bad) {
		Error.Function("new");
		Error.Value("Derivative");
		Error.Malloc();
		goto ExitError;
	}
	for (x = 0; x < size.width * size.height; x++) {
		Derivative[x].x = Image_dx[x];
		Derivative[x].y = Image_dy[x];
	}
	delete[] Image_dx;
	delete[] Image_dy;
	return Derivative;
// Error
ExitError:
	delete[] Derivative;
	delete[] Image_dy;
	delete[] Image_dx;
	return nullptr;
}


double *
Derivation_abs(VECTOR_2D *Derivative_2D, SIZE size)
{
	ERROR Error("Derivation_abs");

	double *Derivative = nullptr;
	int x;

	if (Derivative_2D == nullptr) {
		Error.Value("Derivative_2D");
		Error.PointerNull();
		goto ExitError;
	}
	try {
		Derivative = new double[size.width * size.height];
	}
	catch (const std::bad_alloc &bad) {
		Error.Function("new");
		Error.Value("Derivative");
		Error.Malloc();
		goto ExitError;
	}
	for (x = 0; x < size.width * size.height; x++) {
		Derivative[x] = sqrt(POW2(Derivative_2D[x].x) + POW2(Derivative_2D[x].y));
	}
	return Derivative;
// Error
ExitError:
	return nullptr;
}


double *
Filterer(double *Image, SIZE size, double *Filter, SIZE size_f, int Mirroring)
{
	ERROR Error("Filterer");

	double *Filtered = nullptr;
	SIZE center;
	int x, y;
	int m, n;
	int xx, yy;
	double sum;

	if (Image == nullptr) {
		Error.Value("Image");
		Error.PointerNull();
		goto ExitError;
	} else if (Filter == nullptr) {
		Error.Value("Filter");
		Error.PointerNull();
		goto ExitError;
	} else if (size.width < 0 || size.height < 0) {
		Error.Value("size");
		Error.ValueIncorrect();
		goto ExitError;
	} else if (size_f.width < 0 || size_f.height < 0) {
		Error.Value("size_f");
		Error.ValueIncorrect();
		goto ExitError;
	}

	try {
		Filtered = new double[size.width * size.height];
	}
	catch (const std::bad_alloc &bad) {
		Error.Function("new");
		Error.Value("Filtered");
		Error.Malloc();
		goto ExitError;
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
// Error
ExitError:
	delete[] Filtered;
	return nullptr;
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

