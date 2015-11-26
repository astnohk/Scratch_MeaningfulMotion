#include "../Scratch_MeaningfulMotion.h"

//#define SHOW_GAUSSIAN_FILTER



double
HorizontalMedian(ImgVector<double> *img, int x, int y, int width)
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
	} else if (x + width / 2 >= img->width()) {
		m_s = -(width - 1) / 2;
		m_e = img->width() - 1 - x;
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
		array[m - m_s] = img->get(x + m, y);
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


ImgVector<double> *
EpsilonFilter(ImgVector<double> *img, FILTER_PARAM Param)
{
	ERROR Error("EpsilonFilter");

	ImgVector<double> *Epsilon = nullptr;
	int filter_width_2 = (int)floor(Param.size.width / 2.0);
	int filter_height_2 = (int)floor(Param.size.height / 2.0);
	double Div_coeff = 1.0 / (Param.size.width * Param.size.height);
	int center;
	int x, y;
	int fx, fy;

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
		Epsilon = new ImgVector<double>(img->width(), img->height());
	}
	catch (const std::bad_alloc &bad) {
		Error.Function("new");
		Error.Value("Epsilon");
		Error.Malloc();
		goto ExitError;
	}
#pragma omp parallel for private(x, fx, fy, center)
	for (y = 0; y < img->height(); y++) {
		for (x = 0; x < img->width(); x++) {
			center = img->get(x, y);
			Epsilon->at(x, y) = .0;
			for (fy = -filter_height_2; fy <= filter_height_2; fy++) {
				for (fx = -filter_width_2; fx <= filter_width_2; fx++) {
					if (fabs(center - img->get(x + fx, y + fy)) <= Param.epsilon) {
						Epsilon->at(x, y) += img->get_mirror(x + fx, y + fy);
					} else {
						Epsilon->at(x, y) += center;
					}
				}
			}
			Epsilon->at(x, y) *= Div_coeff;
		}
	}
	return Epsilon;
// Error
ExitError:
	delete Epsilon;
	return nullptr;
}


ImgVector<double> *
Gaussian(ImgVector<double> *img, FILTER_PARAM Param)
{
	ERROR Error("Gaussian");

	ImgVector<double> *Gauss = nullptr;
	ImgVector<double> *blurred = nullptr;
	int filter_width_2 = 0;
	int filter_height_2 = 0;
	int m, n, y, x;
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
		Gauss = new ImgVector<double>(Param.size.width, Param.size.height);
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
					Gauss->at(n, m) =
					    exp(
					    -(POW2(m - filter_height_2) + POW2(n - filter_width_2))
					    / (2.0 * POW2(Param.std_deviation))
					    );
					sum += Gauss->get(n, m);
				}
#if defined(SHOW_GAUSSIAN_FILTER)
				printf(" %.5f", Gauss->get(n, m));
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
				Gauss->at(n, m) =
				    exp(
				    -(POW2(m - filter_height_2) + POW2(n - filter_width_2))
				    / (2.0 * POW2(Param.std_deviation))
				    );
				sum += Gauss->get(n, m);
#if defined(SHOW_GAUSSIAN_FILTER)
				printf(" %.5f", Gauss->get(n, m));
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
		(*Gauss)[m] /= sum;
	}
	try {
		blurred = new ImgVector<double>(img->width(), img->height());
	}
	catch (const std::bad_alloc &bad) {
		Error.Function("new");
		Error.Value("blurred");
		Error.Malloc();
		goto ExitError;
	}
#pragma omp parallel for private(n, x, y)
	for (m = 0; m < blurred->height(); m++) {
		for (n = 0; n < blurred->width(); n++) {
			for (y = -filter_height_2; y <= filter_height_2; y++) {
				for (x = -filter_width_2; x <= filter_width_2; x++) {
					blurred->at(n, m) += img->get(n + x, m + y) * Gauss->get(x + filter_width_2, y + filter_height_2);
				}
			}
		}
	}
	delete Gauss;
	Gauss = nullptr;
	return blurred;
// Errors
ExitError:
	delete blurred;
	delete Gauss;
	return nullptr;
}


ImgVector<double> *
DerivativeAngler(ImgVector<double> *img)
{
	ERROR Error("DerivativeAngler");

	ImgVector<VECTOR_2D<double> > *Derivative = nullptr;
	ImgVector<double> *angles = nullptr;
	int n;
	double dx, dy;
 
	if (img == nullptr) {
		Error.Value("img");
		Error.PointerNull();
		goto ExitError;
	} else if (img->width() < 0 || img->height() < 0) {
		Error.Value("size");
		Error.ValueIncorrect();
		goto ExitError;
	}

	if ((Derivative = Derivator(img, "Sobel")) == nullptr) {
		Error.Function("Derivator");
		Error.Value("Derivative");
		Error.FunctionFail();
		goto ExitError;
	}
	try {
		angles = new ImgVector<double>(img->width(), img->height());
	}
	catch (const std::bad_alloc &bad) {
		Error.Function("new");
		Error.Value("angles");
		Error.Malloc();
		goto ExitError;
	}
	for (n = 0; n < img->size(); n++) {
		dx = (*Derivative)[n].x;
		dy = (*Derivative)[n].y;
		if (fabs(dx) <= DERIVATIVE_MINIMUM && fabs(dy) <= DERIVATIVE_MINIMUM) {
			(*angles)[n] = -2.0 * ANGLE_MAX;
		} else {
			(*angles)[n] = atan2(dy, dx) / M_PI + 0.5; // Add 0.5 for pi/2 rotation
			if ((*angles)[n] > ANGLE_MAX) {
				(*angles)[n] -= ANGLE_MAX;
			} else if ((*angles)[n] < 0.0) {
				(*angles)[n] += ANGLE_MAX;
			}
		}
	}
	delete Derivative;
	return angles;
// Error
ExitError:
	delete angles;
	return nullptr;
}


ImgVector<VECTOR_2D<double> > *
Derivator(ImgVector<double> *Image, const char *Type)
{
	ERROR Error("Derivator");

	ImgVector<VECTOR_2D<double> > *Derivative = nullptr;
	ImgVector<double> *Image_dx = nullptr;
	ImgVector<double> *Image_dy = nullptr;
	/* Note that convolution invert *Filter */
	double DiffFilter_array_x[4] = {-0.5, 0.5, -0.5, 0.5};
	double DiffFilter_array_y[4] = {-0.5, -0.5, 0.5, 0.5};
	double SobelFilter_array_x[9] = {-0.25, 0, 0.25, -0.5, 0, 0.5, -0.25, 0, 0.25};
	double SobelFilter_array_y[9] = {-0.25, -0.5, -0.25, 0, 0, 0, 0.25, 0.5, 0.25};
	ImgVector<double> DiffFilter_x(2, 2, DiffFilter_array_x);
	ImgVector<double> DiffFilter_y(2, 2, DiffFilter_array_y);
	ImgVector<double> SobelFilter_x(3, 3, SobelFilter_array_x);
	ImgVector<double> SobelFilter_y(3, 3, SobelFilter_array_y);
	ImgVector<double> *Dx = nullptr;
	ImgVector<double> *Dy = nullptr;
	int x;
 
	if (Image == nullptr) {
		Error.Value("Image");
		Error.PointerNull();
		goto ExitError;
	} else if (Image->width() < 0 || Image->height() < 0) {
		Error.Value("size");
		Error.ValueIncorrect();
		goto ExitError;
	}

	if (Type == nullptr || strcmp(Type, "Normal") == 0) {
		/* use four pixels differential filter */
		Dx = &DiffFilter_x;
		Dy = &DiffFilter_y;
	} else if (strcmp(Type, "Sobel") == 0) {
		/* use Sobel filter */
		Dx = &SobelFilter_x;
		Dy = &SobelFilter_y;
	} else {
		Error.Others("Filter options incorrect");
		goto ExitError;
	}

	Image_dx = Filterer(Image, Dx, false);
	Image_dy = Filterer(Image, Dy, false);

	try {
		Derivative = new ImgVector<VECTOR_2D<double> >(Image->width(), Image->height());
	}
	catch (const std::bad_alloc &bad) {
		Error.Function("new");
		Error.Value("Derivative");
		Error.Malloc();
		goto ExitError;
	}
	for (x = 0; x < Derivative->size(); x++) {
		(*Derivative)[x].x = Image_dx->get(x);
		(*Derivative)[x].y = Image_dy->get(x);
	}
	delete Image_dx;
	delete Image_dy;
	return Derivative;
// Error
ExitError:
	delete Derivative;
	delete Image_dy;
	delete Image_dx;
	return nullptr;
}


ImgVector<double> *
Derivation_abs(ImgVector<VECTOR_2D<double> > *Derivative_2D)
{
	ERROR Error("Derivation_abs");

	ImgVector<double> *Derivative = nullptr;
	int n;

	if (Derivative_2D == nullptr) {
		Error.Value("Derivative_2D");
		Error.PointerNull();
		goto ExitError;
	}
	try {
		Derivative = new ImgVector<double>(Derivative_2D->width(), Derivative_2D->height());
	}
	catch (const std::bad_alloc &bad) {
		Error.Function("new");
		Error.Value("Derivative");
		Error.Malloc();
		goto ExitError;
	}
	for (n = 0; n < Derivative_2D->size(); n++) {
		(*Derivative)[n] = sqrt(POW2((*Derivative_2D)[n].x) + POW2((*Derivative_2D)[n].y));
	}
	return Derivative;
// Error
ExitError:
	return nullptr;
}


ImgVector<double> *
Filterer(ImgVector<double> *Image, ImgVector<double> *Filter, bool Mirroring)
{
	ERROR Error("Filterer");

	ImgVector<double> *Filtered = nullptr;
	int center_x, center_y;
	int x, y;
	int m, n;
	double sum;

	if (Image == nullptr) {
		Error.Value("Image");
		Error.PointerNull();
		goto ExitError;
	} else if (Filter == nullptr) {
		Error.Value("Filter");
		Error.PointerNull();
		goto ExitError;
	} else if (Image->width() < 0 || Image->height() < 0) {
		Error.Value("size");
		Error.ValueIncorrect();
		goto ExitError;
	} else if (Filter->width() < 0 || Filter->height() < 0) {
		Error.Value("size_f");
		Error.ValueIncorrect();
		goto ExitError;
	}

	try {
		Filtered = new ImgVector<double>(Image->width(), Image->height());
	}
	catch (const std::bad_alloc &bad) {
		Error.Function("new");
		Error.Value("Filtered");
		Error.Malloc();
		goto ExitError;
	}
	center_x = Filter->width() / 2;
	center_y = Filter->height() / 2;
#pragma omp parallel for private(x, m, n, sum)
	for (y = 0; y < Image->height(); y++) {
		for (x = 0; x < Image->width(); x++) {
			sum = .0;
			for (m = 0; m < Filter->height(); m++) {
				for (n = 0; n < Filter->width(); n++) {
					if (Mirroring) {
						sum += Image->get_mirror(x + center_x - n, y + center_y - m) * Filter->get(n, m);
					} else {
						sum += Image->get_zeropad(x + center_x - n, y + center_y - m) * Filter->get(n, m);
					}
				}
			}
			Filtered->at(x, y) = sum;
		}
	}
	return Filtered;
// Error
ExitError:
	delete Filtered;
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

