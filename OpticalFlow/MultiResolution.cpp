#include "MultiResolution.h"
#include "../Scratch_MeaningfulMotion.h"

#define DEBUG_PYRAMID




ImgVector<double> *
Pyramider(const ImgVector<double> *img, const int MaxLevel)
{
#define WEIGHTED_FILTER_SIZE 5
	ERROR Error("Pyramider");

	ImgVector<double> *Pyramid = nullptr;
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
		Pyramid = new ImgVector<double>[MaxLevel + 1];
	}
	catch (const std::bad_alloc &bad) {
		std::cerr << bad.what() << std::endl;
		Error.Value("Pyramid");
		Error.Malloc();
		goto ExitError;
	}
	Pyramid[0].reset(img->width(), img->height(), img->data()); // MaxLevel 0 is equal to original image
	for (l = 1; l <= MaxLevel; l++) {
		size_l.width = int(ceil(img->width() * pow_int(0.5, l)));
		size_l.height = int(ceil(img->height() * pow_int(0.5, l)));
		if (size_l.width <= 0 || size_l.height <= 0) {
			Error.OthersWarning("The image reaches minimum size in Pyramid");
			break;
		}
		Pyramid[l].reset(size_l.width, size_l.height);
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

	// Make image pyramid (l > 0)
	for (l = 1; l <= MaxLevel; l++) {
		size_lm1 = size_l;
		size_l.width = int(ceil(img->width() * pow_int(0.5, l)));
		size_l.height = int(ceil(img->height() * pow_int(0.5, l)));
		if (size_l.width <= 0 || size_l.height <= 0) {
			Error.OthersWarning("The image reaches minimum size in Pyramid");
			break;
		}
		for (x = 0; x < size_l.width; x++) {
			for (y = 0; y < size_l.height; y++) {
				Pyramid[l].at(x, y) = .0;
				for (m = 0; m < WEIGHTED_FILTER_SIZE; m++) {
					ym = 2 * y + m - int(floor(WEIGHTED_FILTER_SIZE / 2.0));
					for (n = 0; n < WEIGHTED_FILTER_SIZE; n++) {
						xn = 2 * x + n - int(floor(WEIGHTED_FILTER_SIZE / 2.0));
						Pyramid[l].at(x, y) += w[m] * w[n] * Pyramid[l - 1].get_mirror(xn, ym);
					}
				}
			}
		}
	}
#if defined(DEBUG_PYRAMID)
	// Test Output
	for (l = 0; l <= MaxLevel; l++) {
		pnm.copy(PORTABLE_GRAYMAP_BINARY, Pyramid[l].width(), Pyramid[l].height(), 255, Pyramid[l].data(), 256.0);
		sprintf(filename, "Pyramid_%04d.pgm", l);
		pnm.write(filename);
		pnm.free();
	}
#endif
	return Pyramid;
// Error
ExitError:
	delete Pyramid;
	return nullptr;
}


ImgVector<VECTOR_2D<double> > *
grad_Pyramid(const ImgVector<double> *img_t_levels, const ImgVector<double> *img_tp1_levels, const int MaxLevel)
{
	ERROR Error("grad_Pyramid");

	// Note that convolution invert *Filter
	ImgVector<VECTOR_2D<double> > *grad_levels = nullptr;
	int x, y;
	int m, n;
	int l;

	// img_tp1_levels is allowed to be NULL
	if (img_t_levels == nullptr) {
		Error.Value("img_t_levels");
		Error.PointerNull();
		goto ExitError;
	}
	try {
		grad_levels = new ImgVector<VECTOR_2D<double> >[MaxLevel + 1];
	}
	catch (const std::bad_alloc &bad) {
		std::cerr << bad.what() << std::endl;
		Error.Value("grad_levels");
		Error.Malloc();
		goto ExitError;
	}
	for (l = 0; l <= MaxLevel; l++) {
		grad_levels[l].reset(img_t_levels[l].width(), img_t_levels[l].height());
		for (m = 0; m < img_t_levels[l].height(); m++) {
			y = SATURATE(m, 0, img_t_levels[l].height() - 2);
			for (n = 0; n < img_t_levels[l].width(); n++) {
				x = SATURATE(n, 0, img_t_levels[l].width() - 2);
				// dx
				grad_levels[l].at(n, m).x =
				    (img_t_levels[l].get(x + 1, y) - img_t_levels[l].get(x, y)
				    + img_t_levels[l].get(x + 1, y + 1) - img_t_levels[l].get(x, y + 1))
				    / 2.0;
				// dy
				grad_levels[l].at(n, m).y =
				    (img_t_levels[l].get(x, y + 1) - img_t_levels[l].get(x, y)
				    + img_t_levels[l].get(x + 1, y + 1) - img_t_levels[l].get(x + 1, y))
				    / 2.0;
				if (img_tp1_levels != NULL) {
					// dx
					grad_levels[l].at(n, m).x +=
					    (img_tp1_levels[l].get(x + 1, y) - img_tp1_levels[l].get(x, y)
					    + img_tp1_levels[l].get(x + 1, y + 1) - img_tp1_levels[l].get(x, y + 1))
					    / 2.0;
					// dy
					grad_levels[l].at(n, m).y +=
					    (img_tp1_levels[l].get(x, y + 1) - img_tp1_levels[l].get(x, y)
					    + img_tp1_levels[l].get(x + 1, y + 1) - img_tp1_levels[l].get(x + 1, y))
					    / 2.0;
				}
			}
		}
	}
	return grad_levels;
// Error
ExitError:
	delete grad_levels;
	return nullptr;
}


ImgVector<double> *
dt_Pyramid(const ImgVector<double> *img_t_levels, const ImgVector<double> *img_tp1_levels, const int MaxLevel)
{
	ERROR Error("grad_Pyramid");

	// Note that convolution invert *Filter
	ImgVector<double> *dt_levels = nullptr;
	int x, y;
	int m, n;
	int l;

	if (img_t_levels == nullptr) {
		Error.Value("img_t_levels");
		Error.PointerNull();
		goto ExitError;
	} else if (img_tp1_levels == nullptr) {
		Error.Value("img_tp1_levels");
		Error.PointerNull();
		goto ExitError;
	}
	try {
		dt_levels = new ImgVector<double>[MaxLevel + 1];
	}
	catch (const std::bad_alloc &bad) {
		std::cerr << bad.what() << std::endl;
		Error.Value("dt_levels");
		Error.Malloc();
		goto ExitError;
	}
	for (l = 0; l <= MaxLevel; l++) {
		dt_levels[l].reset(img_t_levels[l].width(), img_t_levels[l].height());
		for (m = 0; m < img_t_levels[l].height(); m++) {
			y = SATURATE(m, 0, img_t_levels[l].height() - 2);
			for (n = 0; n < img_t_levels[l].width(); n++) {
				x = SATURATE(n, 0, img_t_levels[l].width() - 2);
				// dt
				dt_levels[l].at(n, m) =
				    (img_tp1_levels[l].get(x, y) - img_t_levels[l].get(x, y)
				    + img_tp1_levels[l].get(x + 1, y) - img_t_levels[l].get(x + 1, y)
				    + img_tp1_levels[l].get(x, y + 1) - img_t_levels[l].get(x, y + 1)
				    + img_tp1_levels[l].get(x + 1, y + 1) - img_t_levels[l].get(x + 1, y + 1))
				    / 4.0;
			}
		}
	}
	return dt_levels;
// Error
ExitError:
	delete dt_levels;
	return nullptr;
}

