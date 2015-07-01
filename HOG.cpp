#include "Scratch_MeaningfulMotion.h"
#include "HOG.h"


HOG *
HisotgramOfGradient(const PNM_DOUBLE &Img, SIZE size)
{
	ERROR Error("HistogramOfGradient");
	HOG *hog = nullptr;
	double *orient = nullptr;
	SIZE cell(3, 3);
	int bins = 16;

	orient = Orientation(Img, HOG_ORIENTATION_SIGNED);
	if (orient == nullptr) {
		Error.Function("Orientation");
		Error.Value("orient");
		Error.FunctionFail();
		goto ExitError;
	}
	hog = ComputeHistogramOfGradient(orient, size, cell, bins, HOG_ORIENTATION_SIGNED);
	if (hog == nullptr) {
		Error.Function("ComputeHistogramOfGradient");
		Error.Value("hog");
		Error.FunctionFail();
		goto ExitError;
	}
	return hog;
// Error
ExitError:
	delete[] orient;
	delete[] hog;
	return nullptr;
}


double *
Orientation(const PNM_DOUBLE &Img, bool sign)
{
	ERROR Error("Orientation");
	SIZE size;
	double *orient = nullptr;
	VECTOR_2D *grad = nullptr;
	int x, y;

	size.width = Img.Width();
	size.height = Img.Height();
	try {
		grad = new VECTOR_2D[size.width * size.height];
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("grad");
		Error.Malloc();
		goto ExitError;
	}
	for (y = 0; y < size.height; y++) {
		for (x = 0; x < size.width; x++) {
			grad[size.width * y + x].x = Img.Image(x + 1, y) - Img.Image(x - 1, y);
			grad[size.width * y + x].y = Img.Image(x, y + 1) - Img.Image(x, y - 1);
		}
	}
	try {
		orient = new double[size.width * size.height];
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("orient");
		Error.Malloc();
		goto ExitError;
	}
	for (int i = 0; i < size.width * size.height; i++) {
		switch (sign) {
			case HOG_ORIENTATION_SIGNED:
				orient[i] = atan2(grad[i].y, grad[i].x) / M_PI;
				break;
			case HOG_ORIENTATION_UNSIGNED:
			default:
				orient[i] = atan2(grad[i].y, grad[i].x) / M_PI;
				if (orient[i] < .0) {
					orient[i] += 1.0;
				}
		}
	}
	delete[] grad;
	grad = nullptr;
	return orient;
// Error
ExitError:
	delete[] grad;
	delete[] orient;
	return nullptr;
}

HOG *
ComputeHistogramOfGradient(const double *orient, SIZE size, SIZE cell, int bins, bool sign)
{
	ERROR Error("ComputeHistogramOfGradient");
	SIZE Cells;
	HOG *hog = nullptr;
	int x, y;
	int m, n;

	Cells.width = (int)ceil(size.width / cell.width);
	Cells.height = (int)ceil(size.height / cell.height);
	try {
		hog = new HOG[Cells.width * Cells.height];
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("hog");
		Error.Malloc();
		goto ExitError;
	}
	for (y = 0; y < Cells.height; y++) {
		int Y = cell.height * y;
		for (x = 0; x < Cells.width; x++) {
			int X = cell.width * x;
			hog[Cells.width * y + x].reset(sign, bins); // reset histogram with number of bins
			for (m = 0; m < cell.height; m++) {
				for (n = 0; n < cell.width; n++) {
					int dir;
					dir = (int)floor(bins * orient[size.width * (Y + m) + X + n]);
					if (dir == bins) {
						dir = dir - 1;
					}
					hog[Cells.width * y + x].AddHist(dir);
				}
			}
		}
	}
	return hog;
// Error
ExitError:
	delete[] hog;
	return nullptr;
}


HOG::HOG(void)
{
	orient_signed = false;
	bins = 0;
	hist = nullptr;
}

HOG::HOG(bool init_signed, int num_bins)
{
	ERROR Error("HOG::HOG(int)");

	if (init_signed == false) {
		orient_signed = false;
	} else {
		orient_signed = true;
	}
	bins = num_bins;
	try {
		hist = new int[bins];
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("hist");
		Error.Malloc();
		hist = nullptr;
		bins = 0;
		return;
	}
	for (int i = 0; i < bins; i++) {
		hist[i] = 0;
	}
}

HOG::~HOG(void)
{
	orient_signed = false;
	bins = 0;
	delete[] hist;
	hist = nullptr;
}

bool
HOG::reset(bool init_signed, int num_bins)
{
	ERROR Error("HOG::reset");

	delete[] hist;
	hist = nullptr;
	if (init_signed == false) {
		orient_signed = false;
	} else {
		orient_signed = true;
	}
	bins = num_bins;
	try {
		hist = new int[bins];
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("hist");
		Error.Malloc();
		hist = nullptr;
		bins = 0;
		return false;
	}
	for (int i = 0; i < bins; i++) {
		hist[i] = 0;
	}
	return true;
}

void
HOG::setSign(bool init_signed)
{
	if (init_signed == false) {
		orient_signed = false;
	} else {
		orient_signed = true;
	}
}

bool
HOG::Signed(void) const
{
	return orient_signed;
}

int
HOG::Bins(void) const
{
	return bins;
}

const int *
HOG::Hist(void) const
{
	return hist;
}

bool
HOG::AddHist(int bin)
{
	if (bin < 0 || bins <= bin) {
		return false;
	}
	hist[bin] += 1;
	return true;
}

bool
HOG::SubHist(int bin)
{
	if (bin < 0 || bins <= bin) {
		return false;
	}
	if (hist[bin] > 0) {
		hist[bin] -= 1;
	}
	return true;
}

