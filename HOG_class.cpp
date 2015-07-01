#include "Scratch_MeaningfulMotion.h"
#include "HOG.h"




Histogram::Histogram(void)
{
	bins = 0;
	hist = nullptr;
}

Histogram::Histogram(int init_bins)
{
	ERROR Error("Histogram::Histogram");

	bins = init_bins;
	try {
		hist = new double[bins];
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("hist");
		Error.Malloc();
		bins = 0;
		return;
	}
	for (int i = 0; i < bins; i++) {
		hist[i] = .0;
	}
}

Histogram::~Histogram(void)
{
	bins = 0;
	delete[] hist;
	hist = nullptr;
}

bool
Histogram::reset(int init_bins)
{
	ERROR Error("Histogram::reset");

	bins = init_bins;
	try {
		hist = new double[bins];
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("hist");
		Error.Malloc();
		bins = 0;
		return false;
	}
	for (int i = 0; i < bins; i++) {
		hist[i] = .0;
	}
	return true;
}

int
Histogram::Bins(void) const
{
	return bins;
}

double
Histogram::Hist(int bin) const
{
	if (bin < 0 || bins <= bin) {
		return -2;
	}
	return hist[bin];
}

const double *
Histogram::Data(void) const
{
	return hist;
}

bool
Histogram::Add(int bin, double val)
{
	if (bin < 0 || bins <= bin) {
		return false;
	}
	hist[bin] += val;
	return true;
}




HOG::HOG(void)
{
	orient_signed = false;
	bins = 0;
	width = 0;
	height = 0;
	hist = nullptr;
}

HOG::HOG(bool init_signed, int init_width, int init_height, int init_bins)
{
	ERROR Error("HOG::HOG(int)");

	if (init_signed == false) {
		orient_signed = false;
	} else {
		orient_signed = true;
	}
	bins = init_bins;
	width = init_width;
	height = init_height;
	try {
		hist = new Histogram[width * height];
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("hist");
		Error.Malloc();
		hist = nullptr;
		bins = 0;
		width = 0;
		height = 0;
		return;
	}
	for (int i = 0; i < width * height; i++) {
		if (hist[i].reset(bins) == false) {
			Error.Value("hist[i]");
			Error.Malloc();
			delete[] hist;
			hist = nullptr;
			bins = 0;
			width = 0;
			height = 0;
			return;
		}
	}
}

bool
HOG::reset(bool init_signed, int init_width, int init_height, int init_bins)
{
	ERROR Error("HOG::reset");

	delete[] hist;
	hist = nullptr;
	if (init_signed == false) {
		orient_signed = false;
	} else {
		orient_signed = true;
	}
	bins = init_bins;
	width = init_width;
	height = init_height;
	try {
		hist = new Histogram[width * height];
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("hist");
		Error.Malloc();
		hist = nullptr;
		bins = 0;
		width = 0;
		height = 0;
		return false;
	}
	for (int i = 0; i < width * height; i++) {
		if (hist[i].reset(bins) == false) {
			Error.Value("hist[i]");
			Error.Malloc();
			delete[] hist;
			hist = nullptr;
			bins = 0;
			width = 0;
			height = 0;
			return false;
		}
	}
	return true;
}

HOG::~HOG(void)
{
	orient_signed = false;
	bins = 0;
	width = 0;
	height = 0;
	delete[] hist;
	hist = nullptr;
}

void
HOG::free(void)
{
	orient_signed = false;
	bins = 0;
	width = 0;
	height = 0;
	delete[] hist;
	hist = nullptr;
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

int
HOG::Width(void) const
{
	return width;
}

int
HOG::Height(void) const
{
	return height;
}

double
HOG::Hist(int x, int y, int bin) const
{
	if (x < 0 || width <= x
	    || y < 0 || height <= y) {
		return -1.0;
	}
	return hist[width * y + x].Hist(bin);
}

const Histogram *
HOG::Data(void) const
{
	return hist;
}

bool
HOG::AddHist(int x, int y, int bin, double val)
{
	if (x < 0 || width <= x
	    || y < 0 || height <= y) {
		return false;
	}
	if (hist[width * y + x].Add(bin, val) == false) {
		return false;
	}
	return true;
}

