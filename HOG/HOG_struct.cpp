#include <cstdio>
#include <new>
#include "HOG_struct.h"

#define nullptr 0




HOG::HOG(void)
{
	orient_signed = false;
	bins = 0;
	width = 0;
	height = 0;
	hist = nullptr;
}

HOG::HOG(const HOG &copy)
{
	orient_signed = false;
	bins = 0;
	width = 0;
	height = 0;
	hist = nullptr;
	if (copy.width > 0 && copy.height > 0) {
		Histogram *tmp_hist = nullptr;
		try {
			tmp_hist = new Histogram[copy.width * copy.height];
		}
		catch (const std::bad_alloc &bad) {
			fprintf(stderr, "HOG::HOG(const HOG &) error : Cannot allocate memory\n");
			throw;
		}
		hist = tmp_hist;
		for (int i = 0; i < copy.width * copy.height; i++) {
			hist[i].copy(copy.hist[i]);
		}
		orient_signed = copy.orient_signed;
		bins = copy.bins;
		width = copy.bins;
		height = copy.bins;
	}
}

HOG::HOG(bool init_signed, int init_width, int init_height, int init_bins)
{
	orient_signed = false;
	bins = 0;
	width = 0;
	height = 0;
	hist = nullptr;
	if (init_width > 0 && init_height > 0) {
		Histogram *tmp_hist = nullptr;
		try {
			tmp_hist = new Histogram[init_width * init_height];
		}
		catch (const std::bad_alloc &bad) {
			fprintf(stderr, "HOG::HOG(bool, int, int, int) error : Cannot allocate memory\n");
			throw;
		}
		for (int i = 0; i < init_width * init_height; i++) {
			tmp_hist[i].reset(init_bins);
		}
		orient_signed = init_signed;
		bins = init_bins;
		width = init_width;
		height = init_height;
		hist = tmp_hist;
	}
}

HOG &
HOG::copy(const HOG &copy)
{
	Histogram *tmp_hist = nullptr;

	if (this == &copy) {
		return *this;
	}
	try {
		tmp_hist = new Histogram[copy.width * copy.height];
	}
	catch (const std::bad_alloc &bad) {
		fprintf(stderr, "HOG::HOG(const HOG &) error : Cannot allocate memory\n");
		throw;
		return *this;
	}
	for (int i = 0; i < copy.width * copy.height; i++) {
		tmp_hist[i].copy(copy.hist[i]);
	}
	delete[] hist;
	hist = tmp_hist;
	orient_signed = copy.orient_signed;
	bins = copy.bins;
	width = copy.width;
	height = copy.height;
	return *this;
}

bool
HOG::reset(bool init_signed, int init_width, int init_height, int init_bins)
{
	Histogram *tmp_hist = nullptr;

	try {
		tmp_hist = new Histogram[init_width * init_height];
	}
	catch (const std::bad_alloc &bad) {
		fprintf(stderr, "HOG::reset(bool, int, int, int) error : Cannot allocate memory\n");
		return false;
	}
	for (int i = 0; i < init_width * init_height; i++) {
		tmp_hist[i].reset(init_bins);
	}
	delete[] hist;
	hist = tmp_hist;
	orient_signed = init_signed;
	bins = init_bins;
	width = init_width;
	height = init_height;
	return true;
}

HOG::~HOG(void)
{
	delete[] hist;
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
HOG::setSign(bool set_signed)
{
	orient_signed = set_signed;
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
	return hist[width * y + x].get(bin);
}

const Histogram *
HOG::Data(void) const
{
	return hist;
}

const Histogram *
HOG::Data(int x, int y) const
{
	return &(hist[width * y + x]);
}

bool
HOG::AddHist(int x, int y, int bin, double val)
{
	if (x < 0 || width <= x
	    || y < 0 || height <= y) {
		return false;
	}
	if (hist[width * y + x].add(bin, val) == false) {
		return false;
	}
	return true;
}

