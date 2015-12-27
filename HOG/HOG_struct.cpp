#include <cstdio>
#include <new>
#include "HOG_struct.h"
#include "../Scratch_MeaningfulMotion.h"




HOG::HOG(void)
{
	orient_signed = false;
	bins = 0;
	width = 0;
	height = 0;
	hist = nullptr;
}

HOG::HOG(const HOG& copy)
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
			std::cerr << bad.what() << std::endl;
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

HOG::HOG(const bool init_signed, const int init_width, const int init_height, const int init_bins)
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
			std::cerr << bad.what() << std::endl;
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
HOG::copy(const HOG& copy)
{
	Histogram *tmp_hist = nullptr;

	if (this == &copy) {
		return *this;
	}
	try {
		tmp_hist = new Histogram[copy.width * copy.height];
	}
	catch (const std::bad_alloc &bad) {
		std::cerr << bad.what() << std::endl;
		fprintf(stderr, "HOG::HOG(const HOG &) error : Cannot allocate memory\n");
		throw;
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
HOG::reset(const bool init_signed, const int init_width, const int init_height, const int init_bins)
{
	Histogram *tmp_hist = nullptr;

	try {
		tmp_hist = new Histogram[init_width * init_height];
	}
	catch (const std::bad_alloc &bad) {
		std::cerr << bad.what() << std::endl;
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
HOG::setSign(const bool set_signed)
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
HOG::Hist(const int x, const int y, const int bin) const
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
HOG::Data(const int x, const int y) const
{
	return &(hist[width * y + x]);
}

bool
HOG::AddHist(const int x, const int y, const int bin, const double val)
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



#define HOG_PARAM_Bins 16
#define HOG_PARAM_Dense true
#define HOG_PARAM_SignedOrient true
HOG_PARAM::HOG_PARAM(void)
{
	Bins = HOG_PARAM_Bins;
	Dense = HOG_PARAM_Dense;
	SignedOrient= HOG_PARAM_SignedOrient;
}

void
HOG_PARAM::set_default(const char* name)
{
	if (strcmp(name, "Bins") == 0) {
		Bins = HOG_PARAM_Bins;
	} else if (strcmp(name, "Dense") == 0) {
		Dense = HOG_PARAM_Dense;
	} else if (strcmp(name, "SignedOrient") == 0) {
		SignedOrient = HOG_PARAM_SignedOrient;
	} else {
		fprintf(stderr, "*** HOG_PARAM::set_default() error - There are NOT such a parameter '%s' ***\n", name);
	}
}

void
HOG_PARAM::set_value(const char* name, const void* value)
{
	if (strcmp(name, "Bins") == 0) {
		Bins = *static_cast<const int*>(value);
		if (Bins < 1) {
			Bins = 1;
		}
	} else if (strcmp(name, "Dense") == 0) {
		Dense = *static_cast<const bool*>(value);
	} else if (strcmp(name, "SignedOrient") == 0) {
		SignedOrient = *static_cast<const bool*>(value);
	} else {
		fprintf(stderr, "*** HOG_PARAM::set_default() error - There are NOT such a parameter '%s' ***\n", name);
	}
}

