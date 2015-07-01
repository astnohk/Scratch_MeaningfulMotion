#include "Scratch_MeaningfulMotion.h"
#include "HOG.h"




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
HOG::Data(void) const
{
	return hist;
}

int
HOG::Hist(int bin) const
{
	if (bin < 0 || bins <= bin) {
		return -1;
	}
	return hist[bin];
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

