#include <cstdio>
#include "ImgStruct.h"




SIZE::SIZE(void)
{
	width = 0;
	height = 0;
}

SIZE::SIZE(int w, int h)
{
	width = w;
	height = h;
}

void
SIZE::reset(void)
{
	width = 0;
	height = 0;
}

void
SIZE::set_size(const SIZE &size)
{
	width = size.width;
	height = size.height;
}

void
SIZE::set_size(const SIZE *size)
{
	if (size != nullptr) {
		width = size->width;
		height = size->height;
	}
}



COORDINATE::COORDINATE(void)
{
	x = 0;
	y = 0;
}

COORDINATE::COORDINATE(int init_x, int init_y)
{
	x = init_x;
	y = init_y;
}



COORDINATE_3D::COORDINATE_3D(void)
{
	x = .0;
	y = .0;
	z = .0;
}

COORDINATE_3D::COORDINATE_3D(double init_x, double init_y, double init_z)
{
	x = init_x;
	y = init_y;
	z = init_z;
}

void
COORDINATE_3D::set(double set_x, double set_y, double set_z)
{
	x = set_x;
	y = set_y;
	z = set_z;
}



VECTOR_2D::VECTOR_2D(void)
{
	x = .0;
	y = .0;
}

VECTOR_2D::VECTOR_2D(double init_x, double init_y)
{
	x = init_x;
	y = init_y;
}

void
VECTOR_2D::reset(void)
{
	x = .0;
	y = .0;
}



Histogram::Histogram(void)
{
	_bins = 0;
	_hist = nullptr;
}

Histogram::Histogram(const Histogram &copy) // copy constructor
{
	_bins = 0;
	_hist = nullptr;
	try {
		_hist = new double[_bins];
	}
	catch (const std::bad_alloc &bad) {
		fprintf(stderr, "Histogram::Histogram(const Histogram &) error : Cannot allocate memory\n");
		_hist = nullptr;
		return;
	}
	_bins = copy._bins;
	for (int i = 0; i < _bins; i++) {
		_hist[i] = copy._hist[i];
	}
}

Histogram::Histogram(int init_bins)
{
	_bins = 0;
	try {
		_hist = new double[_bins];
	}
	catch (const std::bad_alloc &bad) {
		fprintf(stderr, "Histogram::copy(const Histogram &) error : Cannot allocate memory\n");
		_hist = nullptr;
		return;
	}
	_bins = init_bins;
	for (int i = 0; i < _bins; i++) {
		_hist[i] = .0;
	}
}

Histogram::~Histogram(void)
{
	delete[] _hist;
}

void
Histogram::free(void)
{
	_bins = 0;
	delete[] _hist;
	_hist = nullptr;
}

Histogram &
Histogram::copy(const Histogram &copy)
{
	if (this != &copy) {
		double *tmp_hist = nullptr;
		try {
			tmp_hist = new double[_bins];
		}
		catch (const std::bad_alloc &bad) {
			fprintf(stderr, "Histogram::copy(const Histogram &) error : Cannot allocate memory\n");
			return *this;
		}
		delete[] _hist;
		_hist = tmp_hist;
		_bins = copy._bins;
		for (int i = 0; i < _bins; i++) {
			_hist[i] = copy._hist[i];
		}
	}
	return *this;
}

Histogram &
Histogram::reset(int init_bins)
{
	if (init_bins >= 0) {
		double *tmp_hist;
		try {
			tmp_hist = new double[_bins];
		}
		catch (const std::bad_alloc &bad) {
			fprintf(stderr, "Histogram::reset(int) error : Cannot allocate memory\n");
			return *this;
		}
		delete[] _hist;
		_hist = tmp_hist;
		_bins = init_bins;
		for (int i = 0; i < _bins; i++) {
			_hist[i] = .0;
		}
	}
	return *this;
}

const double *
Histogram::data(void) const
{
	return _hist;
}

int
Histogram::bins(void) const
{
	return _bins;
}

double
Histogram::get(int bin) const
{
	if (bin < 0 || _bins <= bin) {
		return 0;
	}
	return _hist[bin];
}

bool
Histogram::add(int bin, double val)
{
	if (bin < 0 || _bins <= bin) {
		return false;
	}
	_hist[bin] += val;
	return true;
}

