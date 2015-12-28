#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "Class.h"




ERROR::ERROR(void)
    : FunctionName("[undefined]")
{
}

ERROR::ERROR(const char *name)
    : FunctionName(name)
{
}

const char *
ERROR::OutputFunctionName(void)
{
	return FunctionName.c_str();
}

void
ERROR::Function(const char *name)
{
	ErrorFunctionName = name;
}

void
ERROR::Value(const char *name)
{
	ValueName = name;
}

void
ERROR::File(const char *name)
{
	FileName = name;
}

void
ERROR::Others(const char *error)
{
	fprintf(stderr, "*** %s() error - %s ***\n", FunctionName.c_str(), error);
}

void
ERROR::OthersWarning(const char *error)
{
	fprintf(stderr, "*** %s() warning - %s ***\n", FunctionName.c_str(), error);
}

void
ERROR::Malloc(void)
{
	fprintf(stderr, "*** %s() error - Cannot allocate memory for (*%s) ***\n", FunctionName.c_str(), ValueName.c_str());
}

void
ERROR::FunctionFail(void)
{
	fprintf(stderr, "*** %s() error - %s() failed to compute (%s) ***\n", FunctionName.c_str(), ErrorFunctionName.c_str(), ValueName.c_str());
}

void
ERROR::PointerNull(void)
{
	fprintf(stderr, "*** %s() error - The pointer (*%s) is NULL ***\n", FunctionName.c_str(), ValueName.c_str());
}

void
ERROR::ValueIncorrect(void)
{
	fprintf(stderr, "*** %s() error - The value (%s) is invalid value ***\n", FunctionName.c_str(), ValueName.c_str());
}

void
ERROR::ImageSize(void)
{
	fprintf(stderr, "*** %s() error - The size of image is varied from First Frame ***\n", FunctionName.c_str());
}

void
ERROR::FileRead(void)
{
	fprintf(stderr, "*** %s() error - Failed to read the file \"%s\" by %s() ***\n", FunctionName.c_str(), FileName.c_str(), ErrorFunctionName.c_str());
}

void
ERROR::FileWrite(void)
{
	fprintf(stderr, "*** %s error - Failed to write the file \"%s\" by %s() ***\n", FunctionName.c_str(), FileName.c_str(), ErrorFunctionName.c_str());
}




ATAN2_DIV_PI::ATAN2_DIV_PI(void)
{
	width = 0;
	height = 0;
	table = nullptr;
}

ATAN2_DIV_PI::ATAN2_DIV_PI(const ATAN2_DIV_PI &copy)
{
	const double *data = nullptr;

	width = copy.Width();
	height = copy.Height();
	try {
		table = new double[width * height];
	}
	catch (const std::bad_alloc& bad) {
		std::cerr << bad.what() << std::endl
		    << "error : ATAN2_DIV_PI::ATAN2_DIV_PI(const ATAN2_DIV_PI&) : Cannot allocate memory" << std::endl;
		throw;
	}
	data = copy.Data();
	for (int n = 0; n < width * height; n++) {
		table[n] = data[n];
	}
}

ATAN2_DIV_PI::ATAN2_DIV_PI(int W, int H)
{
	width = W;
	height = H;
	try {
		table = new double[width * height];
	}
	catch (const std::bad_alloc& bad) {
		std::cerr << bad.what() << std::endl
		    << "ATAN2_DIV_PI::ATAN2_DIV_PI(int, int) : Cannot allocate memory" << std::endl;
		throw;
	}
	for (int m = 0; m < height; m++) {
		for (int n = 0; n < width; n++) {
			table[width * m + n] = atan2(double(m), double(n)) / M_PI;
		}
	}
}

ATAN2_DIV_PI::~ATAN2_DIV_PI(void)
{
	width = 0;
	height = 0;
	delete[] table;
	table = nullptr;
}

int
ATAN2_DIV_PI::Width(void) const
{
	return width;
}

int
ATAN2_DIV_PI::Height(void) const
{
	return height;
}

const double *
ATAN2_DIV_PI::Data(void) const
{
	return table;
}

void
ATAN2_DIV_PI::reset(int W, int H)
{
	width = W;
	height = H;
	if (table != nullptr) {
		delete[] table;
	}
	try {
		table = new double[width * height];
	}
	catch (const std::bad_alloc& bad) {
		std::cerr << bad.what() << std::endl
		    << "error : void ATAN2_DIV_PI::reset(int, int) : Cannot allocate memory" << std::endl;
		throw;
	}
	for (int m = 0; m < height; m++) {
		for (int n = 0; n < width; n++) {
			table[width * m + n] = atan2(double(m), double(n)) / M_PI;
		}
	}
}

double
ATAN2_DIV_PI::val(int y, int x) const
{
	if (table == nullptr) {
		std::cerr << "error : double ATAN2_DIV_PI::val(int, int) : table is NULL" << std::endl;
		return 0.0;
	}

	double angle;
	if (abs(x) > width || abs(y) > height) {
		angle = atan2(double(y), double(x)) / M_PI;
	} else if (x == 0 && y == 0) {
		angle = 0.0;
	} else if (x == 0) {
		angle = y > 0 ? 0.5 : -0.5;
	} else if (y == 0) {
		angle = x > 0 ? 0.0 : 1.0;
	} else {
		if (x * y > 0) {
			angle = table[width * abs(y) + abs(x)];
			if (x < 0) {
				angle -= 1.0;
			}
		} else if (x < 0) {
			angle = 1.0 - table[width * abs(y) + abs(x)];
		} else {
			angle = -table[width * abs(y) + abs(x)];
		}
	}
	return angle;
}

