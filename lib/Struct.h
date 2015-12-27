#include "ImgStruct.h"

#ifndef LIB_Struct
#define LIB_Struct

struct FRAGMENT
{
	int start;
	int end;
	double Pr;
	FRAGMENT(void);
	FRAGMENT(const int s, const int e, const double& prob);
};

struct SEGMENT
{
	int n; // Start point x
	int m; // Start point y
	int x; // End point x
	int y; // End point y
	double Pr; // Probability
	SEGMENT(void);
	SEGMENT(const int stx, const int sty, const int endx, const int endy, const double& prob);
};

struct LINEPOLE
{
	double r;
	double theta;
	double cos;
	double sin;
	LINEPOLE(void);
	LINEPOLE(const double& rad, const double& th, const double& icos, const double& isin);
};

#endif

