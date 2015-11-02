#include "ImgStruct.h"



#ifndef LIB_Struct
#define LIB_Struct

struct FRAGMENT
{
	int start;
	int end;
	double Pr;
	FRAGMENT(void);
	FRAGMENT(int s, int e, double prob);
};

struct SEGMENT
{
	int n; // Start point x
	int m; // Start point y
	int x; // End point x
	int y; // End point y
	double Pr; // Probability
	SEGMENT(void);
	SEGMENT(int stx, int sty, int endx, int endy, double prob);
};

struct LINEPOLE
{
	double r;
	double theta;
	double cos;
	double sin;
	LINEPOLE(void);
	LINEPOLE(double rad, double th, double icos, double isin);
};

#endif

