// Definition of The Constructors and The Member functions of struct
#include <cstdio>
#include <cstring>
#include "Class.h"
#include "Struct.h"




FRAGMENT::FRAGMENT(void)
{
	start = 0;
	end = 0;
	Pr = .0;
}

FRAGMENT::FRAGMENT(const int s, const int e, const double& prob)
{
	start = s;
	end = e;
	Pr = prob;
}

SEGMENT::SEGMENT(void)
{
	n = 0;
	m = 0;
	x = 0;
	y = 0;
	Pr = .0;
}

SEGMENT::SEGMENT(const int stx, const int sty, const int endx, const int endy, const double& prob)
{
	n = stx;
	m = sty;
	x = endx;
	y = endy;
	Pr = prob;
}

LINEPOLE::LINEPOLE(void)
{
	r = .0;
	theta = .0;
	cos = .0;
	sin = .0;
}

LINEPOLE::LINEPOLE(const double& rad, const double& th, const double& icos, const double& isin)
{
	r = rad;
	theta = th;
	cos = icos;
	sin = isin;
}

