#include "Plot_X11_Struct.h"



SEGMENT_X11::SEGMENT_X11(void)
{
	start = (XPoint){0, 0};
	end = (XPoint){0, 0};
}


XPLOT::XPLOT(void)
{
	point = (XPoint){0, 0};
	z = .0;
}

void
XPLOT::set(XPoint _point, double _z)
{
	point = _point;
	z = _z;
}

void
XPLOT::set(int point_x, int point_y, double _z)
{
	point.x = point_x;
	point.x = point_y;
	z = _z;
}

