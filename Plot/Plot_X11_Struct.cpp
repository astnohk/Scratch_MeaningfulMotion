#include "Plot_X11_Struct.h"



SEGMENT_X11::SEGMENT_X11(void)
{
	start.x = 0;
	start.y = 0;
	end.x = 0;
	end.y = 0;
}


XPLOT::XPLOT(void)
{
	point.x = 0;
	point.y = 0;
	z = .0;
}

void
XPLOT::set(XPoint _point, double _z)
{
	point = _point;
	z = _z;
}

void
XPLOT::set(short int point_x, short int point_y, double _z)
{
	point.x = point_x;
	point.x = point_y;
	z = _z;
}

