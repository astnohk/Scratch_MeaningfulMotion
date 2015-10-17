#include "../Scratch_MeaningfulMotion.h"
#include "Plot_X11.h"



X11_PARAM::X11_PARAM(void)
{
	Int_interval = 1;
	Latitude = 0;
	Longitude = 0;
	Center_x = .0;
	Center_y = .0;
	Center_z = .0;
	Scale = 1.0;
	Plot_Z_Scale = DEFAULT_PLOT_Z_SCALE;
	RotateSwitch = 0;
	ModeSwitch = 0;
	FillSwitch = 0;
}

SEGMENT_X11::SEGMENT_X11(void)
{
	start = (XPoint){0, 0};
	end = (XPoint){0, 0};
}


COORDINATE_3D::COORDINATE_3D(void)
{
	x = .0;
	y = .0;
	z = .0;
}

COORDINATE_3D::COORDINATE_3D(double ix, double iy, double iz)
{
	x = ix;
	y = iy;
	z = iz;
}

void
COORDINATE_3D::set(double sx, double sy, double sz)
{
	x = sx;
	y = sy;
	z = sz;
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

