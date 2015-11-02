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

