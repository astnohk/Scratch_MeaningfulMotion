#include "Vector.h"




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




VECTOR_2D_W_SCORE::VECTOR_2D_W_SCORE(void)
{
	x = .0;
	y = .0;
	score = .0;
}

VECTOR_2D_W_SCORE::VECTOR_2D_W_SCORE(double ix, double iy, double iscore)
{
	x = ix;
	y = iy;
	score = iscore;
}

void
VECTOR_2D_W_SCORE::reset(void)
{
	x = .0;
	y = .0;
	score = .0;
}

