#include "ExtVector.h"




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



VECTOR_AFFINE::VECTOR_AFFINE(void)
{
	for (int i = 0; i < NUM_AFFINE_PARAMETER; i++) {
		a[i] = .0;
	}
}

void
VECTOR_AFFINE::reset(void)
{
	for (int i = 0; i < NUM_AFFINE_PARAMETER; i++) {
		a[i] = .0;
	}
}

