#include "../ImgClass/Vector.h"


#ifndef LIB_Ext_Vector
#define LIB_Ext_Vector

struct VECTOR_2D_W_SCORE
{
	double x;
	double y;
	double score;
	// member function
	VECTOR_2D_W_SCORE(void);
	VECTOR_2D_W_SCORE(double ix, double iy, double iscore);
	void reset(void);
};

#define NUM_AFFINE_PARAMETER 6
struct VECTOR_AFFINE
{
	double a[NUM_AFFINE_PARAMETER];
	VECTOR_AFFINE(void);
	void reset(void);
};
#endif

