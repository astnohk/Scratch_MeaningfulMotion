#ifndef LIB_Vector
#define LIB_Vector
struct VECTOR_2D
{
	double x;
	double y;
	VECTOR_2D(void);
	VECTOR_2D(double ix, double iy);
	void reset(void);
};

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
#endif

