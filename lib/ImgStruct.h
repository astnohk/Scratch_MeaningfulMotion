#include <new>

/*
#ifndef nullptr
#define nullptr 0
#endif
*/


#ifndef LIB_ImgStruct
#define LIB_ImgStruct

struct SIZE
{
	int width;
	int height;
	SIZE(void);
	SIZE(int w, int h);
	void reset(void);
	void set_size(const SIZE &size);
	void set_size(const SIZE *size);
};

struct COORDINATE
{
	int x;
	int y;
	COORDINATE(void);
	COORDINATE(int ix, int iy);
};

struct COORDINATE_3D
{
	double x;
	double y;
	double z;
	COORDINATE_3D(void);
	COORDINATE_3D(double ix, double iy, double iz);
	void set(double sx, double sy, double sz);
};

#endif

