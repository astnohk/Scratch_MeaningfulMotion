#include <X11/Xlib.h>



#ifndef LIB_PLOT_X11_STRUCT
#define LIB_PLOT_X11_STRUCT

struct SEGMENT_X11
{
	XPoint start;
	XPoint end;
	SEGMENT_X11(void);
};

struct XPLOT
{
	XPoint point;
	double z;
	XPLOT(void);
	void set(XPoint _point, double _z);
	void set(short int point_x, short int point_y, double _z);
};

#endif

