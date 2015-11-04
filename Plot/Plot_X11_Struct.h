#include <X11/Xlib.h>



#ifndef LIB_PLOT_X11_STRUCT
#define LIB_PLOT_X11_STRUCT

#define DEFAULT_INTENSITY_INTERVAL 8
#define DEFAULT_PLOT_Z_SCALE 0.1

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
	void set(int point_x, int point_y, double _z);
};

#endif
