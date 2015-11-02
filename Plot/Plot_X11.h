#include <X11/Xlib.h>
#include "../lib/ImgStruct.h"



#ifndef LIB_PLOT_X11
#define LIB_PLOT_X11

#define WAIT_TIME (1E6 / 20.0)
#define WINDOW_X_DEFAULT 800
#define WINDOW_Y_DEFAULT 800

/* Constants for Events */
#define X_NO_EVENT 0
#define X_ESCAPE 1
/* /Constants for Events */


/* Xlib variables */
extern const char *ProgramName;
extern SIZE Window_size;
extern Display *disp;
extern Window win;
extern Pixmap pix;
extern GC GCmono;
#define RGB_COLOR 3
extern GC GCcol[RGB_COLOR];
extern GC GCcol_dark[RGB_COLOR];
extern Colormap cmap;
/* /Xlib variables */


#define ROTATE_ANGLE_MAX 3600
extern double cos_a[ROTATE_ANGLE_MAX];
extern double sin_a[ROTATE_ANGLE_MAX];

#endif

