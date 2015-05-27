#define WAIT_TIME (1E6 / 20.0)
#define WINDOW_X_DEFAULT 800
#define WINDOW_Y_DEFAULT 800

/* Constants for Events */
#define X_NO_EVENT 0
#define X_ESCAPE 1
/* /Constants for Events */


/* Xlib variables */
extern char *ProgramName;
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

/* Other Parameters for Plot_X11 */
#define DEFAULT_INTENSITY_INTERVAL 8
#define DEFAULT_PLOT_Z_SCALE 0.1
#define NUMBER_OF_MODE 3
/* * Plot Mode */
#define X11_Plot_Point 0
#define X11_Plot_Point_A_Segment 1
#define X11_Plot_Grid_A_Segment 2
#define X11_Plot_Garaxy 3
#define X11_Plot_GravityCorrupt 4
/* * /Plot Mode */

