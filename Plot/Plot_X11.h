#include <X11/Xlib.h>
#include "../lib/ImgStruct.h"



#ifndef LIB_PLOT_X11
#define LIB_PLOT_X11

#define WAIT_TIME (1000000 / 20)
#define WINDOW_X_DEFAULT 800
#define WINDOW_Y_DEFAULT 800


// Plotting parameters
#define DEFAULT_INTENSITY_INTERVAL 8
#define DEFAULT_PLOT_Z_SCALE 0.1


// Constants for Events
#define X_NO_EVENT 0
#define X_ESCAPE 1


// Xlib variables
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


#define ROTATE_ANGLE_MAX 3600
extern double cos_a[ROTATE_ANGLE_MAX];
extern double sin_a[ROTATE_ANGLE_MAX];




void freeXWindow(void);
void Init_X11(X11_PARAM* X11_Param, SIZE Img_size);
int XEventor(X11_PARAM* X11_Param, SIZE Img_size);
void SwitchEventer(X11_PARAM *X11_Param);

void TransRotate_3DSegment(const X11_PARAM& X11_Param, const SEGMENT* segments, SEGMENT_X11* segments_plot, const unsigned int Num_Segments, const SIZE& Img_size, const SIZE& Img_size_resample);
void TransRotate_3DPoint(const X11_PARAM& X11_Param, const ImgVector<int>* Img, const int MaxInt, ImgVector<XPLOT>* Img_plot);
void TransGaraxy_3DPoint(const X11_PARAM& X11_Param, const ImgVector<int>* Img, ImgVector<COORDINATE_3D>* Img_coord, ImgVector<COORDINATE_3D>* Img_vel, const COORDINATE_3D& GaraxyCenter, ImgVector<XPLOT>* Img_plot);
void TransGravity_3DPoint(const X11_PARAM& X11_Param, const ImgVector<int>* Img, ImgVector<COORDINATE_3D>* Img_coord, ImgVector<COORDINATE_3D>* Img_vel, ImgVector<XPLOT>* Img_plot);

void Plot_3DPoints(const X11_PARAM& X11_Param, const ImgVector<int>* Img, ImgVector<XPLOT>* Img_plot, size_t* Img_index);
void Plot_3DGrid(const X11_PARAM& X11_Param, const ImgVector<int>* Img, ImgVector<XPLOT>* Img_plot, size_t* Img_index);
void Plot_3DSegment(const X11_PARAM& X11_Param, const SEGMENT_X11* segments_plot, const unsigned int Num_Segments);
void PlotParameters(const X11_PARAM& X11_Param);

void Set_Pixmap2Window(void);
void reset_index(size_t* Img_index, const size_t N);
void sort_index(const ImgVector<XPLOT>* Img_plot, size_t* Index, size_t* Index_tmp, const size_t N);

#endif

