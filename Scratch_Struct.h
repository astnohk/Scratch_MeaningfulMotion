#include "lib/ImgStruct.h"
#include "HOG/HOG_struct.h"



#ifndef LIB_SCRATCH_STRUCT
#define LIB_SCRATCH_STRUCT

#define NUM_FILTER_TYPE 3
#define FILTER_ID_UNDEFINED 0
#define FILTER_ID_EPSILON 1
#define FILTER_ID_GAUSSIAN 2
struct FILTER_PARAM
{
	int type;
	SIZE size;
	double std_deviation;
	double epsilon;
	FILTER_PARAM(void);
	bool ChangeFilter(const char *name);
	void set_default(const char *name);
	void set_value(const char *name, const void *value);
};

struct MULTIPLE_MOTION_PARAM
{
	int Level;
	int IRLS_Iter_Max;
	double Error_Min_Threshold;
	double lambdaD;
	double lambdaS;
	double sigmaD;
	double sigmaS;
	int BlockMatching_BlockSize;

	MULTIPLE_MOTION_PARAM(void);
	void set_default(const char *name);
	void set_value(const char *name, void *value);
};

// X11 structures
// * X11 Plotting Parameters
struct X11_PARAM
{
	int Int_interval;
	int Latitude;
	int Longitude;
	double Center_x;
	double Center_y;
	double Center_z;
	double Scale;
	double Plot_Z_Scale;
	short RotateSwitch;
	short ModeSwitch;
	short FillSwitch;
	X11_PARAM(void);
};

struct OPTIONS
{
	SIZE ResampleSize;
	int ResampleMethod;
	int mode;
	int Max_Length;
	int Max_Output_Length;
	bool ExclusivePrinciple;
	int Superimpose;
	int PlotOptions;
	int s_med;
	int s_avg;
	double p;
	double ep;
	double Exclusive_Max_Radius;
	MULTIPLE_MOTION_PARAM MultipleMotion_Param;
	HOG_PARAM HOG_Param;
	bool x11_plot;

	OPTIONS(void);
	bool ChangeResampleMethod(const char *name);
	void set_default(const char *name);
	void set_value(const char *name, const void *value);
};
// Mode Options
#define MODE_OUTPUT_FILTERED_IMAGE 0x0010
#define MODE_OUTPUT_BINARY_IMAGE 0x0020
#define MODE_OUTPUT_MULTIPLE_MOTIONS_AFFINE 0x0040
#define MODE_OUTPUT_OPTICALFLOW 0x0080
#define MODE_OUTPUT_AFFINE_BLOCKMATCHING 0x0100
#define MODE_OUTPUT_HISTOGRAMS_OF_ORIENTED_GRADIENTS_RAW_HOG 0x1000
#define MODE_OUTPUT_HISTOGRAMS_OF_ORIENTED_GRADIENTS 0x2000
#define MODE_OUTPUT_HISTOGRAMS_OF_ORIENTED_GRADIENTS_MATCHING_VECTOR 0x4000
// PlotOptions
#define PLOT_NEGATE 0x01
#define PLOT_AS_RESAMPLED 0x02
#define PLOT_RESAMPLED_IMG_ONLY 0x04

#endif

