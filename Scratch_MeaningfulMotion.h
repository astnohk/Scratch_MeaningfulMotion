/*
 * In this program. the coordinate system defined as below.
 * --->x     --->n
 * |         |
 * V         V
 * y         m
*/


/* DEBUG Options */
//#define SHOW_GAUSSIAN_FILTER
/* /DEBUG Options */


#define POW2(x) ((x) * (x))
#define POW3(x) ((x) * (x) * (x))
#define SIGN(x) ((x) >= 0 ? (x) > 0 ? 1 : 0 : -1)
#define SIGN_NOZERO(x) ((x) >= 0 ? 1 : -1)
#define SATURATE(x, min, max) (min <= (x) ? (x) <= max ? (x) : max : min)
#define MIN(x, y) ((x) <= (y) ? (x) : (y))
#define MAX(x, y) ((x) >= (y) ? (x) : (y))



#define _USE_MATH_DEFINES

#define MEANINGFUL_TRUE 1
#define MEANINGFUL_FALSE 0
#define MEANINGFUL_SUCCESS 1
#define MEANINGFUL_FAILURE 0


// C++
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <list>
#include <string>
#include "pnm.h"
#if defined(_OPENMP)
#include <omp.h>
#endif
#include <X11/Xlib.h>
#include <X11/keysym.h>


#define REGEXP_MAX_DIGITS 10
#define NUM_PROGRESS 64

#define PPM_EXTENSION_LENGTH 3
#define OPTION_INSUFFICIENT 1
#define OPTION_INCORRECT 2
#define OPTION_UNKNOWN 4

/* Superimpose Color constants */
#define OVERLAY_COLOR_PATTERNS 3
#define NOT_SUPERIMPOSE 0
#define RED 1
#define GREEN 2
#define BLUE 3

/* ANGLE_MAX is used as ANGLE_MAX * pi */
#define ANGLE_MAX 2.0

/* Filter Default Values */
#define EPSILONFILTER_SIZE (SIZE){21, 21}
#define EPSILONFILTER_EPSILON 20
#define GAUSSIAN_SIZE (SIZE){21, 21}
#define GAUSSIAN_STD_DEVIATION 5.0


/*  Scratch Detection Parameters
 *
 *                       A
 *                      / \
 *                     /   \
 *                    /     \
 *                   /       \
 *                  /         \
 *                 /           \
 *       ----------             ----------
 *       <--------><-----------><-------->
 *         Left       Scratch      Right
 *         Ave.       Width        Ave.
 *       <---------.           .--------->
 *  AVE_MAX_FAR                     AVE_MAX_FAR
*/
//#define SCRATCH_WIDTH 24
//#define AVE_MAX_FAR 40
#define SCRATCH_WIDTH 3
#define AVE_MAX_FAR 5
#define AVE_FAR (SCRATCH_WIDTH / 2 + AVE_MAX_FAR)
#define MEAN_WIDTH SCRATCH_WIDTH
#define SCRATCH_MED_THRESHOLD 3
#define SCRATCH_AVG_THRESHOLD 20


#define DERIVATIVE_MINIMUM 0
#define EPSILON 1
/* DIR_PROBABILITY independent of DIV_ANGLE */
#define DIR_PROBABILITY (1.0 / 16.0)
/* DIV_ANGLE affects on segments searching but not on precision */
#define DIV_ANGLE 40
/* DIV_ANGLE_VERTICAL limits the segments angle between Vertical line and Segments less than or equal to ((pi / 2) / DIV_ANGLE_VERTICAL) */
#define DIV_ANGLE_VERTICAL 9.0

#define EXCLUSIVE_PRINCIPLE_MAX_RADIUS 1.5

#define PLOT_INTENSITY_MAX 255



class ERROR
{
	private:
		std::string FunctionName;
		std::string ErrorFunctionName;
		std::string ValueName;
		std::string FileName;
	public:
		ERROR();
		// Error data set
		void Function(const char *name);
		void ErrorFunction(const char *name);
		void Value(const char *name);
		void File(const char *name);
		// Error output
		void Others(const char *error);
		void Malloc(void);
		void FunctionFail(void);
		void PointerNull(void);
		void ImageSize(void);
		void FileRead(void);
		void FileWrite(void);
};



struct SIZE
{
	int width;
	int height;
	SIZE(void);
	SIZE(int w, int h);
	void reset(void);
};

struct COORDINATE
{
	int x;
	int y;
	COORDINATE(void);
	COORDINATE(int ix, int iy);
};

struct FRAGMENT
{
	int start;
	int end;
	double Pr;
	FRAGMENT *next;
	FRAGMENT(void);
	FRAGMENT(int s, int e, double prob, FRAGMENT *p);
};

struct SEGMENT
{
	int n; /* Start point x */
	int m; /* Start point y */
	int x; /* End point x */
	int y; /* End point y */
	double Pr; /* Probability */
	SEGMENT(void);
	SEGMENT(int stx, int sty, int endx, int endy, double prob);
};

struct LINEPOLE
{
	double r;
	double theta;
	double cos;
	double sin;
	LINEPOLE(void);
	LINEPOLE(double rad, double th, double icos, double isin);
};

struct FILTER_PARAM
{
	int type;
	SIZE size;
	double std_deviation;
	double epsilon;
	FILTER_PARAM(void);
	void ChangeFilter(char type);
};
#define NUM_FILTER_TYPE 3
#define FILTER_ID_UNDEFINED 0
#define FILTER_ID_EPSILON 1
#define FILTER_ID_GAUSSIAN 2

struct VECTOR_2D
{
	double x;
	double y;
	VECTOR_2D(void);
	VECTOR_2D(double ix, double iy);
};

#define TUPLE_VECTOR_SIZE 6
struct TUPLE_VEC_SCALAR{
	double vector[TUPLE_VECTOR_SIZE];
	double scalar;
	TUPLE_VEC_SCALAR(void);
	TUPLE_VEC_SCALAR operator+(TUPLE_VEC_SCALAR tup);
	TUPLE_VEC_SCALAR operator*(double C);
	double operator*(TUPLE_VEC_SCALAR tup);
	double norm(SIZE W_l);
};

struct OPTICALFLOW_PARAM
{
	int Level;
	SIZE WindowSize;
	int IRLS_Iter_Max;
	double IRLS_Convergence_Threshold;
	double IRLS_Min_C;
	OPTICALFLOW_PARAM(void);
};
#define OPTICALFLOW_PARAM_DEFAULT (OPTICALFLOW_PARAM){1, (SIZE){25, 25}, 1600, 1.0E-20, 8}


#define NUM_AFFINE_PARAMETER 6
struct VECTOR_AFFINE
{
	double a[NUM_AFFINE_PARAMETER];
	VECTOR_AFFINE(void);
};

struct MULTIPLE_MOTION_PARAM
{
	int Level;
	SIZE WindowSize;
	int IRLS_Iter_Max;
	double IRLS_Min_Threshold;
	double lambdaD;
	double lambdaS;
	double sigmaD;
	double sigmaS;
	MULTIPLE_MOTION_PARAM(void);
};
#define MULTIPLE_MOTION_PARAM_DEFAULT (MULTIPLE_MOTION_PARAM){4, (SIZE){25, 25}, 16, 0.0001, 5, 1, 12.72, 2.121}
/* /MultipleMotion.c */


struct OPTIONS
{
	SIZE ResampleSize;
	std::string ResampleMethod;
	int mode;
	int Max_Length;
	int Max_Output_Length;
	int ExclusivePrinciple;
	int Superimpose;
	int PlotOptions;
	int s_med;
	int s_avg;
	double p;
	double ep;
	double Exclusive_Max_Radius;
	OPTICALFLOW_PARAM OpticalFlow_Param;
	MULTIPLE_MOTION_PARAM MultipleMotion_Param;
	OPTIONS(void);
};
/* Mode Options */
#define MODE_OUTPUT_FILTERED_IMAGE 0x0010
#define MODE_OUTPUT_BINARY_IMAGE 0x0020
#define MODE_OUTPUT_MULTIPLE_MOTIONS_AFFINE 0x0040
#define MODE_OUTPUT_MULTIPLE_MOTIONS_OPTICALFLOW 0x0080
#define MODE_OUTPUT_OPTICALFLOW_AFFINE_PARAMETER 0x0100
/* PlotOptions */
#define PLOT_NEGATE 0x01
#define PLOT_AS_RESAMPLE 0x02
#define PLOT_RESAMPLED_IMG_ONLY 0x04

#define OPTIONS_DEFAULT (OPTIONS){SIZE_ZERO, "z-hold", 0, 0, 0, MEANINGFUL_FALSE, 0, 0, SCRATCH_MED_THRESHOLD, SCRATCH_AVG_THRESHOLD, DIR_PROBABILITY, EPSILON, EXCLUSIVE_PRINCIPLE_MAX_RADIUS, OPTICALFLOW_PARAM_DEFAULT, MULTIPLE_MOTION_PARAM_DEFAULT}


/* X11 structures */
/* * X11 Plotting Parameters */
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

struct SEGMENT_X11
{
	XPoint start;
	XPoint end;
	SEGMENT_X11(void);
};

struct COORDINATE_3D
{
	float x;
	float y;
	float z;
	COORDINATE_3D(void);
	COORDINATE_3D(float ix, float iy, float iz);
};

struct XPLOT
{
	XPoint point;
	double z;
	XPLOT(void);
};
/* /X11 structures */


extern const std::string Progress[NUM_PROGRESS];
extern const std::string Progress_End;





/* Prototype of functions */
int Scratch_MeaningfulMotion(char *OutputName, char *InputName, unsigned int OutputNameLength, unsigned int InputNameLength, int Start, int End, OPTIONS Options, FILTER_PARAM FilterParam);

/* Memory Libraries */
void list_free(std::list<SEGMENT>* plist);
void segments_free(FRAGMENT *segments);

/* Mathematical Libraries */
double pow_int(double x, int a);
double atan2_div_pi_table(int y, int x, SIZE *size);
int* Calc_k_l(SIZE size, double p, double ep);
double Pr(int k, int l, double p);

/* Image Libraries */
double* Gaussian(double *img, SIZE size, FILTER_PARAM Param);
double* EpsilonFilter(double *img, SIZE size, FILTER_PARAM Param);
double HorizontalMedian(double *img, int N, int x, int y, int width);
double* DerivativeAngler(double *img, SIZE size);
VECTOR_2D* Derivator(double *Image, SIZE size, char *Type);
double* Derivation_abs(VECTOR_2D *Derivative_2D, SIZE size);
double* Filterer(double *Image, SIZE size, double *Filter, SIZE size_f, int Mirroring);
int IndexOfMirroring(int x, int size);

/* Other Libraries */
char* regexp(char *s);


/* Scratch Detection */
double* DetectScratch(PNM *pnm, double s_med, double s_avg, FILTER_PARAM FilterParam, int Do_Detection);
#define DO_DETECTION 1
#define DO_NOT_DETECTION 0

/* Meaningful Alignments */
SEGMENT* AlignedSegment_vertical(double *angles, SIZE size, int *k_list, int l_min, double *Pr_table, int *Num_segments, int Max_Length, int Max_Output_Length);
int AlignedCheck(double *angles, SIZE size, int *k_list, double *Pr_table, std::list<SEGMENT>* list_start, int l_min, int m, int n, int x, int y, int Max_Length, int Max_Output_Length);
SEGMENT* ExclusivePrinciple(double *angles, SIZE size, int *k_list, double *Pr_table, SEGMENT *MaximalSegments, int *Num_Segments, double Exclusive_max_radius);

/* Plotting */
int* PlotSegment(SEGMENT *coord_array, int Num_Segments, SIZE size, SIZE size_out, int Negate);
int Superimposer(PNM *pnm_out, PNM *pnm_in, int *Plot, SIZE size, int Color, int Negate);


/* X11 Plotting */
int ShowSegments_X11(int *Img, SIZE Img_size, SIZE Img_size_resample, int MaxInt, SEGMENT *Segment_Array, unsigned int Num_Segments);
int Init_X11(X11_PARAM *X11_Param, SIZE Img_size);
int XEventor(X11_PARAM *X11_Param, SIZE Img_size);
void SwitchEventer(X11_PARAM *X11_Param);
int TransRotate_3DSegment(X11_PARAM X11_Param, SEGMENT *segments, SEGMENT_X11 *segments_plot, unsigned int Num_Segments, SIZE Img_size, SIZE Img_size_resample);
int TransRotate_3DPoint(X11_PARAM X11_Param, int *Img, SIZE size, int MaxInt, XPLOT *Img_plot);
int TransGaraxy_3DPoint(X11_PARAM X11_Param, int *Img, SIZE Img_size, COORDINATE_3D *Img_coord, COORDINATE_3D *Img_vel, COORDINATE_3D GaraxyCenter, XPLOT *Img_plot);
int TransGravity_3DPoint(X11_PARAM X11_Param, int *Img, SIZE size, COORDINATE_3D *Img_coord, COORDINATE_3D *Img_vel, XPLOT *Img_plot);
int Plot_3DPoint(X11_PARAM X11_Param, int *Img, XPLOT *Img_plot, int *Img_index, SIZE size);
int Plot_3DPointANDSegment(X11_PARAM X11_Param, int *Img, XPLOT *Img_plot, int *Img_index, SIZE size, SEGMENT_X11 *segments_plot, unsigned int Num_Segments);
int Plot_3DGridANDSegment(X11_PARAM X11_Param, int *Img, XPLOT *Img_plot, int *Img_index, SIZE size, SEGMENT_X11 *segments_plot, unsigned int Num_Segments);
void PlotParameters(X11_PARAM X11_Param);
int reset_index(int *Img_index, int N);
int sort_index(XPLOT *Img_plot, int *Index, int *Index_tmp, int N);
/* /X11 Plotting */

