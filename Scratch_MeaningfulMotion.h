/*
 * In this program. the coordinate system defined as below.
 * --->x     --->n
 * |         |
 * V         V
 * y         m
*/

// Options for Old Compiler which is not available to use -std=c++11
#ifndef nullptr
#define nullptr NULL
#endif


// DEBUG Options
//#define SHOW_GAUSSIAN_FILTER
// /DEBUG Options


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
#define MEANINGFUL_SUCCESS true
#define MEANINGFUL_FAILURE false


// C++
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <list>
#include <string>
#if defined(_OPENMP)
#include <omp.h>
#endif
#include <X11/Xlib.h>
#include <X11/keysym.h>
#include <unistd.h>

#include "lib/Class.h"
#include "lib/ImgStruct.h"
#include "lib/Struct.h"
#include "lib/ExtVector.h"
#include "ImgClass/ImgClass.h"
#include "ImgClass/ImgStatistics.h"
#include "Plot/Plot_X11_Struct.h"
#include "HOG/HOG_struct.h"
#include "Scratch_Struct.h"
#include "PNM/pnm.h"


#define REGEXP_MAX_DIGITS 10
#define NUM_PROGRESS 64

#define PPM_EXTENSION_LENGTH 3
#define OPTION_INSUFFICIENT 1
#define OPTION_INCORRECT 2
#define OPTION_UNKNOWN 4

// Superimpose Color constants
#define OVERLAY_COLOR_PATTERNS 3
#define NOT_SUPERIMPOSE 0
#define RED 1
#define GREEN 2
#define BLUE 3

// ANGLE_MAX is used as ANGLE_MAX * pi
#define ANGLE_MAX 2.0

// Filter Default Values
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
// DIR_PROBABILITY independent of DIV_ANGLE
#define DIR_PROBABILITY (1.0 / 16.0)
// DIV_ANGLE affects on segments searching but not on precision
#define DIV_ANGLE 40
// DIV_ANGLE_VERTICAL limits the segments angle between Vertical line and Segments less than or equal to ((pi / 2) / DIV_ANGLE_VERTICAL)
#define DIV_ANGLE_VERTICAL 9.0

#define EXCLUSIVE_PRINCIPLE_MAX_RADIUS 1.5

#define PLOT_INTENSITY_MAX 255




// Multiple Motion Parameters
#define MULTIPLE_MOTION_PARAM_Level 5
#define MULTIPLE_MOTION_PARAM_IRLS_Iter_Max 300
#define MULTIPLE_MOTION_PARAM_Error_Min_Threshold (1.0E-6)
#define MULTIPLE_MOTION_PARAM_lambdaD 5
#define MULTIPLE_MOTION_PARAM_lambdaS 1
#define MULTIPLE_MOTION_PARAM_sigmaD 12.72
#define MULTIPLE_MOTION_PARAM_sigmaS 2.121
#define MULTIPLE_MOTION_PARAM_BlockMatching_BlockSize 17;




extern const std::string Progress[NUM_PROGRESS];
extern const std::string Progress_End;




// Prototype of functions
int Scratch_MeaningfulMotion(char *OutputName, char *InputName, unsigned int OutputNameLength, unsigned int InputNameLength, int Start, int End, OPTIONS Options, FILTER_PARAM FilterParam);

// Mathematical Libraries
double pow_int(double x, int a);
int* Calc_k_l(SIZE &size, double p, double ep);
double Pr(int k, int l, double p);

// Image Libraries
ImgVector<double>* Gaussian(ImgVector<double> *img, FILTER_PARAM Param);
ImgVector<double>* EpsilonFilter(ImgVector<double> *img, FILTER_PARAM Param);
double HorizontalMedian(ImgVector<double> *img, int x, int y, int width);
ImgVector<double>* DerivativeAngler(ImgVector<double> *img);
ImgVector<VECTOR_2D<double> >* Derivator(ImgVector<double> *Image, const char *Type);
ImgVector<double>* Derivation_abs(ImgVector<VECTOR_2D<double> > *Derivative_2D);
ImgVector<double>* Filterer(ImgVector<double> *Image, ImgVector<double> *Filter, bool Mirroring);
int IndexOfMirroring(int x, int size);

// Other Libraries
char* regexp(char *s);


// Scratch Detection
ImgVector<double>* DetectScratch(const PNM &pnm, double s_med, double s_avg, FILTER_PARAM FilterParam, bool Do_Detection);
#define DO_DETECTION true
#define DO_NOT_DETECTION false

// Meaningful Alignments
SEGMENT* AlignedSegment_vertical(ImgVector<double> *angles, int *k_list, int l_min, ImgVector<double> *Pr_table, int *Num_segments, int Max_Length, int Max_Output_Length);
std::list<FRAGMENT>* AlignedCheck(ImgVector<double> *angles, int *k_list, ImgVector<double> *Pr_table, int l_min, int m, int n, int x, int y, int Max_Length);
bool MaximalMeaningfulness(std::list<SEGMENT>* list_segment, std::list<FRAGMENT>* list_fragment, int m, int n, int x, int y, int Max_Output_Length);
SEGMENT* ExclusivePrinciple(ImgVector<double> *angles, int *k_list, ImgVector<double> *Pr_table, SEGMENT *MaximalSegments, int *Num_Segments, double Exclusive_max_radius);
ImgVector<int>* ExclusiveIndexMap(SIZE size, SEGMENT *MaximalSegments, int *Num_Segments, double Exclusive_max_radius);
SEGMENT* ExclusiveSegments(ImgVector<int> *IndexMap, ImgVector<double> *angles, SEGMENT *MaximalSegments, int *Num_Segments, int *k_list, ImgVector<double> *Pr_table);

// Plotting
int* PlotSegment(SEGMENT *coord_array, int Num_Segments, SIZE size, SIZE size_out, int Negate);
bool Superimposer(PNM *pnm_out, const PNM &pnm_in, int *Plot, SIZE size, int Color, int Negate);


// X11 Plotting
bool ShowSegments_X11(ImgVector<int> *Img, SIZE Img_size_resample, int MaxInt, SEGMENT *Segment_Array, unsigned int Num_Segments);
bool Init_X11(X11_PARAM *X11_Param, SIZE Img_size);
int XEventor(X11_PARAM *X11_Param, SIZE Img_size);
void SwitchEventer(X11_PARAM *X11_Param);
bool TransRotate_3DSegment(X11_PARAM X11_Param, SEGMENT *segments, SEGMENT_X11 *segments_plot, unsigned int Num_Segments, SIZE Img_size, SIZE Img_size_resample);
bool TransRotate_3DPoint(X11_PARAM X11_Param, ImgVector<int> *Img, int MaxInt, ImgVector<XPLOT> *Img_plot);
bool TransGaraxy_3DPoint(X11_PARAM X11_Param, ImgVector<int> *Img, ImgVector<COORDINATE_3D> *Img_coord, ImgVector<COORDINATE_3D> *Img_vel, COORDINATE_3D GaraxyCenter, ImgVector<XPLOT> *Img_plot);
bool TransGravity_3DPoint(X11_PARAM X11_Param, ImgVector<int> *Img, ImgVector<COORDINATE_3D> *Img_coord, ImgVector<COORDINATE_3D> *Img_vel, ImgVector<XPLOT> *Img_plot);
bool Plot_3DPoints(X11_PARAM X11_Param, ImgVector<int> *Img, ImgVector<XPLOT> *Img_plot, int *Img_index);
bool Plot_3DGrid(X11_PARAM X11_Param, ImgVector<int> *Img, ImgVector<XPLOT> *Img_plot, int *Img_index);
bool Plot_3DSegment(X11_PARAM X11_Param, SEGMENT_X11 *segments_plot, unsigned int Num_Segments);
void PlotParameters(X11_PARAM X11_Param);
bool Set_Pixmap2Window(void);
bool reset_index(int *Img_index, int N);
bool sort_index(ImgVector<XPLOT> *Img_plot, int *Index, int *Index_tmp, int N);
// /X11 Plotting

