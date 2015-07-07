/* Definition of The Constructors and The Member functions of struct */
#include "Scratch_MeaningfulMotion.h"
#include "Plot_X11.h"



SIZE::SIZE(void)
{
	width = 0;
	height = 0;
}

SIZE::SIZE(int w, int h)
{
	width = w;
	height = h;
}

void
SIZE::reset(void)
{
	width = 0;
	height = 0;
}

COORDINATE::COORDINATE(void)
{
	x = 0;
	y = 0;
}

COORDINATE::COORDINATE(int ix, int iy)
{
	x = ix;
	y = iy;
}

FRAGMENT::FRAGMENT(void)
{
	start = 0;
	end = 0;
	Pr = .0;
}

FRAGMENT::FRAGMENT(int s, int e, double prob)
{
	start = s;
	end = e;
	Pr = prob;
}

SEGMENT::SEGMENT(void)
{
	n = 0;
	m = 0;
	x = 0;
	y = 0;
	Pr = .0;
}

SEGMENT::SEGMENT(int stx, int sty, int endx, int endy, double prob)
{
	n = stx;
	m = sty;
	x = endx;
	y = endy;
	Pr = prob;
}

LINEPOLE::LINEPOLE(void)
{
	r = .0;
	theta = .0;
	cos = .0;
	sin = .0;
}

LINEPOLE::LINEPOLE(double rad, double th, double icos, double isin)
{
	r = rad;
	theta = th;
	cos = icos;
	sin = isin;
}

FILTER_PARAM::FILTER_PARAM(void)
{
	type = FILTER_ID_UNDEFINED;
	std_deviation = .0;
	epsilon = .0;
}

void
FILTER_PARAM::ChangeFilter(char newtype)
{
	switch (newtype) {
		case 'e': // Epsilon filter
			type = FILTER_ID_EPSILON;
			size = EPSILONFILTER_SIZE;
			std_deviation = .0;
			epsilon = EPSILONFILTER_EPSILON;
			break;
		case 'g': // Gaussian lowpass filter
			type = FILTER_ID_GAUSSIAN;
			size = GAUSSIAN_SIZE;
			std_deviation = GAUSSIAN_STD_DEVIATION;
			epsilon = 0;
	}
}

VECTOR_2D::VECTOR_2D(void)
{
	x = .0;
	y = .0;
}

VECTOR_2D::VECTOR_2D(double ix, double iy)
{
	x = ix;
	y = iy;
}

void
VECTOR_2D::reset(void)
{
	x = .0;
	y = .0;
}

VECTOR_2D_W_SCORE::VECTOR_2D_W_SCORE(void)
{
	x = .0;
	y = .0;
	score = .0;
}

VECTOR_2D_W_SCORE::VECTOR_2D_W_SCORE(double ix, double iy, double iscore)
{
	x = ix;
	y = iy;
	score = iscore;
}

void
VECTOR_2D_W_SCORE::reset(void)
{
	x = .0;
	y = .0;
	score = .0;
}


VECTOR_AFFINE::VECTOR_AFFINE(void)
{
	for (int i = 0; i < NUM_AFFINE_PARAMETER; i++) {
		a[i] = .0;
	}
}

void
VECTOR_AFFINE::reset(void)
{
	for (int i = 0; i < NUM_AFFINE_PARAMETER; i++) {
		a[i] = .0;
	}
}


// Multiple Motion Parameters
#define MULTIPLE_MOTION_PARAM_Level 5
#define MULTIPLE_MOTION_PARAM_IRLS_Iter_Max 300
#define MULTIPLE_MOTION_PARAM_Error_Min_Threshold (1.0E-6)
#define MULTIPLE_MOTION_PARAM_lambdaD 5
#define MULTIPLE_MOTION_PARAM_lambdaS 1
#define MULTIPLE_MOTION_PARAM_sigmaD 12.72
#define MULTIPLE_MOTION_PARAM_sigmaS 2.121
MULTIPLE_MOTION_PARAM::MULTIPLE_MOTION_PARAM(void)
{
	Level = MULTIPLE_MOTION_PARAM_Level;
	IRLS_Iter_Max = MULTIPLE_MOTION_PARAM_IRLS_Iter_Max;
	Error_Min_Threshold = MULTIPLE_MOTION_PARAM_Error_Min_Threshold;
	lambdaD = MULTIPLE_MOTION_PARAM_lambdaD;
	lambdaS = MULTIPLE_MOTION_PARAM_lambdaS;
	sigmaD = MULTIPLE_MOTION_PARAM_sigmaD;
	sigmaS = MULTIPLE_MOTION_PARAM_sigmaS;
}

void
MULTIPLE_MOTION_PARAM::set_default(const char *val)
{
	if (strcmp(val, "Level") == 0) {
		Level = MULTIPLE_MOTION_PARAM_Level;
	} else if (strcmp(val, "IRLS_Iter_Max") == 0) {
		IRLS_Iter_Max = MULTIPLE_MOTION_PARAM_IRLS_Iter_Max;
	} else if (strcmp(val, "Error_Min_Threshold") == 0) {
		Error_Min_Threshold = MULTIPLE_MOTION_PARAM_Error_Min_Threshold;
	} else if (strcmp(val, "lambdaD") == 0) {
		lambdaD = MULTIPLE_MOTION_PARAM_lambdaD;
	} else if (strcmp(val, "lambdaS") == 0) {
		lambdaS = MULTIPLE_MOTION_PARAM_lambdaS;
	} else if (strcmp(val, "sigmaD") == 0) {
		sigmaD = MULTIPLE_MOTION_PARAM_sigmaD;
	} else if (strcmp(val, "sigmaS") == 0) {
		sigmaS = MULTIPLE_MOTION_PARAM_sigmaS;
	} else {
		fprintf(stderr, "*** MULTIPLE_MOTION_PARAM::set_default() error - There are NOT such a parameter '%s' ***\n", val);
	}
}


OPTIONS::OPTIONS(void)
{
	ResampleMethod = 0;
	mode = 0;
	Max_Length = 0;
	Max_Output_Length = 0;
	ExclusivePrinciple = MEANINGFUL_FALSE;
	Superimpose = 0;
	PlotOptions = 0;
	s_med = SCRATCH_MED_THRESHOLD;
	s_avg = SCRATCH_AVG_THRESHOLD;
	p = DIR_PROBABILITY;
	ep = EPSILON;
	Exclusive_Max_Radius = EXCLUSIVE_PRINCIPLE_MAX_RADIUS;
}


X11_PARAM::X11_PARAM(void)
{
	Int_interval = 1;
	Latitude = 0;
	Longitude = 0;
	Center_x = .0;
	Center_y = .0;
	Center_z = .0;
	Scale = 1.0;
	Plot_Z_Scale = DEFAULT_PLOT_Z_SCALE;
	RotateSwitch = 0;
	ModeSwitch = 0;
	FillSwitch = 0;
}

SEGMENT_X11::SEGMENT_X11(void)
{
	start = (XPoint){0, 0};
	end = (XPoint){0, 0};
}

COORDINATE_3D::COORDINATE_3D(void)
{
	x = .0;
	y = .0;
	z = .0;
}

COORDINATE_3D::COORDINATE_3D(double ix, double iy, double iz)
{
	x = ix;
	y = iy;
	z = iz;
}

void
COORDINATE_3D::set(double sx, double sy, double sz)
{
	x = sx;
	y = sy;
	z = sz;
}

XPLOT::XPLOT(void)
{
	point = (XPoint){0, 0};
	z = .0;
}

