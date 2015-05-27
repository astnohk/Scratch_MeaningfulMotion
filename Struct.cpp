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
	next = NULL;
}

FRAGMENT::FRAGMENT(int s, int e, double prob, FRAGMENT *p)
{
	start = s;
	end = e;
	Pr = prob;
	next = p;
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
FILTER_PARAM::ChangeFilter(char type)
{
	switch (type) {
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

TUPLE_VEC_SCALAR::TUPLE_VEC_SCALAR(void)
{
	for (int i = 0; i < TUPLE_VECTOR_SIZE; i++) {
		vector[i] = .0;
	}
	scalar = .0;
}

OPTICALFLOW_PARAM::OPTICALFLOW_PARAM(void) : WindowSize(25, 25)
{
	Level = 4;
	IRLS_Iter_Max = 16;
	IRLS_Convergence_Threshold = 1.0E-4;
	IRLS_Min_C = 1.0E-20;
}

OPTIONS::OPTIONS(void) : ResampleMethod("z-hold")
{
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
	x = .0f;
	y = .0f;
	z = .0f;
}

COORDINATE_3D::COORDINATE_3D(float ix, float iy, float iz)
{
	x = ix;
	y = iy;
	z = iz;
}

XPLOT::XPLOT(void)
{
	point = (XPoint){0, 0};
	z = .0;
}

