// Definition of The Constructors and The Member functions of struct
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

void
SIZE::set_size(const SIZE *size)
{
	if (size != nullptr) {
		width = size->width;
		height = size->height;
	}
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


// Filter Parameters
FILTER_PARAM::FILTER_PARAM(void)
{
	type = FILTER_ID_UNDEFINED;
	std_deviation = .0;
	epsilon = .0;
}

void
FILTER_PARAM::set_default(const char *name)
{
	if (strcmp(name, "type") == 0) {
		type = FILTER_ID_UNDEFINED;
	} else if (strcmp(name, "std_deviation") == 0) {
		std_deviation = .0;
	} else if (strcmp(name, "epsilon") == 0) {
		epsilon = .0;
	} else {
		fprintf(stderr, "*** FILTER_PARAM::set_default() error - There are NOT such a parameter '%s' ***\n", name);
	}
}

void
FILTER_PARAM::set_value(const char *name, const void *value)
{
	if (strcmp(name, "type") == 0) {
		type = *(int *)value;
		if (type < 0 || type >= NUM_FILTER_TYPE) {
			type = 0;
		}
	} else if (strcmp(name, "std_deviation") == 0) {
		std_deviation = *(double *)value;
		if (std_deviation < .0) {
			std_deviation = .0;
		}
	} else if (strcmp(name, "epsilon") == 0) {
		epsilon = *(double *)value;
		if (epsilon < .0) {
			epsilon = .0;
		}
	} else {
		fprintf(stderr, "*** FILTER_PARAM::set_default() error - There are NOT such a parameter '%s' ***\n", name);
	}
}

bool
FILTER_PARAM::ChangeFilter(const char *name)
{
	ERROR Error("FILTER_PARAM::ChangeFilter");
	char *lower = nullptr;

	if (strlen(name) > 1) {
		try {
			lower = new char[strlen(name) + 1];
		}
		catch (const std::bad_alloc &bad) {
			Error.Value("lower");
			Error.Malloc();
			return false;
		}
		for (size_t i = 0; i < strlen(name); i++) {
			lower[i] = tolower(name[i]);
		}
		lower[strlen(lower) + 1] = '\0';
		if (strcmp(lower, "epsilon") == 0) {
			type = FILTER_ID_EPSILON;
			size = EPSILONFILTER_SIZE;
			std_deviation = .0;
			epsilon = EPSILONFILTER_EPSILON;
		} else if (strcmp(lower, "gaussian") == 0) {
			type = FILTER_ID_GAUSSIAN;
			size = GAUSSIAN_SIZE;
			std_deviation = GAUSSIAN_STD_DEVIATION;
			epsilon = 0;
		}
		delete[] lower;
		lower = nullptr;
	} else {
		switch (*name) {
			case 'e': // Epsilon filter
			case 'E':
				type = FILTER_ID_EPSILON;
				size = EPSILONFILTER_SIZE;
				std_deviation = .0;
				epsilon = EPSILONFILTER_EPSILON;
				break;
			case 'g': // Gaussian lowpass filter
			case 'G':
				type = FILTER_ID_GAUSSIAN;
				size = GAUSSIAN_SIZE;
				std_deviation = GAUSSIAN_STD_DEVIATION;
				epsilon = 0;
		}
	}
	return true;
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
MULTIPLE_MOTION_PARAM::set_default(const char *name)
{
	if (strcmp(name, "Level") == 0) {
		Level = MULTIPLE_MOTION_PARAM_Level;
	} else if (strcmp(name, "IRLS_Iter_Max") == 0) {
		IRLS_Iter_Max = MULTIPLE_MOTION_PARAM_IRLS_Iter_Max;
	} else if (strcmp(name, "Error_Min_Threshold") == 0) {
		Error_Min_Threshold = MULTIPLE_MOTION_PARAM_Error_Min_Threshold;
	} else if (strcmp(name, "lambdaD") == 0) {
		lambdaD = MULTIPLE_MOTION_PARAM_lambdaD;
	} else if (strcmp(name, "lambdaS") == 0) {
		lambdaS = MULTIPLE_MOTION_PARAM_lambdaS;
	} else if (strcmp(name, "sigmaD") == 0) {
		sigmaD = MULTIPLE_MOTION_PARAM_sigmaD;
	} else if (strcmp(name, "sigmaS") == 0) {
		sigmaS = MULTIPLE_MOTION_PARAM_sigmaS;
	} else {
		fprintf(stderr, "*** MULTIPLE_MOTION_PARAM::set_default() error - There are NOT such a parameter '%s' ***\n", name);
	}
}

void
MULTIPLE_MOTION_PARAM::set_value(const char *name, void *value)
{
	if (strcmp(name, "Level") == 0) {
		Level = *(int *)value;
		if (Level < 1) {
			Level = 1;
		}
	} else if (strcmp(name, "IRLS_Iter_Max") == 0) {
		IRLS_Iter_Max = *(int *)value;
		if (IRLS_Iter_Max < 0) {
			IRLS_Iter_Max = 0;
		}
	} else if (strcmp(name, "Error_Min_Threshold") == 0) {
		Error_Min_Threshold = *(double *)value;
		if (Error_Min_Threshold < .0) {
			Error_Min_Threshold = .0;
		}
	} else if (strcmp(name, "lambdaD") == 0) {
		lambdaD = *(double *)value;
	} else if (strcmp(name, "lambdaS") == 0) {
		lambdaS = *(double *)value;
	} else if (strcmp(name, "sigmaD") == 0) {
		sigmaD = *(double *)value;
	} else if (strcmp(name, "sigmaS") == 0) {
		sigmaS = *(double *)value;
	} else {
		fprintf(stderr, "*** MULTIPLE_MOTION_PARAM::set_default() error - There are NOT such a parameter '%s' ***\n", name);
	}
}

#define HOG_PARAM_Bins 16
#define HOG_PARAM_Dense true
#define HOG_PARAM_SignedOrient true
HOG_PARAM::HOG_PARAM(void)
{
	Bins = HOG_PARAM_Bins;
	Dense = HOG_PARAM_Dense;
	SignedOrient= HOG_PARAM_SignedOrient;
}

void
HOG_PARAM::set_default(const char *name)
{
	if (strcmp(name, "Bins") == 0) {
		Bins = HOG_PARAM_Bins;
	} else if (strcmp(name, "Dense") == 0) {
		Dense = HOG_PARAM_Dense;
	} else if (strcmp(name, "SignedOrient") == 0) {
		SignedOrient = HOG_PARAM_SignedOrient;
	} else {
		fprintf(stderr, "*** HOG_PARAM::set_default() error - There are NOT such a parameter '%s' ***\n", name);
	}
}

void
HOG_PARAM::set_value(const char *name, const void *value)
{
	if (strcmp(name, "Bins") == 0) {
		Bins = *(int *)value;
		if (Bins < 1) {
			Bins = 1;
		}
	} else if (strcmp(name, "Dense") == 0) {
		Dense = *(bool *)value;
	} else if (strcmp(name, "SignedOrient") == 0) {
		SignedOrient = *(bool *)value;
	} else {
		fprintf(stderr, "*** HOG_PARAM::set_default() error - There are NOT such a parameter '%s' ***\n", name);
	}
}


OPTIONS::OPTIONS(void)
{
	ResampleMethod = 0;
	mode = 0;
	Max_Length = 0;
	Max_Output_Length = 0;
	ExclusivePrinciple = false;
	Superimpose = 0;
	PlotOptions = 0;
	s_med = SCRATCH_MED_THRESHOLD;
	s_avg = SCRATCH_AVG_THRESHOLD;
	p = DIR_PROBABILITY;
	ep = EPSILON;
	Exclusive_Max_Radius = EXCLUSIVE_PRINCIPLE_MAX_RADIUS;
}

bool
OPTIONS::ChangeResampleMethod(const char *name)
{
	ERROR Error("OPTIONS::ChangeResampleMethod");
	char *lower = nullptr;

	try {
		lower = new char[strlen(name) + 1];
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("lower");
		Error.Malloc();
		return false;
	}
	for (size_t i = 0; i < strlen(name); i++) {
		lower[i] = tolower(name[i]);
	}
	lower[strlen(lower) + 1] = '\0';
	if (strcmp(lower, "z-hold") == 0) {
		ResampleMethod = PNM_Resize_ZeroOrderHold;
	} else if (strcmp(lower, "bicubic") == 0) {
		ResampleMethod = PNM_Resize_Bicubic;
	}
	delete[] lower;
	return true;
}

void
OPTIONS::set_default(const char *name)
{
	if (strcmp(name, "ResampleSize") == 0) {
		ResampleSize.reset();
	} else if (strcmp(name, "") == 0) {
		mode = 0;
	} else if (strcmp(name, "") == 0) {
		Max_Length = 0;
	} else if (strcmp(name, "") == 0) {
		Max_Output_Length = 0;
	} else if (strcmp(name, "") == 0) {
		ExclusivePrinciple = false;
	} else if (strcmp(name, "") == 0) {
		Superimpose = 0;
	} else if (strcmp(name, "") == 0) {
		PlotOptions = 0;
	} else if (strcmp(name, "") == 0) {
		s_med = SCRATCH_MED_THRESHOLD;
	} else if (strcmp(name, "") == 0) {
		s_avg = SCRATCH_AVG_THRESHOLD;
	} else if (strcmp(name, "") == 0) {
		p = DIR_PROBABILITY;
	} else if (strcmp(name, "") == 0) {
		ep = EPSILON;
	} else if (strcmp(name, "") == 0) {
		Exclusive_Max_Radius = EXCLUSIVE_PRINCIPLE_MAX_RADIUS;
	} else {
		fprintf(stderr, "*** OPTIONS::set_default() error - There are NOT such a parameter '%s' ***\n", name);
	}
}

void
OPTIONS::set_value(const char *name, const void *value)
{
	if (strcmp(name, "ResampleSize") == 0) {
		ResampleSize.set_size((SIZE *)value);
	} else if (strcmp(name, "") == 0) {
		mode = *(int *)value;
	} else if (strcmp(name, "") == 0) {
		Max_Length = *(int *)value;
	} else if (strcmp(name, "") == 0) {
		Max_Output_Length = *(int *)value;
	} else if (strcmp(name, "") == 0) {
		ExclusivePrinciple = *(bool *)value;
	} else if (strcmp(name, "") == 0) {
		Superimpose = *(int *)value;
	} else if (strcmp(name, "") == 0) {
		PlotOptions = *(int *)value;
	} else if (strcmp(name, "") == 0) {
		s_med = *(int *)value;
	} else if (strcmp(name, "") == 0) {
		s_avg = *(int *)value;
	} else if (strcmp(name, "") == 0) {
		p = *(double *)value;
	} else if (strcmp(name, "") == 0) {
		ep = *(double *)value;
	} else if (strcmp(name, "") == 0) {
		Exclusive_Max_Radius = *(double *)value;
	} else {
		fprintf(stderr, "*** OPTIONS::set_default() error - There are NOT such a parameter '%s' ***\n", name);
	}
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

