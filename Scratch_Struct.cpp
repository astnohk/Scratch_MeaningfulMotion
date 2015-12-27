#include "Scratch_MeaningfulMotion.h"
#include "Plot/Plot_X11.h"



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
		type = *static_cast<const int*>(value);
		if (type < 0 || type >= NUM_FILTER_TYPE) {
			type = 0;
		}
	} else if (strcmp(name, "std_deviation") == 0) {
		std_deviation = *static_cast<const double*>(value);
		if (std_deviation < .0) {
			std_deviation = .0;
		}
	} else if (strcmp(name, "epsilon") == 0) {
		epsilon = *static_cast<const double*>(value);
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

	if (strlen(name) > 1) {
		char *lower = nullptr;
		try {
			lower = new char[strlen(name) + 1];
		}
		catch (const std::bad_alloc &bad) {
			std::cerr << bad.what() << std::endl;
			Error.Value("lower");
			Error.Malloc();
			return false;
		}
		for (size_t i = 0; i < strlen(name); i++) {
			lower[i] = char(tolower(name[i]));
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




MULTIPLE_MOTION_PARAM::MULTIPLE_MOTION_PARAM(void)
{
	Level = MULTIPLE_MOTION_PARAM_Level;
	IRLS_Iter_Max = MULTIPLE_MOTION_PARAM_IRLS_Iter_Max;
	Error_Min_Threshold = MULTIPLE_MOTION_PARAM_Error_Min_Threshold;
	lambdaD = MULTIPLE_MOTION_PARAM_lambdaD;
	lambdaS = MULTIPLE_MOTION_PARAM_lambdaS;
	sigmaD = MULTIPLE_MOTION_PARAM_sigmaD;
	sigmaS = MULTIPLE_MOTION_PARAM_sigmaS;
	BlockMatching_BlockSize = MULTIPLE_MOTION_PARAM_BlockMatching_BlockSize;
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
	} else if (strcmp(name, "BlockMatching_BlockSize") == 0) {
		BlockMatching_BlockSize = MULTIPLE_MOTION_PARAM_BlockMatching_BlockSize;
	} else {
		fprintf(stderr, "*** MULTIPLE_MOTION_PARAM::set_default() error - There are NOT such a parameter '%s' ***\n", name);
	}
}

void
MULTIPLE_MOTION_PARAM::set_value(const char *name, void *value)
{
	if (strcmp(name, "Level") == 0) {
		Level = *static_cast<const int*>(value);
	} else if (strcmp(name, "IRLS_Iter_Max") == 0) {
		IRLS_Iter_Max = *static_cast<const int*>(value);
		if (IRLS_Iter_Max < 0) {
			IRLS_Iter_Max = 0;
		}
	} else if (strcmp(name, "Error_Min_Threshold") == 0) {
		Error_Min_Threshold = *static_cast<const double*>(value);
		if (Error_Min_Threshold < .0) {
			Error_Min_Threshold = .0;
		}
	} else if (strcmp(name, "lambdaD") == 0) {
		lambdaD = *static_cast<const double*>(value);
	} else if (strcmp(name, "lambdaS") == 0) {
		lambdaS = *static_cast<const double*>(value);
	} else if (strcmp(name, "sigmaD") == 0) {
		sigmaD = *static_cast<const double*>(value);
	} else if (strcmp(name, "sigmaS") == 0) {
		sigmaS = *static_cast<const double*>(value);
	} else if (strcmp(name, "BlockMatching_BlockSize") == 0) {
		BlockMatching_BlockSize = *static_cast<const int*>(value);
	} else {
		fprintf(stderr, "*** MULTIPLE_MOTION_PARAM::set_default() error - There are NOT such a parameter '%s' ***\n", name);
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
	x11_plot = false;
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
		std::cerr << bad.what() << std::endl;
		Error.Value("lower");
		Error.Malloc();
		return false;
	}
	for (size_t i = 0; i < strlen(name); i++) {
		lower[i] = char(tolower(name[i]));
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
	} else if (strcmp(name, "mode") == 0) {
		mode = 0;
	} else if (strcmp(name, "Max_Length") == 0) {
		Max_Length = 0;
	} else if (strcmp(name, "Max_Output_Length") == 0) {
		Max_Output_Length = 0;
	} else if (strcmp(name, "ExclusivePrinciple") == 0) {
		ExclusivePrinciple = false;
	} else if (strcmp(name, "Superimpose") == 0) {
		Superimpose = 0;
	} else if (strcmp(name, "PlotOptions") == 0) {
		PlotOptions = 0;
	} else if (strcmp(name, "s_med") == 0) {
		s_med = SCRATCH_MED_THRESHOLD;
	} else if (strcmp(name, "s_avg") == 0) {
		s_avg = SCRATCH_AVG_THRESHOLD;
	} else if (strcmp(name, "p") == 0) {
		p = DIR_PROBABILITY;
	} else if (strcmp(name, "ep") == 0) {
		ep = EPSILON;
	} else if (strcmp(name, "Exclusive_Max_Radius") == 0) {
		Exclusive_Max_Radius = EXCLUSIVE_PRINCIPLE_MAX_RADIUS;
	} else if (strcmp(name, "x11_plot") == 0) {
		x11_plot = false;
	} else {
		fprintf(stderr, "*** OPTIONS::set_default() error - There are NOT such a parameter '%s' ***\n", name);
	}
}

void
OPTIONS::set_value(const char *name, const void *value)
{
	if (strcmp(name, "ResampleSize") == 0) {
		ResampleSize.set_size(static_cast<const SIZE*>(value));
	} else if (strcmp(name, "mode") == 0) {
		mode = *static_cast<const int*>(value);
	} else if (strcmp(name, "Max_Length") == 0) {
		Max_Length = *static_cast<const int*>(value);
	} else if (strcmp(name, "Max_Output_Length") == 0) {
		Max_Output_Length = *static_cast<const int*>(value);
	} else if (strcmp(name, "ExclusivePrinciple") == 0) {
		ExclusivePrinciple = *static_cast<const bool*>(value);
	} else if (strcmp(name, "Superimpose") == 0) {
		Superimpose = *static_cast<const int*>(value);
	} else if (strcmp(name, "PlotOptions") == 0) {
		PlotOptions = *static_cast<const int*>(value);
	} else if (strcmp(name, "s_med") == 0) {
		s_med = *static_cast<const int*>(value);
	} else if (strcmp(name, "s_avg") == 0) {
		s_avg = *static_cast<const int*>(value);
	} else if (strcmp(name, "p") == 0) {
		p = *static_cast<const double*>(value);
	} else if (strcmp(name, "ep") == 0) {
		ep = *static_cast<const double*>(value);
	} else if (strcmp(name, "Exclusive_Max_Radius") == 0) {
		Exclusive_Max_Radius = *static_cast<const double*>(value);
	} else if (strcmp(name, "x11_plot") == 0) {
		x11_plot = *static_cast<const bool*>(value);
	} else {
		fprintf(stderr, "*** OPTIONS::set_default() error - There are NOT such a parameter '%s' ***\n", name);
	}
}

