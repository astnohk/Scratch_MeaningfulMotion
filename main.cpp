#include "Scratch_MeaningfulMotion.h"

const char Progress[NUM_PROGRESS][6] = {
    "\0\0\0\0\0", "\x1b[1C\0", "\x1b[2C\0", "\x1b[3C\0", "\x1b[4C\0",
    "\x1b[5C\0", "\x1b[6C\0", "\x1b[7C\0", "\x1b[8C\0", "\x1b[9C\0",
    "\x1b[10C", "\x1b[11C", "\x1b[12C", "\x1b[13C", "\x1b[14C",
    "\x1b[15C", "\x1b[16C", "\x1b[17C", "\x1b[18C", "\x1b[19C",
    "\x1b[20C", "\x1b[21C", "\x1b[22C", "\x1b[23C", "\x1b[24C",
    "\x1b[25C", "\x1b[26C", "\x1b[27C", "\x1b[28C", "\x1b[29C",
    "\x1b[30C", "\x1b[31C", "\x1b[32C", "\x1b[33C", "\x1b[34C",
    "\x1b[35C", "\x1b[36C", "\x1b[37C", "\x1b[38C", "\x1b[39C",
    "\x1b[40C", "\x1b[41C", "\x1b[42C", "\x1b[43C", "\x1b[44C",
    "\x1b[45C", "\x1b[46C", "\x1b[47C", "\x1b[48C", "\x1b[49C",
    "\x1b[50C", "\x1b[51C", "\x1b[52C", "\x1b[53C", "\x1b[54C",
    "\x1b[55C", "\x1b[56C", "\x1b[57C", "\x1b[58C", "\x1b[59C",
    "\x1b[60C", "\x1b[61C", "\x1b[62C", "\x1b[63C"};
const char Progress_End[] = "\x1b[64C|";




int
main(int argc, char *argv[])
{
	ERROR Error("main");
	const char help[] =
	    "\n"
	    "     - Line Scratch Detection by Meaningful Alignments -\n"
	    "\n"
	    "  Detect Line Scratch from The Picture using Meaningful Alignments.\n"
	    "  It can read tiff and pnm image format.\n"
	    "\n"
	    "\n"
	    "  format : Scratch_MeaningfulA -i files_####.pgm -o files_out_####.pgm -s [start_num] -e [end_num] [option [argument]]\n"
	    "                                        ^Sharp means some numbers.\n"
	    "\n"
	    "\n"
	    "    General option:\n"
	    "      -h, --help                       : show this manual\n"
	    "      -i [input]                       : set input filename\n"
	    "      -o [output]                      : set output filename\n"
	    "      -s [start_num]                   : set output filename\n"
	    "      -e [end_num]                     : set output filename\n"
	    "      --filtered                       : output first filtered image\n"
	    "      --binary                         : output middle data at Line Scratch detection\n"
	    "      --multiple_affine                : output multiple motions' affine parameters estimated by method of M.J.Black\n"
	    "      --affine_blockmatching           : output multiple motions' optical flow estimated by method of M.J.Black\n"
	    "      --opticalflow_blockmatching      : output multiple motions' optical flow estimated by method of M.J.Black\n"
	    "      --mm_level [value]               : set maximum level of multi resolution gaussian pyramid (default : 4)\n"
	    "      --HOG                            : output block normalized Histograms of Oriented Gradients\n"
	    "      --HOG_raw                        : output raw Histograms of Oriented Gradients\n"
	    "      --HOG_matching_vector            : output vector of matching by Histograms of Oriented Gradients\n"
	    "      --HOG_bins [value]               : set the number of bins of histogram of gradients (default : 16)\n"
	    "      --HOG_densely                    : compute HOG feature at each pixels (default)\n"
	    "      --HOG_less_densely               : compute HOG feature at NOT overlapped cells (NOT default)\n"
	    "      --HOG_signed                     : use signed orientation for HOG feature (default)\n"
	    "      --HOG_unsigned                   : neglect sign of orientation for HOG feature (NOT default)\n"
	    "      --resample [(W)x(H)]             : resampling the input image to size of [WxH] before processing (use as : --resample 128x128)\n"
	    "      --resample_method [method]       : set resampling method (z-hold, bicubic)\n"
	    "      --plot_as_resampled              : output size is same as resampled image\n"
	    "      --plot_resampled_only            : output the resampled image without any other processing\n"
	    "      --x11_plot                       : show image and detected lines in 3D world by X11 library\n"
	    "\n"
	    "    Line Scratch Detection option:\n"
	    "      --filter_size [(width)x(height)] : set Filter size [width]x[height] (default value : 21x21)\n"
	    "      --filter_type [name]             : set which Filter will used (Epsilon, Gaussian or None)\n"
	    "      --gauss_var   [value]            : set Gaussian filter's Standard Deviation (default value : 5.0)\n"
	    "      --filter_ep   [value]            : set Epsilon Filter Threshold (default value : 20)\n"
	    "      --s_med       [value]            : set s_med the threshold of difference on center from median (default value : 3)\n"
	    "      --s_avg       [value]            : set s_avg the threshold of difference of averages on both sides of scratch (default value : 20)\n"
	    "\n"
	    "    Meaningful Alignments option:\n"
	    "      -l [value]                       : set maximum length limits of lines when detecting\n"
	    "      -L [value]                       : set maximum length limits of lines when output\n"
	    "      -n                               : set output with negative (fg:black, bg:white)\n"
	    "      --epsilon       [value]          : set epsilon for threshold of NFA(Number of False Alarm)\n"
	    "      --exclusive_rad [value]          : set Exclusive Principle Max Radius (default value : 1.5)\n"
	    "      --exclusive                      : use Exclusive Principle to remove redundant lines\n"
	    "      --superimpose   [color]          : superimpose the segments on original image with colors (red, green, blue)\n"
	    "\n";
	int errors = 0;

	char *delimiter = nullptr;
	char c_tmp = 0;
	std::string InputName;
	std::string OutputName;
	int Start = 0;
	int End = 0;
	OPTIONS Options;
	FILTER_PARAM FilterParam;
	char Superimpose_Color[8];
	const std::string Superimpose_Color_List[OVERLAY_COLOR_PATTERNS] = {"red", "green", "blue"};

	int inf = 0, outf = 0;
	char *strtmp = nullptr;
	char *strdivp = nullptr;
	bool bval;
	int ival;
	double dval;
	SIZE tmpsize;
	int i, k;

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			if (strlen(argv[i]) > 2) {
				if (strcmp(argv[i], "--help") == 0) {
					printf("%s\n", help);
					return EXIT_SUCCESS;
				} else if (strcmp(argv[i], "--binary") == 0) {
					Options.mode |= MODE_OUTPUT_BINARY_IMAGE;
				} else if (strcmp(argv[i], "--epsilon") == 0) {
					if (i + 1 >= argc) {
						fprintf(stderr, "*** Please input value after '--epsilon' option ***\n");
						errors |= OPTION_INSUFFICIENT;
					} else {
						i++;
						sscanf(argv[i], "%10lf", &dval);
						Options.set_value("ep", &dval);
					}
				} else if (strcmp(argv[i], "--exclusive") == 0) {
					Options.ExclusivePrinciple = true;
				} else if (strcmp(argv[i], "--exclusive_rad") == 0) {
					if (i + 1 >= argc) {
						fprintf(stderr, "*** Please input value after '--exclusive_rad' option ***\n");
						errors |= OPTION_INSUFFICIENT;
					} else {
						i++;
						if (sscanf(argv[i], "%10lf", &dval) != 1) {
							fprintf(stderr, "*** Cannot read value for Exclusive Principle Max Radius correctly ***\n*** Use default value instead ***\n");
							errors |= OPTION_INCORRECT;
						}
						Options.set_value("Exclusive_Max_Radius", &dval);
					}
				} else if (strcmp(argv[i], "--filtered") == 0) {
					Options.mode |= MODE_OUTPUT_FILTERED_IMAGE;
				} else if (strcmp(argv[i], "--filter_ep") == 0) {
					if (i + 1 >= argc) {
						fprintf(stderr, "*** Please input value after '--filter_ep' option ***\n");
						errors |= OPTION_INSUFFICIENT;
					} else {
						i++;
						if (sscanf(argv[i], "%10lf", &dval) != 1) {
							fprintf(stderr, "*** Cannot read value for Epsilon Filter's epsilon correctly ***\n");
							errors |= OPTION_INCORRECT;
						}
						FilterParam.set_value("epsilon", &dval);
					}
				} else if (strcmp(argv[i], "--filter_size") == 0) {
					if (i + 1 >= argc) {
						fprintf(stderr, "*** Please input value after '--filter_size' option ***\n");
						errors |= OPTION_INSUFFICIENT;
					} else {
						i++;
						try {
							strtmp = new char[strlen(argv[i])];
						}
						catch (const std::bad_alloc &bad) {
							std::cerr << bad.what() << std::endl;
							Error.Value("strtmp");
							Error.Malloc();
							goto ExitError;
						}
						strcpy(strtmp, argv[i]);
						strdivp = strchr(strtmp, 'x');
						if (strdivp != nullptr) {
							if (sscanf(strdivp + 1, "%7d", &tmpsize.height) != 1) {
								Error.Function("sscanf");
								Error.Value("tmpsize.height");
								Error.FunctionFail();
								goto ExitError;
							}
							*strdivp = '\0';
							strdivp = nullptr;
							if (sscanf(strtmp, "%7d", &tmpsize.width) != 1) {
								Error.Function("sscanf");
								Error.Value("tmpsize.width");
								Error.FunctionFail();
								goto ExitError;
							}
						} else {
							if (sscanf(argv[i], "%7d", &tmpsize.width) != 1) {
								fprintf(stderr, "*** Cannot read value for Filter Size correctly ***\n*** Use default value instead ***\n");
								errors |= OPTION_INCORRECT;
							}
							tmpsize.height = tmpsize.width;
						}
						delete[] strtmp;
						strtmp = nullptr;
						FilterParam.set_value("size", &tmpsize);
					}
				} else if (strcmp(argv[i], "--filter_type") == 0) {
					if (i + 1 >= argc) {
						fprintf(stderr, "*** Please input value after '--filter_type' option ***\n");
						errors |= OPTION_INSUFFICIENT;
					} else {
						i++;
						FilterParam.ChangeFilter(argv[i]);
					}
				} else if (strcmp(argv[i], "--gauss_stddev") == 0) {
					if (i + 1 >= argc) {
						fprintf(stderr, "*** Please input value after '--gauss_stddev' option ***\n");
						errors |= OPTION_INSUFFICIENT;
					} else {
						i++;
						if (sscanf(argv[i], "%10lf", &dval) != 1) {
							fprintf(stderr, "*** Cannot read value for Gaussian Filter Variance correctly ***\n*** Use default value instead ***\n");
							errors |= OPTION_INCORRECT;
						}
						FilterParam.set_value("std_deviation", &dval);
					}
				} else if (strcmp(argv[i], "--HOG") == 0) {
					Options.mode = MODE_OUTPUT_HISTOGRAMS_OF_ORIENTED_GRADIENTS;
				} else if (strcmp(argv[i], "--HOG_raw") == 0) {
					Options.mode = MODE_OUTPUT_HISTOGRAMS_OF_ORIENTED_GRADIENTS_RAW_HOG;
				} else if (strcmp(argv[i], "--HOG_matching_vector") == 0) {
					Options.mode = MODE_OUTPUT_HISTOGRAMS_OF_ORIENTED_GRADIENTS_MATCHING_VECTOR;
				} else if (strcmp(argv[i], "--HOG_bins") == 0) {
					if (i + 1 >= argc) {
						fprintf(stderr, "*** Please input value after '--HOG_bins' option ***\n");
						errors |= OPTION_INSUFFICIENT;
					} else {
						i++;
						if (sscanf(argv[i], "%7d", &ival) != 1) {
							fprintf(stderr, "*** Cannot read value for the number of Bins of histograms ***\n*** Use default value instead ***\n");
							errors |= OPTION_INCORRECT;
						}
						Options.HOG_Param.set_value("Bins", &ival);
					}
				} else if (strcmp(argv[i], "--HOG_densely") == 0) {
					Options.HOG_Param.Dense = true;
				} else if (strcmp(argv[i], "--HOG_less_densely") == 0) {
					Options.HOG_Param.Dense = false;
				} else if (strcmp(argv[i], "--HOG_signed") == 0) {
					Options.HOG_Param.SignedOrient = true;
				} else if (strcmp(argv[i], "--HOG_unsigned") == 0) {
					Options.HOG_Param.SignedOrient = false;
				} else if (strcmp(argv[i], "--multiple_affine") == 0) {
					Options.mode = MODE_OUTPUT_MULTIPLE_MOTIONS_AFFINE;
				} else if (strcmp(argv[i], "--affine_blockmatching") == 0) {
					Options.mode = MODE_OUTPUT_AFFINE_BLOCKMATCHING;
				} else if (strcmp(argv[i], "--opticalflow_blockmatching") == 0) {
					Options.mode = MODE_OUTPUT_OPTICALFLOW;
				} else if (strcmp(argv[i], "--mm_level") == 0) {
					if (i + 1 >= argc) {
						fprintf(stderr, "*** Please input value after '--mm_level' option ***\n");
						errors |= OPTION_INSUFFICIENT;
					} else {
						i++;
						if (sscanf(argv[i], "%7d", &ival) != 1) {
							fprintf(stderr, "*** Cannot read value for the maximum level of Multi-Resolution Pyramid ***\n*** Use default value instead ***\n");
							errors |= OPTION_INCORRECT;
						}
						Options.MultipleMotion_Param.set_value("Level", &ival);
					}
				} else if (strcmp(argv[i], "--plot_as_resampled") == 0) {
					Options.PlotOptions |= PLOT_AS_RESAMPLED;
				} else if (strcmp(argv[i], "--plot_resampled_only") == 0) {
					Options.PlotOptions |= PLOT_RESAMPLED_IMG_ONLY;
				} else if (strcmp(argv[i], "--resample") == 0) {
					if (i + 1 >= argc) {
						fprintf(stderr, "*** Please input Width, Height or both after '--resample' option (e.g. : --resample 128x128) ***\n");
						errors |= OPTION_INSUFFICIENT;
					} else {
						i++;
						if ((delimiter = strchr(argv[i], 'x')) == nullptr) {
							if (sscanf(argv[i], "%7d", &tmpsize.height) != 1) {
								fprintf(stderr, "*** Cannot read value for ResampleSize.height ***\n");
								errors |= OPTION_INCORRECT;
							}
							tmpsize.width = tmpsize.height;
						} else {
							c_tmp = *delimiter;
							*delimiter = '\0'; // set end point of width
							if (sscanf(argv[i], "%7d", &tmpsize.width) != 1) {
								fprintf(stderr, "*** Cannot read value for ResampleSize.width ***\n");
								errors |= OPTION_INCORRECT;
							}
							if (sscanf(delimiter + 1, "%7d", &tmpsize.height) != 1) {
								fprintf(stderr, "*** Cannot read value for ResampleSize.height ***\n");
								errors |= OPTION_INCORRECT;
							}
							*delimiter = c_tmp; // recover original character
						}
						Options.set_value("ResampleSize", &tmpsize);
					}
				} else if (strcmp(argv[i], "--resample_method") == 0) {
					if (i + 1 >= argc) {
						fprintf(stderr, "*** Please input method name after '--resample_method' option ***\n");
						errors |= OPTION_INSUFFICIENT;
					} else {
						i++;
						if (strlen(argv[i]) < 1) {
							fprintf(stderr, "*** '%s' is NOT method name ***\n", argv[i]);
							errors |= OPTION_INCORRECT;
						} else {
							Options.ChangeResampleMethod(argv[i]);
						}
					}
				} else if (strcmp(argv[i], "--s_avg") == 0) { // set s_avg
					if (i + 1 >= argc) {
						fprintf(stderr, "*** Please input color name after '--s_avg' option ***\n");
						errors |= OPTION_INSUFFICIENT;
					} else {
						i++;
						if (sscanf(argv[i], "%7d", &ival) != 1) {
							fprintf(stderr, "*** Cannot read the value of 's_avg' correctly ***\n");
							errors |= OPTION_INCORRECT;
						}
						Options.set_value("s_avg", &ival);
					}
				} else if (strcmp(argv[i], "--s_med") == 0) { // set s_med
					if (i + 1 >= argc) {
						fprintf(stderr, "*** Please input color name after '--s_med' option ***\n");
						errors |= OPTION_INSUFFICIENT;
					} else {
						i++;
						if (sscanf(argv[i], "%7d", &ival) != 1) {
							fprintf(stderr, "*** Cannot read the value of 's_med' correctly ***\n");
							errors |= OPTION_INCORRECT;
						}
						Options.set_value("s_med", &ival);
					}
				} else if (strcmp(argv[i], "--superimpose") == 0) { // Color output superimpose
					if (i + 1 >= argc) {
						fprintf(stderr, "*** Please input color name after '--superimpose' option ***\n       - e.g. --superimpose green\n");
						errors |= OPTION_INSUFFICIENT;
					} else {
						i++;
						if (sscanf(argv[i], "%7s", Superimpose_Color) != 1) {
							fprintf(stderr, "*** Cannot read color name for superimpose output color correctly ***\n");
							errors |= OPTION_INCORRECT;
						}
						for (k = 0; k < OVERLAY_COLOR_PATTERNS; k++){
							if (strcmp(argv[i], Superimpose_Color_List[k].c_str()) == 0) {
								Options.Superimpose = k + 1;
								break;
							}
						}
						if (k >= OVERLAY_COLOR_PATTERNS) {
							fprintf(stderr, "*** Color name is incorrect (%s) ***\n", argv[i]);
							errors |= OPTION_INCORRECT;
						}
					}
				} else if (strcmp(argv[i], "--x11_plot") == 0) { // X11 plotting mode
					bval = true;
					Options.set_value("x11_plot", &bval);
				} else {
					fprintf(stderr, " *** Unknown option \"%s\" ***\n      - Please begin with \"--\" for long name options\n", argv[i]);
					errors |= OPTION_UNKNOWN;
				}
			} else {
				switch (argv[i][1]) {
					case 'e': // End Number
						if (i + 1 >= argc) {
							fprintf(stderr, "*** Please input value after '-e' option ***\n");
							errors |= OPTION_INSUFFICIENT;
						} else {
							i++;
							sscanf(argv[i], "%10d", &End);
							if (End <= 0) {
								End = 0;
							}
						}
						break;
					case 'h':
						printf("%s\n", help);
						return EXIT_SUCCESS;
					case 'i':
						if (i + 1 >= argc) {
							fprintf(stderr, "*** Please input INPUT FILENAME after '-i' option ***\n");
							errors |= OPTION_INSUFFICIENT;
						} else {
							i++;
							inf = i;
						}
						break;
					case 'L':
						if ((i + 1) >= argc) {
							fprintf(stderr, "*** Please input value after '-l' option ***\n");
							errors |= OPTION_INSUFFICIENT;
						} else {
							i++;
							if (sscanf(argv[i], "%7d", &ival) != 1) {
								fprintf(stderr, "*** Cannot read value for Maximum Length Limits of Output Segments correctly ***\n*** Use default value instead ***\n");
								errors |= OPTION_INCORRECT;
							}
							Options.set_value("Max_Output_Length", &ival);
						}
						break;
					case 'l':
						if ((i + 1) >= argc) {
							fprintf(stderr, "*** Please input value after '-l' option ***\n");
							errors |= OPTION_INSUFFICIENT;
						} else {
							i++;
							if (sscanf(argv[i], "%7d", &ival) != 1) {
								fprintf(stderr, "*** Cannot read value for Maximum Length Limits of Searching Segments correctly ***\n*** Use default value instead ***\n");
								errors |= OPTION_INCORRECT;
							}
							Options.set_value("Max_Length", &ival);
						}
						break;
					case 'n':
						Options.PlotOptions |= PLOT_NEGATE;
						break;
					case 'o':
						if (i + 1 >= argc) {
							fprintf(stderr, "*** Please input OUTPUT FILENAME after '-o' option ***\n");
							errors |= OPTION_INSUFFICIENT;
						} else {
							i++;
							outf = i;
						}
						break;
					case 's': // Start Number
						if (i + 1 >= argc) {
							fprintf(stderr, "*** Please input value after '-s' option ***\n");
							errors |= OPTION_INSUFFICIENT;
						} else {
							i++;
							sscanf(argv[i], "%10d", &Start);
							if (Start <= 0) {
								Start = 0;
							}
						}
						break;
					default:
						fprintf(stderr, "*** Unknown option \"%s\" ***\n", argv[i]);
						errors |= OPTION_UNKNOWN;
				}
			}
		}
	}
	if (End == 0) {
		End = Start;
	}
	printf("\n");
	if ((errors & OPTION_INSUFFICIENT) != 0) {
		fprintf(stderr, "*** FATAL main error - the last option needs argument ***\n");
		exit(EXIT_FAILURE);
	} else if ((inf == 0) || (outf == 0)) {
	 	fprintf(stderr, "*** FATAL main error - Cannot find INPUT or OUTPUT Filename ***\n");
		exit(EXIT_FAILURE);
	} else if (Start > End) {
		fprintf(stderr, "*** FATAL main error - Start number exceeds End number (Start : %d, End : %d) ***\n", Start, End);
		exit(EXIT_FAILURE);
	}
	if ((errors & (OPTION_INCORRECT | OPTION_UNKNOWN)) != 0) {
		if ((errors & OPTION_INCORRECT) != 0) {
			fprintf(stderr, "*** main error - Some options are used incorrectly ***\n");
		}
		if ((errors & OPTION_UNKNOWN) != 0) {
			fprintf(stderr, "*** main error - Found some unknown options ***\n");
		}
		exit(EXIT_FAILURE);
	}

	if (Options.Superimpose != 0) {
		strcpy(argv[outf] + strlen(argv[outf]) - PPM_EXTENSION_LENGTH, "ppm");
		printf("Output as %s line Superimpose on Original image\nChange output filename to \"%s\" due to --superimpose option\n", Superimpose_Color_List[(Options.Superimpose - 1) % OVERLAY_COLOR_PATTERNS].c_str(), argv[outf]);
	}
	// Regular expression
	InputName = argv[inf];
	OutputName = argv[outf];
	// main routine
	try {
		Scratch_MeaningfulMotion(OutputName, InputName, Start, End, Options, FilterParam);
	}
	catch (const std::logic_error& logic) {
		std::cerr << logic.what() << std::endl;
		goto ExitError;
	}
	catch (const std::runtime_error& runtime) {
		std::cerr << runtime.what() << std::endl;
		goto ExitError;
	}
	return EXIT_SUCCESS;
// Error
ExitError:
	std::cout << std::endl
	    << "        *** FATAL main error ***" << std::endl;
	return EXIT_FAILURE;
}
