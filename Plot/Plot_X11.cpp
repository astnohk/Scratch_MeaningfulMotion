#include "../Scratch_MeaningfulMotion.h"
#include "Plot_X11.h"

/* Coordinate System
 *
 *    Z o----> X
 *      |
 *      |
 *      V
 *      Y
 */


const char *ProgramName = "Scratch_MeaningfulMotion";
SIZE Window_size(WINDOW_X_DEFAULT, WINDOW_Y_DEFAULT);
Display *disp = nullptr;
Window win;
static Atom atom1, atom2;
Pixmap pix;
GC GCmono;
GC GCcol[RGB_COLOR];
GC GCcol_dark[RGB_COLOR];
Colormap cmap;

double cos_a[ROTATE_ANGLE_MAX];
double sin_a[ROTATE_ANGLE_MAX];


// Plot Modes
#define NUMBER_OF_MODE 6
#define NUMBER_OF_MODE_SELECTED_BY_SWITCH 4

#define X11_Plot_Points 0
#define X11_Plot_Points_And_Segments 1
#define X11_Plot_Grid 2
#define X11_Plot_Grid_And_Segments 3
#define X11_Plot_Garaxy 4
#define X11_Plot_GravityCorrupt 5
static const char* Plot_Mode[NUMBER_OF_MODE] = {
    "Points",
    "Points and Segments",
    "Grid",
    "Grid and Segments",
    "Garaxy",
    "Gravity Corruption"
};




void
freeXWindow(void)
{
	if (disp == nullptr) {
		return;
	}
	XFreeGC(disp, GCmono);
	for (int k = 0; k < RGB_COLOR; k++) {
		XFreeGC(disp, GCcol[k]);
		XFreeGC(disp, GCcol_dark[k]);
	}
	XFreePixmap(disp, pix);
	XFreeColormap(disp, cmap);
	XDestroyWindow(disp, win);
	XCloseDisplay(disp);
}


void
ShowSegments_X11(const ImgVector<pnm_img>* Img, const SIZE& Img_size_resample, const int MaxInt, const SEGMENT* segments, const unsigned int Num_Segments)
{
	ERROR Error("ShowBounds_X11");

	X11_PARAM X11_Param;
	ImgVector<COORDINATE_3D> *Img_coord = nullptr;
	ImgVector<COORDINATE_3D> *Img_vel = nullptr;
	ImgVector<XPLOT> *Img_plot = nullptr;
	SEGMENT_X11 *segments_plot = nullptr;
	size_t *Img_index = nullptr;
	size_t *Img_index_tmp = nullptr;
	SIZE Img_size;

	COORDINATE_3D GaraxyCenter;
	int loop = 0;
	int event_state = 0;
	int cur_mode = 0;
	double r;
	int m, n;
	double x, y, z;
	double X, Y;
	unsigned int k;

	if (MaxInt < 1) {
		Error.Value("MaxInt");
		Error.ValueIncorrect();
		freeXWindow();
		throw std::invalid_argument("const int MaxInt");
	}
	Img_size.width = Img->width();
	Img_size.height = Img->height();

	// Initialize X11 Window
	try {
		Init_X11(&X11_Param, Img_size);
	}
	catch (const std::invalid_argument& arg) {
		std::cerr << arg.what() << std::endl;
		Error.Function("Init_X11");
		Error.FunctionFail();
		freeXWindow();
		throw std::invalid_argument("const int MaxInt");
	}

	// Memory Allocation
	try {
		Img_plot = new ImgVector<XPLOT>(Img->width(), Img->height());
	}
	catch (const std::bad_alloc& bad) {
		std::cerr << bad.what() << std::endl;
		Error.Function("new");
		Error.Value("Img_plot");
		Error.Malloc();
		freeXWindow();
		throw;
	}
	try {
		segments_plot = new SEGMENT_X11[Num_Segments];
	}
	catch (const std::bad_alloc& bad) {
		std::cerr << bad.what() << std::endl;
		Error.Function("new");
		Error.Value("segments_plot");
		Error.Malloc();
		delete Img_plot;
		freeXWindow();
		throw;
	}
	try {
		Img_coord = new ImgVector<COORDINATE_3D>(Img->width(), Img->height());
	}
	catch (const std::bad_alloc& bad) {
		std::cerr << bad.what() << std::endl;
		Error.Function("new");
		Error.Value("Img_coord");
		Error.Malloc();
		delete Img_plot;
		delete[] segments_plot;
		freeXWindow();
		throw;
	}
	try {
		Img_vel = new ImgVector<COORDINATE_3D>(Img->width(), Img->height());
	}
	catch (const std::bad_alloc& bad) {
		std::cerr << bad.what() << std::endl;
		Error.Function("new");
		Error.Value("Img_vel");
		Error.Malloc();
		delete Img_plot;
		delete[] segments_plot;
		delete Img_coord;
		freeXWindow();
		throw;
	}
	try {
		Img_index = new size_t[Img->size()];
	}
	catch (const std::bad_alloc& bad) {
		std::cerr << bad.what() << std::endl;
		Error.Function("new");
		Error.Value("Img_index");
		Error.Malloc();
		delete Img_plot;
		delete[] segments_plot;
		delete Img_coord;
		delete Img_vel;
		freeXWindow();
		throw;
	}
	try {
		Img_index_tmp = new size_t[Img->size()];
	}
	catch (const std::bad_alloc& bad) {
		std::cerr << bad.what() << std::endl;
		Error.Function("new");
		Error.Value("Img_index_tmp");
		Error.Malloc();
		delete Img_plot;
		delete[] segments_plot;
		delete Img_coord;
		delete Img_vel;
		delete[] Img_index;
		freeXWindow();
		throw;
	}

	GaraxyCenter.set(X11_Param.Center_x, X11_Param.Center_y, X11_Param.Center_z);
	// Infinite Loop
	loop = 1;
	while (loop != 0) {
		event_state = XEventor(&X11_Param, Img_size);
		if (event_state == X_ESCAPE) {
			printf("EXIT X11 Plotting\n");
			loop = 0;
		}
		SwitchEventer(&X11_Param);
		if (cur_mode != X11_Param.ModeSwitch) {
			cur_mode = X11_Param.ModeSwitch;
			// Set initial velocity and position
			switch (cur_mode) {
				case X11_Plot_Garaxy:
					GaraxyCenter.x = X11_Param.Center_x;
					GaraxyCenter.x = X11_Param.Center_y;
					GaraxyCenter.x = X11_Param.Center_z;
					for (m = 0; m < Img->height(); m++) {
						y = m - GaraxyCenter.y;
						for (n = 0; n < Img->width(); n++) {
							x = n - GaraxyCenter.x;
							z = (Img->get(n, m) - MaxInt / 2.0) * X11_Param.Plot_Z_Scale * 2.0;
							Img_coord->at(n, m).x = n;
							Img_coord->at(n, m).y = m;
							Img_coord->at(n, m).z = z;
							r = sqrt(POW2(x) + POW2(y) + POW2(z));
							if (r <= 1.0E-3) {
								r = 1.0E-3;
							}
							if (fabs(z) > 1.0E-6) {
								X = -(y / z + z * y);
								Y = z * x + x / z;
								Img_vel->at(n, m).x = X / sqrt(POW2(X) + POW2(Y)) / sqrt(r);
								Img_vel->at(n, m).y = Y / sqrt(POW2(X) + POW2(Y)) / sqrt(r);
							} else {
								Img_vel->at(n, m).x = -y / r / sqrt(r);
								Img_vel->at(n, m).y = -x / r / sqrt(r);
							}
							if (z < -1.0E-6) {
								Img_vel->at(n, m).x *= -1.0;
								Img_vel->at(n, m).y *= -1.0;
							}
							Img_vel->at(n, m).z = 0.0;
						}
					}
					break;
				case X11_Plot_GravityCorrupt:
					for (m = 0; m < Img_size.height; m++) {
						for (n = 0; n < Img_size.width; n++) {
							Img_coord->at(n, m).x = n;
							Img_coord->at(n, m).y = m;
							Img_coord->at(n, m).z = (Img->get(n, m) - MaxInt / 2.0) * X11_Param.Plot_Z_Scale * 2.0;
							Img_vel->at(n, m).x = 0.0;
							Img_vel->at(n, m).y = 0.0;
							Img_vel->at(n, m).z = 0.0;
						}
					}
			}
		}
		try {
			switch (X11_Param.ModeSwitch) {
				case X11_Plot_Points: // Plot Image Intensity with Points
					TransRotate_3DSegment(X11_Param, segments, segments_plot, Num_Segments, Img_size, Img_size_resample);
					TransRotate_3DPoint(X11_Param, Img, MaxInt, Img_plot);
					reset_index(Img_index, Img->size());
					sort_index(Img_plot, Img_index, Img_index_tmp, Img->size());
					Plot_3DPoints(X11_Param, Img, Img_plot, Img_index);
					break;
				case X11_Plot_Points_And_Segments: // Plot Image Intensity with Points
					TransRotate_3DSegment(X11_Param, segments, segments_plot, Num_Segments, Img_size, Img_size_resample);
					TransRotate_3DPoint(X11_Param, Img, MaxInt, Img_plot);
					reset_index(Img_index, Img->size());
					sort_index(Img_plot, Img_index, Img_index_tmp, Img->size());
					Plot_3DPoints(X11_Param, Img, Img_plot, Img_index);
					Plot_3DSegment(X11_Param, segments_plot, Num_Segments);
					break;
				case X11_Plot_Grid: // Plot Image Intensity on Grid
					TransRotate_3DPoint(X11_Param, Img, MaxInt, Img_plot);
					reset_index(Img_index, Img->size());
					sort_index(Img_plot, Img_index, Img_index_tmp, Img->size());
					Plot_3DGrid(X11_Param, Img, Img_plot, Img_index);
					break;
				case X11_Plot_Grid_And_Segments: // Plot Image Intensity on Grid
					TransRotate_3DSegment(X11_Param, segments, segments_plot, Num_Segments, Img_size, Img_size_resample);
					TransRotate_3DPoint(X11_Param, Img, MaxInt, Img_plot);
					reset_index(Img_index, Img->size());
					sort_index(Img_plot, Img_index, Img_index_tmp, Img->size());
					Plot_3DGrid(X11_Param, Img, Img_plot, Img_index);
					Plot_3DSegment(X11_Param, segments_plot, Num_Segments);
					break;
				case X11_Plot_Garaxy: // Plot Image Gravity Motion (Garaxy)
					TransGaraxy_3DPoint(X11_Param, Img, Img_coord, Img_vel, GaraxyCenter, Img_plot);
					reset_index(Img_index, Img->size());
					sort_index(Img_plot, Img_index, Img_index_tmp, Img->size());
					Plot_3DPoints(X11_Param, Img, Img_plot, Img_index);
					break;
				case X11_Plot_GravityCorrupt: // Plot Image Gravity Motion
					TransGravity_3DPoint(X11_Param, Img, Img_coord, Img_vel, Img_plot);
					reset_index(Img_index, Img->size());
					sort_index(Img_plot, Img_index, Img_index_tmp, Img->size());
					Plot_3DPoints(X11_Param, Img, Img_plot, Img_index);
			}
		}
		catch (const std::invalid_argument& arg) {
			std::cerr << arg.what() << std::endl;
			delete Img_plot;
			delete[] segments_plot;
			delete Img_coord;
			delete Img_vel;
			delete[] Img_index;
			delete[] Img_index_tmp;
			throw;
		}
		PlotParameters(X11_Param);
		Set_Pixmap2Window();
		usleep(WAIT_TIME);
	}

	// Freeing Resources
	XFreeGC(disp, GCmono);
	for (k = 0; k < RGB_COLOR; k++) {
		XFreeGC(disp, GCcol[k]);
		XFreeGC(disp, GCcol_dark[k]);
	}
	XFreePixmap(disp, pix);
	XFreeColormap(disp, cmap);
	XDestroyWindow(disp, win);
	XCloseDisplay(disp);
	delete[] Img_index_tmp;
	delete[] Img_index;
	delete Img_vel;
	delete Img_coord;
	delete Img_plot;
	delete[] segments_plot;
}


void
Init_X11(X11_PARAM *X11_Param, SIZE Img_size)
{
	ERROR Error("Init_X11()");
	XColor col, exact;
	XEvent noev;

	X11_Param->Int_interval = DEFAULT_INTENSITY_INTERVAL;
	X11_Param->Latitude = 0;
	X11_Param->Longitude = 0;
	X11_Param->Center_x = Img_size.width / 2;
	X11_Param->Center_y = Img_size.height / 2;
	X11_Param->Center_z = .0;
	X11_Param->Scale = Window_size.width * 0.8 / Img_size.width;
	X11_Param->Plot_Z_Scale = DEFAULT_PLOT_Z_SCALE;
	X11_Param->RotateSwitch = 0;
	X11_Param->ModeSwitch = 0;
	X11_Param->FillSwitch = 0;

	if ((disp = XOpenDisplay(nullptr)) == nullptr) {
		Error.Others("Cannot open display. Abort X11 Plotting");
		throw std::runtime_error("XOpenDisplay(nullptr)");
	}
	win = XCreateSimpleWindow(disp, RootWindow(disp, 0), 0, 0, static_cast<unsigned int>(Window_size.width), static_cast<unsigned int>(Window_size.height), 0, BlackPixel(disp, 0), WhitePixel(disp, 0));
	XSelectInput(disp, win, ExposureMask | StructureNotifyMask | KeyPressMask | ButtonPressMask | ButtonMotionMask);
	XStoreName(disp, win, ProgramName);
	XMapWindow(disp, win);

	// Set behavior when the window close button is pressed
	atom1 = XInternAtom(disp, "WM_PROTOCOLS", False);
	atom2 = XInternAtom(disp, "WM_DELETE_WINDOW", False);
	XSetWMProtocols(disp, win, &atom2, 1);

	// Graphic Context
	GCmono = XCreateGC(disp, win, 0, 0);
	for (int i = 0; i < RGB_COLOR; i++) {
		GCcol[i] = XCreateGC(disp, win, 0, 0);
		GCcol_dark[i] = XCreateGC(disp, win, 0, 0);
	}
	// set GC monochrome
	XSetForeground(disp, GCmono, BlackPixel(disp, 0));
	// Colormap
	cmap = DefaultColormap(disp, 0);
	// set GC color
	if (!(XAllocNamedColor(disp, cmap, "Red", &col, &exact))) {
		Error.Function("XAllocNamedColor");
		Error.Value("cmap (Red)");
		Error.FunctionFail();
		throw std::runtime_error("XAllocNamedColor(disp, cmap, \"Red\", &col, &exact)");
	}
	XSetForeground(disp, GCcol[0], col.pixel);
	if (!(XAllocNamedColor(disp, cmap, "Green", &col, &exact))) {
		Error.Function("XAllocNamedColor");
		Error.Value("cmap (Green)");
		Error.FunctionFail();
		throw std::runtime_error("XAllocNamedColor(disp, cmap, \"Green\", &col, &exact)");
	}
	XSetForeground(disp, GCcol[1], col.pixel);
	if (!(XAllocNamedColor(disp, cmap, "Blue", &col, &exact))) {
		Error.Function("XAllocNamedColor");
		Error.Value("cmap (Blue)");
		Error.FunctionFail();
		throw std::runtime_error("XAllocNamedColor(disp, cmap, \"Blue\", &col, &exact)");
	}
	XSetForeground(disp, GCcol[2], col.pixel);
	// set GC color dark
	if (!(XAllocNamedColor(disp, cmap, "Red", &col, &exact))) {
		Error.Function("XAllocNamedColor");
		Error.Value("cmap (Red)");
		Error.FunctionFail();
		throw std::runtime_error("XAllocNamedColor(disp, cmap, \"Red\", &col, &exact)");
	}
	XSetForeground(disp, GCcol_dark[0], col.pixel);
	if (!(XAllocNamedColor(disp, cmap, "Green", &col, &exact))) {
		Error.Function("XAllocNamedColor");
		Error.Value("cmap (Green)");
		Error.FunctionFail();
		throw std::runtime_error("XAllocNamedColor(disp, cmap, \"Green\", &col, &exact)");
	}
	XSetForeground(disp, GCcol_dark[1], col.pixel);
	if (!(XAllocNamedColor(disp, cmap, "Blue", &col, &exact))) {
		Error.Function("XAllocNamedColor");
		Error.Value("cmap (Blue)");
		Error.FunctionFail();
		throw std::runtime_error("XAllocNamedColor(disp, cmap, \"Blue\", &col, &exact)");
	}
	XSetForeground(disp, GCcol_dark[2], col.pixel);
	// Pixmap
	pix = XCreatePixmap(disp, win, static_cast<unsigned int>(Window_size.width), static_cast<unsigned int>(Window_size.height), DefaultDepth(disp, 0));
	// Pass through XMapWindow() event
	XMaskEvent(disp, ExposureMask, &noev);
	// Initialize cos_a[] and sin_a[]
	for (int i = 0; i < ROTATE_ANGLE_MAX; i++) {
		cos_a[i] = cos(2.0 * M_PI * double(i) / double(ROTATE_ANGLE_MAX));
		sin_a[i] = sin(2.0 * M_PI * double(i) / double(ROTATE_ANGLE_MAX));
	}
	X11_Param->RotateSwitch = 0;
}


int
XEventor(X11_PARAM *X11_Param, SIZE Img_size)
{
	static COORDINATE Current_Mice_Pos;
	XEvent event;
	KeySym key;
	double dx, dy;

	while (XPending(disp) != 0) {
		XNextEvent(disp, &event);
		if (event.xclient.message_type == atom1
		    && static_cast<unsigned long>(event.xclient.data.l[0]) == atom2) {
			// The window close button is pressed
			return X_ESCAPE;
		}
		switch (event.type) {
			case ConfigureNotify: // Window State Changes (Resize, etc.)
				Window_size.width = event.xconfigure.width;
				Window_size.height = event.xconfigure.height;
				XFreePixmap(disp, pix);
				pix = XCreatePixmap(disp, win, static_cast<unsigned int>(Window_size.width), static_cast<unsigned int>(Window_size.height), DefaultDepth(disp, 0));
				break;
			case ButtonPress: // Mice button
				if (event.xbutton.button == Button1
				    || event.xbutton.button == Button2
				    || event.xbutton.button == Button3) {
					Current_Mice_Pos.x = event.xbutton.x;
					Current_Mice_Pos.y = event.xbutton.y;
				} else if (event.xbutton.button == Button4) {
					X11_Param->Scale *= 1.125;
				} else if (event.xbutton.button == Button5) {
					if (X11_Param->Scale > 1E-4) {
						X11_Param->Scale /= 1.125;
					}
				}
				break;
			case MotionNotify: // Mice motion
				if (event.xbutton.state & Button1Mask) {
					X11_Param->Longitude -= event.xbutton.x - Current_Mice_Pos.x;
					if (X11_Param->Longitude >= ROTATE_ANGLE_MAX) {
						X11_Param->Longitude = int(fmod(X11_Param->Longitude, ROTATE_ANGLE_MAX));
					} else if (X11_Param->Longitude < 0) {
						X11_Param->Longitude = int(ROTATE_ANGLE_MAX + fmod(X11_Param->Longitude, ROTATE_ANGLE_MAX));
					}
					X11_Param->Latitude += event.xbutton.y - Current_Mice_Pos.y;
					if (X11_Param->Latitude >= ROTATE_ANGLE_MAX) {
						X11_Param->Latitude = int(fmod(X11_Param->Latitude, ROTATE_ANGLE_MAX));
					} else if (X11_Param->Latitude < 0) {
						X11_Param->Latitude = int(ROTATE_ANGLE_MAX + fmod(X11_Param->Latitude, ROTATE_ANGLE_MAX));
					}
				} else if (event.xbutton.state & Button2Mask) {
					dx = event.xbutton.x - Current_Mice_Pos.x;
					dy = event.xbutton.y - Current_Mice_Pos.y;
					X11_Param->Center_x -= (dx * cos_a[X11_Param->Longitude] + dy * sin_a[X11_Param->Longitude] * cos_a[X11_Param->Latitude]) / X11_Param->Scale;
					X11_Param->Center_y -= (dy * cos_a[X11_Param->Longitude] * cos_a[X11_Param->Latitude] - dx * sin_a[X11_Param->Longitude]) / X11_Param->Scale;
					X11_Param->Center_z -= -dy * sin_a[X11_Param->Latitude] / X11_Param->Scale;
				} else if (event.xbutton.state & Button3Mask) {
					dy = event.xbutton.y - Current_Mice_Pos.y;
					if (dy > 0) {
						X11_Param->Plot_Z_Scale /= pow_int(1.01, int(dy));
					} else {
						X11_Param->Plot_Z_Scale *= pow_int(1.01, int(-dy));
					}
				}
				Current_Mice_Pos.x = event.xbutton.x;
				Current_Mice_Pos.y = event.xbutton.y;
				break;
			case KeyPress:
				if ((key = XLookupKeysym(&event.xkey, 0)) == XK_Escape
				    || key == XK_q) {
					return X_ESCAPE;
				}
				switch (key) {
					case XK_0: // Reset vision
						X11_Param->Longitude = X11_Param->Latitude = 0;
						X11_Param->Center_x = Img_size.width / 2;
						X11_Param->Center_y = Img_size.height / 2;
						X11_Param->Center_z = .0;
						X11_Param->Scale = Window_size.width * 0.8 / Img_size.width;
						X11_Param->Plot_Z_Scale = DEFAULT_PLOT_Z_SCALE;
						break;
					case XK_y: // Reset all
						X11_Param->Longitude = X11_Param->Latitude = 0;
						X11_Param->Center_x = Img_size.width / 2;
						X11_Param->Center_y = Img_size.height / 2;
						X11_Param->Center_z = .0;
						X11_Param->Scale = Window_size.width * 0.8 / Img_size.width;
						X11_Param->Plot_Z_Scale = DEFAULT_PLOT_Z_SCALE;
						X11_Param->ModeSwitch = 0;
						X11_Param->FillSwitch = 0;
						X11_Param->RotateSwitch = 0;
						break;
					case XK_c:
						X11_Param->ModeSwitch = X11_Param->ModeSwitch != X11_Plot_GravityCorrupt ? X11_Plot_GravityCorrupt : 0;
						break;
					case XK_f:
						X11_Param->FillSwitch = !X11_Param->FillSwitch;
						break;
					case XK_g:
						X11_Param->ModeSwitch = X11_Param->ModeSwitch != X11_Plot_Garaxy ? X11_Plot_Garaxy : 0;
						break;
					case XK_m:
						X11_Param->ModeSwitch = (X11_Param->ModeSwitch + 1) % NUMBER_OF_MODE_SELECTED_BY_SWITCH;
						break;
					case XK_r:
						X11_Param->RotateSwitch = (X11_Param->RotateSwitch + 1) % 5;
						break;
				}
		}
	}
	return X_NO_EVENT;
}


void
SwitchEventer(X11_PARAM *X11_Param)
{
	if (X11_Param->RotateSwitch != 0) {
		switch (X11_Param->RotateSwitch) {
			case 1: X11_Param->Longitude += 2; break;
			case 2: X11_Param->Longitude += 6; break;
			case 3: X11_Param->Longitude -= 2; break;
			case 4: X11_Param->Longitude -= 6; break;
		}
		if (X11_Param->Longitude >= ROTATE_ANGLE_MAX) {
			X11_Param->Longitude = int(fmod(X11_Param->Longitude, ROTATE_ANGLE_MAX));
		} else if (X11_Param->Longitude < 0) {
			X11_Param->Longitude = int(ROTATE_ANGLE_MAX + fmod(X11_Param->Longitude, ROTATE_ANGLE_MAX));
		}
	}
}


void
TransRotate_3DSegment(const X11_PARAM& X11_Param, const SEGMENT* segments, SEGMENT_X11* segments_plot, const unsigned int Num_Segments, const SIZE& Img_size, const SIZE& Img_size_resample)
{
	ERROR Error("TransRotate_3DSegment");
	double Scale_x = 1.0;
	double Scale_y = 1.0;

	if (Num_Segments < 1) {
		// Do NOT anything
		segments_plot = nullptr;
		return;
	} else if (segments_plot == nullptr) {
		Error.Value("segments_plot");
		Error.PointerNull();
		throw std::invalid_argument("error : TransRotate_3DSegment(const X11_PARAM&, const SEGMENT*, SEGMENT_X11*, const unsigned int, const SIZE&, const SIZE&) : SEGMENT_X11* segments_plot");
	}

	if (Img_size_resample.width > 0) {
		Scale_x = Img_size.width / Img_size_resample.width;
	}
	if (Img_size_resample.height > 0) {
		Scale_y = Img_size.height / Img_size_resample.height;
	}
	// Pixel Coordinate
	for (unsigned int n = 0; n < Num_Segments; n++) {
		// Start
		{
			double x = (Scale_x * segments[n].n - X11_Param.Center_x) * X11_Param.Scale;
			double y = (Scale_y * segments[n].m - X11_Param.Center_y) * X11_Param.Scale;
			double z = 0.0;
			segments_plot[n].start.x = short(Window_size.width / 2.0 + round(x * cos_a[X11_Param.Longitude] - y * sin_a[X11_Param.Longitude]));
			segments_plot[n].start.y = short(Window_size.height / 2.0 + round((y * cos_a[X11_Param.Longitude] + x * sin_a[X11_Param.Longitude]) * cos_a[X11_Param.Latitude] - z * sin_a[X11_Param.Latitude]));
		}
		// End
		{
			double x = (Scale_x * segments[n].x - X11_Param.Center_x) * X11_Param.Scale;
			double y = (Scale_y * segments[n].y - X11_Param.Center_y) * X11_Param.Scale;
			double z = 0.0;
			segments_plot[n].end.x = short(Window_size.width / 2.0 + round(x * cos_a[X11_Param.Longitude] - y * sin_a[X11_Param.Longitude]));
			segments_plot[n].end.y = short(Window_size.height / 2.0 + round((y * cos_a[X11_Param.Longitude] + x * sin_a[X11_Param.Longitude]) * cos_a[X11_Param.Latitude] - z * sin_a[X11_Param.Latitude]));
		}
	}
}


void
TransRotate_3DPoint(const X11_PARAM& X11_Param, const ImgVector<int>* Img, const int MaxInt, ImgVector<XPLOT>* Img_plot)
{
	ERROR Error("TransRotate_3DPoint");
	int m;

	if (Img == nullptr) {
		Error.Value("Img");
		Error.PointerNull();
		throw std::invalid_argument("error : void TransRotate_3DPoint(const X11_PARAM&, const ImgVector<int>*, const int, ImgVector<XPLOT>*) : ImgVector<int>* Img");
	} else if (Img_plot == nullptr) {
		Error.Value("Img_plot");
		Error.PointerNull();
		throw std::invalid_argument("error : void TransRotate_3DPoint(const X11_PARAM&, const ImgVector<int>*, const int, ImgVector<XPLOT>*) : ImgVector<XPLOT>* Img_plot");
	}
	// Pixel Coordinate
#ifdef _OPENMP
#pragma omp parallel for num_threads(8)
#endif
	for (m = 0; m < Img->height(); m++) {
		double y = (m - X11_Param.Center_y) * X11_Param.Scale;
		for (int n = 0; n < Img->width(); n++) {
			double x = (n - X11_Param.Center_x) * X11_Param.Scale;
			double z = ((-Img->get(n, m) + MaxInt / 2.0) - X11_Param.Center_z) * X11_Param.Plot_Z_Scale * X11_Param.Scale;
			Img_plot->at(n, m).point.x = short(Window_size.width / 2.0 + round(x * cos_a[X11_Param.Longitude] - y * sin_a[X11_Param.Longitude]));
			Img_plot->at(n, m).point.y = short(Window_size.height / 2.0 + round((y * cos_a[X11_Param.Longitude] + x * sin_a[X11_Param.Longitude]) * cos_a[X11_Param.Latitude] - z * sin_a[X11_Param.Latitude]));
			Img_plot->at(n, m).z = round(z * cos_a[X11_Param.Latitude] + (y * cos_a[X11_Param.Longitude] + x * sin_a[X11_Param.Longitude]) * sin_a[X11_Param.Latitude]);
		}
	}
}


void
TransGaraxy_3DPoint(const X11_PARAM& X11_Param, const ImgVector<int>* Img, ImgVector<COORDINATE_3D>* Img_coord, ImgVector<COORDINATE_3D>* Img_vel, const COORDINATE_3D& GaraxyCenter, ImgVector<XPLOT>* Img_plot)
{
	ERROR Error("TransGaraxy_3DPoint");
	const double Radius_Minimum = 0.01;
	const double dt = 0.5;

	if (Img == nullptr) {
		Error.Value("Img");
		Error.PointerNull();
		throw std::invalid_argument("error : void TransGaraxy_3DPoint(X11_PARAM&, ImgVector<int>*, ImgVector<COORDINATE_3D>*, ImgVector<COORDINATE_3D>*, COORDINATE_3D&, ImgVector<XPLOT>*) : ImgVector<int>* Img");
	} else if (Img_coord == nullptr) {
		Error.Value("Img_coord");
		Error.PointerNull();
		throw std::invalid_argument("error : void TransGaraxy_3DPoint(X11_PARAM&, ImgVector<int>*, ImgVector<COORDINATE_3D>*, ImgVector<COORDINATE_3D>*, COORDINATE_3D&, ImgVector<XPLOT>*) : ImgVector<COORDINATE_3D>* Img_coord");
	} else if (Img_vel == nullptr) {
		Error.Value("Img_vel");
		Error.PointerNull();
		throw std::invalid_argument("error : void TransGaraxy_3DPoint(X11_PARAM&, ImgVector<int>*, ImgVector<COORDINATE_3D>*, ImgVector<COORDINATE_3D>*, COORDINATE_3D&, ImgVector<XPLOT>*) : ImgVector<COORDINATE_3D>* Img_vel");
	} else if (Img_plot == nullptr) {
		Error.Value("Img_plot");
		Error.PointerNull();
		throw std::invalid_argument("error : void TransGaraxy_3DPoint(X11_PARAM&, ImgVector<int>*, ImgVector<COORDINATE_3D>*, ImgVector<COORDINATE_3D>*, COORDINATE_3D&, ImgVector<XPLOT>*) : ImgVector<COORDINATE_3D>* Img_plot");
	}
	// Gravity Motion
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (size_t i = 0; i < Img->size(); i++) {
		double r = sqrt(POW2(GaraxyCenter.x - Img_coord->at(i).x)
		    + POW2(GaraxyCenter.y - Img_coord->at(i).y)
		    + POW2(GaraxyCenter.z - Img_coord->at(i).z));
		if (r < Radius_Minimum) {
			r = Radius_Minimum;
		}
		Img_vel->at(i).x += dt * (GaraxyCenter.x - Img_coord->get(i).x) / pow_int(r, 3);
		Img_vel->at(i).y += dt * (GaraxyCenter.y - Img_coord->get(i).y) / pow_int(r, 3);
		Img_vel->at(i).z += dt * (GaraxyCenter.z - Img_coord->get(i).z) / pow_int(r, 3);
		Img_coord->at(i).x += Img_vel->get(i).x * dt;
		Img_coord->at(i).y += Img_vel->get(i).y * dt;
		Img_coord->at(i).z += Img_vel->get(i).z * dt;
	}
	// Pixel Coordinate
#ifdef _OPENMP
#pragma omp parallel for num_threads(8)
#endif
	for (size_t i = 0; i < Img->size(); i++) {
		double x = (Img_coord->get(i).x - X11_Param.Center_x) * X11_Param.Scale;
		double y = (Img_coord->get(i).y - X11_Param.Center_y) * X11_Param.Scale;
		double z = (-Img_coord->get(i).z - X11_Param.Center_z) * X11_Param.Scale;
		Img_plot->at(i).point.x = short(Window_size.width / 2.0 + round(x * cos_a[X11_Param.Longitude] - y * sin_a[X11_Param.Longitude]));
		Img_plot->at(i).point.y = short(Window_size.height / 2.0 + round((y * cos_a[X11_Param.Longitude] + x * sin_a[X11_Param.Longitude]) * cos_a[X11_Param.Latitude] - z * sin_a[X11_Param.Latitude]));
		Img_plot->at(i).z = round(z * cos_a[X11_Param.Latitude] + (y * cos_a[X11_Param.Longitude] + x * sin_a[X11_Param.Longitude]) * sin_a[X11_Param.Latitude]);
	}
}


void
TransGravity_3DPoint(const X11_PARAM& X11_Param, const ImgVector<int>* Img, ImgVector<COORDINATE_3D>* Img_coord, ImgVector<COORDINATE_3D>* Img_vel, ImgVector<XPLOT>* Img_plot)
{
	ERROR Error("TransGravity_3DPoint");
	const double Radius_Minimum = 0.01;
	const double dt = 0.5;
	size_t *core = nullptr;
	size_t Num_Cores = 0;
	int maxint = 0;

	if (Img == nullptr) {
		Error.Value("Img");
		Error.PointerNull();
		throw std::invalid_argument("error : void TransGravity_3DPoint(X11_PARAM&, ImgVector<int>*, ImgVector<COORDINATE_3D>*, ImgVector<COORDINATE_3D>*, ImgVector<XPLOT>*) : ImgVector<int>* Img");
	} else if (Img_coord == nullptr) {
		Error.Value("Img_coord");
		Error.PointerNull();
		throw std::invalid_argument("error : void TransGravity_3DPoint(X11_PARAM&, ImgVector<int>*, ImgVector<COORDINATE_3D>*, ImgVector<COORDINATE_3D>*, ImgVector<XPLOT>*) : ImgVector<COORDINATE_3D>* Img_coord");
	} else if (Img_vel == nullptr) {
		Error.Value("Img_vel");
		Error.PointerNull();
		throw std::invalid_argument("error : void TransGravity_3DPoint(X11_PARAM&, ImgVector<int>*, ImgVector<COORDINATE_3D>*, ImgVector<COORDINATE_3D>*, ImgVector<XPLOT>*) : ImgVector<COORDINATE_3D>* Img_vel");
	} else if (Img_plot == nullptr) {
		Error.Value("Img_plot");
		Error.PointerNull();
		throw std::invalid_argument("error : void TransGravity_3DPoint(X11_PARAM&, ImgVector<int>*, ImgVector<COORDINATE_3D>*, ImgVector<COORDINATE_3D>*, ImgVector<XPLOT>*) : ImgVector<XPLOT>* Img_plot");
	}
	try {
		core = new size_t[Img->size()];
	}
	catch (const std::bad_alloc& bad) {
		std::cerr << bad.what() << std::endl;
		Error.Value("core");
		Error.Malloc();
	}
	// List Cores
	for (size_t i = 0; i < Img->size(); i++) {
		if (maxint < Img->get(i)) {
			maxint = Img->get(i);
		}
	}
	for (size_t i = 0; i < Img->size(); i++) {
		if (Img->get(i) > maxint * 0.95) {
			core[Num_Cores] = i;
			Num_Cores++;
		}
	}
	// Gravity Motion
	{
		size_t i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (i = 0; i < Img->size(); i++) {
			for (size_t j = 0; j < Num_Cores; j++) {
				double M = double(Img->get(core[j])) / double(maxint);
				double r = sqrt(POW2(Img_coord->get(core[j]).x - Img_coord->get(i).x)
				    + POW2(Img_coord->get(core[j]).y - Img_coord->get(i).y)
				    + POW2(Img_coord->get(core[j]).z - Img_coord->get(i).z));
				if (r < Radius_Minimum) {
					r = Radius_Minimum;
				}
				Img_vel->at(i).x += dt * M * (Img_coord->get(core[j]).x - Img_coord->get(i).x) / pow_int(r, 3);
				Img_vel->at(i).y += dt * M * (Img_coord->get(core[j]).y - Img_coord->get(i).y) / pow_int(r, 3);
				Img_vel->at(i).z += dt * M * (Img_coord->get(core[j]).z - Img_coord->get(i).z) / pow_int(r, 3);
			}
			Img_coord->at(i).x += Img_vel->get(i).x * dt;
			Img_coord->at(i).y += Img_vel->get(i).y * dt;
			Img_coord->at(i).z += Img_vel->get(i).z * dt;
		}
	}
	delete[] core;
	// Pixel Coordinate
	{
		size_t i;
#ifdef _OPENMP
#pragma omp parallel for num_threads(8)
#endif
		for (i = 0; i < Img->size(); i++) {
			double x = (Img_coord->get(i).x - X11_Param.Center_x) * X11_Param.Scale;
			double y = (Img_coord->get(i).y - X11_Param.Center_y) * X11_Param.Scale;
			double z = (-Img_coord->get(i).z - X11_Param.Center_z) * X11_Param.Scale;
			Img_plot->at(i).point.x = short(Window_size.width / 2.0 + round(x * cos_a[X11_Param.Longitude] - y * sin_a[X11_Param.Longitude]));
			Img_plot->at(i).point.y = short(Window_size.height / 2.0 + round((y * cos_a[X11_Param.Longitude] + x * sin_a[X11_Param.Longitude]) * cos_a[X11_Param.Latitude] - z * sin_a[X11_Param.Latitude]));
			Img_plot->at(i).z = round(z * cos_a[X11_Param.Latitude] + (y * cos_a[X11_Param.Longitude] + x * sin_a[X11_Param.Longitude]) * sin_a[X11_Param.Latitude]);
		}
	}
}


void
Plot_3DPoints(const X11_PARAM& X11_Param, const ImgVector<int>* Img, ImgVector<XPLOT>* Img_plot, size_t* Img_index)
{
	ERROR Error("Plot_3DPoint");
	int Min_Intensity, Max_Intensity;
	SIZE rectsize;

	if (Img == nullptr) {
		Error.Value("Img");
		Error.PointerNull();
		throw std::invalid_argument("void Plot_3DPoints(const X11_PARAM&, const ImgVector<int>*, ImgVector<XPLOT>*, int*) : ImgVector<int>* Img");
	} else if (Img_plot == nullptr) {
		Error.Value("Img_plot");
		Error.PointerNull();
		throw std::invalid_argument("void Plot_3DPoints(const X11_PARAM&, const ImgVector<int>*, ImgVector<XPLOT>*, int*) : ImgVector<XPLOT>* Img_plot");
	} else if (Img_index == nullptr) {
		Error.Value("Img_index");
		Error.PointerNull();
		throw std::invalid_argument("void Plot_3DPoints(const X11_PARAM&, const ImgVector<int>*, ImgVector<XPLOT>*, int*) : int* Img_index");
	}

	rectsize.width = MAX(1, int(round(X11_Param.Scale * 0.25)));
	rectsize.height = MAX(1, int(floor(X11_Param.Scale * 0.25)));
	// Scan Intensity MIN and MAX
	Min_Intensity = Max_Intensity = Img->get(0);
	for (size_t i = 1; i < Img->size(); i++) {
		if (Img->get(i) < Min_Intensity) {
			Min_Intensity = Img->get(i);
		} else if (Img->get(i) > Max_Intensity) {
			Max_Intensity = Img->get(i);
		}
	}
	// Fill the window with Black
	XSetForeground(disp, GCmono, BlackPixel(disp, 0));
	XFillRectangle(disp, pix, GCmono, 0, 0, static_cast<unsigned int>(Window_size.width), static_cast<unsigned int>(Window_size.height));
	// Draw The Points
	for (size_t i = 0; i < Img->size(); i++) {
		size_t index = Img_index[i];
		if (0 <= Img_plot->get(index).point.x && Img_plot->get(index).point.x < Window_size.width
		    && 0 <= Img_plot->get(index).point.y && Img_plot->get(index).point.y < Window_size.height) {
			if (Img->get(index) > (Max_Intensity - Min_Intensity) * 0.75 + Min_Intensity) {
				if (rectsize.width == 1) {
					XDrawPoint(disp, pix, GCcol[0], Img_plot->get(index).point.x, Img_plot->get(index).point.y);
				} else {
					XFillRectangle(disp, pix, GCcol[0], Img_plot->get(index).point.x, Img_plot->get(index).point.y, static_cast<unsigned int>(rectsize.width), static_cast<unsigned int>(rectsize.height));
				}
			} else if (Img->get(index) > (Max_Intensity - Min_Intensity) * 0.25 + Min_Intensity) {
				if (rectsize.width == 1) {
					XDrawPoint(disp, pix, GCcol[1], Img_plot->get(index).point.x, Img_plot->get(index).point.y);
				} else {
					XFillRectangle(disp, pix, GCcol[1], Img_plot->get(index).point.x, Img_plot->get(index).point.y, static_cast<unsigned int>(rectsize.width), static_cast<unsigned int>(rectsize.height));
				}
			} else {                                                                         
				if (rectsize.width == 1) {
					XDrawPoint(disp, pix, GCcol[2], Img_plot->get(index).point.x, Img_plot->get(index).point.y);
				} else {
					XFillRectangle(disp, pix, GCcol[2], Img_plot->get(index).point.x, Img_plot->get(index).point.y, static_cast<unsigned int>(rectsize.width), static_cast<unsigned int>(rectsize.height));
				}
			}
		}
	}
}


void
Plot_3DGrid(const X11_PARAM& X11_Param, const ImgVector<int>* Img, ImgVector<XPLOT>* Img_plot, size_t* Img_index)
{
	ERROR Error("Plot_3DGridANDSegment");
	XPoint Triplet[8];
	int Min_Intensity, Max_Intensity;
	int local_Min, local_Max;
	int local_Min2, local_Max2;
	int count1, count2;

	if (Img == nullptr) {
		Error.Value("Img");
		Error.PointerNull();
		throw std::invalid_argument("void Plot_3DGrid(const X11_PARAM&, const ImgVector<int>*, ImgVector<XPLOT>*, int*) : ImgVector<int>* Img");
	} else if (Img_plot == nullptr) {
		Error.Value("Img_plot");
		Error.PointerNull();
		throw std::invalid_argument("void Plot_3DGrid(const X11_PARAM&, const ImgVector<int>*, ImgVector<XPLOT>*, int*) : ImgVector<XPLOT>* Img_plot");
	} else if (Img_index == nullptr) {
		Error.Value("Img_index");
		Error.PointerNull();
		throw std::invalid_argument("void Plot_3DGrid(const X11_PARAM&, const ImgVector<int>*, ImgVector<XPLOT>*, int*) : int* Img_index");
	}

	// Scan Intensity MIN and MAX
	Min_Intensity = Max_Intensity = Img->get(0);
	for (size_t n = 1; n < Img->size(); n++) {
		if (Img->get(n) < Min_Intensity) {
			Min_Intensity = Img->get(n);
		} else if (Img->get(n) > Max_Intensity) {
			Max_Intensity = Img->get(n);
		}
	}
	// Fill the window with Black
	XSetForeground(disp, GCmono, BlackPixel(disp, 0));
	XFillRectangle(disp, pix, GCmono, 0, 0, static_cast<unsigned int>(Window_size.width), static_cast<unsigned int>(Window_size.height));
	// Draw The Grid
	for (size_t n = 0; n < Img->size(); n++) {
		short x = short(Img_index[n] % size_t(Img->width()));
		short y = short(Img_index[n] / size_t(Img->width()));
		if (x == Img->width() - 1 || y == Img->height() - 1) {
			continue;
		}
		// Set Triplet
		if (Img_plot->get(x, y).z > Img_plot->get(x + 1, y).z) {
			Triplet[0].x = x;
			Triplet[0].y = y;
			Triplet[1].x = x;
			Triplet[1].y = y + 1;
			Triplet[2].x = x + 1;
			Triplet[2].y = y;
			Triplet[4].x = x + 1;
			Triplet[4].y = y;
			Triplet[5].x = x;
			Triplet[5].y = y + 1;
			Triplet[6].x = x + 1;
			Triplet[6].y = y + 1;
		} else {
			Triplet[0].x = x + 1;
			Triplet[0].y = y;
			Triplet[1].x = x;
			Triplet[1].y = y + 1;
			Triplet[2].x = x + 1;
			Triplet[2].y = y + 1;
			Triplet[4].x = x;
			Triplet[4].y = y;
			Triplet[5].x = x;
			Triplet[5].y = y + 1;
			Triplet[6].x = x + 1;
			Triplet[6].y = y;
		}
		Triplet[3] = Triplet[0];
		Triplet[7] = Triplet[4];
		local_Min = local_Max = Img->get(Triplet[0].x, Triplet[0].y);
		// Local Maximum
		if (Img->get(Triplet[1].x, Triplet[1].y) > local_Max) {
			local_Max = Img->get(Triplet[1].x, Triplet[1].y);
		}
		if (Img->get(Triplet[2].x, Triplet[2].y) > local_Max) {
			local_Max = Img->get(Triplet[2].x, Triplet[2].y);
		}
		// Local Minimum
		if (Img->get(Triplet[1].x, Triplet[1].y) < local_Min) {
			local_Min = Img->get(Triplet[1].x, Triplet[1].y);
		}
		if (Img->get(Triplet[2].x, Triplet[2].y) < local_Min) {
			local_Min = Img->get(Triplet[2].x, Triplet[2].y);
		}
		local_Min2 = local_Max2 = Img->get(Triplet[4].x, Triplet[4].y);
		// Local Maximum
		if (Img->get(Triplet[5].x, Triplet[5].y) > local_Max2) {
			local_Max2 = Img->get(Triplet[5].x, Triplet[5].y);
		}
		if (Img->get(Triplet[6].x, Triplet[6].y) > local_Max2) {
			local_Max2 = Img->get(Triplet[6].x, Triplet[6].y);
		}
		// Local Minimum
		if (Img->get(Triplet[5].x, Triplet[5].y) < local_Min2) {
			local_Min2 = Img->get(Triplet[5].x, Triplet[5].y);
		}
		if (Img->get(Triplet[6].x, Triplet[6].y) < local_Min2) {
			local_Min2 = Img->get(Triplet[6].x, Triplet[6].y);
		}
		// Convert Triplet coordinate to Real point coordinate
		count1 = count2 = 0;
		for (int i = 0; i < 4; i++) {
			Triplet[i] = Img_plot->get(Triplet[i].x, Triplet[i].y).point;
			if (Triplet[i].x < 0 || Window_size.width <= Triplet[i].x
			    || Triplet[i].y < 0 || Window_size.height <= Triplet[i].y) {
				count1--;
			}
			Triplet[4 + i] = Img_plot->get(Triplet[4 + i].x, Triplet[4 + i].y).point;
			if (Triplet[4 + i].x < 0 || Window_size.width <= Triplet[4 + i].x
			    || Triplet[4 + i].y < 0 || Window_size.height <= Triplet[4 + i].y) {
				count2--;
			}
		}
		// Triplet 1
		if (count1 > -4) {
			if (X11_Param.FillSwitch != 0) {
				XFillPolygon(disp, pix, GCmono, Triplet, 4, Convex, CoordModeOrigin);
			}
			if (local_Max > (Max_Intensity - Min_Intensity) * 0.75 + Min_Intensity) {
				XDrawLines(disp, pix, GCcol[0], Triplet, 4, CoordModeOrigin);
			} else if (local_Min < (Max_Intensity - Min_Intensity) * 0.25 + Min_Intensity) {
				XDrawLines(disp, pix, GCcol[2], Triplet, 4, CoordModeOrigin);
			} else {
				XDrawLines(disp, pix, GCcol[1], Triplet, 4, CoordModeOrigin);
			}
		}
		// Triplet 2
		if (count2 > -4) {
			if (X11_Param.FillSwitch != 0) {
				XFillPolygon(disp, pix, GCmono, Triplet + 4, 4, Convex, CoordModeOrigin);
			}
			if (local_Max2 > (Max_Intensity - Min_Intensity) * 0.75 + Min_Intensity) {
				XDrawLines(disp, pix, GCcol[0], Triplet + 4, 4, CoordModeOrigin);
			} else if (local_Min2 < (Max_Intensity - Min_Intensity) * 0.25 + Min_Intensity) {
				XDrawLines(disp, pix, GCcol[2], Triplet + 4, 4, CoordModeOrigin);
			} else {
				XDrawLines(disp, pix, GCcol[1], Triplet + 4, 4, CoordModeOrigin);
			}
		}
	}
}


void
Plot_3DSegment(const X11_PARAM& X11_Param, const SEGMENT_X11* segments_plot, const unsigned int Num_Segments)
{
	unsigned int num;

	if (Num_Segments < 1) {
		return;
	}
	// Draw The Segments
	// * Set Line width
	XSetLineAttributes(disp, GCcol[0], MAX(1u, static_cast<unsigned int>(round(X11_Param.Scale * 0.5))), LineSolid, CapNotLast, JoinMiter);
	for (num = 0; num < Num_Segments; num++) {
		XDrawLine(disp, pix, GCcol[0], segments_plot[num].start.x, segments_plot[num].start.y, segments_plot[num].end.x, segments_plot[num].end.y);
	}
	// * Reset Line width
	XSetLineAttributes(disp, GCcol[0], 0, LineSolid, CapNotLast, JoinMiter);
}


void
PlotParameters(const X11_PARAM& X11_Param)
{
#define X_STRING_LENGTH 256
	char Str[X_STRING_LENGTH];

	XSetForeground(disp, GCmono, WhitePixel(disp, 0));
	sprintf(Str, "Center (%.2f, %.2f), Latitude %d, Longitude %d, Scale %.2f, WinSize (%d, %d)", X11_Param.Center_x, X11_Param.Center_y, X11_Param.Latitude, X11_Param.Longitude, X11_Param.Scale, Window_size.width, Window_size.height);
	XDrawString(disp, pix, GCmono, 4, 12, Str, int(strlen(Str)));
	sprintf(Str, "[0] Reset Vision, [y] Reset all");
	XDrawString(disp, pix, GCmono, 4, 24, Str, int(strlen(Str)));
	// Switches mode
	sprintf(Str, "[m]ode :");
	XDrawString(disp, pix, GCmono, 4, Window_size.height - 27, Str, int(strlen(Str)));
	XDrawString(disp, pix, GCmono, 60, Window_size.height - 27, Plot_Mode[X11_Param.ModeSwitch], int(strlen(Plot_Mode[X11_Param.ModeSwitch])));
	// Switches
	sprintf(Str, "[r]otate");
	if (X11_Param.RotateSwitch != 0) {
		if (X11_Param.RotateSwitch <= 2) {
			XDrawString(disp, pix, GCcol[1], 4, Window_size.height - 8, Str, int(strlen(Str)));
		} else {
			XDrawString(disp, pix, GCcol[0], 4, Window_size.height - 8, Str, int(strlen(Str)));
		}
	} else {
		XDrawString(disp, pix, GCmono, 4, Window_size.height - 8, Str, int(strlen(Str)));
	}
	sprintf(Str, "[f]ill_grid");
	if (X11_Param.FillSwitch != 0) {
		XDrawString(disp, pix, GCcol[0], 64, Window_size.height - 8, Str, int(strlen(Str)));
	} else {
		XDrawString(disp, pix, GCmono, 64, Window_size.height - 8, Str, int(strlen(Str)));
	}
	sprintf(Str, "[g]araxy");
	if (X11_Param.ModeSwitch == X11_Plot_Garaxy) {
		XDrawString(disp, pix, GCcol[0], 142, Window_size.height - 8, Str, int(strlen(Str)));
	} else {
		XDrawString(disp, pix, GCmono, 142, Window_size.height - 8, Str, int(strlen(Str)));
	}
	sprintf(Str, "[c]orrupt");
	if (X11_Param.ModeSwitch == X11_Plot_GravityCorrupt) {
		XDrawString(disp, pix, GCcol[0], 202, Window_size.height - 8, Str, int(strlen(Str)));
	} else {
		XDrawString(disp, pix, GCmono, 202, Window_size.height - 8, Str, int(strlen(Str)));
	}
}


void
Set_Pixmap2Window(void)
{
	XEvent noev;
	// Copy pixmap to the window
	XCopyArea(disp, pix, win, GCmono, 0, 0, static_cast<unsigned int>(Window_size.width), static_cast<unsigned int>(Window_size.height), 0, 0);
	XMaskEvent(disp, ExposureMask, &noev);
}


void
reset_index(size_t* Img_index, const size_t N)
{
	if (Img_index == nullptr) {
		std::cerr << "error : void reset_index(int*, const int) : int* Img_index" << std::endl;
		throw std::invalid_argument("int* Img_index");
	}
	for (size_t n = 0; n < N; n++) {
		Img_index[n] = n;
	}
}


void
sort_index(const ImgVector<XPLOT>* Img_plot, size_t* Index, size_t* Index_tmp, const size_t N)
{
	if (Img_plot == nullptr) {
		std::cerr << "error : void sort_index(const ImgVector<XPLOT>*, int*, int*, const int) : const ImgVector<XPLOT>* Img_plot" << std::endl;
		throw std::invalid_argument("const ImgVector<XPLOT>* Img_plot");
	} else if (Index == nullptr) {
		std::cerr << "error : void sort_index(const ImgVector<XPLOT>*, int*, int*, const int) : int* index" << std::endl;
		throw std::invalid_argument("int* Index");
	} else if (Index_tmp == nullptr) {
		std::cerr << "error : void sort_index(const ImgVector<XPLOT>*, int*, int*, const int) : int* Index_tmp" << std::endl;
		throw std::invalid_argument("int* Index_tmp");
	}
	for (size_t n = 0; n < N; n++) {
		if (Index[n] >= N) {
			std::cerr << "error : void sort_index(const ImgVector<XPLOT>*, int*, int*, const int) : int* Index[n] out of range" << std::endl;
			throw std::out_of_range("int* Index[n] out of range");
		}
	}
	size_t step = 2;
	for (int div = 0; div < ceil(log(N) / log(2.0)); div++) {
		for (size_t n = 0; n < N; n += step) {
			size_t l = 0;
			size_t r = step / 2;
			for (size_t k = n; k < n + step && k < N; k++) {
				if (l < step / 2 && n + l < N
				    && ((r >= step || n + r >= N) || Img_plot->get(Index[n + l]).z > Img_plot->get(Index[n + r]).z)) {
					Index_tmp[k] = Index[n + l];
					l++;
				} else if (n + r < N) {
					Index_tmp[k] = Index[n + r];
					r++;
				}
			}
		}
		for (size_t n = 0; n < N; n++) {
			Index[n] = Index_tmp[n];
		}
		step *= 2;
	}
}

