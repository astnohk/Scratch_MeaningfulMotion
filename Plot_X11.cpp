#include "Scratch_MeaningfulMotion.h"
#include "Plot_X11.h"

/* Coordinate System
 *
 *    Z o----> X
 *      |
 *      |
 *      V
 *      Y
 */


char *ProgramName = "Scratch_MeaningfulA";
SIZE Window_size = (SIZE){WINDOW_X_DEFAULT, WINDOW_Y_DEFAULT};
Display *disp = NULL;
Window win;
Pixmap pix;
GC GCmono;
GC GCcol[RGB_COLOR];
GC GCcol_dark[RGB_COLOR];
Colormap cmap;

double cos_a[ROTATE_ANGLE_MAX];
double sin_a[ROTATE_ANGLE_MAX];



int
ShowSegments_X11(int *Img, SIZE Img_size, SIZE Img_size_resample, int MaxInt, SEGMENT *segments, unsigned int Num_Segments)
{
	char *FunctionName = "ShowBounds_X11()";
	char *ErrorFunctionName = "";
	char *ErrorValueName = "";

	X11_PARAM X11_Param;
	COORDINATE_3D *Img_coord = NULL;
	COORDINATE_3D *Img_vel = NULL;
	XPLOT *Img_plot = NULL;
	SEGMENT_X11 *segments_plot = NULL;
	int *Img_index = NULL;
	int *Img_index_tmp = NULL;
	COORDINATE_3D GaraxyCenter;
	int loop = 0;
	int event_state = 0;
	int cur_mode = 0;
	double r;
	int m, n;
	double x, y, z;
	double X, Y;
	unsigned int k;

	if (segments == NULL) {
		ErrorValueName = "segments";
		goto ErrorPointerNull;
	} else if (MaxInt < 1) {
		ErrorValueName = "MaxInt";
		goto ErrorValueIncorrect;
	}

	/* Initialize X11 Window */
	if (Init_X11(&X11_Param, Img_size) == MEANINGFUL_FAILURE) {
		ErrorFunctionName = "Init_X11()";
		goto ErrorFunctionFailed;
	}

	/* Memory Allocation */
	if ((Img_plot = (XPLOT *)calloc((size_t)Img_size.width * Img_size.height, sizeof(XPLOT))) == NULL) {
		ErrorFunctionName = "calloc()";
		ErrorValueName = "Img_plot";
		goto ErrorMalloc;
	}
	if ((segments_plot = (SEGMENT_X11 *)calloc((size_t)Num_Segments, sizeof(SEGMENT_X11))) == NULL) {
		ErrorFunctionName = "calloc()";
		ErrorValueName = "segments_plot";
		goto ErrorMalloc;
	}
	if ((Img_coord = (COORDINATE_3D *)calloc((size_t)Img_size.width * Img_size.height, sizeof(COORDINATE_3D))) == NULL) {
		ErrorFunctionName = "calloc()";
		ErrorValueName = "Img_coord";
		goto ErrorMalloc;
	}
	if ((Img_vel = (COORDINATE_3D *)calloc((size_t)Img_size.width * Img_size.height, sizeof(COORDINATE_3D))) == NULL) {
		ErrorFunctionName = "calloc()";
		ErrorValueName = "Img_vel";
		goto ErrorMalloc;
	}
	if ((Img_index = (int *)calloc((size_t)Img_size.width * Img_size.height, sizeof(int))) == NULL) {
		ErrorFunctionName = "calloc()";
		ErrorValueName = "Img_index";
		goto ErrorMalloc;
	}
	if ((Img_index_tmp = (int *)calloc((size_t)Img_size.width * Img_size.height, sizeof(int))) == NULL) {
		ErrorFunctionName = "calloc()";
		ErrorValueName = "Img_index_tmp";
		goto ErrorMalloc;
	}

	GaraxyCenter = (COORDINATE_3D){X11_Param.Center_x, X11_Param.Center_y, X11_Param.Center_z};
	/* Infinite Loop */
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
			/* Set initial velocity and position */
			switch (cur_mode) {
				case X11_Plot_Garaxy:
					GaraxyCenter = (COORDINATE_3D){X11_Param.Center_x, X11_Param.Center_y, X11_Param.Center_z};
					for (m = 0; m < Img_size.height; m++) {
						y = m - GaraxyCenter.y;
						for (n = 0; n < Img_size.width; n++) {
							x = n - GaraxyCenter.x;
							z = (Img[Img_size.width * m + n] - MaxInt / 2.0) * X11_Param.Plot_Z_Scale * 2.0;
							Img_coord[Img_size.width * m + n].x = n;
							Img_coord[Img_size.width * m + n].y = m;
							Img_coord[Img_size.width * m + n].z = z;
							r = sqrt(POW2(x) + POW2(y) + POW2(z));
							if (r <= 1.0E-3) {
								r = 1.0E-3;
							}
							if (fabs(z) > 1.0E-6) {
								X = -(y / z + z * y);
								Y = z * x + x / z;
								Img_vel[Img_size.width * m + n].x = X / sqrt(POW2(X) + POW2(Y)) / sqrt(r);
								Img_vel[Img_size.width * m + n].y = Y / sqrt(POW2(X) + POW2(Y)) / sqrt(r);
							} else {
								Img_vel[Img_size.width * m + n].x = -y / r / sqrt(r);
								Img_vel[Img_size.width * m + n].y = -x / r / sqrt(r);
							}
							if (z < -1.0E-6) {
								Img_vel[Img_size.width * m + n].x *= -1.0;
								Img_vel[Img_size.width * m + n].y *= -1.0;
							}
							Img_vel[Img_size.width * m + n].z = 0.0;
						}
					}
					break;
				case X11_Plot_GravityCorrupt:
					for (m = 0; m < Img_size.height; m++) {
						for (n = 0; n < Img_size.width; n++) {
							Img_coord[Img_size.width * m + n].x = n;
							Img_coord[Img_size.width * m + n].y = m;
							Img_coord[Img_size.width * m + n].z = (Img[Img_size.width * m + n] - MaxInt / 2.0) * X11_Param.Plot_Z_Scale * 2.0;
							Img_vel[Img_size.width * m + n].x = 0.0;
							Img_vel[Img_size.width * m + n].y = 0.0;
							Img_vel[Img_size.width * m + n].z = 0.0;
						}
					}
			}
		}
		switch (X11_Param.ModeSwitch) {
			case X11_Plot_Point: /* Plot Image Intensity with Points */
				TransRotate_3DSegment(X11_Param, segments, segments_plot, Num_Segments, Img_size, Img_size_resample);
				TransRotate_3DPoint(X11_Param, Img, Img_size, MaxInt, Img_plot);
				if (reset_index(Img_index, Img_size.width * Img_size.height) == MEANINGFUL_FAILURE) {
					ErrorFunctionName = "reset_index";
					ErrorValueName = "Img_index";
					goto ErrorFunctionFailed;
				}
				sort_index(Img_plot, Img_index, Img_index_tmp, Img_size.width * Img_size.height);
				Plot_3DPoint(X11_Param, Img, Img_plot, Img_index, Img_size);
				break;
			case X11_Plot_Point_A_Segment: /* Plot Image Intensity with Points */
				TransRotate_3DSegment(X11_Param, segments, segments_plot, Num_Segments, Img_size, Img_size_resample);
				TransRotate_3DPoint(X11_Param, Img, Img_size, MaxInt, Img_plot);
				if (reset_index(Img_index, Img_size.width * Img_size.height) == MEANINGFUL_FAILURE) {
					ErrorFunctionName = "reset_index";
					ErrorValueName = "Img_index";
					goto ErrorFunctionFailed;
				}
				sort_index(Img_plot, Img_index, Img_index_tmp, Img_size.width * Img_size.height);
				Plot_3DPointANDSegment(X11_Param, Img, Img_plot, Img_index, Img_size, segments_plot, Num_Segments);
				break;
			case X11_Plot_Grid_A_Segment: /* Plot Image Intensity on Grid */
				TransRotate_3DSegment(X11_Param, segments, segments_plot, Num_Segments, Img_size, Img_size_resample);
				TransRotate_3DPoint(X11_Param, Img, Img_size, MaxInt, Img_plot);
				if (reset_index(Img_index, Img_size.width * Img_size.height) == MEANINGFUL_FAILURE) {
					ErrorFunctionName = "reset_index";
					ErrorValueName = "Img_index";
					goto ErrorFunctionFailed;
				}
				sort_index(Img_plot, Img_index, Img_index_tmp, Img_size.width * Img_size.height);
				Plot_3DGridANDSegment(X11_Param, Img, Img_plot, Img_index, Img_size, segments_plot, Num_Segments);
				break;
			case X11_Plot_Garaxy: /* Plot Image Gravity Motion (Garaxy) */
				TransGaraxy_3DPoint(X11_Param, Img, Img_size, Img_coord, Img_vel, GaraxyCenter, Img_plot);
				if (reset_index(Img_index, Img_size.width * Img_size.height) == MEANINGFUL_FAILURE) {
					ErrorFunctionName = "reset_index";
					ErrorValueName = "Img_index";
					goto ErrorFunctionFailed;
				}
				sort_index(Img_plot, Img_index, Img_index_tmp, Img_size.width * Img_size.height);
				Plot_3DPoint(X11_Param, Img, Img_plot, Img_index, Img_size);
				break;
			case X11_Plot_GravityCorrupt: /* Plot Image Gravity Motion */
				TransGravity_3DPoint(X11_Param, Img, Img_size, Img_coord, Img_vel, Img_plot);
				if (reset_index(Img_index, Img_size.width * Img_size.height) == MEANINGFUL_FAILURE) {
					ErrorFunctionName = "reset_index";
					ErrorValueName = "Img_index";
					goto ErrorFunctionFailed;
				}
				sort_index(Img_plot, Img_index, Img_index_tmp, Img_size.width * Img_size.height);
				Plot_3DPoint(X11_Param, Img, Img_plot, Img_index, Img_size);
		}
		usleep(WAIT_TIME);
	}

	/* Freeing Resources */
	XFreeGC(disp, GCmono);
	for (k = 0; k < RGB_COLOR; k++) {
		XFreeGC(disp, GCcol[k]);
		XFreeGC(disp, GCcol_dark[k]);
	}
	XFreePixmap(disp, pix);
	XFreeColormap(disp, cmap);
	XDestroyWindow(disp, win);
	XCloseDisplay(disp);
	free(Img_index_tmp);
	free(Img_index);
	free(Img_vel);
	free(Img_coord);
	free(segments_plot);
	free(Img_plot);
	return MEANINGFUL_SUCCESS;
/* Error */
ErrorMalloc:
	fprintf(stderr, "*** %s error - Cannot allocate memory for (*%s) by %s ***\n", FunctionName, ErrorValueName, ErrorFunctionName);
	goto ErrorReturn;
ErrorPointerNull:
	fprintf(stderr, "*** %s error - The pointer (*%s) is NULL ***\n", FunctionName, ErrorValueName);
	goto ErrorReturn;
ErrorValueIncorrect:
	fprintf(stderr, "*** %s error - The variable (%s) has incorrect value ***\n", FunctionName, ErrorValueName);
	goto ErrorReturn;
ErrorFunctionFailed:
	fprintf(stderr, "*** %s error - %s failed to compute (%s) ***\n", FunctionName, ErrorFunctionName, ErrorValueName);
ErrorReturn:
	if (disp != NULL) {
		XFreeGC(disp, GCmono);
		for (k = 0; k < RGB_COLOR; k++) {
			XFreeGC(disp, GCcol[k]);
			XFreeGC(disp, GCcol_dark[k]);
		}
		XFreePixmap(disp, pix);
		XFreeColormap(disp, cmap);
		XDestroyWindow(disp, win);
		XCloseDisplay(disp);
	}
	free(Img_index_tmp);
	free(Img_index);
	free(Img_vel);
	free(Img_coord);
	free(segments_plot);
	free(Img_plot);
	return MEANINGFUL_FAILURE;
}


int
Init_X11(X11_PARAM *X11_Param, SIZE Img_size)
{
	char *FunctionName = "Init_X11()";
	char *ErrorFunctionName = "";
	char *ErrorValueName = "";

	XColor col, exact;
	XEvent noev;
	int i;

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

	if ((disp = XOpenDisplay(NULL)) == NULL) {
		fprintf(stderr, "*** %s error - Cannot open display. Abort X11 Plotting ***\n", FunctionName);
		return MEANINGFUL_FAILURE;
	}
	win = XCreateSimpleWindow(disp, RootWindow(disp, 0), 0, 0, Window_size.width, Window_size.height, 0, BlackPixel(disp, 0), WhitePixel(disp, 0));
	XSelectInput(disp, win, ExposureMask | StructureNotifyMask | KeyPressMask | ButtonPressMask | ButtonMotionMask);
	XStoreName(disp, win, ProgramName);
	XMapWindow(disp, win);

	/* Graphic Context */
	GCmono = XCreateGC(disp, win, 0, 0);
	for (i = 0; i < RGB_COLOR; i++) {
		GCcol[i] = XCreateGC(disp, win, 0, 0);
		GCcol_dark[i] = XCreateGC(disp, win, 0, 0);
	}
	/* set GC monochrome */
	XSetForeground(disp, GCmono, BlackPixel(disp, 0));
	/* Colormap */
	cmap = DefaultColormap(disp, 0);
	/* set GC color */
	if (!(XAllocNamedColor(disp, cmap, "Red", &col, &exact))) {
		ErrorFunctionName = "XAllocNamedColor";
		ErrorValueName = "cmap (Red)";
		goto ErrorFunctionFailed;
	}
	XSetForeground(disp, GCcol[0], col.pixel);
	if (!(XAllocNamedColor(disp, cmap, "Green", &col, &exact))) {
		ErrorFunctionName = "XAllocNamedColor";
		ErrorValueName = "cmap (Green)";
		goto ErrorFunctionFailed;
	}
	XSetForeground(disp, GCcol[1], col.pixel);
	if (!(XAllocNamedColor(disp, cmap, "Blue", &col, &exact))) {
		ErrorFunctionName = "XAllocNamedColor";
		ErrorValueName = "cmap (Blue)";
		goto ErrorFunctionFailed;
	}
	XSetForeground(disp, GCcol[2], col.pixel);
	/* set GC color dark */
	if (!(XAllocNamedColor(disp, cmap, "Red", &col, &exact))) {
		ErrorFunctionName = "XAllocNamedColor";
		ErrorValueName = "cmap (Red)";
		goto ErrorFunctionFailed;
	}
	XSetForeground(disp, GCcol_dark[0], col.pixel);
	if (!(XAllocNamedColor(disp, cmap, "Green", &col, &exact))) {
		ErrorFunctionName = "XAllocNamedColor";
		ErrorValueName = "cmap (Green)";
		goto ErrorFunctionFailed;
	}
	XSetForeground(disp, GCcol_dark[1], col.pixel);
	if (!(XAllocNamedColor(disp, cmap, "Blue", &col, &exact))) {
		ErrorFunctionName = "XAllocNamedColor";
		ErrorValueName = "cmap (Blue)";
		goto ErrorFunctionFailed;
	}
	XSetForeground(disp, GCcol_dark[2], col.pixel);
	/* Pixmap */
	pix = XCreatePixmap(disp, win, Window_size.width, Window_size.height, DefaultDepth(disp, 0));
	/* Pass through XMapWindow() event */
	XMaskEvent(disp, ExposureMask, &noev);
	/* Initialize cos_a[] and sin_a[] */
	for (i = 0; i < ROTATE_ANGLE_MAX; i++) {
		cos_a[i] = cos(2.0 * M_PI * (double)i / ROTATE_ANGLE_MAX);
		sin_a[i] = sin(2.0 * M_PI * (double)i / ROTATE_ANGLE_MAX);
	}
	X11_Param->RotateSwitch = 0;

	return MEANINGFUL_SUCCESS;
/* Error */
ErrorFunctionFailed:
	fprintf(stderr, "*** %s error - %s() failed to compute (%s) ***\n", FunctionName, ErrorFunctionName, ErrorValueName);
	return MEANINGFUL_FAILURE;
}


int
XEventor(X11_PARAM *X11_Param, SIZE Img_size)
{
	static COORDINATE Current_Mice_Pos = COORDINATE_ZERO;
	XEvent event;
	int key;
	double dx, dy;

	while (XPending(disp) != 0) {
		XNextEvent(disp, &event);
		switch (event.type) {
			case ConfigureNotify: /* Window State Changes (Resize, etc.)*/
				Window_size.width = event.xconfigure.width;
				Window_size.height = event.xconfigure.height;
				XFreePixmap(disp, pix);
				pix = XCreatePixmap(disp, win, Window_size.width, Window_size.height, DefaultDepth(disp, 0));
				break;
			case ButtonPress: /* Mice button */
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
			case MotionNotify: /* Mice motion */
				if (event.xbutton.state & Button1Mask) {
					X11_Param->Longitude -= event.xbutton.x - Current_Mice_Pos.x;
					if (X11_Param->Longitude >= ROTATE_ANGLE_MAX) {
						X11_Param->Longitude = fmod(X11_Param->Longitude, ROTATE_ANGLE_MAX);
					} else if (X11_Param->Longitude < 0) {
						X11_Param->Longitude = ROTATE_ANGLE_MAX + fmod(X11_Param->Longitude, ROTATE_ANGLE_MAX);
					}
					X11_Param->Latitude += event.xbutton.y - Current_Mice_Pos.y;
					if (X11_Param->Latitude >= ROTATE_ANGLE_MAX) {
						X11_Param->Latitude = fmod(X11_Param->Latitude, ROTATE_ANGLE_MAX);
					} else if (X11_Param->Latitude < 0) {
						X11_Param->Latitude = ROTATE_ANGLE_MAX + fmod(X11_Param->Latitude, ROTATE_ANGLE_MAX);
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
						X11_Param->Plot_Z_Scale /= pow_int(1.01, dy);
					} else {
						X11_Param->Plot_Z_Scale *= pow_int(1.01, -dy);
					}
				}
				Current_Mice_Pos.x = event.xbutton.x;
				Current_Mice_Pos.y = event.xbutton.y;
				break;
			case KeyPress:
				if ((key = XLookupKeysym(&event.xkey, 0)) == XK_Escape) {
					return X_ESCAPE;
				}
				switch (key) {
					case XK_0:
						X11_Param->Longitude = X11_Param->Latitude = 0;
						X11_Param->Center_x = Img_size.width / 2;
						X11_Param->Center_y = Img_size.height / 2;
						X11_Param->Center_z = .0;
						X11_Param->Scale = Window_size.width * 0.8 / Img_size.width;
						X11_Param->Plot_Z_Scale = 1.0;
						break;
					case XK_c:
						X11_Param->ModeSwitch = X11_Plot_GravityCorrupt;
						break;
					case XK_f:
						X11_Param->FillSwitch = !X11_Param->FillSwitch;
						break;
					case XK_g:
						X11_Param->ModeSwitch = X11_Plot_Garaxy;
						break;
					case XK_m:
						X11_Param->ModeSwitch = (X11_Param->ModeSwitch + 1) % NUMBER_OF_MODE;
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
			X11_Param->Longitude = fmod(X11_Param->Longitude, ROTATE_ANGLE_MAX);
		} else if (X11_Param->Longitude < 0) {
			X11_Param->Longitude = ROTATE_ANGLE_MAX + fmod(X11_Param->Longitude, ROTATE_ANGLE_MAX);
		}
	}
}


int
TransRotate_3DSegment(X11_PARAM X11_Param, SEGMENT *segments, SEGMENT_X11 *segments_plot, unsigned int Num_Segments, SIZE Img_size, SIZE Img_size_resample)
{
	const char *FunctionName = "TransRotate_3DSegment()";
	char *ErrorValueName = "";
	unsigned int n;
	double x, y, z;
	double Scale_x = 1.0;
	double Scale_y = 1.0;

	if (segments == NULL) {
		ErrorValueName = "segments";
		goto ErrorPointerNull;
	} else if (segments_plot == NULL) {
		ErrorValueName = "segments_plot";
		goto ErrorPointerNull;
	}

	if (Img_size_resample.width > 0) {
		Scale_x = Img_size.width / Img_size_resample.width;
	}
	if (Img_size_resample.height > 0) {
		Scale_y = Img_size.height / Img_size_resample.height;
	}
	/* Pixel Coordinate */
	for (n = 0; n < Num_Segments; n++) {
		/* Start */
		x = (Scale_x * segments[n].n - X11_Param.Center_x) * X11_Param.Scale;
		y = (Scale_y * segments[n].m - X11_Param.Center_y) * X11_Param.Scale;
		z = 0.0;
		segments_plot[n].start.x = Window_size.width / 2.0 + round(x * cos_a[X11_Param.Longitude] - y * sin_a[X11_Param.Longitude]);
		segments_plot[n].start.y = Window_size.height / 2.0 + round((y * cos_a[X11_Param.Longitude] + x * sin_a[X11_Param.Longitude]) * cos_a[X11_Param.Latitude] - z * sin_a[X11_Param.Latitude]);
		/* End */
		x = (Scale_x * segments[n].x - X11_Param.Center_x) * X11_Param.Scale;
		y = (Scale_y * segments[n].y - X11_Param.Center_y) * X11_Param.Scale;
		z = 0.0;
		segments_plot[n].end.x = Window_size.width / 2.0 + round(x * cos_a[X11_Param.Longitude] - y * sin_a[X11_Param.Longitude]);
		segments_plot[n].end.y = Window_size.height / 2.0 + round((y * cos_a[X11_Param.Longitude] + x * sin_a[X11_Param.Longitude]) * cos_a[X11_Param.Latitude] - z * sin_a[X11_Param.Latitude]);
	}
	return MEANINGFUL_SUCCESS;
/* Error */
ErrorPointerNull:
	fprintf(stderr, "*** %s error - The pointer (*%s) is NULL ***\n", FunctionName, ErrorValueName);
	return MEANINGFUL_FAILURE;
}


int
TransRotate_3DPoint(X11_PARAM X11_Param, int *Img, SIZE size, int MaxInt, XPLOT *Img_plot)
{
	const char *FunctionName = "TransRotate_3DPoint()";
	char *ErrorValueName = "";
	int m, n;
	double x, y, z;

	if (Img == NULL) {
		ErrorValueName = "Img";
		goto ErrorPointerNull;
	} else if (Img_plot == NULL) {
		ErrorValueName = "Img_plot";
		goto ErrorPointerNull;
	}

	/* Pixel Coordinate */
#pragma omp parallel for private(n, x, y, z) num_threads(8)
	for (m = 0; m < size.height; m++) {
		y = (m - X11_Param.Center_y) * X11_Param.Scale;
		for (n = 0; n < size.width; n++) {
			x = (n - X11_Param.Center_x) * X11_Param.Scale;
			z = ((-Img[size.width * m + n] + MaxInt / 2.0) - X11_Param.Center_z) * X11_Param.Plot_Z_Scale * X11_Param.Scale;
			Img_plot[size.width * m + n].point.x = Window_size.width / 2.0 + round(x * cos_a[X11_Param.Longitude] - y * sin_a[X11_Param.Longitude]);
			Img_plot[size.width * m + n].point.y = Window_size.height / 2.0 + round((y * cos_a[X11_Param.Longitude] + x * sin_a[X11_Param.Longitude]) * cos_a[X11_Param.Latitude] - z * sin_a[X11_Param.Latitude]);
			Img_plot[size.width * m + n].z = round(z * cos_a[X11_Param.Latitude] + (y * cos_a[X11_Param.Longitude] + x * sin_a[X11_Param.Longitude]) * sin_a[X11_Param.Latitude]);
		}
	}
	return MEANINGFUL_SUCCESS;
/* Error */
ErrorPointerNull:
	fprintf(stderr, "*** %s error - The pointer (*%s) is NULL ***\n", FunctionName, ErrorValueName);
	return MEANINGFUL_FAILURE;
}


int
TransGaraxy_3DPoint(X11_PARAM X11_Param, int *Img, SIZE size, COORDINATE_3D *Img_coord, COORDINATE_3D *Img_vel, COORDINATE_3D GaraxyCenter, XPLOT *Img_plot)
{
	const char *FunctionName = "TransGaraxy_3DPoint()";
	char *ErrorValueName = "";

	const double Radius_Minimum = 0.01;
	const double dt = 0.5;
	double r;
	double x, y, z;
	int i;

	if (Img == NULL) {
		ErrorValueName = "Img";
		goto ErrorPointerNull;
	} else if (Img_coord == NULL) {
		ErrorValueName = "Img_coord";
		goto ErrorPointerNull;
	} else if (Img_vel == NULL) {
		ErrorValueName = "Img_vel";
		goto ErrorPointerNull;
	} else if (Img_plot == NULL) {
		ErrorValueName = "Img_plot";
		goto ErrorPointerNull;
	}

	/* Gravity Motion */
#pragma omp parallel for private(r)
	for (i = 0; i < size.width * size.height; i++) {
		r = sqrt(POW2(GaraxyCenter.x - Img_coord[i].x)
		    + POW2(GaraxyCenter.y - Img_coord[i].y)
		    + POW2(GaraxyCenter.z - Img_coord[i].z));
		if (r < Radius_Minimum) {
			r = Radius_Minimum;
		}
		Img_vel[i].x += dt * (GaraxyCenter.x - Img_coord[i].x) / pow_int(r, 3);
		Img_vel[i].y += dt * (GaraxyCenter.y - Img_coord[i].y) / pow_int(r, 3);
		Img_vel[i].z += dt * (GaraxyCenter.z - Img_coord[i].z) / pow_int(r, 3);
		Img_coord[i].x += Img_vel[i].x * dt;
		Img_coord[i].y += Img_vel[i].y * dt;
		Img_coord[i].z += Img_vel[i].z * dt;
	}
	/* Pixel Coordinate */
#pragma omp parallel for private(x, y, z) num_threads(8)
	for (i = 0; i < size.width * size.height; i++) {
		x = (Img_coord[i].x - X11_Param.Center_x) * X11_Param.Scale;
		y = (Img_coord[i].y - X11_Param.Center_y) * X11_Param.Scale;
		z = (-Img_coord[i].z - X11_Param.Center_z) * X11_Param.Scale;
		Img_plot[i].point.x = Window_size.width / 2.0 + round(x * cos_a[X11_Param.Longitude] - y * sin_a[X11_Param.Longitude]);
		Img_plot[i].point.y = Window_size.height / 2.0 + round((y * cos_a[X11_Param.Longitude] + x * sin_a[X11_Param.Longitude]) * cos_a[X11_Param.Latitude] - z * sin_a[X11_Param.Latitude]);
		Img_plot[i].z = round(z * cos_a[X11_Param.Latitude] + (y * cos_a[X11_Param.Longitude] + x * sin_a[X11_Param.Longitude]) * sin_a[X11_Param.Latitude]);
	}
	return MEANINGFUL_SUCCESS;
/* Error */
ErrorPointerNull:
	fprintf(stderr, "*** %s error - The pointer (*%s) is NULL ***\n", FunctionName, ErrorValueName);
	return MEANINGFUL_FAILURE;
}


int
TransGravity_3DPoint(X11_PARAM X11_Param, int *Img, SIZE size, COORDINATE_3D *Img_coord, COORDINATE_3D *Img_vel, XPLOT *Img_plot)
{
	const char *FunctionName = "TransGravity_3DPoint()";
	char *ErrorValueName = "";

	const double Radius_Minimum = 0.01;
	const double dt = 0.5;
	int *core = NULL;
	int Num_Cores = 0;
	double M, r;
	int maxint = 0;
	double x, y, z;
	int i, j;

	if (Img == NULL) {
		ErrorValueName = "Img";
		goto ErrorPointerNull;
	} else if (Img_coord == NULL) {
		ErrorValueName = "Img_coord";
		goto ErrorPointerNull;
	} else if (Img_vel == NULL) {
		ErrorValueName = "Img_vel";
		goto ErrorPointerNull;
	} else if (Img_plot == NULL) {
		ErrorValueName = "Img_plot";
		goto ErrorPointerNull;
	}
	if ((core = (int *)calloc((size_t)size.width * size.height, sizeof(int))) == NULL) {
		ErrorValueName = "core";
		goto ErrorMalloc;
	}

	/* List Cores */
	for (i = 0; i <size.width * size.height; i++) {
		if (maxint < Img[i]) {
			maxint = Img[i];
		}
	}
	for (i = 0; i <size.width * size.height; i++) {
		if (Img[i] > maxint * 0.95) {
			core[Num_Cores] = i;
			Num_Cores++;
		}
	}
	/* Gravity Motion */
#pragma omp parallel for private(j, M, r)
	for (i = 0; i < size.width * size.height; i++) {
		for (j = 0; j < Num_Cores; j++) {
			M = (double)Img[core[j]] / maxint;
			r = sqrt(POW2(Img_coord[core[j]].x - Img_coord[i].x)
			    + POW2(Img_coord[core[j]].y - Img_coord[i].y)
			    + POW2(Img_coord[core[j]].z - Img_coord[i].z));
			if (r < Radius_Minimum) {
				r = Radius_Minimum;
			}
			Img_vel[i].x += dt * M * (Img_coord[core[j]].x - Img_coord[i].x) / pow_int(r, 3);
			Img_vel[i].y += dt * M * (Img_coord[core[j]].y - Img_coord[i].y) / pow_int(r, 3);
			Img_vel[i].z += dt * M * (Img_coord[core[j]].z - Img_coord[i].z) / pow_int(r, 3);
		}
		Img_coord[i].x += Img_vel[i].x * dt;
		Img_coord[i].y += Img_vel[i].y * dt;
		Img_coord[i].z += Img_vel[i].z * dt;
	}
	/* Pixel Coordinate */
#pragma omp parallel for private(x, y, z) num_threads(8)
	for (i = 0; i < size.width * size.height; i++) {
		x = (Img_coord[i].x - X11_Param.Center_x) * X11_Param.Scale;
		y = (Img_coord[i].y - X11_Param.Center_y) * X11_Param.Scale;
		z = (-Img_coord[i].z - X11_Param.Center_z) * X11_Param.Scale;
		Img_plot[i].point.x = Window_size.width / 2.0 + round(x * cos_a[X11_Param.Longitude] - y * sin_a[X11_Param.Longitude]);
		Img_plot[i].point.y = Window_size.height / 2.0 + round((y * cos_a[X11_Param.Longitude] + x * sin_a[X11_Param.Longitude]) * cos_a[X11_Param.Latitude] - z * sin_a[X11_Param.Latitude]);
		Img_plot[i].z = round(z * cos_a[X11_Param.Latitude] + (y * cos_a[X11_Param.Longitude] + x * sin_a[X11_Param.Longitude]) * sin_a[X11_Param.Latitude]);
	}
	free(core);
	return MEANINGFUL_SUCCESS;
/* Error */
ErrorMalloc:
	fprintf(stderr, "*** %s error - Cannot allocate memory for (*%s) ***\n", FunctionName, ErrorValueName);
	goto ErrorReturn;
ErrorPointerNull:
	fprintf(stderr, "*** %s error - The pointer (*%s) is NULL ***\n", FunctionName, ErrorValueName);
ErrorReturn:
	free(core);
	return MEANINGFUL_FAILURE;
}


int
Plot_3DPoint(X11_PARAM X11_Param, int *Img, XPLOT *Img_plot, int *Img_index, SIZE size)
{
	const char *FunctionName = "Plot_3DPoint()";
	char *ErrorValueName = "";

	XEvent noev;
	int Min_Intensity, Max_Intensity;
	SIZE rectsize;
	int index;
	int i;

	if (Img == NULL) {
		ErrorValueName = "Img";
		goto ErrorPointerNull;
	} else if (Img_plot == NULL) {
		ErrorValueName = "Img_plot";
		goto ErrorPointerNull;
	}

	rectsize.width = MAX(1, round(X11_Param.Scale * 0.25));
	rectsize.height = MAX(1, floor(X11_Param.Scale * 0.25));
	/* Scan Intensity MIN and MAX */
	Min_Intensity = Max_Intensity = Img[0];
	for (i = 1; i < size.width * size.height; i++) {
		if (Img[i] < Min_Intensity) {
			Min_Intensity = Img[i];
		} else if (Img[i] > Max_Intensity) {
			Max_Intensity = Img[i];
		}
	}
	/* Fill the window with Black */
	XSetForeground(disp, GCmono, BlackPixel(disp, 0));
	XFillRectangle(disp, pix, GCmono, 0, 0, Window_size.width, Window_size.height);
	/* Draw The Points */
	for (i = 0; i < size.width * size.height; i++) {
		index = Img_index[i];
		if (0 <= Img_plot[index].point.x && Img_plot[index].point.x < Window_size.width
		    && 0 <= Img_plot[index].point.y && Img_plot[index].point.y < Window_size.height) {
			if (Img[index] > (Max_Intensity - Min_Intensity) * 0.75 + Min_Intensity) {
				if (rectsize.width == 1) {
					XDrawPoint(disp, pix, GCcol[0], Img_plot[index].point.x, Img_plot[index].point.y);
				} else {
					XFillRectangle(disp, pix, GCcol[0], Img_plot[index].point.x, Img_plot[index].point.y, rectsize.width, rectsize.height);
				}
			} else if (Img[index] > (Max_Intensity - Min_Intensity) * 0.25 + Min_Intensity) {
				if (rectsize.width == 1) {
					XDrawPoint(disp, pix, GCcol[1], Img_plot[index].point.x, Img_plot[index].point.y);
				} else {
					XFillRectangle(disp, pix, GCcol[1], Img_plot[index].point.x, Img_plot[index].point.y, rectsize.width, rectsize.height);
				}
			} else {                                                                         
				if (rectsize.width == 1) {
					XDrawPoint(disp, pix, GCcol[2], Img_plot[index].point.x, Img_plot[index].point.y);
				} else {
					XFillRectangle(disp, pix, GCcol[2], Img_plot[index].point.x, Img_plot[index].point.y, rectsize.width, rectsize.height);
				}
			}
		}
	}
	/* Show Parameters */
	PlotParameters(X11_Param);
	/* Copy pixmap to the window */
	XCopyArea(disp, pix, win, GCmono, 0, 0, Window_size.width, Window_size.height, 0, 0);
	XMaskEvent(disp, ExposureMask, &noev);
	return MEANINGFUL_SUCCESS;
/* Error */
ErrorPointerNull:
	fprintf(stderr, "*** %s error - The pointer (*%s) is NULL ***\n", FunctionName, ErrorValueName);
	return MEANINGFUL_FAILURE;
}


int
Plot_3DPointANDSegment(X11_PARAM X11_Param, int *Img, XPLOT *Img_plot, int *Img_index, SIZE size, SEGMENT_X11 *segments_plot, unsigned int Num_Segments)
{
	const char *FunctionName = "Plot_3DPointANDSegment()";
	char *ErrorValueName = "";

	XEvent noev;
	int Min_Intensity, Max_Intensity;
	SIZE rectsize;
	int index;
	int i;
	unsigned int n;

	if (Img == NULL) {
		ErrorValueName = "Img";
		goto ErrorPointerNull;
	} else if (Img_plot == NULL) {
		ErrorValueName = "Img_plot";
		goto ErrorPointerNull;
	} else if (segments_plot == NULL) {
		ErrorValueName = "segments_plot";
		goto ErrorPointerNull;
	}

	rectsize.width = MAX(1, round(X11_Param.Scale * 0.25));
	rectsize.height = MAX(1, floor(X11_Param.Scale * 0.25));
	/* Scan Intensity MIN and MAX */
	Min_Intensity = Max_Intensity = Img[0];
	for (i = 1; i < size.width * size.height; i++) {
		if (Img[i] < Min_Intensity) {
			Min_Intensity = Img[i];
		} else if (Img[i] > Max_Intensity) {
			Max_Intensity = Img[i];
		}
	}
	/* Fill the window with Black */
	XSetForeground(disp, GCmono, BlackPixel(disp, 0));
	XFillRectangle(disp, pix, GCmono, 0, 0, Window_size.width, Window_size.height);
	/* Draw The Points */
	for (i = 0; i < size.width * size.height; i++) {
		index = Img_index[i];
		if (0 <= Img_plot[index].point.x && Img_plot[index].point.x < Window_size.width
		    && 0 <= Img_plot[index].point.y && Img_plot[index].point.y < Window_size.height) {
			if (Img[index] > (Max_Intensity - Min_Intensity) * 0.75 + Min_Intensity) {
				if (rectsize.width == 1) {
					XDrawPoint(disp, pix, GCcol[0], Img_plot[index].point.x, Img_plot[index].point.y);
				} else {
					XFillRectangle(disp, pix, GCcol[0], Img_plot[index].point.x, Img_plot[index].point.y, rectsize.width, rectsize.height);
				}
			} else if (Img[index] > (Max_Intensity - Min_Intensity) * 0.25 + Min_Intensity) {
				if (rectsize.width == 1) {
					XDrawPoint(disp, pix, GCcol[1], Img_plot[index].point.x, Img_plot[index].point.y);
				} else {
					XFillRectangle(disp, pix, GCcol[1], Img_plot[index].point.x, Img_plot[index].point.y, rectsize.width, rectsize.height);
				}
			} else {                                                                         
				if (rectsize.width == 1) {
					XDrawPoint(disp, pix, GCcol[2], Img_plot[index].point.x, Img_plot[index].point.y);
				} else {
					XFillRectangle(disp, pix, GCcol[2], Img_plot[index].point.x, Img_plot[index].point.y, rectsize.width, rectsize.height);
				}
			}
		}
		if (X11_Param.Scale < 1.0) {
			i += round(1.0 / X11_Param.Scale);
			i = i % size.width;
			i--;
		}
	}
	/* Draw The Segments */
	/* * Set Line width */
	XSetLineAttributes(disp, GCcol[0], MAX(1, (int)round(X11_Param.Scale * 0.5)), LineSolid, CapNotLast, JoinMiter);
	for (n = 0; n < Num_Segments; n++) {
		XDrawLine(disp, pix, GCcol[0], segments_plot[n].start.x, segments_plot[n].start.y, segments_plot[n].end.x, segments_plot[n].end.y);
	}
	/* * Reset Line width */
	XSetLineAttributes(disp, GCcol[0], 0, LineSolid, CapNotLast, JoinMiter);
	/* Show Parameters */
	PlotParameters(X11_Param);
	/* Copy pixmap to the window */
	XCopyArea(disp, pix, win, GCmono, 0, 0, Window_size.width, Window_size.height, 0, 0);
	XMaskEvent(disp, ExposureMask, &noev);
	return MEANINGFUL_SUCCESS;
/* Error */
ErrorPointerNull:
	fprintf(stderr, "*** %s error - The pointer (*%s) is NULL ***\n", FunctionName, ErrorValueName);
	return MEANINGFUL_FAILURE;
}


int
Plot_3DGridANDSegment(X11_PARAM X11_Param, int *Img, XPLOT *Img_plot, int *Img_index, SIZE size, SEGMENT_X11 *segments_plot, unsigned int Num_Segments)
{
	const char *FunctionName = "Plot_3DGridANDSegment()";
	char *ErrorValueName = "";

	XPoint Triplet[8];
	XEvent noev;
	int Min_Intensity, Max_Intensity;
	int local_Min, local_Max;
	int local_Min2, local_Max2;
	int x, y;
	int n, i;
	unsigned int num;
	int count1, count2;

	if (Img == NULL) {
		ErrorValueName = "Img";
		goto ErrorPointerNull;
	} else if (Img_plot == NULL) {
		ErrorValueName = "Img_plot";
		goto ErrorPointerNull;
	} else if (Img_index == NULL) {
		ErrorValueName = "Img_index";
		goto ErrorPointerNull;
	} else if (segments_plot == NULL) {
		ErrorValueName = "segments_plot";
		goto ErrorPointerNull;
	}

	/* Scan Intensity MIN and MAX */
	Min_Intensity = Max_Intensity = Img[0];
	for (n = 1; n < size.width * size.height; n++) {
		if (Img[n] < Min_Intensity) {
			Min_Intensity = Img[n];
		} else if (Img[n] > Max_Intensity) {
			Max_Intensity = Img[n];
		}
	}
	/* Fill the window with Black */
	XSetForeground(disp, GCmono, BlackPixel(disp, 0));
	XFillRectangle(disp, pix, GCmono, 0, 0, Window_size.width, Window_size.height);
	/* Draw The Grid */
	for (n = 0; n < size.width * size.height; n++) {
		x = n % size.width;
		y = (int)floor(n / size.width);
		if (x == size.width - 1 || y == size.height - 1) {
			continue;
		}
		/* Set Triplet */
		if (Img_plot[size.width * y + x].z > Img_plot[size.width * y + x + 1].z) {
			Triplet[0] = (XPoint){x, y};
			Triplet[1] = (XPoint){x, y + 1};
			Triplet[2] = (XPoint){x + 1, y};
			Triplet[4] = (XPoint){x + 1, y};
			Triplet[5] = (XPoint){x, y + 1};
			Triplet[6] = (XPoint){x + 1, y + 1};
		} else {
			Triplet[0] = (XPoint){x + 1, y};
			Triplet[1] = (XPoint){x, y + 1};
			Triplet[2] = (XPoint){x + 1, y + 1};
			Triplet[4] = (XPoint){x, y};
			Triplet[5] = (XPoint){x, y + 1};
			Triplet[6] = (XPoint){x + 1, y};
		}
		Triplet[3] = Triplet[0];
		Triplet[7] = Triplet[4];
		local_Min = local_Max = Img[size.width * Triplet[0].y + Triplet[0].x];
		/* Local Maximum */
		if (Img[size.width * Triplet[1].y + Triplet[1].x] > local_Max) {
			local_Max = Img[size.width * Triplet[1].y + Triplet[1].x];
		}
		if (Img[size.width * Triplet[2].y + Triplet[2].x] > local_Max) {
			local_Max = Img[size.width * Triplet[2].y + Triplet[2].x];
		}
		/* Local Minimum */
		if (Img[size.width * Triplet[1].y + Triplet[1].x] < local_Min) {
			local_Min = Img[size.width * Triplet[1].y + Triplet[1].x];
		}
		if (Img[size.width * Triplet[2].y + Triplet[2].x] < local_Min) {
			local_Min = Img[size.width * Triplet[2].y + Triplet[2].x];
		}
		local_Min2 = local_Max2 = Img[size.width * Triplet[4].y + Triplet[4].x];
		/* Local Maximum */
		if (Img[size.width * Triplet[5].y + Triplet[5].x] > local_Max2) {
			local_Max2 = Img[size.width * Triplet[5].y + Triplet[5].x];
		}
		if (Img[size.width * Triplet[6].y + Triplet[6].x] > local_Max2) {
			local_Max2 = Img[size.width * Triplet[6].y + Triplet[6].x];
		}
		/* Local Minimum */
		if (Img[size.width * Triplet[5].y + Triplet[5].x] < local_Min2) {
			local_Min2 = Img[size.width * Triplet[5].y + Triplet[5].x];
		}
		if (Img[size.width * Triplet[6].y + Triplet[6].x] < local_Min2) {
			local_Min2 = Img[size.width * Triplet[6].y + Triplet[6].x];
		}
		/* Convert Triplet coordinate to Real point coordinate */
		count1 = count2 = 0;
		for (i = 0; i < 4; i++) {
			Triplet[i] = Img_plot[size.width * Triplet[i].y + Triplet[i].x].point;
			if (Triplet[i].x < 0 || Window_size.width <= Triplet[i].x
			    || Triplet[i].y < 0 || Window_size.height <= Triplet[i].y) {
				count1--;
			}
			Triplet[4 + i] = Img_plot[size.width * Triplet[4 + i].y + Triplet[4 + i].x].point;
			if (Triplet[4 + i].x < 0 || Window_size.width <= Triplet[4 + i].x
			    || Triplet[4 + i].y < 0 || Window_size.height <= Triplet[4 + i].y) {
				count2--;
			}
		}
		/* Triplet 1 */
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
		/* Triplet 2 */
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
	/* Draw The Segments */
	/* * Set Line width */
	XSetLineAttributes(disp, GCcol[0], MAX(1, (int)round(X11_Param.Scale * 0.5)), LineSolid, CapNotLast, JoinMiter);
	for (num = 0; num < Num_Segments; num++) {
		XDrawLine(disp, pix, GCcol[0], segments_plot[num].start.x, segments_plot[num].start.y, segments_plot[num].end.x, segments_plot[num].end.y);
	}
	/* * Reset Line width */
	XSetLineAttributes(disp, GCcol[0], 0, LineSolid, CapNotLast, JoinMiter);
	/* Show Parameters */
	PlotParameters(X11_Param);
	/* Copy pixmap to the window */
	XCopyArea(disp, pix, win, GCmono, 0, 0, Window_size.width, Window_size.height, 0, 0);
	XMaskEvent(disp, ExposureMask, &noev);
	return MEANINGFUL_SUCCESS;
/* Error */
ErrorPointerNull:
	fprintf(stderr, "*** %s error - The pointer (*%s) is NULL ***\n", FunctionName, ErrorValueName);
	return MEANINGFUL_FAILURE;
}


void
PlotParameters(X11_PARAM X11_Param)
{
#define X_STRING_LENGTH 256
	char Str[X_STRING_LENGTH];

	XSetForeground(disp, GCmono, WhitePixel(disp, 0));
	sprintf(Str, "Center (%.2f, %.2f), Latitude %d, Longitude %d, Scale %.2f, WinSize (%d, %d)", X11_Param.Center_x, X11_Param.Center_y, X11_Param.Latitude, X11_Param.Longitude, X11_Param.Scale, Window_size.width, Window_size.height);
	XDrawString(disp, pix, GCmono, 3, 12, Str, strlen(Str));
	//sprintf(Str, "[0]reset, [r]otate, [m]ode, [f]ill grid, [g]araxy, [c]orrupt");
	//XDrawString(disp, pix, GCmono, 3, Window_size.height - 12, Str, strlen(Str));
	sprintf(Str, "[r]otate");
	if (X11_Param.RotateSwitch != 0) {
		if (X11_Param.RotateSwitch <= 2) {
			XDrawString(disp, pix, GCcol[1], 3, Window_size.height - 12, Str, strlen(Str));
		} else {
			XDrawString(disp, pix, GCcol[0], 3, Window_size.height - 12, Str, strlen(Str));
		}
	} else {
		XDrawString(disp, pix, GCmono, 3, Window_size.height - 12, Str, strlen(Str));
	}
	sprintf(Str, "[f]ill_grid");
	if (X11_Param.FillSwitch != 0) {
		XDrawString(disp, pix, GCcol[0], 60, Window_size.height - 12, Str, strlen(Str));
	} else {
		XDrawString(disp, pix, GCmono, 60, Window_size.height - 12, Str, strlen(Str));
	}
	sprintf(Str, "[g]araxy");
	if (X11_Param.ModeSwitch == X11_Plot_Garaxy) {
		XDrawString(disp, pix, GCcol[0], 135, Window_size.height - 12, Str, strlen(Str));
	} else {
		XDrawString(disp, pix, GCmono, 135, Window_size.height - 12, Str, strlen(Str));
	}
	sprintf(Str, "[c]orrupt");
	if (X11_Param.ModeSwitch == X11_Plot_GravityCorrupt) {
		XDrawString(disp, pix, GCcol[0], 192, Window_size.height - 12, Str, strlen(Str));
	} else {
		XDrawString(disp, pix, GCmono, 192, Window_size.height - 12, Str, strlen(Str));
	}
}


int
reset_index(int *Img_index, int N)
{
	char *FunctionName = "reset_index()";
	char *ErrorValueName = "";
	int n;

	if (Img_index == NULL) {
		ErrorValueName = "Img_index";
		goto ErrorPointerNull;
	}
	for (n = 0; n < N; n++) {
		Img_index[n] = n;
	}
	return MEANINGFUL_SUCCESS;
/* Error */
ErrorPointerNull:
	fprintf(stderr, "*** %s error - The pointer (*%s) is NULL ***\n", FunctionName, ErrorValueName);
	return MEANINGFUL_FAILURE;
}


int
sort_index(XPLOT *Img_plot, int *Index, int *Index_tmp, int N)
{
	const char *FunctionName = "sort_index()";
	char *ErrorValueName = "";

	int div;
	int step;
	int n, k;
	int l, r;

	if (Img_plot == NULL) {
		ErrorValueName = "Img_plot";
		goto ErrorPointerNull;
	} else if (Index == NULL) {
		ErrorValueName = "Index";
		goto ErrorPointerNull;
	} else if (Index_tmp == NULL) {
		ErrorValueName = "Index_tmp";
		goto ErrorPointerNull;
	}
	for (n = 0; n < N; n++) {
		if (Index[n] < 0 || Index[n] >= N) {
			ErrorValueName = "Index[]";
			goto ErrorDataCorrupted;
		}
	}
	step = 2;
	for (div = 0; div < ceil(log(N) / log(2.0)); div++) {
		for (n = 0; n < N; n += step) {
			l = 0;
			r = step / 2;
			for (k = n; k < n + step && k < N; k++) {
				if (l < step / 2 && n + l < N
				    && ((r >= step || n + r >= N) || Img_plot[Index[n + l]].z > Img_plot[Index[n + r]].z)) {
					Index_tmp[k] = Index[n + l];
					l++;
				} else if (n + r < N) {
					Index_tmp[k] = Index[n + r];
					r++;
				}
			}
		}
		for (n = 0; n < N; n++) {
			Index[n] = Index_tmp[n];
		}
		step *= 2;
	}
	return MEANINGFUL_SUCCESS;
/* Error */
ErrorPointerNull:
	fprintf(stderr, "*** %s error - The pointer (*%s) is NULL ***\n", FunctionName, ErrorValueName);
	return MEANINGFUL_FAILURE;
ErrorDataCorrupted:
	fprintf(stderr, "*** %s error - The data (%s) is corrupted ***\n", FunctionName, ErrorValueName);
	return MEANINGFUL_FAILURE;
}

