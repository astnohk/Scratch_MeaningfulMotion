#include "Scratch_MeaningfulA.h"



int*
PlotSegment(SEGMENT *coord_array, int Num_Segments, SIZE size, SIZE size_out, int Negate)
{
	int *segments = NULL;
	int Foreground = PLOT_INTENSITY_MAX;
	double scale_x, scale_y;
	int lx, ly, L;
	double dx, dy;
	int m, n, x, y;
	int tmpx, tmpy;
	int i, t;

	if ((segments = (int *)calloc((size_t)size_out.height * size_out.width, sizeof(int))) == NULL) {
		fprintf(stderr, "calloc error on PlotSegment() (*segments)\n");
		return NULL;
	}
	scale_x = (double)size_out.width / size.width;
	scale_y = (double)size_out.height / size.height;
	if (Negate != 0) {
		for (i = 0; i < size_out.height * size_out.width; i++) {
			segments[i] = PLOT_INTENSITY_MAX;
		}
		Foreground = 0;
	}
	for (i = 0; i < Num_Segments; i++) {
		n = (int)round(coord_array[i].n * scale_x);
		m = (int)round(coord_array[i].m * scale_y);
		x = (int)round(coord_array[i].x * scale_x);
		y = (int)round(coord_array[i].y * scale_y);
		lx = abs(x - n);
		ly = abs(y - m);
		L = (lx > ly) ? lx : ly;
		dx = (double)(x - n) / L;
		dy = (double)(y - m) / L;
		for (t = 0; t <= L; t++) {
			tmpx = (int)round(n + dx * t);
			tmpx = (tmpx >= 0) ? (tmpx < size_out.width) ? tmpx : size_out.width - 1 : 0;
			tmpy = (int)round(m + dy * t);
			tmpy = (tmpy >= 0) ? (tmpy < size_out.height) ? tmpy : size_out.height - 1 : 0;
			segments[size_out.width * tmpy + tmpx] = Foreground;
		}
	}
	return segments;
}


int
Superimposer(PNM *pnm_out, PNM *pnm_in, int *Plot, SIZE size, int Color, int Negate)
{
	char *FunctionName = "Superimposer";
	char *ErrorFunctionName = "";
	char *ErrorValueName = "";

	int i;
	int tmp;

	if (pnm_out == NULL) {
		ErrorValueName = "pnm_out";
		goto ErrorPointerNull;
	} else if (pnm_in == NULL) {
		ErrorValueName = "pnm_in";
		goto ErrorPointerNull;
	} else if (pnm_isNULL(pnm_in) == PNM_TRUE) {
		ErrorValueName = "pnm_in";
		goto ErrorValueIncorrect;
	} else if (Plot == NULL) {
		ErrorValueName = "Plot";
		goto ErrorPointerNull;
	} else if (size.width <= 0 || size.height <= 0) {
		ErrorValueName = "size";
		goto ErrorValueIncorrect;
	}

	if (pnm_isNULL(pnm_out) == PNM_FALSE) {
		fprintf(stderr, "*** %s() warning - Output (struct PNM) pnm_out is NOT NULL or NOT initialized ***\n", FunctionName);
		pnmfree(pnm_out);
	}
	tmp = pnm_isRGB(pnm_in);
	if (tmp == PNM_FUNCTION_ERROR) {
		fprintf(stderr, "*** %s() error - Cannot check whether the input image is RGB or not. (Possibly (*pnm_in) is NULL) ***\n", FunctionName);
		goto ErrorReturn;
	} else if (tmp == PNM_TRUE) {
		if (pnmcp(pnm_out, pnm_in) != PNM_FUNCTION_SUCCESS) {
			ErrorFunctionName = "pnmcp";
			ErrorValueName = "(pnm_in -> pnm_out)";
			goto ErrorFunctionFailed;
		}
	} else {
		if (pnm_Gray2RGB(pnm_out, pnm_in) != PNM_FUNCTION_SUCCESS) {
			ErrorFunctionName = "pnmcp_Gray2RGB";
			ErrorValueName = "(pnm_in -> pnm_out)";
			goto ErrorFunctionFailed;
		}
	}
	pnm_out->desc = PORTABLE_PIXMAP_BINARY;

	if (pnm_in->maxint > PLOT_INTENSITY_MAX) {
		for (i = 0; i < size.width * size.height; i++) {
			if (Plot[i] > 0) {
				Plot[i] = (int)round(Plot[i] * ((double)pnm_in->maxint / PLOT_INTENSITY_MAX));
			}
		}
	}
	switch (Color) {
		default: /* DEFAULT is Red */
		case RED: /* RED */
			if (Negate == 0) {
				for (i = 0; i < size.height * size.width; i++) {
					if (Plot[i] > 0) {
						pnm_out->img[i] += Plot[i];
						if (pnm_out->img[i] > (int)pnm_out->maxint) {
							pnm_out->img[i] = pnm_out->maxint;
						}
						pnm_out->img[size.height * size.width + i] /= 2;
						pnm_out->img[2 * size.height * size.width + i] /= 2;
					}
				}
			} else { /* Negative Superimpose */
				for (i = 0; i < size.height * size.width; i++) {
					pnm_out->img[i] = (double)Plot[i] / pnm_out->maxint;
				}
			}
			break;
		case GREEN: /* GREEN */
			if (Negate == 0) {
				for (i = 0; i < size.height * size.width; i++) {
					if (Plot[i] > 0) {
						pnm_out->img[i] /= 2;
						pnm_out->img[size.height * size.width + i] += Plot[i];
						if (pnm_out->img[size.height * size.width + i] > (int)pnm_out->maxint) {
							pnm_out->img[size.height * size.width + i] = pnm_out->maxint;
						}
						pnm_out->img[2 * size.height * size.width + i] /= 2;
					}
				}
			} else { /* Negative Superimpose */
				for (i = 0; i < size.height * size.width; i++) {
					pnm_out->img[size.height * size.width + i] = (double)Plot[i] / pnm_out->maxint;
				}
			}
			break;
		case BLUE: /* BLUE */
			if (Negate == 0) {
				for (i = 0; i < size.height * size.width; i++) {
					if (Plot[i] > 0) {
						pnm_out->img[i] /= 2;
						pnm_out->img[size.height * size.width + i] /= 2;
						pnm_out->img[2 * size.height * size.width + i] += Plot[i];
						if (pnm_out->img[2 * size.height * size.width + i] > (int)pnm_out->maxint) {
							pnm_out->img[2 * size.height * size.width + i] = pnm_out->maxint;
						}
					}
				}
			} else { /* Negative Superimpose */
				for (i = 0; i < size.height * size.width; i++) {
					pnm_out->img[2 * size.height * size.width + i] = (double)Plot[i] / pnm_out->maxint;
				}
			}
	}
	return MEANINGFUL_SUCCESS;
/* Error */
ErrorFunctionFailed:
	fprintf(stderr, "*** %s() error - %s() failed to compute (%s) ***\n", FunctionName, ErrorFunctionName, ErrorValueName);
	goto ErrorReturn;
ErrorPointerNull:
	fprintf(stderr, "*** %s() error - The pointer (*%s) is NULL ***\n", FunctionName, ErrorValueName);
	goto ErrorReturn;
ErrorValueIncorrect:
	fprintf(stderr, "*** %s() error - The variable (%s) has incorrect value ***\n", FunctionName, ErrorValueName);
ErrorReturn:
	pnmfree(pnm_out);
	return MEANINGFUL_FAILURE;
}

