#include "Scratch_MeaningfulMotion.h"



int*
PlotSegment(SEGMENT *coord_array, int Num_Segments, SIZE size, SIZE size_out, int Negate)
{
	ERROR Error("PlotSegment");

	int *segments = NULL;
	int Foreground = PLOT_INTENSITY_MAX;
	double scale_x, scale_y;
	int lx, ly, L;
	double dx, dy;
	int m, n, x, y;
	int tmpx, tmpy;
	int i, t;

	if ((segments = (int *)calloc((size_t)size_out.height * size_out.width, sizeof(int))) == NULL) {
		Error.Others("calloc error on PlotSegment() (*segments)");
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
	ERROR Error("Superimposer");
	int i;
	int tmp;

	if (pnm_out == NULL) {
		Error.Value("pnm_out");
		Error.PointerNull();
		goto ExitError;
	} else if (pnm_in == NULL) {
		Error.Value("pnm_in");
		Error.PointerNull();
		goto ExitError;
	} else if (pnm_isNULL(pnm_in) == PNM_TRUE) {
		Error.Value("pnm_in");
		Error.ValueIncorrect();
		goto ExitError;
	} else if (Plot == NULL) {
		Error.Value("Plot");
		Error.PointerNull();
		goto ExitError;
	} else if (size.width <= 0 || size.height <= 0) {
		Error.Value("size");
		Error.ValueIncorrect();
		goto ExitError;
	}

	if (pnm_isNULL(pnm_out) == PNM_FALSE) {
		Error.OthersWarning("Output (struct PNM) pnm_out is NOT NULL or NOT initialized");
		pnmfree(pnm_out);
	}
	tmp = pnm_isRGB(pnm_in);
	if (tmp == PNM_FUNCTION_ERROR) {
		Error.Others("Cannot check whether the input image is RGB or not. (Possibly (*pnm_in) is NULL)");
		goto ExitError;
	} else if (tmp == PNM_TRUE) {
		if (pnmcp(pnm_out, pnm_in) != PNM_FUNCTION_SUCCESS) {
			Error.Function("pnmcp");
			Error.Value("(pnm_in -> pnm_out)");
			Error.FunctionFail();
			goto ExitError;
		}
	} else {
		if (pnm_Gray2RGB(pnm_out, pnm_in) != PNM_FUNCTION_SUCCESS) {
			Error.Function("pnmcp_Gray2RGB");
			Error.Value("(pnm_in -> pnm_out)");
			Error.FunctionFail();
			goto ExitError;
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
ExitError:
	pnmfree(pnm_out);
	return MEANINGFUL_FAILURE;
}

