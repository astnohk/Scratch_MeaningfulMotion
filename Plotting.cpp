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
Superimposer(PNM *pnm_out, PNM &pnm_in, int *Plot, SIZE size, int Color, int Negate)
{
	ERROR Error("Superimposer");
	int i;
	pnm_img *img_out = nullptr;

	if (pnm_out == nullptr) {
		Error.Value("pnm_out");
		Error.PointerNull();
		goto ExitError;
	} else if (pnm_in.isNULL() != false) {
		Error.Value("pnm_in");
		Error.ValueIncorrect();
		goto ExitError;
	} else if (Plot == nullptr) {
		Error.Value("Plot");
		Error.PointerNull();
		goto ExitError;
	} else if (size.width <= 0 || size.height <= 0) {
		Error.Value("size");
		Error.ValueIncorrect();
		goto ExitError;
	}

	if (pnm_in.isRGB() != false) {
		if (pnm_out->copy(pnm_in) == PNM_FUNCTION_ERROR) {
			Error.Function("pnm_out->copy");
			Error.Value("pnm_out <- pnm_in");
			Error.FunctionFail();
			goto ExitError;
		}
	} else {
		if (pnm_out->Gray2RGB(pnm_in) == PNM_FUNCTION_ERROR) {
			Error.Function("pnm_out->Gray2RGB");
			Error.Value("pnm_out <- pnm_in");
			Error.FunctionFail();
			goto ExitError;
		}
	}

	if (pnm_in.MaxInt() > PLOT_INTENSITY_MAX) {
		for (i = 0; i < size.width * size.height; i++) {
			if (Plot[i] > 0) {
				Plot[i] = (int)round(Plot[i] * ((double)pnm_in.MaxInt() / PLOT_INTENSITY_MAX));
			}
		}
	}
	img_out = pnm_out->Data();
	switch (Color) {
		default: // DEFAULT is Red
		case RED: // RED
			if (Negate == 0) {
				for (i = 0; i < size.height * size.width; i++) {
					if (Plot[i] > 0) {
						img_out[i] += Plot[i];
						if (img_out[i] > (int)pnm_out->MaxInt()) {
							img_out[i] = pnm_out->MaxInt();
						}
						img_out[size.height * size.width + i] /= 2;
						img_out[2 * size.height * size.width + i] /= 2;
					}
				}
			} else { // Negative Superimpose
				for (i = 0; i < size.height * size.width; i++) {
					img_out[i] = (double)Plot[i] / pnm_out->MaxInt();
				}
			}
			break;
		case GREEN: // GREEN
			if (Negate == 0) {
				for (i = 0; i < size.height * size.width; i++) {
					if (Plot[i] > 0) {
						img_out[i] /= 2;
						img_out[size.height * size.width + i] += Plot[i];
						if (img_out[size.height * size.width + i] > (int)pnm_out->MaxInt()) {
							img_out[size.height * size.width + i] = pnm_out->MaxInt();
						}
						img_out[2 * size.height * size.width + i] /= 2;
					}
				}
			} else { // Negative Superimpose
				for (i = 0; i < size.height * size.width; i++) {
					img_out[size.height * size.width + i] = (double)Plot[i] / pnm_out->MaxInt();
				}
			}
			break;
		case BLUE: // BLUE
			if (Negate == 0) {
				for (i = 0; i < size.height * size.width; i++) {
					if (Plot[i] > 0) {
						img_out[i] /= 2;
						img_out[size.height * size.width + i] /= 2;
						img_out[2 * size.height * size.width + i] += Plot[i];
						if (img_out[2 * size.height * size.width + i] > (int)pnm_out->MaxInt()) {
							img_out[2 * size.height * size.width + i] = pnm_out->MaxInt();
						}
					}
				}
			} else { // Negative Superimpose
				for (i = 0; i < size.height * size.width; i++) {
					img_out[2 * size.height * size.width + i] = (double)Plot[i] / pnm_out->MaxInt();
				}
			}
	}
	img_out = nullptr;
	return MEANINGFUL_SUCCESS;
// Error
ExitError:
	pnm_out->free();
	return MEANINGFUL_FAILURE;
}

