#include "../Scratch_MeaningfulMotion.h"



int *
PlotSegment(const SEGMENT* coord_array, const unsigned int Num_Segments, const SIZE& size, const SIZE& size_out, const bool Negate)
{
	ERROR Error("PlotSegment");

	int *segments = nullptr;
	int Foreground = PLOT_INTENSITY_MAX;
	double scale_x, scale_y;

	try {
		segments = new int[size_out.width * size_out.height];
	}
	catch (const std::bad_alloc& bad) {
		std::cerr << bad.what() << std::endl;
		Error.Value("segments");
		Error.Malloc();
		return nullptr;
	}
	scale_x = double(size_out.width) / double(size.width);
	scale_y = double(size_out.height) / double(size.height);
	if (Negate) {
		for (int i = 0; i < size_out.width * size_out.height; i++) {
			segments[i] = PLOT_INTENSITY_MAX;
		}
		Foreground = 0;
	} else {
		for (int i = 0; i < size_out.width * size_out.height; i++) {
			segments[i] = 0;
		}
	}
	for (unsigned int i = 0; i < Num_Segments; i++) {
		int lx, ly, L;
		double dx, dy;
		int m, n, x, y;
		n = int(round(coord_array[i].n * scale_x));
		m = int(round(coord_array[i].m * scale_y));
		x = int(round(coord_array[i].x * scale_x));
		y = int(round(coord_array[i].y * scale_y));
		lx = abs(x - n);
		ly = abs(y - m);
		L = (lx > ly) ? lx : ly;
		dx = double(x - n) / L;
		dy = double(y - m) / L;
		for (int t = 0; t <= L; t++) {
			int tmpx, tmpy;
			tmpx = int(round(n + dx * t));
			tmpx = (tmpx >= 0) ? (tmpx < size_out.width) ? tmpx : size_out.width - 1 : 0;
			tmpy = int(round(m + dy * t));
			tmpy = (tmpy >= 0) ? (tmpy < size_out.height) ? tmpy : size_out.height - 1 : 0;
			segments[size_out.width * tmpy + tmpx] = Foreground;
		}
	}
	return segments;
}


void
Superimposer(PNM* pnm_out, const PNM& pnm_in, int* Plot, const SIZE& size, const int Color, const bool Negate)
{
	ERROR Error("Superimposer");
	pnm_img *img_out = nullptr;

	if (pnm_out == nullptr) {
		Error.Value("pnm_out");
		Error.PointerNull();
		throw std::invalid_argument("void Superimposer(PNM*, const PNM&, int*, const SIZE&, const int, bool) : PNM* pnm_out");
	} else if (pnm_in.isNULL() != false) {
		Error.Value("pnm_in");
		Error.ValueIncorrect();
		throw std::invalid_argument("void Superimposer(PNM*, const PNM&, int*, const SIZE&, const int, bool) : PNM& pnm_in");
	} else if (Plot == nullptr) {
		Error.Value("Plot");
		Error.PointerNull();
		throw std::invalid_argument("void Superimposer(PNM*, const PNM&, int*, const SIZE&, const int, bool) : int* Plot");
	} else if (size.width <= 0 || size.height <= 0) {
		Error.Value("size");
		Error.ValueIncorrect();
		throw std::invalid_argument("void Superimposer(PNM*, const PNM&, int*, const SIZE&, const int, bool) : SIZE& size");
	}

	if (pnm_in.isRGB() != false) {
		if (pnm_out->copy(pnm_in) == PNM_FUNCTION_ERROR) {
			Error.Function("pnm_out->copy");
			Error.Value("pnm_out <- pnm_in");
			Error.FunctionFail();
			throw std::runtime_error("void Superimposer(PNM*, const PNM&, int*, const SIZE&, const int, bool) : pnm_out->copy(PNM &)");
		}
	} else {
		if (pnm_out->Gray2RGB(pnm_in) == PNM_FUNCTION_ERROR) {
			Error.Function("pnm_out->Gray2RGB");
			Error.Value("pnm_out <- pnm_in");
			Error.FunctionFail();
			throw std::runtime_error("void Superimposer(PNM*, const PNM&, int*, const SIZE&, const int, bool) : pnm_out->Gray2RGB(PNM &)");
		}
	}

	if (pnm_in.MaxInt() > PLOT_INTENSITY_MAX) {
		for (int i = 0; i < size.width * size.height; i++) {
			if (Plot[i] > 0) {
				Plot[i] = int(round(Plot[i] * (double(pnm_in.MaxInt()) / PLOT_INTENSITY_MAX)));
			}
		}
	}
	img_out = pnm_out->Data();
	switch (Color) {
		default: // DEFAULT is Red
		case RED: // RED
			if (Negate) { // negative Superimpose
				for (int i = 0; i < size.height * size.width; i++) {
					img_out[i] = pnm_img(double(Plot[i]) / pnm_out->MaxInt());
				}
			} else {
				for (int i = 0; i < size.height * size.width; i++) {
					if (Plot[i] > 0) {
						img_out[i] += Plot[i];
						if (img_out[i] > int(pnm_out->MaxInt())) {
							img_out[i] = pnm_out->MaxInt();
						}
						img_out[size.height * size.width + i] /= 2;
						img_out[2 * size.height * size.width + i] /= 2;
					}
				}
			}
			break;
		case GREEN: // GREEN
			if (Negate) { // Negative Superimpose
				for (int i = 0; i < size.height * size.width; i++) {
					img_out[size.height * size.width + i] = pnm_img(double(Plot[i]) / pnm_out->MaxInt());
				}
			} else {
				for (int i = 0; i < size.height * size.width; i++) {
					if (Plot[i] > 0) {
						img_out[i] /= 2;
						img_out[size.height * size.width + i] += Plot[i];
						if (img_out[size.height * size.width + i] > int(pnm_out->MaxInt())) {
							img_out[size.height * size.width + i] = pnm_out->MaxInt();
						}
						img_out[2 * size.height * size.width + i] /= 2;
					}
				}
			}
			break;
		case BLUE: // BLUE
			if (Negate) { // Negative Superimpose
				for (int i = 0; i < size.height * size.width; i++) {
					img_out[2 * size.height * size.width + i] = pnm_img(double(Plot[i]) / pnm_out->MaxInt());
				}
			} else {
				for (int i = 0; i < size.height * size.width; i++) {
					if (Plot[i] > 0) {
						img_out[i] /= 2;
						img_out[size.height * size.width + i] /= 2;
						img_out[2 * size.height * size.width + i] += Plot[i];
						if (img_out[2 * size.height * size.width + i] > int(pnm_out->MaxInt())) {
							img_out[2 * size.height * size.width + i] = pnm_out->MaxInt();
						}
					}
				}
			}
	}
}

