#include "Scratch_MeaningfulMotion.h"
#include "HOG.h"


bool
HistogramsOfOrientedGradients(HOG *hog, const PNM_DOUBLE &Img)
{
	ERROR Error("HistogramOfOrientedGradients");
	double *orient = nullptr;
	SIZE size;
	SIZE cell(3, 3);
	int bins = 16;

	size.width = Img.Width();
	size.height = Img.Height();
	orient = Orientation(Img, HOG_ORIENTATION_SIGNED);
	if (orient == nullptr) {
		Error.Function("Orientation");
		Error.Value("orient");
		Error.FunctionFail();
		goto ExitError;
	}
	ComputeHistogramsOfOrientedGradients(hog, orient, size, cell, bins, HOG_ORIENTATION_SIGNED);
	delete[] orient;
	return true;
// Error
ExitError:
	delete[] orient;
	return false;
}


double *
Orientation(const PNM_DOUBLE &Img, bool sign)
{
	ERROR Error("Orientation");
	SIZE size;
	double *orient = nullptr;
	VECTOR_2D *grad = nullptr;
	int x, y;

	size.width = Img.Width();
	size.height = Img.Height();
	try {
		grad = new VECTOR_2D[size.width * size.height];
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("grad");
		Error.Malloc();
		goto ExitError;
	}
	for (y = 0; y < size.height; y++) {
		for (x = 0; x < size.width; x++) {
			grad[size.width * y + x].x = Img.Image(x + 1, y) - Img.Image(x - 1, y);
			grad[size.width * y + x].y = Img.Image(x, y + 1) - Img.Image(x, y - 1);
		}
	}
	try {
		orient = new double[size.width * size.height];
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("orient");
		Error.Malloc();
		goto ExitError;
	}
	for (int i = 0; i < size.width * size.height; i++) {
		switch (sign) {
			case HOG_ORIENTATION_SIGNED:
				orient[i] = atan2(grad[i].y, grad[i].x) / M_PI;
				break;
			case HOG_ORIENTATION_UNSIGNED:
			default:
				orient[i] = atan2(grad[i].y, grad[i].x) / M_PI;
				if (orient[i] < .0) {
					orient[i] += 1.0;
				}
		}
	}
	delete[] grad;
	grad = nullptr;
	return orient;
// Error
ExitError:
	delete[] grad;
	delete[] orient;
	return nullptr;
}

bool
ComputeHistogramsOfOrientedGradients(HOG *hog, const double *orient, SIZE size, SIZE cell, int bins, bool sign)
{
	ERROR Error("ComputeHistogramsOfOrientedGradients");
	SIZE Cells;
	int x, y;
	int m, n;

	Cells.width = (int)ceil(size.width / cell.width);
	Cells.height = (int)ceil(size.height / cell.height);
	if (hog->reset(sign, Cells.width, Cells.height, bins) == false) {
		Error.Value("hog");
		Error.Malloc();
		goto ExitError;
	}
	for (y = 0; y < Cells.height; y++) {
		int Y = cell.height * y;
		for (x = 0; x < Cells.width; x++) {
			int X = cell.width * x;
			for (m = 0; m < cell.height; m++) {
				for (n = 0; n < cell.width; n++) {
					int dir;
					dir = (int)floor(bins * orient[size.width * (Y + m) + X + n]);
					if (dir == bins) {
						dir = dir - 1;
					}
					hog->AddHist(x, y, dir);
				}
			}
		}
	}
	return true;
// Error
ExitError:
	hog->free();
	return false;
}


bool
HOG_write(const HOG &hog, const char *filename)
{
	ERROR Error("HOG_write");
	FILE *fp = nullptr;

	fp = fopen(filename, "bw");
	if (fp == nullptr) {
		Error.Function("fopen");
		Error.File("fp");
		Error.FileRead();
		goto ExitError;
	}
	fprintf(fp, "%d\n", (int)hog.Signed());
	fprintf(fp, "%d %d\n", hog.Width(), hog.Height());
	fprintf(fp, "%d\n", hog.Bins());
	for (int y = 0; y < hog.Height(); y++) {
		for (int x = 0; x < hog.Height(); x++) {
			fprintf(fp, "%d", hog.Hist(x, y, 0));
			for (int bin = 1; bin < hog.Bins(); bin++) {
				fprintf(fp, " %d", hog.Hist(x, y, bin));
			}
			fprintf(fp, "\n");
		}
	}
	fclose(fp);
	return true;
// Error
ExitError:
	return false;
}
