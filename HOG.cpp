#include "Scratch_MeaningfulMotion.h"
#include "HOG.h"


bool
HistogramsOfOrientedGradients(HOG *hog, const PNM_DOUBLE &Img)
{
	ERROR Error("HistogramOfOrientedGradients");
	double *magnitude = nullptr;
	double *orient = nullptr;
	SIZE size;
	SIZE cell(5, 5);
	SIZE block(3, 3);
	int bins_unsigned = 16;

	size.width = Img.Width();
	size.height = Img.Height();
	try {
		magnitude = new double[size.width * size.height];
		orient = new double[size.width * size.height];
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("magnitude, orient");
		Error.Malloc();
		goto ExitError;
	}
	if (Orientation(magnitude, orient, Img, HOG_ORIENTATION_SIGNED) == false) {
		Error.Function("Orientation");
		Error.Value("orient");
		Error.FunctionFail();
		goto ExitError;
	}
	ComputeHistogramsOfOrientedGradients(hog, magnitude, orient, size, cell, bins_unsigned, HOG_ORIENTATION_SIGNED);
	delete[] magnitude;
	delete[] orient;
	return true;
// Error
ExitError:
	delete[] magnitude;
	delete[] orient;
	return false;
}


bool
Orientation(double *magnitude, double *orient, const PNM_DOUBLE &Img, bool sign)
{
	ERROR Error("Orientation");
	SIZE size;
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
	for (int i = 0; i < size.width * size.height; i++) {
		switch (sign) {
			case HOG_ORIENTATION_SIGNED:
				magnitude[i] = sqrt(POW2(grad[i].x) + POW2(grad[i].y));
				orient[i] = atan2(grad[i].y, grad[i].x) / M_PI + 1.0;
				break;
			case HOG_ORIENTATION_UNSIGNED:
			default:
				magnitude[i] = sqrt(POW2(grad[i].x) + POW2(grad[i].y));
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
	return nullptr;
}

bool
ComputeHistogramsOfOrientedGradients(HOG *hog, const double *magnitude, const double *orient, SIZE size, SIZE cell, int bins_unsigned, bool sign)
{
	ERROR Error("ComputeHistogramsOfOrientedGradients");
	SIZE Cells;
	int bins;
	int x, y;
	int m, n;

	if (sign == false) {
		bins = bins_unsigned;
	} else {
		bins = 2 * bins_unsigned;
	}
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
					dir = (int)floor(bins_unsigned * orient[size.width * (Y + m) + X + n]);
					if (dir == bins) {
						dir = 0;
					}
					if (hog->AddHist(x, y, dir, magnitude[size.width * (Y + m) + X + n]) == false) {
						Error.OthersWarning("The bin is out of bounds");
					}
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

	fp = fopen(filename, "wb");
	if (fp == nullptr) {
		Error.Function("fopen");
		Error.File(filename);
		Error.FileWrite();
		goto ExitError;
	}
	printf("* Output HOG to '%s'\n", filename);
	fprintf(fp, "%d\n", (int)hog.Signed());
	fprintf(fp, "%d %d\n", hog.Width(), hog.Height());
	fprintf(fp, "%d\n", hog.Bins());
	for (int y = 0; y < hog.Height(); y++) {
		for (int x = 0; x < hog.Height(); x++) {
			/*
			double tmp = hog.Hist(x, y, 0);
			if (fwrite(&tmp, sizeof(double), 1, fp) < 1) {
				Error.Function("fwrite");
				Error.Value("hist(x, y, bin)");
				Error.FunctionFail();
				goto ExitError;
			}
			*/
			fprintf(fp, "%f", hog.Hist(x, y, 0));
			for (int bin = 1; bin < hog.Bins(); bin++) {
				fprintf(fp, " %f", hog.Hist(x, y, bin));
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

