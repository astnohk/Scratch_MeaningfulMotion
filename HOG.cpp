#include "Scratch_MeaningfulMotion.h"
#include "HOG.h"


bool
HistogramsOfOrientedGradients(HOG *hog, HOG *block, const PNM_DOUBLE &Img)
{
	ERROR Error("HistogramOfOrientedGradients");
	double *magnitude = nullptr;
	double *orient = nullptr;
	SIZE size;
	SIZE cellsize(5, 5);
	SIZE blocksize(3, 3);
	bool orient_signed = HOG_ORIENTATION_UNSIGNED;
	int bins = 9;

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
	if (Orientation(magnitude, orient, Img) == false) {
		Error.Function("Orientation");
		Error.Value("orient");
		Error.FunctionFail();
		goto ExitError;
	}
	if (ComputeHistogramsOfOrientedGradients(hog, magnitude, orient, size, cellsize, bins, orient_signed)
	    == false) {
		Error.Function("ComputeHistogramOfOrientedGradients");
		Error.Value("hog");
		Error.FunctionFail();
		goto ExitError;
	}
	delete[] magnitude;
	delete[] orient;
	if (HOG_BlockNormalize(block, hog, blocksize)
	    == false) {
		Error.Function("HOG_BlockNormalize");
		Error.Value("block");
		Error.FunctionFail();
		goto ExitError;
	}
	return true;
// Error
ExitError:
	delete[] magnitude;
	delete[] orient;
	return false;
}


bool
Orientation(double *magnitude, double *orient, const PNM_DOUBLE &Img)
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
		magnitude[i] = sqrt(POW2(grad[i].x) + POW2(grad[i].y));
		orient[i] = atan2(grad[i].y, grad[i].x) / M_PI;
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
ComputeHistogramsOfOrientedGradients(HOG *hog, const double *magnitude, const double *orient, SIZE size, SIZE cell, int bins, bool sign)
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
					double angle;

					if (sign == false) {
						if (orient[size.width * (Y + m) + X + n] < .0) {
							angle = 1.0 + orient[size.width * (Y + m) + X + n];
						} else {
							angle = orient[size.width * (Y + m) + X + n];
						}
					} else {
						angle = orient[size.width * (Y + m) + X + n] / 2.0 + 1.0;
					}
					dir = (int)floor(bins * angle);
					if (dir == bins) {
						dir = 0;
					}
					if (hog->AddHist(x, y, dir, magnitude[size.width * (Y + m) + X + n]) == false) {
						Error.OthersWarning("The bin is out of bounds");
						printf("orient = %f, angle = %f, bin = %d\n", orient[size.width * (Y + m) + X + n], angle, dir);
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
HOG_BlockNormalize(HOG *block, const HOG *hog, SIZE blocksize)
{
	ERROR Error("HOG_BlockNormalize");
	SIZE size;
	double *orig_hist = nullptr;
	const double ep = 1E-6;
	int x, y;

	size.width = hog->Width() - (blocksize.width - 1);
	size.height = hog->Height() - (blocksize.height - 1);
	if (block->reset(hog->Signed(), size.width, size.height, blocksize.width * blocksize.height * hog->Bins())
	    == false) {
		Error.Function("block->reset");
		Error.Value("NULL");
		Error.FunctionFail();
		goto ExitError;
	}
	try {
		orig_hist = new double[size.width * size.height * hog->Bins()];
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("orig_hist");
		Error.Malloc();
		goto ExitError;
	}
	for (y = 0; y < size.height; y++) {
		for (x = 0; x < size.width; x++) {
			int num = 0;
			double norm = .0;
			for (int m = 0; m < size.height; m++) {
				for (int n = 0; n < size.width; n++) {
					for (int bin = 0; bin < hog->Bins(); bin++) {
						orig_hist[num] = hog->Hist(x + n, y + m, bin);
						norm += POW2(orig_hist[num]);
						num++;
					}
				}
			}
			for (int bin = 0; bin < block->Bins(); bin++) {
				double coeff = 1.0 / sqrt(norm + POW2(ep));
				block->AddHist(x, y, bin, orig_hist[bin] * coeff);
			}
		}
	}
	delete[] orig_hist;
	return true;
// Error
ExitError:
	delete[] orig_hist;
	block->free();
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
		for (int x = 0; x < hog.Width(); x++) {
			double tmp;

			for (int bin = 0; bin < hog.Bins(); bin ++) {
				tmp = hog.Hist(x, y, bin);
				if (fwrite(&tmp, sizeof(double), 1, fp) < 1) {
					Error.Function("fwrite");
					Error.Value("hog.Hist(x, y, bin)");
					Error.FunctionFail();
					goto ExitError;
				}
			}
			//fprintf(fp, "%f", hog.Hist(x, y, 0));
			//for (int bin = 1; bin < hog.Bins(); bin++) {
			//	fprintf(fp, " %f", hog.Hist(x, y, bin));
			//}
			//fprintf(fp, "\n");
		}
	}
	fclose(fp);
	return true;
// Error
ExitError:
	return false;
}

