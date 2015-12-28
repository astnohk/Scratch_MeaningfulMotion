#include "../Scratch_MeaningfulMotion.h"
#include "HOG.h"


bool
HistogramsOfOrientedGradients(HOG* hog, HOG* block, const PNM_DOUBLE& Img, const HOG_PARAM& HOG_Param)
{
	ERROR Error("HistogramOfOrientedGradients");
	double *magnitude = nullptr;
	int *orient = nullptr;
	SIZE size;
	SIZE cellsize(7, 7);
	SIZE blocksize(3, 3);
	SIZE distance(4, 4);

	// Show Parameters
	printf("The number of Bins  : %d\n", HOG_Param.Bins);
	printf("Compute HOG densely : %s\n", HOG_Param.Dense ? "Yes" : "No");
	printf("Orientation signed  : %s\n\n", HOG_Param.SignedOrient ? "Yes" : "No");

	size.width = Img.Width();
	size.height = Img.Height();
	try {
		magnitude = new double[size.width * size.height];
		orient = new int[size.width * size.height];
	}
	catch (const std::bad_alloc &bad) {
		std::cerr << bad.what() << std::endl;
		Error.Value("magnitude, orient");
		Error.Malloc();
		goto ExitError;
	}
	printf("* * Compute Orientation\n");
	if (Orientation(magnitude, orient, Img, HOG_Param.Bins, HOG_Param.SignedOrient) == false) {
		Error.Function("Orientation");
		Error.Value("orient");
		Error.FunctionFail();
		goto ExitError;
	}
	printf("* * Compute Histogram\n");
	if (ComputeHistogramsOfOrientedGradients(hog, magnitude, orient, size, cellsize, HOG_Param.Bins, HOG_Param.SignedOrient, HOG_Param.Dense)
	    == false) {
		Error.Function("ComputeHistogramOfOrientedGradients");
		Error.Value("hog");
		Error.FunctionFail();
		goto ExitError;
	}
	delete[] magnitude;
	delete[] orient;
	printf("* * Normalize the block\n");
	if (HOG_BlockNormalize(block, hog, blocksize, distance) == false) {
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
Orientation(double* magnitude, int* orient, const PNM_DOUBLE& Img, const int bins, const bool sign)
{
	ERROR Error("Orientation");
	SIZE size;
	VECTOR_2D<double> *grad = nullptr;
	int x, y;

	size.width = Img.Width();
	size.height = Img.Height();
	try {
		grad = new VECTOR_2D<double>[size.width * size.height];
	}
	catch (const std::bad_alloc &bad) {
		std::cerr << bad.what() << std::endl;
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
		double tmp;
		double angle;

		magnitude[i] = sqrt(POW2(grad[i].x) + POW2(grad[i].y));
		tmp = atan2(grad[i].y, grad[i].x) / M_PI;
		if (sign == false) {
			if (tmp < .0) {
				angle = 1.0 + tmp;
			} else {
				angle = tmp;
			}
		} else {
			angle = (tmp + 1.0) / 2.0;
		}
		orient[i] = int(floor(bins * angle));
		if (orient[i] == bins) {
			orient[i] = 0;
		}
	}
	delete[] grad;
	grad = nullptr;
	return true;
// Error
ExitError:
	delete[] grad;
	return false;
}

bool
ComputeHistogramsOfOrientedGradients(HOG* hog, const double* magnitude, const int* orient, const SIZE& size, const SIZE& cell, const int bins, const bool sign, const bool denseHOG)
{
	ERROR Error("ComputeHistogramsOfOrientedGradients");
	SIZE Cells;

	if (denseHOG == false) {
		Cells.width = int(ceil(size.width / cell.width));
		Cells.height = int(ceil(size.height / cell.height));
	} else {
		Cells.width = size.width - (cell.width - 1);
		Cells.height = size.height - (cell.height - 1);
	}
	if (hog->reset(sign, Cells.width, Cells.height, bins) == false) {
		Error.Value("hog");
		Error.Malloc();
		goto ExitError;
	}
	int y;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (y = 0; y < Cells.height; y++) {
		for (int x = 0; x < Cells.width; x++) {
			int X, Y;
			if (denseHOG == false) {
				X = cell.width * x;
				Y = cell.height * y;
			} else {
				X = x;
				Y = y;
			}
			for (int m = 0; m < cell.height; m++) {
				for (int n = 0; n < cell.width; n++) {
					int dir = orient[size.width * (Y + m) + X + n];
					if (hog->AddHist(x, y, dir, magnitude[size.width * (Y + m) + X + n]) == false) {
						Error.OthersWarning("The bin is out of bounds");
						printf("orient = %.0f, bin = %d\n", double(orient[size.width * (Y + m) + X + n]) / bins * 360.0, dir);
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
HOG_BlockNormalize(HOG* block, const HOG* hog, const SIZE& blocksize)
{
	ERROR Error("HOG_BlockNormalize");
	SIZE size;
	double *integral_hist_norm = nullptr;
	const double ep = 1E-6;

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
		integral_hist_norm = new double[(size.width + 1) * (size.height + 1)];
	}
	catch (const std::bad_alloc &bad) {
		std::cerr << bad.what() << std::endl;
		Error.Value("integral_hist_norm");
		Error.Malloc();
		goto ExitError;
	}
	for (int y = 0; y < size.height; y++) {
		double norm = 0.0;
		for (int x = 0; x < size.width; x++) {
			for (int bin = 0; bin < hog->Bins(); bin++) {
				norm += POW2(hog->Hist(x, y, bin));
			}
			integral_hist_norm[size.width * (y + 1) + x + 1] = norm + integral_hist_norm[size.width * y + x + 1];
		}
	}
	int y;
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (y = 0; y < size.height; y++) {
		for (int x = 0; x < size.width; x++) {
			double norm = integral_hist_norm[size.width * (y + blocksize.height) + x + blocksize.width];
			norm -= integral_hist_norm[size.width * (y + blocksize.height) + x];
			norm -= integral_hist_norm[size.width * y + x + blocksize.width];
			norm += integral_hist_norm[size.width * y + x];
			double coeff = 1.0 / sqrt(norm + POW2(ep));
			for (int m = 0; m < blocksize.height; m++) {
				for (int n = 0; n < blocksize.width; n++) {
					for (int bin = 0; bin < hog->Bins(); bin++) {
						block->AddHist(x, y, (m * blocksize.width + n) * hog->Bins() + bin, hog->Hist(x + n, y + m, bin) * coeff);
					}
				}
			}
		}
	}
	delete[] integral_hist_norm;
	return true;
// Error
ExitError:
	delete[] integral_hist_norm;
	block->free();
	return false;
}

bool
HOG_BlockNormalize(HOG* block, const HOG* hog, const SIZE& blocksize, const SIZE& distance)
{
	// For dense trajectories
	ERROR Error("HOG_BlockNormalize");
	SIZE size;
	SIZE margin;
	const double ep = 1E-6;
	double norm;
	double coeff;
	int x, y;
	int m, n;
	int M, N;
	int bin;

	margin.width = (blocksize.width - 1) / 2 * distance.width;
	margin.height = (blocksize.height - 1) / 2 * distance.height;
	size.width = hog->Width() - 2 * margin.width;
	size.height = hog->Height() - 2 * margin.height;
	if (block->reset(hog->Signed(), size.width, size.height, blocksize.width * blocksize.height * hog->Bins())
	    == false) {
		Error.Function("block->reset");
		Error.Value("NULL");
		Error.FunctionFail();
		goto ExitError;
	}
#ifdef _OPENMP
#pragma omp parallel for private(x, norm, m, M, n, N, bin, coeff)
#endif
	for (y = 0; y < size.height; y++) {
		for (x = 0; x < size.width; x++) {
			norm = 0.0;
			for (m = 0; m < blocksize.height; m++) {
				M = m * distance.height;
				for (n = 0; n < blocksize.width; n++) {
					N = n * distance.width;
					for (bin = 0; bin < hog->Bins(); bin++) {
						norm += POW2(hog->Hist(x + N, y + M, bin));
					}
				}
			}
			coeff = 1.0 / sqrt(norm + POW2(ep));
			for (m = 0; m < blocksize.height; m++) {
				M = m * distance.height;
				for (n = 0; n < blocksize.width; n++) {
					N = n * distance.width;
					for (bin = 0; bin < hog->Bins(); bin++) {
						block->AddHist(x, y, (m * blocksize.width + n) * hog->Bins() + bin, hog->Hist(x + N, y + M, bin) * coeff);
					}
				}
			}
		}
	}
	return true;
// Error
ExitError:
	block->free();
	return false;
}


bool
HOG_write(const HOG& hog, const std::string& filename)
{
	ERROR Error("HOG_write");
	FILE *fp = nullptr;

	fp = fopen(filename.c_str(), "wb");
	if (fp == nullptr) {
		Error.Function("fopen");
		Error.File(filename.c_str());
		Error.FileWrite();
		goto ExitError;
	}
	printf("* Output HOG to '%s'\n", filename.c_str());
	fprintf(fp, "%d\n", int(hog.Signed()));
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
		}
	}
	fclose(fp);
	return true;
// Error
ExitError:
	return false;
}

