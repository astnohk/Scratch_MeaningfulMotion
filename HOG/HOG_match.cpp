#include "../Scratch_MeaningfulMotion.h"
#include "HOG.h"




VECTOR_2D_W_SCORE *
HOG_Matching(const HOG *hog_prv, const HOG *hog_cur)
{
	ERROR Error("HOG_Matching");
	const double ep = 1.0E-6;
	VECTOR_2D_W_SCORE *vector = nullptr;
	SIZE search_region(65, 65); // recommends using odd number for the size of search region
	double d1;
	double d2;
	int W, H;
	double d;
	int x, y;
	int xc, yc;

	W = hog_cur->Width();
	H = hog_cur->Height();
	try {
		vector = new VECTOR_2D_W_SCORE[W * H];
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("vector");
		Error.Malloc();
		goto ExitError;
	}
#pragma omp parallel for private(x, xc, yc, d, d1, d2)
	for (y = 0; y < H; y++) {
		for (x = 0; x < W; x++) {
			d1 = 1.0E10;
			d2 = 1.0E10;
			for (yc = -search_region.height / 2; yc < search_region.height / 2; yc++) {
				if (y + yc < 0 || y + yc >= H) {
					continue;
				}
				for (xc = -search_region.width / 2; xc < search_region.width / 2; xc++) {
					if (x + xc < 0 || x + xc >= W) {
						continue;
					}
					d = HOG_Distance(hog_prv->Data(x, y), hog_cur->Data(x + xc, y + yc));
					if (d < d1) {
						d2 = d1;
						d1 = d;
						vector[W * y + x].x = xc;
						vector[W * y + x].y = yc;
					} else if (d < d2) {
						d2 = d;
					}
				}
			}
			vector[W * y + x].score = (d2 - d1) / (d1 + ep);
		}
	}
	return vector;
// Error
ExitError:
	delete[] vector;
	return nullptr;
}

double
HOG_Distance(const Histogram *Hist1, const Histogram *Hist2)
{
	double sum = 0.0;

	for (int i = 0; i < Hist1->bins(); i++) {
		sum += POW2(Hist1->get(i) - Hist2->get(i));
	}
	return sqrt(sum);
}


bool
HOG_vector_write(const VECTOR_2D_W_SCORE *vector, int width, int height, const char *filename)
{
	ERROR Error("HOG_vector_write");
	FILE *fp = nullptr;

	fp = fopen(filename, "wb");
	if (fp == nullptr) {
		Error.Function("fopen");
		Error.File(filename);
		Error.FileWrite();
		goto ExitError;
	}
	printf("* Output HOG matching vector to '%s'\n", filename);
	fprintf(fp, "%d %d\n", width, height);
	for (int i = 0; i < width * height; i++) {
		double tmp;
		tmp = vector[i].x;
		if (fwrite(&tmp, sizeof(double), 1, fp) < 1) {
			Error.Function("fwrite");
			Error.Value("vector[].x");
			Error.FunctionFail();
			goto ExitError;
		}
		tmp = vector[i].y;
		if (fwrite(&tmp, sizeof(double), 1, fp) < 1) {
			Error.Function("fwrite");
			Error.Value("vector[].y");
			Error.FunctionFail();
			goto ExitError;
		}
		tmp = vector[i].score;
		if (fwrite(&tmp, sizeof(double), 1, fp) < 1) {
			Error.Function("fwrite");
			Error.Value("vector[].score");
			Error.FunctionFail();
			goto ExitError;
		}
	}
	fclose(fp);
	return true;
// Error
ExitError:
	return false;
}

