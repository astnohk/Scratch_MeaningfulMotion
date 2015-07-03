#include "Scratch_MeaningfulMotion.h"
#include "HOG.h"




VECTOR_2D *
HOG_Matching(const HOG *hog_prv, const HOG *hog_cur)
{
	ERROR Error("HOG_Matching");
	int *match = nullptr;
	VECTOR_2D *vector = nullptr;
	double d1 = 1.0E5;
	double d2 = 1.0E5;
	int index;
	int W, H;
	double d;
	int x, y;
	int m, n;

	W = hog_cur->Width();
	H = hog_cur->Height();
	try {
		match = new int[W * H];
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("match");
		Error.Malloc();
		goto ExitError;
	}
#pragma omp parallel for private(x, m, n, d, d1, d2, index)
	for (y = 0; y < hog_prv->Height(); y++) {
		for (x = 0; x < hog_prv->Width(); x++) {
			for (m = 0; m < hog_cur->Height(); m++) {
				for (n = 0; n < hog_cur->Width(); n++) {
					d = HOG_Distance(hog_prv->Data(x, y), hog_cur->Data(n, m));
					if (d < d1) {
						d2 = d1;
						d1 = d;
						index = m * hog_cur->Width() + n;
					} else if (d < d2) {
						d2 = d;
					}
				}
			}
			match[W * y + x] = index;
		}
	}
	try {
		vector = new VECTOR_2D[W * H];
	}
	catch (const std::bad_alloc &bad) {
		Error.Value("vector");
		Error.Malloc();
		goto ExitError;
	}
	for (int i = 0; i < W * H; i++) {
		vector[i].x = match[i] % W;
		vector[i].y = match[i] / W;
	}
	delete[] match;
	return vector;
// Error
ExitError:
	delete[] match;
	delete[] vector;
	return nullptr;
}

double
HOG_Distance(const Histogram *Hist1, const Histogram *Hist2)
{
	double sum = 0.0;

	for (int i = 0; i < Hist1->Bins(); i++) {
		sum += POW2(Hist1->Hist(i) - Hist2->Hist(i));
	}
	return sqrt(sum);
}


bool
HOG_vector_write(const VECTOR_2D *vector, int width, int height, const char *filename)
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
	}
	fclose(fp);
	return true;
// Error
ExitError:
	return false;
}

