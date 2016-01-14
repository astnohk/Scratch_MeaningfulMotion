#include "../ImgClass/MotionCompensation.h"
#include "HOG.h"

#include "../Scratch_MeaningfulMotion.h"




VECTOR_2D_W_SCORE *
HOG_Matching(const HOG* hog_prv, const HOG* hog_cur)
{
	ERROR Error("HOG_Matching");
	const double ep = 1.0E-6;
	VECTOR_2D_W_SCORE *vector = nullptr;
	SIZE search_region(65, 65); // recommends using odd number for the size of search region
	int W, H;
	int y;

	W = hog_cur->Width();
	H = hog_cur->Height();
	try {
		vector = new VECTOR_2D_W_SCORE[W * H];
	}
	catch (const std::bad_alloc &bad) {
		std::cerr << bad.what() << std::endl;
		Error.Value("vector");
		Error.Malloc();
		goto ExitError;
	}
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (y = 0; y < H; y++) {
		for (int x = 0; x < W; x++) {
			double d1 = 1.0E10;
			double d2 = 1.0E10;
			for (int yc = -search_region.height / 2; yc < search_region.height / 2; yc++) {
				if (y + yc < 0 || y + yc >= H) {
					continue;
				}
				for (int xc = -search_region.width / 2; xc < search_region.width / 2; xc++) {
					if (x + xc < 0 || x + xc >= W) {
						continue;
					}
					double d = HOG_Distance(hog_prv->Data(x, y), hog_cur->Data(x + xc, y + yc));
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
HOG_Distance(const Histogram* Hist1, const Histogram* Hist2)
{
	double sum = 0.0;

	for (int i = 0; i < Hist1->bins(); i++) {
		sum += POW2(Hist1->get(i) - Hist2->get(i));
	}
	return sqrt(sum);
}


bool
HOG_vector_write(const VECTOR_2D_W_SCORE* vector, const int width, const int height, const std::string& filename)
{
	ERROR Error("HOG_vector_write");
	FILE *fp = nullptr;

	fp = fopen(filename.c_str(), "wb");
	if (fp == nullptr) {
		Error.Function("fopen");
		Error.File(filename.c_str());
		Error.FileWrite();
		goto ExitError;
	}
	printf("* Output HOG matching vector to '%s'\n", filename.c_str());
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


void
HOG_vector_compensated_write(const ImgVector<double>& img_prev, const ImgVector<double>& img_next, const VECTOR_2D_W_SCORE* vector, const int width, const int height, const std::string& filename)
{
	ERROR Error("HOG_vector_write");
	ImgVector<VECTOR_2D<double> > vector2d(width, height);
	std::string filename_compensated;
	MotionCompensation<double> compensated;
	PNM pnm;

	for (size_t i = 0; i < vector2d.size(); i++) {
		vector2d[i].x = vector[i].x;
		vector2d[i].y = vector[i].y;
	}
	compensated.set(img_prev, img_next, vector2d); // initialize
	compensated.create_image_compensated(); // Make compensated image
	filename_compensated = filename.substr(0, filename.length() - 4) + "compensated" + filename.substr(filename.length() - 4);
	printf("* Output The compensated image by HOG matching vector to '%s'(binary)\n", filename_compensated.c_str());
	pnm.copy(PORTABLE_GRAYMAP_BINARY, compensated.width(), compensated.height(), 255, compensated.ref_image_compensated().data(), 1.0);
	pnm.write(filename_compensated.c_str());
	pnm.free();
}

