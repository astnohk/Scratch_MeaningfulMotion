#include <string>
#include "../ImgClass/ImgClass.h"
#include "../lib/ExtVector.h"
#include "../pnm_lib_cpp/pnm.h"
#include "HOG_struct.h"



// HOG
bool HistogramsOfOrientedGradients(HOG* hog, HOG* block, const PNM_DOUBLE& Img, const HOG_PARAM& HOG_Param);
bool Orientation(double* magnitude, int* orient, const PNM_DOUBLE& Img, const int bins, const bool sign);
bool ComputeHistogramsOfOrientedGradients(HOG* hog, const double* magnitude, const int* orient, const SIZE& size, const SIZE& cell, const int bins, const bool sign, const bool denseHOG);
bool HOG_BlockNormalize(HOG* block, const HOG* hog, const SIZE& blocksize);
bool HOG_BlockNormalize(HOG* block, const HOG* hog, const SIZE& blocksize, const SIZE& distance);

bool HOG_write(const HOG& hog, const std::string& filename);

// Matching
VECTOR_2D_W_SCORE* HOG_Matching(const HOG* hog_prv, const HOG* hog_cur);
double HOG_Distance(const Histogram* Hist1, const Histogram* Hist2);

bool HOG_vector_write(const VECTOR_2D_W_SCORE* vector, const int width, const int height, const std::string& filename);
void HOG_vector_compensated_write(const ImgVector<double>& img_prev, const ImgVector<double>& img_next, const VECTOR_2D_W_SCORE* vector, const int width, const int height, const std::string& filename);

