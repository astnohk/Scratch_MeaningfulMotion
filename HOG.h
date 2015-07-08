// HOG
bool HistogramsOfOrientedGradients(HOG *hog, HOG *block, const PNM_DOUBLE &Img, HOG_PARAM &HOG_Param);
bool Orientation(double *magnitude, int *orient, const PNM_DOUBLE &Img, int bins, bool sign);
bool ComputeHistogramsOfOrientedGradients(HOG *hog, const double *magnitude, const int *orient, SIZE size, SIZE cell, int bins, bool sign, bool denseHOG);
bool HOG_BlockNormalize(HOG *block, const HOG *hog, SIZE blocksize);
bool HOG_BlockNormalize(HOG *block, const HOG *hog, SIZE blocksize, SIZE distance);

bool HOG_write(const HOG &hog, const char *filename);

// Matching
VECTOR_2D_W_SCORE* HOG_Matching(const HOG *hog_prv, const HOG *hog_cur);
double HOG_Distance(const Histogram *Hist1, const Histogram *Hist2);

bool HOG_vector_write(const VECTOR_2D_W_SCORE *vector, int width, int height, const char *filename);

