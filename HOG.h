#define HOG_ORIENTATION_UNSIGNED false
#define HOG_ORIENTATION_SIGNED true



bool HistogramsOfOrientedGradients(HOG *hog, HOG *block, const PNM_DOUBLE &Img);
bool Orientation(double *magnitude, double *orient, const PNM_DOUBLE &Img);
bool ComputeHistogramsOfOrientedGradients(HOG *hog, const double *magnitude, const double *orient, SIZE size, SIZE cell, int bins, bool sign);
bool HOG_BlockNormalize(HOG *block, const HOG *hog, SIZE blocksize);
bool HOG_write(const HOG &hog, const char *filename);

