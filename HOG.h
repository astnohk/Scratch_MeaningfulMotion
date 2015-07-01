#define HOG_ORIENTATION_UNSIGNED false
#define HOG_ORIENTATION_SIGNED true



HOG* HisotgramOfGradient(const PNM_DOUBLE &Img, SIZE size);
double* Orientation(const PNM_DOUBLE &Img, bool sign);
HOG* ComputeHistogramOfGradient(const double *orient, SIZE size, SIZE cell, int bins, bool sign);
bool HOG_write(const HOG *hog, SIZE size, const char *filename);

