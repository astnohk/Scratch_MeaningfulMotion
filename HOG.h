#define HOG_ORIENTATION_UNSIGNED false
#define HOG_ORIENTATION_SIGNED true



bool HistogramsOfOrientedGradients(HOG *hog, const PNM_DOUBLE &Img);
double* Orientation(const PNM_DOUBLE &Img, bool sign);
bool ComputeHistogramsOfOrientedGradients(HOG *hog, const double *orient, SIZE size, SIZE cell, int bins, bool sign);
bool HOG_write(const HOG &hog, const char *filename);

