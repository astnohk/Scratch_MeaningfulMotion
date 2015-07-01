#define HOG_ORIENTATION_UNSIGNED false
#define HOG_ORIENTATION_SIGNED true



bool HistogramsOfOrientedGradients(HOG *hog, const PNM_DOUBLE &Img);
bool Orientation(double *magnitude, double *orient, const PNM_DOUBLE &Img, bool sign);
bool ComputeHistogramsOfOrientedGradients(HOG *hog, const double *magnitude, const double *orient, SIZE size, SIZE cell, int bins, bool sign);
bool HOG_write(const HOG &hog, const char *filename);

