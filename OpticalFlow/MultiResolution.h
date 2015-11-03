#include "../ImgClass/ImgClass.h"
#include "../lib/ExtVector.h"



// Multi Resolution functions
ImgVector<double>* Pyramider(ImgVector<double> *img, int Level);
ImgVector<VECTOR_2D<double> >* grad_Pyramid(ImgVector<double> *img_t_levels, ImgVector<double> *img_tp1_levels, int Level);
ImgVector<double>* dt_Pyramid(ImgVector<double> *img_t_levels, ImgVector<double> *img_tp1_levels, int Level);

