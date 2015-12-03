#include "../ImgClass/ImgClass.h"
#include "../lib/ExtVector.h"



// Multi Resolution functions
ImgVector<double>* Pyramider(const ImgVector<double> *img, const int MaxLevel);

ImgVector<VECTOR_2D<double> >* grad_Pyramid(const ImgVector<double> *img_t_levels, const ImgVector<double> *img_tp1_levels, const int MaxLevel);
ImgVector<double>* dt_Pyramid(const ImgVector<double> *img_t_levels, const ImgVector<double> *img_tp1_levels, int MaxLevel);

