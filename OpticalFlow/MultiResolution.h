#include "../ImgClass/ImgClass.h"
#include "../lib/Vector.h"



// Multi Resolution functions
ImgVector<double>* Pyramider(ImgVector<double> *img, int Level);
ImgVector<VECTOR_2D>* grad_Pyramid(ImgVector<double> *img_t_levels, ImgVector<double> *img_tp1_levels, int Level);
ImgVector<double>* dt_Pyramid(ImgVector<double> *img_t_levels, ImgVector<double> *img_tp1_levels, int Level);

