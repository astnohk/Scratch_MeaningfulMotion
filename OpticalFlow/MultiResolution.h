// Multi Resolution functions
ImgVector<double>* Pyramider(ImgVector<double> *img, int Level);
ImgVector<VECTOR_2D>* grad_Pyramid(ImgVector<double> *img_t, ImgVector<double> *img_tp1, int Level);
ImgVector<double>* dt_Pyramid(ImgVector<double> *img_t, ImgVector<double> *img_tp1, int Level);

