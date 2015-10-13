/* Multi Resolution functions */
double** Pyramider(double *img, SIZE size, int Level);
VECTOR_2D** grad_Pyramid(double **img_t, double **img_tp1, SIZE size, int Level);
double** dt_Pyramid(double **img_t, double **img_tp1, SIZE size, int Level);

