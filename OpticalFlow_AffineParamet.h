/* Robust MultiResolution Estimation of Parametric Motion
 *
 * This motion estimation based on
 * J.M.Odobez and P.Bouthemy, "Robust Multiresolution Estimation of Parametric Motion Models," Visual Communication And Image Representation, Vol.6, No.4, pp.348-365, 1995.
 */




TUPLE_VEC_SCALAR* OpticalFlow_AffineParamet(double *I_t, double *I_tp1, SIZE size, OPTICALFLOW_PARAM OpticalFlow_Param);
TUPLE_VEC_SCALAR IRLS_OpticalFlow_Affine(double *I_t, double *I_tp1, VECTOR_2D *grad_I_tp1, SIZE size, SIZE W, TUPLE_VEC_SCALAR Theta, TUPLE_VEC_SCALAR DTheta, COORDINATE X_c, double C);

TUPLE_VEC_SCALAR Chi(COORDINATE X_i, double Ix, double Iy);

TUPLE_VEC_SCALAR add_tuple(TUPLE_VEC_SCALAR Theta1, TUPLE_VEC_SCALAR Theta2);
double mult_tuple(TUPLE_VEC_SCALAR Theta1, TUPLE_VEC_SCALAR Theta2);
TUPLE_VEC_SCALAR coeff_tuple(TUPLE_VEC_SCALAR Theta, double C);
double norm_tuple(TUPLE_VEC_SCALAR Theta, SIZE W_l);

double Optical_Bilinear(double *Img, SIZE size, double x, double y);

int OpticalFlow_write(TUPLE_VEC_SCALAR *Theta, SIZE size, char *filename);

