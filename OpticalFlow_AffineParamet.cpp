#include "Scratch_MeaningfulMotion.h"
#include "OpticalFlow_AffineParamet.h"
#include "MultiResolution.h"



TUPLE_VEC_SCALAR *
OpticalFlow_AffineParamet(double *I_t, double *I_tp1, SIZE size, OPTICALFLOW_PARAM OpticalFlow_Param)
{
	char *FunctionName = "OpticalFlow_RMR";
	char *ErrorFunctionName = "";
	char *ErrorValueName = "";

	int Level = OpticalFlow_Param.Level;
	SIZE W = OpticalFlow_Param.WindowSize;
	int iter_max = OpticalFlow_Param.IRLS_Iter_Max;
	double d_th = OpticalFlow_Param.IRLS_Convergence_Threshold;
	double C_p = OpticalFlow_Param.IRLS_Min_C;

	double **I_t_levels = NULL;
	double **I_tp1_levels = NULL;
	TUPLE_VEC_SCALAR *Theta = NULL;
	TUPLE_VEC_SCALAR *DTheta = NULL;
	VECTOR_2D **grad_I_t_levels = NULL;
	VECTOR_2D **grad_I_tp1_levels = NULL;
	VECTOR_2D *grad_I_tp1 = NULL;
	SIZE size_l;
	SIZE W_l;
	double C;
	double C_L;
	COORDINATE X_c;
	int scale;
	int i, l;
	int index;
	int x, y;

	if (I_t == NULL) {
		ErrorValueName = "I_t";
		goto ErrorPointerNull;
	} else if (I_tp1 == NULL) {
		ErrorValueName = "I_tp1";
		goto ErrorPointerNull;
	} else if (size.width <= 0 || size.height <= 0) {
		ErrorValueName = "size";
		goto ErrorValueIncorrect;
	} else if (Level < 1) {
		ErrorValueName = "Level";
		goto ErrorValueIncorrect;
	} else if (W.width <= 0 || W.height <= 0) {
		ErrorValueName = "W";
		goto ErrorValueIncorrect;
	}

	/* Image Pyramid */
	if ((I_t_levels = Pyramider(I_t, size, Level)) == NULL) {
		ErrorFunctionName = "Pyramider";
		ErrorValueName = "I_t_levels";
		goto ErrorFunctionFailed;
	}
	if ((I_tp1_levels = Pyramider(I_tp1, size, Level)) == NULL) {
		ErrorFunctionName = "Pyramider";
		ErrorValueName = "I_tp1_levels";
		goto ErrorFunctionFailed;
	}
	if ((grad_I_t_levels = grad_Pyramid(I_t_levels, NULL, size, Level)) == NULL) {
		ErrorFunctionName = "grad_Pyramid";
		ErrorValueName = "grad_I_t_levels";
		goto ErrorFunctionFailed;
	}
	if ((grad_I_tp1_levels = grad_Pyramid(I_tp1_levels, NULL, size, Level)) == NULL) {
		ErrorFunctionName = "grad_Pyramid";
		ErrorValueName = "grad_I_tp1_levels";
		goto ErrorFunctionFailed;
	}

	/* Parameters */
	if ((Theta = (TUPLE_VEC_SCALAR *)calloc(size.width * size.height, sizeof(TUPLE_VEC_SCALAR))) == NULL) {
		ErrorFunctionName = "calloc";
		ErrorValueName = "Theta";
		goto ErrorMalloc;
	}
	if ((DTheta = (TUPLE_VEC_SCALAR *)calloc(size.width * size.height, sizeof(TUPLE_VEC_SCALAR))) == NULL) {
		ErrorFunctionName = "calloc";
		ErrorValueName = "DTheta";
		goto ErrorMalloc;
	}

	/* const C at L - 1 */
	size_l.width = (int)floor((double)size.width / pow_int(2.0, Level - 1));
	size_l.height = (int)floor((double)size.height / pow_int(2.0, Level - 1));
	C_L = .0;
	for (i = 0; i < size_l.width * size_l.height; i++) {
		if (fabs(I_t_levels[Level - 1][i]) > C_L) {
			C_L = fabs(I_t_levels[Level - 1][i]);
		}
	}
#pragma omp parallel for private(x, l, i, index, scale, C, X_c, size_l, W_l, grad_I_tp1)
	for (y = 0; y < size.width; y++) {
		for (x = 0; x < size.height; x++) {
			index = size.width * y + x;
			/* Level = L - 1 */
			C = C_L;
			DTheta[index] = TUPLE_VEC_SCALAR_NULL;
			Theta[index] = TUPLE_VEC_SCALAR_NULL;
			scale = (int)pow_int(2.0, Level - 1);
			size_l.width = (int)floor((double)size.width / scale);
			size_l.height = (int)floor((double)size.height / scale);
			W_l = W;
			X_c.x = (int)floor((double)x / scale);
			X_c.y = (int)floor((double)y / scale);
			/* Rest levels */
			for (l = Level - 1; l >= 0; l--) {
				scale = (int)pow_int(2.0, l);
				size_l.width = (int)floor((double)size.width / scale);
				size_l.height = (int)floor((double)size.height / scale);
				X_c.x = (int)floor((double)x / scale);
				X_c.y = (int)floor((double)y / scale);
				W_l.width = (int)ceil((double)W.width / scale);
				W_l.height = (int)ceil((double)W.height / scale);
				W_l = W;
				grad_I_tp1 = grad_I_tp1_levels[l];
				i = 0;
				do {
					C = 0.9 * C;
					if (C < C_p) {
						C = C_p;
					}
					DTheta[index] = IRLS_OpticalFlow_Affine(I_t_levels[l], I_tp1_levels[l], grad_I_tp1, size_l, W_l, Theta[index], DTheta[index], X_c, C);
					Theta[index] = add_tuple(Theta[index], DTheta[index]);
					i++;
				} while (i < iter_max && norm_tuple(DTheta[index], W_l) > d_th / scale);
				if (l != 0) {
					Theta[index].vector[0] = 2.0 * Theta[index].vector[0];
					Theta[index].vector[3] = 2.0 * Theta[index].vector[3];
				}
			}
		}
	}

	free(DTheta);
	for (l = 0; l < Level; l++) {
		free(grad_I_t_levels[l]);
		free(grad_I_tp1_levels[l]);
		free(I_t_levels[l]);
		free(I_tp1_levels[l]);
	}
	free(grad_I_t_levels);
	free(grad_I_tp1_levels);
	free(I_t_levels);
	free(I_tp1_levels);
	return Theta;
/* Error */
ErrorMalloc:
	fprintf(stderr, "*** %s() error - Cannot allocate memory for (*%s) by %s() ***\n", FunctionName, ErrorValueName, ErrorFunctionName);
	goto ErrorReturn;
ErrorFunctionFailed:
	fprintf(stderr, "*** %s() error - %s() failed to compute (%s) ***\n", FunctionName, ErrorFunctionName, ErrorValueName);
	goto ErrorReturn;
ErrorPointerNull:
	fprintf(stderr, "*** %s() error - The pointer (*%s) is NULL ***\n", FunctionName, ErrorValueName);
	goto ErrorReturn;
ErrorValueIncorrect:
	fprintf(stderr, "*** %s() error - The variable (%s) has incorrect value ***\n", FunctionName, ErrorValueName);
ErrorReturn:
	free(DTheta);
	free(Theta);
	if (grad_I_tp1_levels != NULL) {
		for (l = 0; l < Level; l++) {
			free(grad_I_tp1_levels[l]);
		}
	}
	free(grad_I_tp1_levels);
	if (I_t_levels != NULL) {
		for (l = 0; l < Level; l++) {
			free(I_t_levels[l]);
		}
	}
	if (I_tp1_levels != NULL) {
		for (l = 0; l < Level; l++) {
			free(I_tp1_levels[l]);
		}
	}
	free(I_t_levels);
	free(I_tp1_levels);
	return NULL;
}


TUPLE_VEC_SCALAR
IRLS_OpticalFlow_Affine(double *I_t, double *I_tp1, VECTOR_2D *grad_I_tp1, SIZE size, SIZE W, TUPLE_VEC_SCALAR Theta, TUPLE_VEC_SCALAR DTheta, COORDINATE X_c, double C)
{
	char *FunctionName = "IRLS_OpticalFlow_Affine";
	char *ErrorValueName = "";

	/* IRLS just compute the minimizing of
	 * E_3(DTheta) = sum_{Chi_i in W}(0.5 * w * (Chi_i * DTheta - Y_i)^2)
	 * DTheta = (sum(w * Chi_i * Chi_i))^(-1) * sum(w * Chi_i * Y_i) */
	VECTOR_2D grad_I_tp1_X_i_BA;
	int iter_max = 8;
	double residual;
	double sum_residual;
	double sum_residual_prev = .0;
	double Y_i;
	TUPLE_VEC_SCALAR Chi_i;
	double w;
	double div = .0;
	TUPLE_VEC_SCALAR sum = TUPLE_VEC_SCALAR_NULL;
	COORDINATE X_i, X_i_BA;
	int i, m, n;

	if (I_t == NULL) {
		ErrorValueName = "I_t";
		goto ErrorPointerNull;
	} else if (I_tp1 == NULL) {
		ErrorValueName = "I_tp1";
		goto ErrorPointerNull;
	} else if (size.width <= 0 || size.height <= 0) {
		ErrorValueName = "size";
		goto ErrorValueIncorrect;
	} else if (W.width <= 0 || W.height <= 0) {
		ErrorValueName = "W";
		goto ErrorValueIncorrect;
	}

	for (i = 0; i < iter_max; i++) {
		sum_residual = .0;
		for (m = 0; m < W.height; m++) {
			X_i.y = X_c.y + m - (int)floor(W.height / 2.0);
			for (n = 0; n < W.width; n++) {
				X_i.x = X_c.x + n - (int)floor(W.width / 2.0);
				X_i_BA.x = X_i.x + Theta.vector[0] + X_i.x * Theta.vector[1] + X_i.y * Theta.vector[2];
				X_i_BA.y = X_i.y + Theta.vector[3] + X_i.x * Theta.vector[4] + X_i.y * Theta.vector[5];
				if (0 <= X_i.x && X_i.x < size.width && 0 <= X_i.y && X_i.y < size.height) {
					Y_i = I_t[size.width * X_i.y + X_i.x];
				} else {
					Y_i = .0;
				}
				if (0 <= X_i_BA.x && X_i_BA.x < size.width && 0 <= X_i_BA.y && X_i_BA.y < size.height) {
					Y_i -= Optical_Bilinear(I_tp1, size, X_i_BA.x, X_i_BA.y);
				}
				Y_i -= Theta.scalar;

				if (0 <= X_i_BA.x && X_i_BA.x < size.width && 0 <= X_i_BA.y && X_i_BA.y < size.height) {
					grad_I_tp1_X_i_BA = grad_I_tp1[size.width * X_i_BA.y + X_i_BA.x];
				} else {
					grad_I_tp1_X_i_BA = VECTOR_2D_ZERO;
				}
				Chi_i = Chi(
				    X_i,
				    grad_I_tp1_X_i_BA.x,
				    grad_I_tp1_X_i_BA.y);
				residual = mult_tuple(Chi_i, DTheta) - Y_i;
				sum_residual += POW2(residual);
				if (residual < C) {
					w = POW2(1.0 - POW2(residual / C));
				} else {
					w = .0;
				}
				div += w * mult_tuple(Chi_i, Chi_i);
				sum = add_tuple(sum, coeff_tuple(Chi_i, w * Y_i));
			}
		}
		if (fabs(div) < 1E-50) {
			div = div < .0 ? -1E-50 : 1E-50;
		}
		DTheta = coeff_tuple(sum, 1.0 / div);
#ifdef SHOW_RESIDUAL
		printf("sum residual[%2d] = %.10f, C = %f\n", i, sum_residual, C);
#endif
		if (fabs(sum_residual - sum_residual_prev) < 1E-6) {
			break;
		}
		sum_residual_prev = sum_residual;
	}
	return DTheta;
/* Error */
ErrorPointerNull:
	fprintf(stderr, "*** %s() error - The pointer (*%s) is NULL ***\n", FunctionName, ErrorValueName);
	return TUPLE_VEC_SCALAR_NULL;
ErrorValueIncorrect:
	fprintf(stderr, "*** %s() error - The variable (%s) has incorrect value ***\n", FunctionName, ErrorValueName);
	return TUPLE_VEC_SCALAR_NULL;
}


TUPLE_VEC_SCALAR
Chi(COORDINATE X_i, double Ix, double Iy)
{
	TUPLE_VEC_SCALAR Chi_i;

	Chi_i.vector[0] = Ix;
	Chi_i.vector[1] = Ix * X_i.x;
	Chi_i.vector[2] = Ix * X_i.y;
	Chi_i.vector[3] = Iy;
	Chi_i.vector[4] = Iy * X_i.x;
	Chi_i.vector[5] = Iy * X_i.y;
	Chi_i.scalar = 1.0;
	return Chi_i;
}


TUPLE_VEC_SCALAR
add_tuple(TUPLE_VEC_SCALAR Theta1, TUPLE_VEC_SCALAR Theta2)
{
	/* (A, b) + (C, d) = (A + C, b + d) */
	TUPLE_VEC_SCALAR ans;
	int i;

	for (i = 0; i < TUPLE_VECTOR_SIZE; i++) {
		ans.vector[i] = Theta1.vector[i] + Theta2.vector[i];
	}
	ans.scalar = Theta1.scalar + Theta2.scalar;
	return ans;
}


double
mult_tuple(TUPLE_VEC_SCALAR Theta1, TUPLE_VEC_SCALAR Theta2)
{
	/* (A, b) * (C, d) = A * C + b * d */
	double ans = 0.0;
	int i;

	for (i = 0; i < TUPLE_VECTOR_SIZE; i++) {
		ans += Theta1.vector[i] + Theta2.vector[i];
	}
	ans += Theta1.scalar + Theta2.scalar;
	return ans;
}


TUPLE_VEC_SCALAR
coeff_tuple(TUPLE_VEC_SCALAR Theta, double C)
{
	/* (A, b) * C = (A * C, b * C) */
	TUPLE_VEC_SCALAR ans;
	int i;

	for (i = 0; i < TUPLE_VECTOR_SIZE; i++) {
		ans.vector[i] = Theta.vector[i] * C;
	}
	ans.scalar = Theta.scalar * C;
	return ans;
}


double
norm_tuple(TUPLE_VEC_SCALAR Theta, SIZE W_l)
{
	double s_j[3] = {1.0, .0, .0};
	double sum = .0;
	int i;

	s_j[1] = floor(W_l.width / 2.0) * (floor(W_l.width / 2.0) + 1.0) / W_l.width;
	s_j[2] = floor(W_l.height / 2.0) * (floor(W_l.width / 2.0) + 1.0) / W_l.height;
	for (i = 0; i < TUPLE_VECTOR_SIZE; i++) {
		sum += s_j[i % 3] * fabs(Theta.vector[i]);
	}
	return sum;
}


double
Optical_Bilinear(double *Img, SIZE size, double x, double y)
{
	char *FunctionName = "Bilinear()";

	int x_f, y_f;
	double tmp_f, tmp_c;

	if (Img == NULL) {
		fprintf(stderr, "*** %s error - The pointer (*Img) is NULL ***\n", FunctionName);
		return .0;
	}
	if (x <= -1.0 || size.width <= x) {
		return .0;
	} else if (y <= -1.0 || size.height <= y) {
		return .0;
	}

	x_f = floor(fabs(x));
	y_f = floor(fabs(y));
	if (x < 0) {
		tmp_f = Img[size.width * y_f] * x
		    + Img[size.width * y_f];
		if (0 <= y && y < size.height - 1) {
			tmp_c = Img[size.width * (y_f + 1)] * x + Img[size.width * (y_f + 1)];
		} else {
			tmp_c = .0;
		}
	} else if (size.width - 1 <= x) {
		tmp_f = Img[size.width * y_f + size.width - 1] * (x_f - x)
		    + Img[size.width * y_f + size.width - 1];
		if (0 <= y && y < size.height - 1) {
			tmp_c = Img[size.width * (y_f + 1) + size.width - 1] * (x_f - x)
			    + Img[size.width * (y_f + 1) + size.width - 1];
		} else {
			tmp_c = .0;
		}
	} else {
		tmp_f = (Img[size.width * y_f + x_f + 1] - Img[size.width * y_f + x_f]) * (x - x_f)
		    + Img[size.width * y_f + x_f];
		if (0 <= y && y < size.height - 1) {
			tmp_c = (Img[size.width * (y_f + 1) + x_f + 1] - Img[size.width * (y_f + 1) + x_f]) * (x - x_f)
			    + Img[size.width * (y_f + 1) + x_f];
		} else {
			tmp_c = .0;
		}
	}
	return (tmp_c - tmp_f) * (y - floor(y)) + tmp_f;
}


int
OpticalFlow_write(TUPLE_VEC_SCALAR *Theta, SIZE size, char *filename)
{
	char *FunctionName = "OpticalFlow_write()";
	char *ErrorFunctionName = "";
	char *ErrorValueName = "";

	FILE *fp = NULL;
	int n, i;

	if (Theta == NULL) {
		ErrorValueName = "Theta";
		goto ErrorPointerNull;
	} else if (size.width <= 0 || size.height <= 0) {
		ErrorValueName = "size";
		goto ErrorValueIncorrect;
	} else if (filename == NULL) {
		ErrorValueName = "filename";
		goto ErrorPointerNull;
	}

	printf("* Output Affine to '%s'\n", filename);
	if ((fp = fopen(filename, "w")) == NULL) {
		ErrorFunctionName = "fopen";
		ErrorValueName = filename;
		goto ErrorFileOpenFailed;
	}
	if (fprintf(fp, "%d %d\n", size.width, size.height) < 0) {
		ErrorFunctionName = "fprintf";
		ErrorValueName = "(size -> *fp)";
		goto ErrorFileOutputFailed;
	}
	for (n = 0; n < size.width * size.height; n++) {
		for (i = 0; i < TUPLE_VECTOR_SIZE; i++) {
			if (fprintf(fp, "%0.16e ", Theta[n].vector[i]) < 0) {
				ErrorFunctionName = "fprintf";
				ErrorValueName = "Theta[n].vector[i]";
				goto ErrorFileOutputFailed;
			}
		}
		if (fprintf(fp, "\n") < 0) {
			ErrorFunctionName = "fprintf";
			ErrorValueName = "\\n";
			goto ErrorFileOutputFailed;
		}
	}
	fclose(fp);
	return MEANINGFUL_SUCCESS;
/* Error */
ErrorFileOutputFailed:
	fprintf(stderr, "*** %s error - %s() failed to write out (%s) ***\n", FunctionName, ErrorFunctionName, ErrorValueName);
	goto ErrorReturn;
ErrorPointerNull:
	fprintf(stderr, "*** %s error - The pointer (*%s) is NULL ***\n", FunctionName, ErrorValueName);
	goto ErrorReturn;
ErrorValueIncorrect:
	fprintf(stderr, "*** %s error - The variable (%s) has incorrect value ***\n", FunctionName, ErrorValueName);
	goto ErrorReturn;
ErrorFileOpenFailed:
	fprintf(stderr, "*** %s error - Cannot open the file '%s' ***\n", FunctionName, ErrorValueName);
ErrorReturn:
	fclose(fp);
	return MEANINGFUL_FAILURE;
}

