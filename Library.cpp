#include "Scratch_MeaningfulMotion.h"




/* Strings Library */
char *
regexp(char *string)
{
	ERROR Error("regexp");

	char *s = string;
	char *start_after = NULL;
	char tmp[4 + REGEXP_MAX_DIGITS];
	char *out = NULL;
	unsigned int len_before = 0;
	unsigned int len_after = 0;
	unsigned int num_sharp = 0;
	unsigned int length = 0;
	unsigned int i;

	while (*s != '\0' && *s != '#') {
		len_before++;
		s++;
	}
	while (*s != '\0' && *s == '#') {
		num_sharp++;
		s++;
	}
	start_after = s;
	while (*s != '\0') {
		len_after++;
		s++;
	}
	if (num_sharp > REGEXP_MAX_DIGITS) {
		Error.Others("Expression digits over REGEXP_MAX_DIGITS");
		exit(EXIT_FAILURE);
	} else if (num_sharp == 0) {
		fprintf(stderr, "*** regexp - do nothing because the filename '%s' does not include '#' ***\n", string);
		try {
			out = new char[strlen(string) + 1u];
		}
		catch (std::bad_alloc bad) {
			Error.Function("new");
			Error.Value("out");
			Error.Malloc();
			goto ExitError;
		}
		for (i = 0; i < strlen(string); i++) {
			out[i] = string[i];
		}
		string[i] = '\0';
		return out;
	}
	sprintf(tmp, "%%0%ud", num_sharp);
	length = len_before + strlen(tmp) + len_after;
	try {
		out = new char[length + 1];
	}
	catch (std::bad_alloc bad) {
		Error.Function("new");
		Error.Value("out");
		Error.Malloc();
		goto ExitError;
	}
	s = out;
	for (i = 0; i < len_before; i++) {
		*s = string[i];
		s++;
	}
	for (i = 0; i < strlen(tmp); i++) {
		*s = tmp[i];
		s++;
	}
	for (i = 0; i < len_after; i++) {
		*s = start_after[i];
		s++;
	}
	*s = '\0';
	return out;
// Error
ExitError:
	free(out);
	out = NULL;
	exit(EXIT_FAILURE);
}


/* Mathematical Library */
double
pow_int(double x, int a)
{
	double ans = 1.0;

	while (a != 0) {
		if ((abs(a) & 1) == 0) {
			x *= x;
			a /= 2;
		} else {
			if (a > 0) {
				ans *= x;
				a--;
			} else {
				ans /= x;
				a++;
			}
		}
	}
	return ans;
}


double
atan2_div_pi_table(int y, int x, SIZE *size)
{
	static double *table = NULL;
	static int M = 0;
	static int N = 0;
	int m, n;
	double angle;

	if (size != NULL) {
		if (size->width == 0 || size->height == 0) {
			free(table);
			table = NULL;
			return 0.0;
		} else {
			M = size->height;
			N = size->width;
			if ((table = (double *)calloc((size_t)size->width * size->height, sizeof(double))) == NULL) {
				fprintf(stderr, "*** atan2_table() error - Cannot allocate memory for (*table) by calloc() ***\n");
				return 0.0;
			}
			for (m = 0; m < size->height; m++) {
				for (n = 0; n < size->width; n++) {
					table[size->width * m + n] = atan2((double)m, (double)n) / M_PI;
				}
			}
		}
	}

	if (abs(x) > N || abs(y) > M) {
		angle = atan2((double)y, (double)x) / M_PI;
	} else if (x == 0 && y == 0) {
		angle = 0.0;
	} else if (x == 0) {
		angle = y > 0 ? 0.5 : -0.5;
	} else if (y == 0) {
		angle = x > 0 ? 0.0 : 1.0;
	} else {
		if (table == NULL) {
			fprintf(stderr, "*** atan2_table() error - This function is NOT initialized ***\n");
			return 0.0;
		}
		if (x * y > 0) {
			angle = table[N * abs(y) + abs(x)];
			if (x < 0) {
				angle -= 1.0;
			}
		} else if (x < 0) {
			angle = 1.0 - table[N * abs(y) + abs(x)];
		} else {
			angle = -table[N * abs(y) + abs(x)];
		}
	}
	return angle;
}


int*
Calc_k_l(SIZE size, double p, double ep)
{
	int *k_list = NULL;
	double C = 2 * p * (1 - p);
	int l, L;
	int k0, k_start;
	int count = 0;
	double progress = 0;

	L = (size.height > size.width) ? size.height : size.width;
	if ((k_list = (int *)calloc((size_t)(L + 1), sizeof(int))) == NULL) {
		fprintf(stderr, "calloc error on Calc_k_l\n");
		return NULL;
	}
	printf("[L =     0]   0.0%% |%s\x1b[1A\n", Progress_End.c_str());
#pragma omp parallel for schedule(dynamic) private(k_start, k0)
	for (l = 1; l <= L; l++) {
		k_start = floor(p * l + sqrt(C * l * (log(DIV_ANGLE) + log(size.height) + 2.0 * log(size.width) - log(ep))));
		if (k_start < 0) {
			k_start = 0;
		}
		for (k0 = k_start; k0 <= l; k0++) {
			k_list[l] = k0;
			if (Pr(k0, l, p) * size.width * size.width * DIV_ANGLE * size.height <= ep) {
				break;
			}
		}
#pragma omp critical
		{
			count++;
			if (round((double)count / L * 1000.0) > progress) {
				progress = round((double)count / L * 1000.0); // Take account of Overflow
				printf("\r[L = %5d] %5.1f%% |%s#\x1b[1A\n", count, progress * 0.1, Progress[NUM_PROGRESS * count / (1 + L)].c_str());
			}
		}
	}
	printf("\nComplete!\n");
	return k_list;
}


double
Pr(int k, int l, double p)
{
	double sum = 0.0;
	double Combination = 1.0;
	double p_i;
	double one_p;
	int i;

	// Initialize
	for (i = 1; i < k; i++) { // l! / ((l - (k - 1))! * (k - 1)!)
		Combination *= (l - i + 1.0) / i;
	}
	p_i = pow_int(p, k - 1); // p^(k - 1)
	one_p = pow_int(1.0 - p, l - k + 1); // (1 - p)^(l - k + 1)
	for (i = k; i <= l; i++) {
		Combination *= (l - i + 1.0) / i;
		p_i *= p;
		one_p /= (1.0 - p);
		sum += Combination * p_i * one_p;
	}
	return sum;
}

