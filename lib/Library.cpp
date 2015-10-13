#include "../Scratch_MeaningfulMotion.h"




// Strings Library
char *
regexp(char *string)
{
	ERROR Error("regexp");

	char *s = nullptr;
	char *start_after = nullptr;
	char tmp[4 + REGEXP_MAX_DIGITS];
	char *out = nullptr;
	unsigned int len_before = 0;
	unsigned int len_after = 0;
	unsigned int num_sharp = 0;
	unsigned int length = 0;
	unsigned int i;

	s = string;
	while (*s != '\0' && *s != '#') {
		len_before++;
		s++;
	}
	if (*s == '\0') { // There are NO # in string
		try {
			out = new char[strlen(string) + 1u];
		}
		catch (const std::bad_alloc &bad) {
			Error.Value("out");
			Error.Malloc();
			goto ExitError;
		}
		strcpy(out, string);
		return out;
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
	}
	sprintf(tmp, "%%0%ud", num_sharp);
	length = len_before + strlen(tmp) + len_after;
	try {
		out = new char[length + 1];
	}
	catch (const std::bad_alloc &bad) {
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
	delete[] out;
	out = nullptr;
	exit(EXIT_FAILURE);
}


// Mathematical Library
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


int*
Calc_k_l(SIZE &size, double p, double ep)
{
	const char *FunctionName = "Calc_k_l";
	const int L = (size.height > size.width) ? size.height : size.width;
	const double C = 2.0 * p * (1.0 - p);
	int *k_list = nullptr;
	int k_start, k0;
	int count = 0;
	double progress = 0.0;

	try {
		k_list = new int[L + 1];
	}
	catch (const std::bad_alloc &bad) {
		fprintf(stderr, "*** %s() error - Cannot allocate memory for (*k_list) ***\n", FunctionName);
		return nullptr;
	}
	printf("[L =     0]   0.0%% |%s\x1b[1A\n", Progress_End.c_str());
#pragma omp parallel for schedule(dynamic) private(k_start, k0)
	for (int l = 1; l <= L; l++) {
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

