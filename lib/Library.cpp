#include "../Scratch_MeaningfulMotion.h"




// Strings Library
char *
regexp(const char* str)
{
	ERROR Error("regexp");

	const char* s = str;
	const char* start_after = nullptr;
	char tmp[4 + REGEXP_MAX_DIGITS];
	char* out = nullptr;
	char* sout = nullptr;
	size_t len_before = 0;
	size_t len_after = 0;

	s = str;
	while (*s != '\0' && *s != '#') {
		len_before++;
		s++;
	}
	if (*s == '\0') { // There are NO # in str
		try {
			out = new char[strlen(str) + 1u];
		}
		catch (const std::bad_alloc& bad) {
			std::cerr << bad.what() << std::endl;
			Error.Value("out");
			Error.Malloc();
			exit(EXIT_FAILURE);
		}
		strcpy(out, str);
		return out;
	}
	size_t num_sharp = 0;
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
	sprintf(tmp, "%%0%lud", num_sharp);
	size_t length = len_before + strlen(tmp) + len_after;
	try {
		out = new char[length + 1];
	}
	catch (const std::bad_alloc& bad) {
		std::cerr << bad.what() << std::endl;
		Error.Value("out");
		Error.Malloc();
		exit(EXIT_FAILURE);
	}
	sout = out;
	for (size_t i = 0; i < len_before; i++) {
		*sout = str[i];
		sout++;
	}
	for (size_t i = 0; i < strlen(tmp); i++) {
		*sout = tmp[i];
		sout++;
	}
	for (size_t i = 0; i < len_after; i++) {
		*sout = start_after[i];
		sout++;
	}
	*sout = '\0';
	return out;
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


int *
Calc_k_l(const SIZE& size, const double& p, const double& ep)
{
	const int L = (size.height > size.width) ? size.height : size.width;
	const double C = 2.0 * p * (1.0 - p);
	int *k_list = nullptr;

	try {
		k_list = new int[L + 1];
	}
	catch (const std::bad_alloc& bad) {
		std::cerr << bad.what() << std::endl
		    << "error : int* Calc_k_l(SIZE&, double, double) Cannot allocate memory" << std::endl;
		return nullptr;
	}

	int count = 0;
	double progress = 0.0;
	printf("[L =     0]   0.0%% |%s\x1b[1A\n", Progress_End.c_str());
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for (int l = 1; l <= L; l++) {
		int k_start = int(floor(p * l + sqrt(C * l * (log(DIV_ANGLE) + log(size.height) + 2.0 * log(size.width) - log(ep)))));
		if (k_start < 0) {
			k_start = 0;
		}
		for (int k0 = k_start; k0 <= l; k0++) {
			k_list[l] = k0;
			if (Pr(k0, l, p) * size.width * size.width * DIV_ANGLE * size.height <= ep) {
				break;
			}
		}
#ifdef _OPENMP
#pragma omp critical
#endif
		{
			count++;
			if (round(double(count) / L * 1000.0) > progress) {
				progress = round(double(count) / L * 1000.0); // Take account of Overflow
				printf("\r[L = %5d] %5.1f%% |%s#\x1b[1A\n", count, progress * 0.1, Progress[NUM_PROGRESS * count / (1 + L)].c_str());
			}
		}
	}
	printf("\nComplete!\n");
	return k_list;
}


double
Pr(const int k, const int l, const double& p)
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

