#include <cassert>
#include <cmath>
#include <cstdio>
#include <new>




template <typename T>
ImgVector<T>::ImgVector(void)
{
	_data = nullptr;
	_width = 0;
	_height = 0;
}


// Copy Constructor
template <typename T>
ImgVector<T>::ImgVector(const ImgVector<T> &target)
{
	_data = nullptr;
	_width = 0;
	_height = 0;
	if (target._width > 0 && target._height > 0) {
		try {
			_data = new T[target._width * target._height]();
		}
		catch (const std::bad_alloc &bad) {
			fprintf(stderr, "ImgVector::ImgVector(T *, int, int) : Cannot Allocate Memory");
			_data = nullptr;
			throw;
			return;
		}
		_width = target._width;
		_height = target._height;
		for (int i = 0; i < _width * _height; i++) {
			_data[i] = target._data[i];
		}
	}
}


template <typename T>
ImgVector<T>::ImgVector(int W, int H)
{
	_data = nullptr;
	_width = 0;
	_height = 0;
	if (W > 0 && H > 0) {
		try {
			_data = new T[W * H]();
		}
		catch (const std::bad_alloc &bad) {
			fprintf(stderr, "ImgVector::ImgVector(T *, int, int) : Cannot Allocate Memory");
			_data = nullptr;
			throw;
			return;
		}
		_width = W;
		_height = H;
	}
}


template <typename T>
ImgVector<T>::ImgVector(int W, int H, const T &value)
{
	_data = nullptr;
	_width = 0;
	_height = 0;
	if (W > 0 && H > 0) {
		try {
			_data = new T[W * H]();
		}
		catch (const std::bad_alloc &bad) {
			fprintf(stderr, "ImgVector::ImgVector(T *, int, int) : Cannot Allocate Memory");
			_data = nullptr;
			throw;
			return;
		}
		_width = W;
		_height = H;
		for (int i = 0; i < _width * _height; i++) {
			_data[i] = value;
		}
	}
}


template <typename T>
ImgVector<T>::ImgVector(int W, int H, const T *array)
{
	_data = nullptr;
	_width = 0;
	_height = 0;
	if (W > 0 && H > 0) {
		try {
			_data = new T[W * H]();
		}
		catch (const std::bad_alloc &bad) {
			fprintf(stderr, "ImgVector::ImgVector(T *, int, int) : Cannot Allocate Memory");
			_data = nullptr;
			throw;
			return;
		}
		_width = W;
		_height = H;
		if (array != nullptr) {
			for (int i = 0; i < _width * _height; i++) {
				_data[i] = array[i];
			}
		}
	}
}


template <typename T>
ImgVector<T>::~ImgVector(void)
{
	delete[] _data;
}


template <typename T>
void
ImgVector<T>::reset(int W, int H)
{
	delete[] _data;
	_data = nullptr;
	_width = 0;
	_height = 0;
	if (W > 0 && H > 0) {
		try {
			_data = new T[W * H]();
		}
		catch (const std::bad_alloc &bad) {
			fprintf(stderr, "ImgVector::ImgVector(T *, int, int) : Cannot Allocate Memory");
			_data = nullptr;
			throw;
			return;
		}
		_width = W;
		_height = H;
	}
}


template <typename T>
void
ImgVector<T>::reset(int W, int H, const T &value)
{
	delete[] _data;
	_data = nullptr;
	_width = 0;
	_height = 0;
	if (W > 0 && H > 0) {
		try {
			_data = new T[W * H]();
		}
		catch (const std::bad_alloc &bad) {
			fprintf(stderr, "ImgVector::ImgVector(T *, int, int) : Cannot Allocate Memory");
			_data = nullptr;
			throw;
			return;
		}
		_width = W;
		_height = H;
		for (int i = 0; i < _width * _height; i++) {
			_data[i] = value;
		}
	}
}


template <typename T>
void
ImgVector<T>::reset(int W, int H, const T *array)
{
	delete[] _data;
	_data = nullptr;
	_width = 0;
	_height = 0;
	if (W > 0 && H > 0) {
		try {
			_data = new T[W * H]();
		}
		catch (const std::bad_alloc &bad) {
			fprintf(stderr, "ImgVector::ImgVector(T *, int, int) : Cannot Allocate Memory");
			_data = nullptr;
			throw;
			return;
		}
		_width = W;
		_height = H;
		if (array != nullptr) {
			for (int i = 0; i < _width * _height; i++) {
				_data[i] = array[i];
			}
		}
	}
}


template <typename T>
void
ImgVector<T>::copy(const ImgVector<T> &target)
{
	if (target._width > 0 && target._height > 0) {
		T *tmp_data = nullptr;
		try {
			tmp_data = new T[target._width * target._height]();
		}
		catch (const std::bad_alloc &bad) {
			fprintf(stderr, "ImgVector::ImgVector(T *, int, int) : Cannot Allocate Memory\n");
			throw;
			return;
		}
		_width = target._width;
		_height = target._height;
		delete[] _data;
		_data = tmp_data;
		for (int i = 0; i < _width * _height; i++) {
			_data[i] = target._data[i];
		}
	}
}


template <typename T>
void
ImgVector<T>::copy(const ImgVector<T> *target)
{
	if (target != nullptr && this != target
	    && target->_width > 0 && target->_height > 0) {
		T *tmp_data = nullptr;
		try {
			tmp_data = new T[target->_width * target->_height]();
		}
		catch (const std::bad_alloc &bad) {
			fprintf(stderr, "ImgVector::ImgVector(T *, int, int) : Cannot Allocate Memory\n");
			throw;
			return;
		}
		_width = target->_width;
		_height = target->_height;
		delete _data;
		_data = tmp_data;
		for (int i = 0; i < _width * _height; i++) {
			_data[i] = target->_data[i];
		}
	}
}


template <typename T>
void
ImgVector<T>::set(int x, int y, const T &value)
{
	if (x < 0 || _width <= x) {
		throw std::out_of_range("int x");
		return;
	} else if (y < 0 || _height <= y) {
		throw std::out_of_range("int y");
		return;
	}
	_data[_width * y + x] = value;
}


template <typename T>
ImgVector<T> &
ImgVector<T>::operator=(const ImgVector<T> &copy)
{
	if (this == &copy) {
		return *this;
	}
	delete[] _data;
	_data = nullptr;
	_width = 0;
	_height = 0;
	if (copy._width > 0 && copy._height > 0) {
		try {
			_data = new T[copy._width * copy._height]();
		}
		catch (const std::bad_alloc &bad) {
			fprintf(stderr, "ImgVector::operator=(ImgVector<T> &) : Cannot Allocate Memory\n");
			_data = nullptr;
			throw;
			return *this;
		}
		_width = copy._width;
		_height = copy._height;
		for (int i = 0; i < _width * _height; i++) {
			_data[i] = copy._data[i];
		}
	}
	return *this;
}


template <typename T>
T *
ImgVector<T>::data(void) const
{
	return _data;
}


template <typename T>
T &
ImgVector<T>::operator[](int n)
{
	assert(0 <= n && n < _width * _height);
	return _data[n];
}


template <typename T>
T &
ImgVector<T>::ref(int x, int y)
{
	assert(0 <= x && x < _width && 0 <= y && y < _height);
	return _data[_width * y + x];
}


template <typename T>
T &
ImgVector<T>::ref_repeat(int x, int y)
{
	int x_repeat, y_repeat;

	if (x >= 0) {
		x_repeat = x % _width;
	} else {
		x_repeat = _width - ((int)std::abs((double)x + 1.0) % _width);
	}
	if (y >= 0) {
		y_repeat = y % _height;
	} else {
		y_repeat = _height - ((int)std::abs((double)y + 1.0) % _height);
	}
	return _data[_width * y_repeat + x_repeat];
}


template <typename T>
T &
ImgVector<T>::ref_mirror(int x, int y)
{
	int x_mirror, y_mirror;

	if (x < 0) {
		x = -x - 1; // Mirroring over negative has offset
	}
	if (y < 0) {
		y = -y - 1;
	}
	x_mirror = (int)round(_width - 0.5 - std::fabs(_width - 0.5 - (x % (2 * _width))));
	y_mirror = (int)round(_height - 0.5 - std::fabs(_height - 0.5 - (y % (2 * _height))));
	return _data[_width * y_mirror + x_mirror];
}


template <typename T>
T
ImgVector<T>::get(int n) const
{
	assert(0 <= n && n < _width * _height);
	return _data[n];
}


template <typename T>
T
ImgVector<T>::get(int x, int y) const
{
	assert(0 <= x && x < _width
	    && 0 <= y && y < _height);
	return _data[_width * y + x];
}


template <typename T>
T
ImgVector<T>::get_zeropad(int x, int y) const
{
	T zero = T();

	if (x < 0 || _width <= x
	    || y < 0 || _height <= y) {
		return zero;
	} else {
		return _data[_width * y + x];
	}
}


template <typename T>
T
ImgVector<T>::get_repeat(int x, int y) const
{
	int x_repeat, y_repeat;

	if (x >= 0) {
		x_repeat = x % _width;
	} else {
		x_repeat = _width - ((int)std::abs((double)x + 1.0) % _width);
	}
	if (y >= 0) {
		y_repeat = y % _height;
	} else {
		y_repeat = _height - ((int)std::abs((double)y + 1.0) % _height);
	}
	return _data[_width * y_repeat + x_repeat];
}


template <typename T>
T
ImgVector<T>::get_mirror(int x, int y) const
{
	int x_mirror, y_mirror;

	if (x < 0) {
		x = -x - 1; // Mirroring over negative has offset
	}
	if (y < 0) {
		y = -y - 1;
	}
	x_mirror = (int)round(_width - 0.5 - std::fabs(_width - 0.5 - (x % (2 * _width))));
	y_mirror = (int)round(_height - 0.5 - std::fabs(_height - 0.5 - (y % (2 * _height))));
	return _data[_width * y_mirror + x_mirror];
}


template <typename T>
int
ImgVector<T>::width(void) const
{
	return _width;
}


template <typename T>
int
ImgVector<T>::height(void) const
{
	return _height;
}


template <typename T>
int
ImgVector<T>::size(void) const
{
	return _width * _height;
}


template <typename T>
bool
ImgVector<T>::isNULL(void) const
{
	if (_data == nullptr) {
		return true;
	} else {
		return false;
	}
}


template <typename T>
void
ImgVector<T>::resize_zerohold(int W, int H)
{
	T *resized = nullptr;
	T additive_identity = T();
	double scale_x = .0;
	double scale_y = .0;
	int area_x;
	int area_y;
	int m, n;
	int x, y;
	T sum;

	if (W <= 0 || H <= 0) {
		throw std::out_of_range("int W, int H");
		return;
	}
	scale_x = (double)W / _width;
	scale_y = (double)H / _height;
	try {
		resized = new T[W * H]();
	}
	catch (const std::bad_alloc &bad) {
		fprintf(stderr, "ImgVector<T>::resize_zerohold(int, int) error : Cannot allocate memory\n");
		throw;
		return;
	}
	area_x = ceil((double)_width / W);
	area_y = ceil((double)_height / H);
	for (y = 0; y < H; y++) {
		for (x = 0; x < W; x++) {
			sum = additive_identity;
			for (m = 0; m < area_y; m++) {
				for (n = 0; n < area_x; n++) {
					sum += this->get((int)floor(x / scale_x) + n, (int)floor(y / scale_y) + m);
				}
			}
			resized[W * y + x] = sum / (area_x * area_y);
		}
	}
	delete[] _data;
	_data = resized;
	_width = W;
	_height = H;
}


/*
    bool ImgVector<T>::resize_bicubic(int W, int H, double min, double max, T (*Nearest_Integer_Method)(double &d), double A)
    int W, int H : width and height of resized image
    double min : minimum value of saturated value
    double max : maximum value of saturated value
    T (*Nearest_Integer_Method)(double &d) : round method (e.g. floor(), round(), etc.)
    A : cubic method's parameter (default A = -0.5 which correspond to Hermite)
*/
template <typename T>
void
ImgVector<T>::resize_bicubic(int W, int H, double min, double max, T (*Nearest_Integer_Method)(double &d), double B, double C)
{
	T *resized = nullptr;
	double *conv = nullptr;
	ImgVector<double> Tmp;
	double scale_x, scale_y;
	double scale_conv;
	int L, L_center;
	double dx, dy;
	int x, y;
	int m, n;
	double sum;

	if (W <= 0 || H <= 0) {
		throw std::out_of_range("int W, int H");
		return;
	}
	scale_x = (double)W / _width;
	scale_y = (double)H / _height;
	Tmp.reset(W, _height);
	// The length of cubic convolution coefficient
	scale_conv = 1.0;
	if (scale_x < 1.0 || scale_y < 1.0) {
		scale_conv = ceil(1.0 / (scale_x < scale_y ? scale_x : scale_y));
	}
	try {
		resized = new T[W * H]();
		conv = new double[(int)scale_conv * 4]();
	}
	catch (const std::bad_alloc &bad) {
		fprintf(stderr, "ImgVector<double>::resize_bicubic(int, int) error : Cannot allocate memory\n");
		delete[] resized;
		delete[] conv;
		throw;
		return;
	}
	// Horizontal convolution
	if (scale_x >= 1.0) {
		scale_conv = 1.0;
		L = 4;
		L_center = floor((L - 1.0) / 2);
	} else {
		scale_conv = 1.0 / scale_x;
		L = 4 * (int)ceil(scale_conv);
		L_center = floor((L - 1.0) / 2);
	}
	for (x = 0; x < W; x++) {
		if (scale_x >= 1.0) {
			dx = (x - (scale_x - 1.0) / 2.0) / scale_x;
			for (n = 0; n < L; n++) {
				conv[n] = ImgVector<T>::cubic((double)(n - L_center) - (dx - floor(dx)), B, C);
			}
		} else {
			dx = x / scale_x + (1.0 / scale_x - 1.0) / 2.0;
			for (n = 0; n < L; n++) {
				conv[n] = ImgVector<T>::cubic(((double)(n - L_center) - (dx - floor(dx))) * scale_x, B, C) / scale_conv;
			}
		}
		for (y = 0; y < _height; y++) {
			sum = 0.0;
			for (n = 0; n < L; n++) {
				sum += conv[n] * this->get_mirror((int)floor(dx) + n - L_center, y);
			}
			Tmp.ref(x, y) = sum;
		}
	}
	// Vertical convolution
	if (scale_y >= 1.0) {
		scale_conv = 1.0;
		L = 4;
		L_center = floor((L - 1.0) / 2);
	} else {
		scale_conv = 1.0 / scale_y;
		L = 4 * (int)ceil(scale_conv);
		L_center = floor((L - 1.0) / 2);
	}
	for (y = 0; y < H; y++) {
		if (scale_y >= 1.0) {
			dy = (y - (scale_y - 1.0) / 2.0) / scale_y;
			for (m = 0; m < L; m++) {
				conv[m] = ImgVector<T>::cubic((double)(m - L_center) - (dy - floor(dy)), B, C);
			}
		} else {
			dy = y / scale_y + (1.0 / scale_y - 1.0) / 2.0;
			for (m = 0; m < L; m++) {
				conv[m] = ImgVector<T>::cubic(((double)(m - L_center) - (dy - floor(dy))) / scale_conv, B, C) / scale_conv;
			}
		}
		for (x = 0; x < W; x++) {
			sum = 0.0;
			for (m = 0; m < L; m++) {
				sum += conv[m] * Tmp.get_mirror(x, (int)floor(dy) + m - L_center);
			}
			if (min != max) {
				sum = sum >= min ? sum <= max ? sum : max : min;
			}
			if (Nearest_Integer_Method != nullptr) {
				resized[W * y + x] = Nearest_Integer_Method(sum);
			} else {
				resized[W * y + x] = sum;
			}
		}
	}
	delete[] conv;
	delete[] _data;
	_data = resized;
	_width = W;
	_height = H;
	return true;
}


template <typename T>
double
ImgVector<T>::cubic(double x, double B, double C)
{
	double x_abs = fabs(x);

	if (x_abs <= 1.0) {
		return ((2.0 - 1.5 * B - C) * x_abs + (-3.0 + 2.0 * B + C)) * x_abs * x_abs + 1.0 - B / 3.0;
	} else if (x_abs < 2.0) {
		return (((-B / 6.0 - C) * x_abs + B + 5.0 * C) * x_abs - 2.0 * B - 8.0 * C) * x_abs + 8.0 / 6.0 * B + 4.0 * C;
	} else {
		return 0.0;
	}
}


template <typename T>
void
ImgVector<T>::map(T (*func)(T &value))
{
	if (func != nullptr) {
		for (int i = 0 ; i < _width * _height; i++) {
			_data[i] = func(_data[i]);
		}
	}
}

