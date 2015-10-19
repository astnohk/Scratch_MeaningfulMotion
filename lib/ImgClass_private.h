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
ImgVector<T>::ImgVector(ImgVector<T> &target)
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
			return;
		}
		_width = W;
		_height = H;
	}
}


template <typename T>
ImgVector<T>::ImgVector(int W, int H, T value)
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
ImgVector<T>::ImgVector(int W, int H, T *array)
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
	_data = nullptr;
	_width = 0;
	_height = 0;
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
			return;
		}
		_width = W;
		_height = H;
	}
}


template <typename T>
void
ImgVector<T>::reset(int W, int H, T value)
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
ImgVector<T>::reset(int W, int H, T *array)
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
ImgVector<T>::copy(ImgVector<T> &target)
{
	delete[] _data;
	_data = nullptr;
	_width = 0;
	_height = 0;

	if (target._width <= 0 || target._height <= 0) {
		return;
	}

	try {
		_data = new T[target._width * target._height]();
	}
	catch (const std::bad_alloc &bad) {
		fprintf(stderr, "ImgVector::ImgVector(T *, int, int) : Cannot Allocate Memory\n");
		_data = nullptr;
		return;
	}
	_width = target._width;
	_height = target._height;
	for (int i = 0; i < _width * _height; i++) {
		_data[i] = target._data[i];
	}
}


template <typename T>
void
ImgVector<T>::copy(ImgVector<T> *target)
{
	delete[] _data;
	_data = nullptr;
	_width = 0;
	_height = 0;

	if (target == nullptr
	    || target->_width <= 0 || target->_height <= 0) {
		return;
	}

	try {
		_data = new T[target->_width * target->_height]();
	}
	catch (const std::bad_alloc &bad) {
		fprintf(stderr, "ImgVector::ImgVector(T *, int, int) : Cannot Allocate Memory\n");
		_data = nullptr;
		return;
	}
	_width = target->_width;
	_height = target->_height;
	for (int i = 0; i < _width * _height; i++) {
		_data[i] = target->_data[i];
	}
}


template <typename T>
void
ImgVector<T>::set(int x, int y, T &value)
{
	if (x < 0 || _width <= x
	    || y < 0 || _height <= y) {
		return;
	}
	_data[_width * y + x] = value;
}


template <typename T>
ImgVector<T> &
ImgVector<T>::operator=(ImgVector<T> &copy)
{
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
		x_repeat = _width - (std::abs(x + 1) % _width);
	}
	if (y >= 0) {
		y_repeat = y % _height;
	} else {
		y_repeat = _height - (std::abs(y + 1) % _height);
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
		x_repeat = _width - (std::abs(x + 1) % _width);
	}
	if (y >= 0) {
		y_repeat = y % _height;
	} else {
		y_repeat = _height - (std::abs(y + 1) % _height);
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

