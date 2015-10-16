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


template <typename T>
ImgVector<T>::ImgVector(ImgVector<T> &target)
{
	try {
		_data = new T[target._width * target._height];
	}
	catch (const std::bad_alloc &bad) {
		fprintf(stderr, "ImgVector::ImgVector(T *, int, int) : Cannot Allocate Memory");
		return;
	}
	_width = target._width;
	_height = target._height;
	for (int i = 0; i < _width * _height; i++) {
		_data[i] = target._data[i];
	}
}


template <typename T>
ImgVector<T>::ImgVector(int W, int H)
{
	try {
		_data = new T[W * H];
	}
	catch (const std::bad_alloc &bad) {
		fprintf(stderr, "ImgVector::ImgVector(T *, int, int) : Cannot Allocate Memory");
		return;
	}
	_width = W;
	_height = H;
}


template <typename T>
ImgVector<T>::ImgVector(int W, int H, T *array)
{
	try {
		_data = new T[W * H];
	}
	catch (const std::bad_alloc &bad) {
		fprintf(stderr, "ImgVector::ImgVector(T *, int, int) : Cannot Allocate Memory");
		return;
	}
	_width = W;
	_height = H;
	for (int i = 0; i < _width * _height; i++) {
		_data[i] = array[i];
	}
}


template <typename T>
ImgVector<T>::~ImgVector(void)
{
	if (_data != nullptr) {
		delete[] _data;
	}
	_data = nullptr;
	_width = 0;
	_height = 0;
}


template <typename T>
void
ImgVector<T>::reset(int W, int H)
{
	if (_data != nullptr) {
		delete[] _data;
	}
	try {
		_data = new T[W * H];
	}
	catch (const std::bad_alloc &bad) {
		fprintf(stderr, "ImgVector::ImgVector(T *, int, int) : Cannot Allocate Memory");
		return;
	}
	_width = W;
	_height = H;
}


template <typename T>
void
ImgVector<T>::reset(int W, int H, T *array)
{
	if (_data != nullptr) {
		delete[] _data;
	}
	try {
		_data = new T[W * H];
	}
	catch (const std::bad_alloc &bad) {
		fprintf(stderr, "ImgVector::ImgVector(T *, int, int) : Cannot Allocate Memory");
		return;
	}
	_width = W;
	_height = H;
	for (int i = 0; i < _width * _height; i++) {
		_data[i] = array[i];
	}
}


template <typename T>
void
ImgVector<T>::copy(ImgVector<T> &target)
{
	try {
		_data = new T[target._width * target._height];
	}
	catch (const std::bad_alloc &bad) {
		fprintf(stderr, "ImgVector::ImgVector(T *, int, int) : Cannot Allocate Memory");
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
ImgVector<T>::set(int x, int y, T &value)
{
	if (x < 0 || _width <= x
	    || y < 0 || _height <= y) {
		return;
	}
	_data[_width * y + x] = value;
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
	//return &(_data[_width * y + x]);
	return _data[_width * y + x];
}


template <typename T>
T &
ImgVector<T>::ref_repeat(int x, int y)
{
	return _data[_width * (y % height) + (x % width)];
}


template <typename T>
T &
ImgVector<T>::ref_mirror(int x, int y)
{
	int x_mirror, y_mirror;

	x_mirror = width - std::abs(x % (2 * width));
	y_mirror = height - std::abs(y % (2 * height));
	return _data[_width * y_mirror + x_mirror];
}


template <typename T>
T
ImgVector<T>::get(int n) const
{
	if (n < 0 || _width * _height <= n) {
		return 0;
	}
	return _data[_width * y + x];
}


template <typename T>
T
ImgVector<T>::get(int x, int y) const
{
	if (x < 0 || _width <= x
	    || y < 0 || _height <= y) {
		return 0;
	}
	return _data[_width * y + x];
}


template <typename T>
T
ImgVector<T>::get_repeat(int x, int y) const
{
	int x_mirror, y_mirror;

	x_mirror = _width - std::abs(_width - (x % (2 * _width)));
	y_mirror = _height - std::abs(_height - (y % (2 * _width)));
	return _data[_width * y_mirror + x_mirror];
}


template <typename T>
T
ImgVector<T>::get_mirror(int x, int y) const
{
	int x_mirror, y_mirror;

	x_mirror = _width - std::abs(_width - (x % (2 * _width)));
	y_mirror = _height - std::abs(_height - (y % (2 * _width)));
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

