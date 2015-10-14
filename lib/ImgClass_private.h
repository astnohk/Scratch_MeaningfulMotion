#include <cstdio>
#include <new>




template <typename T>
ImgVector<T>::ImgVector(void)
{
	Data = nullptr;
	width = 0;
	height = 0;
}


template <typename T>
ImgVector<T>::ImgVector(ImgVector<T> &copy)
{
	try {
		Data = new T[copy.width * copy.height];
	}
	catch (const std::bad_alloc &bad) {
		fprintf(stderr, "ImgVector::ImgVector(T *, int, int) : Cannot Allocate Memory");
		return;
	}
	width = copy.width;
	height = copy.height;
	for (int i = 0; i < width * height; i++) {
		Data[i] = 0.0;
	}
}


template <typename T>
ImgVector<T>::ImgVector(int W, int H)
{
	try {
		Data = new T[W * H];
	}
	catch (const std::bad_alloc &bad) {
		fprintf(stderr, "ImgVector::ImgVector(T *, int, int) : Cannot Allocate Memory");
		return;
	}
	width = W;
	height = H;
	for (int i = 0; i < width * height; i++) {
		Data[i] = 0.0;
	}
}


template <typename T>
ImgVector<T>::ImgVector(int W, int H, T *array)
{
	try {
		Data = new T[W * H];
	}
	catch (const std::bad_alloc &bad) {
		fprintf(stderr, "ImgVector::ImgVector(T *, int, int) : Cannot Allocate Memory");
		return;
	}
	width = W;
	height = H;
	for (int i = 0; i < width * height; i++) {
		Data[i] = array[i];
	}
}


template <typename T>
ImgVector<T>::~ImgVector(void)
{
	if (Data != nullptr) {
		delete[] Data;
	}
	Data = nullptr;
	width = 0;
	height = 0;
}


template <typename T>
void
ImgVector<T>::reset(int W, int H)
{
	if (Data != nullptr) {
		delete[] Data;
	}
	try {
		Data = new T[W * H];
	}
	catch (const std::bad_alloc &bad) {
		fprintf(stderr, "ImgVector::ImgVector(T *, int, int) : Cannot Allocate Memory");
		return;
	}
	width = W;
	height = H;
	for (int i = 0; i < width * height; i++) {
		Data[i] = 0.0;
	}
}


template <typename T>
void
ImgVector<T>::reset(int W, int H, T *array)
{
	if (Data != nullptr) {
		delete[] Data;
	}
	try {
		Data = new T[W * H];
	}
	catch (const std::bad_alloc &bad) {
		fprintf(stderr, "ImgVector::ImgVector(T *, int, int) : Cannot Allocate Memory");
		return;
	}
	width = W;
	height = H;
	for (int i = 0; i < width * height; i++) {
		Data[i] = array[i];
	}
}


template <typename T>
void
ImgVector<T>::set(int x, int y, T value)
{
	if (x < 0 || width <= x
	    || y < 0 || height <= y) {
		return;
	}
	Data[width * y + x] = value;
}


template <typename T>
T *
ImgVector<T>::data(void) const
{
	return Data;
}


template <typename T>
T
ImgVector<T>::get(int x, int y) const
{
	if (x < 0 || width <= x
	    || y < 0 || height <= y) {
		return 0;
	}
	return Data[width * y + x];
}


template <typename T>
T
ImgVector<T>::operator[](int n) const
{
	if (n < 0 || width * height <= n) {
		return 0;
	}
	return Data[n];
}

