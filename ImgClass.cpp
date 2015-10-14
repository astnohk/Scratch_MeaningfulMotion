#include <cstdio>
#include <new>
#include "ImgClass.h"

// Define for old compiler implementation
#define nullptr ((void *)0)




template <typename T>
ImgVector<T>::ImgVector(void)
{
	data = nullptr;
	width = 0;
	height = 0;
}


template <typename T>
ImgVector<T>::ImgVector(int W, int H, T *array)
{
	if (data != nullptr) {
		delete[] data;
	}
	try {
		data = new T[W * H];
	}
	catch (const std::bad_alloc &bad) {
		fprintf(stderr, "ImgVector::ImgVector(T *, int, int) : Cannot Allocate Memory");
		return;
	}
	width = W;
	height = H;
	for (int i = 0; i < width * height; i++) {
		data[i] = array[i];
	}
}


template <typename T>
ImgVector<T>::~ImgVector(void)
{
	if (data != nullptr) {
		delete[] data;
	}
	data = nullptr;
	width = 0;
	height = 0;
}


template <typename T>
void
ImgVector<T>::reset(int W, int H, T *array)
{
	if (data != nullptr) {
		delete[] data;
	}
	try {
		data = new T[W * H];
	}
	catch (const std::bad_alloc &bad) {
		fprintf(stderr, "ImgVector::ImgVector(T *, int, int) : Cannot Allocate Memory");
		return;
	}
	width = W;
	height = H;
	for (int i = 0; i < width * height; i++) {
		data[i] = array[i];
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
	data[width * y + x] = value;
}


template <typename T>
T
ImgVector<T>::get(int x, int y) const
{
	if (x < 0 || width <= x
	    || y < 0 || height <= y) {
		return 0;
	}
	return data[width * y + x];
}


template <typename T>
T
ImgVector<T>::operator[](int n) const
{
	if (n < 0 || width * height <= n) {
		return 0;
	}
	return data[n];
}

