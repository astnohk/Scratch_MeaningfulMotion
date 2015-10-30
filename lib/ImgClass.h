#include <stdexcept>


#ifndef nullptr
#define nullptr 0
#endif



#ifndef LIB_ImgClass
#define LIB_ImgClass

template <typename T>
class ImgVector
{
	private:
		T *_data;
		int _width;
		int _height;
	public:
		ImgVector(void);
		explicit ImgVector(const ImgVector<T> &copy);
		ImgVector(int W, int H);
		ImgVector(int W, int H, const T &value);
		ImgVector(int W, int H, const T *array);
		~ImgVector(void);
		void reset(int W, int H);
		void reset(int W, int H, const T &value);
		void reset(int W, int H, const T *array);
		ImgVector<T>& copy(const ImgVector<T> &copy);
		ImgVector<T>& copy(const ImgVector<T> *copy);
		ImgVector<T>& operator=(const ImgVector<T> &copy);
		void set(int x, int y, const T &value);

		// Data access
		T* data(void) const;
		int width(void) const;
		int height(void) const;
		int size(void) const;
		bool isNULL(void) const;

		T& operator[](int n);
		T& ref(int x, int y);
		T& ref_repeat(int x, int y);
		T& ref_mirror(int x, int y);
		T get(int n) const;
		T get(int x, int y) const;
		T get_zeropad(int x, int y) const;
		T get_repeat(int x, int y) const;
		T get_mirror(int x, int y) const;

		// Resampling
		void resize_zerohold(int W, int H);
		//bool resize_bicubic(int W, int H, double min = 0.0, double max = 0.0, T (*Nearest_Integer_Method)(double &d) = nullptr, double B = (1.0 / 3.0), double C = (1.0 / 3.0));
		void resize_bicubic(int W, int H, double min = 0.0, double max = 0.0, T (*Nearest_Integer_Method)(double &d) = nullptr, double B = (0.0 / 3.0), double C = (1.0 / 2.0));
		double cubic(double x, double B, double C);
		void map(T (*func)(T &value));
};

#include "ImgClass_private.h"

#endif

