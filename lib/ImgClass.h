template <typename T>
class ImgVector
{
	private:
		T *_data;
		int _width;
		int _height;
	public:
		ImgVector(void);
		ImgVector(ImgVector<T> &copy);
		ImgVector(int W, int H);
		ImgVector(int W, int H, T value);
		ImgVector(int W, int H, T *array);
		~ImgVector(void);
		void reset(int W, int H);
		void reset(int W, int H, T value);
		void reset(int W, int H, T *array);
		void copy(ImgVector<T> &copy);
		void copy(ImgVector<T> *copy);
		ImgVector<T>& operator=(ImgVector<T> &copy);
		void set(int x, int y, T &value);

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
		bool resize_zerohold(int W, int H, T (*adder)(T &x1, T &x2) = nullptr, T (*multiplier)(T &value, double &d) = nullptr);
		bool resize_bicubic(int W, int H, double alpha, T (*adder)(T &x1, T &x2) = nullptr, T (*multiplier)(T &value, double &d) = nullptr);
		double bicubic(double x, double a);
		void map(T (*map_def)(T &value));
};


#include "ImgClass_private.h"

