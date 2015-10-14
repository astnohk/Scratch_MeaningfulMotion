template <typename T>
class ImgVector
{
	private:
		T *Data;
		int width;
		int height;
	public:
		ImgVector(void);
		ImgVector(ImgVector<T> &copy);
		ImgVector(int W, int H);
		ImgVector(int W, int H, T *array);
		~ImgVector(void);
		void reset(int W, int H);
		void reset(int W, int H, T *array);
		void set(int x, int y, T value);
		T* data(void) const;
		T get(int x, int y) const;
		T operator[](int n) const;
};


#include "ImgClass_private.h"

