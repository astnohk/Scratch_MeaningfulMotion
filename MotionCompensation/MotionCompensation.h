#include <cstdio>
#include <new>




class MotionCompensation
{
	private:
		ImgVector<double> _image;
		ImgVector<VECTOR_2D> _vector;
	public:
		MotionCompensation(void);
		MotionCompensation(const MotionCompensation &copy);
		MotionCompensation(int width, int height, double *image);
		MotionCompensation(ImgVector<double> &image);
		~MotionCompensation(void);
		copy(MotionCompensation &copy);
		get_image(int n);
		get_image(int x, int y);
		get_vector(int n);
		get_vector(int x, int y);
		operator[](int n); // Get Motion vector
}

