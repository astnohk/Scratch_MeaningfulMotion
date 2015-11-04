#include <cstdio>
#include <new>
#include "../ImgClass/ImgClass.h"
#include "../lib/ExtVector.h"



#ifndef LIB_MotionCompensation
#define LIB_MotionCompensation

class MotionCompensation
{
	private:
		bool motion_compensated;
		int _width;
		int _height;
		ImgVector<double> _image_prev;
		ImgVector<double> _image_next;
		ImgVector<double> _image_compensated;
		ImgVector<VECTOR_2D<double> > _vector;
	public:
		MotionCompensation(void);
		MotionCompensation(const MotionCompensation &copy); // copy constructor
		MotionCompensation(int width, int height, const double *image_prev, const double *image_next, const VECTOR_2D<double> *vector);
		MotionCompensation(int width, int height, const double *image_prev, const double *image_next, int v_width, int v_height, const VECTOR_2D<double> *vector);
		MotionCompensation(const ImgVector<double> &image_prev, const ImgVector<double> &image_next, const ImgVector<VECTOR_2D<double> > &vector);
		MotionCompensation(const ImgVector<double> *image_prev, const ImgVector<double> *image_next, const ImgVector<VECTOR_2D<double> > *vector);
		~MotionCompensation(void);

		MotionCompensation& copy(const MotionCompensation &copy);
		MotionCompensation& set(int width, int height, const double *image_prev, const double *image_next, const VECTOR_2D<double> *vector);
		MotionCompensation& set(int width, int height, const double *image_prev, const double *image_next, int v_width, int v_height, const VECTOR_2D<double> *vector);
		MotionCompensation& set(const ImgVector<double> &image_prev, const ImgVector<double> &image_next, const ImgVector<VECTOR_2D<double> > &vector);
		MotionCompensation& set(const ImgVector<double> *image_prev, const ImgVector<double> *image_next, const ImgVector<VECTOR_2D<double> > *vector);
		ImgVector<double>& image_prev(int n);
		ImgVector<double>& image_next(int n);
		ImgVector<VECTOR_2D<double> >& vector(int n);

		int width(void) const;
		int height(void) const;

		double get_image_prev(int n) const;
		double get_image_prev(int x, int y) const;
		double get_image_next(int n) const;
		double get_image_next(int x, int y) const;
		VECTOR_2D<double> get_vector(int n) const;
		VECTOR_2D<double> get_vector(int x, int y) const;
		double get_image_compensated(int n);
		double get_image_compensated(int x, int y);

		ImgVector<VECTOR_2D<double> >& ref_vector(void); // Get reference to ImgVector<VECTOR_2D<double>>
		VECTOR_2D<double>& ref_vector(int x, int y); // Get reference to motion vector[y][x]

		ImgVector<double>& ref_image_compensated(void); // Get reference to the compensated image
		double& ref_image_compensated(int x, int y); // Get reference to the compensated image
		double& operator[](int n); // Get reference to compensated_image

		void create_image_compensated(void);
		void create_image_masked_compensated(ImgVector<bool> *mask);
};

#endif

