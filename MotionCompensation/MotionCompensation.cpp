#include "MotionCompensation.h"



MotionCompensation::MotionCompensation(void)
{
	motion_estimated = false;
}

MotionCompensation::MotionCompensation(const MotionCompensation &copy) // copy constructor
{
	motion_estimated = copy.motion_estimated;
	_image_prev.copy(copy._image_prev);
	_image_next.copy(copy._image_next);
	_vector.copy(copy._vector);
}

MotionCompensation::MotionCompensation(int width, int height, const double *image_prev, const double *image_next)
{
	motion_estimated = false;
	if (width > 0 && height > 0 && image_prev != nullptr && image_next != nullptr) {
		_image_prev.reset(width, height, image_prev);
		_image_next.reset(width, height, image_next);
		_vector.reset(width, height);
	}
}

MotionCompensation::MotionCompensation(const ImgVector<double> &image_prev, const ImgVector<double> &image_next)
{
	motion_estimated = false;
	if (image_prev.isNULL() != false && image_next.isNULL() != false) {
		_image_prev.copy(image_prev);
		_image_next.copy(image_next);
		_vector.reset(image_prev.width(), image_prev.height());
	}
}

MotionCompensation::~MotionCompensation(void) // Destructor
{
}

MotionCompensation &
MotionCompensation::copy(const MotionCompensation &copy)
{
	motion_estimated = copy.motion_estimated;
	_image_prev.copy(copy._image_prev);
	_image_next.copy(copy._image_next);
	_vector.copy(copy._vector);
	return *this;
}


double
MotionCompensation::get_image_prev(int n) const
{
	return _image_prev.get(n);
}

double
MotionCompensation::get_image_prev(int x, int y) const
{
	return _image_prev.get(x, y);
}

double
MotionCompensation::get_image_next(int n) const
{
	return _image_next.get(n);
}

double
MotionCompensation::get_image_next(int x, int y) const
{
	return _image_next.get(x, y);
}

VECTOR_2D
MotionCompensation::get_vector(int n) const
{
	return _vector.get(n);
}

VECTOR_2D
MotionCompensation::get_vector(int x, int y) const
{
	return _vector.get(x, y);
}


VECTOR_2D &
MotionCompensation::operator[](int n) // Get reference to motion vector[n]
{
	return _vector[n];
}

ImgVector<VECTOR_2D> &
MotionCompensation::ref_vector(void) // Get reference to ImgVector<VECTOR_2D>
{
	return _vector;
}

VECTOR_2D &
MotionCompensation::ref_vector(int x, int y) // Get reference to motion vector[y][x]
{
	return _vector.ref(x, y);
}

