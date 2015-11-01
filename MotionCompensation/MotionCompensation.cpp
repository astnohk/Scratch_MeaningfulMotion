#include <cassert>
#include "MotionCompensation.h"



MotionCompensation::MotionCompensation(void)
{
	motion_compensated = false;
	_width = 0;
	_height = 0;
}

MotionCompensation::MotionCompensation(const MotionCompensation &copy) // copy constructor
{
	motion_compensated = copy.motion_compensated;
	_width = copy._width;
	_height = copy._height;
	_image_prev.copy(copy._image_prev);
	_image_next.copy(copy._image_next);
	_image_compensated.copy(copy._image_compensated);
	_vector.copy(copy._vector);
}

MotionCompensation::MotionCompensation(int width, int height, const double *image_prev, const double *image_next, const VECTOR_2D *vector)
{
	motion_compensated = false;
	_width = 0;
	_height = 0;
	if (width > 0 && height > 0) {
		if (image_prev == nullptr || image_next == nullptr || vector == nullptr) {
			throw std::invalid_argument("pointer is NULL");
		}
		_width = width;
		_height = height;
		_image_prev.reset(width, height, image_prev);
		_image_next.reset(width, height, image_next);
		_vector.reset(width, height, vector);
		_image_compensated.reset(width, height);
	}
}

MotionCompensation::MotionCompensation(int width, int height, const double *image_prev, const double *image_next, int v_width, int v_height, const VECTOR_2D *vector)
{
	motion_compensated = false;
	_width = 0;
	_height = 0;
	if (width > 0 && height > 0) {
		if (image_prev == nullptr || image_next == nullptr || vector == nullptr) {
			throw std::invalid_argument("pointer is NULL");
		}
		_width = width;
		_height = height;
		_image_prev.reset(width, height, image_prev);
		_image_next.reset(width, height, image_next);
		_vector.reset(width, height);
		// Projection of vectors to the scaled plane which has same range of images
		for (int y = 0; y < height; y++) {
			int Y = (int)floor(y * v_height / height);
			for (int x = 0; x < width; x++) {
				int X = (int)floor(x * v_width / width);
				_vector.ref(x, y) = vector[v_width * Y + X];
			}
		}
		_image_compensated.reset(width, height);
	}
}

MotionCompensation::MotionCompensation(const ImgVector<double> &image_prev, const ImgVector<double> &image_next, const ImgVector<VECTOR_2D> &vector)
{
	motion_compensated = false;
	_width = 0;
	_height = 0;
	if (image_prev.isNULL() || image_next.isNULL() || vector.isNULL()) {
		throw std::invalid_argument("variable is empty");
	}
	_width = image_prev.width();
	_height = image_prev.height();
	_image_prev.copy(image_prev);
	_image_next.copy(image_next);
	if (vector.width() == _width && vector.height() == _height) {
		_vector.copy(vector);
	} else {
		// Projection of small vector field to the scaled plane which has same range of images
		_vector.reset(_width, _height);
		for (int y = 0; y < _height; y++) {
			int Y = (int)floor(y * vector.height() / _height);
			for (int x = 0; x < _width; x++) {
				int X = (int)floor(x * vector.width() / _width);
				_vector.ref(x, y) = vector.get(X, Y);
			}
		}
	}
	_image_compensated.reset(image_prev.width(), image_prev.height());
}

MotionCompensation::MotionCompensation(const ImgVector<double> *image_prev, const ImgVector<double> *image_next, const ImgVector<VECTOR_2D> *vector)
{
	motion_compensated = false;
	_width = 0;
	_height = 0;
	if (image_prev == nullptr || image_next == nullptr || vector == nullptr) {
		throw std::invalid_argument("pointer is NULL");
	} else if (image_prev->isNULL() || image_next->isNULL() || vector->isNULL()) {
		throw std::invalid_argument("variable is empty");
	}
	_width = image_prev->width();
	_height = image_prev->height();
	_image_prev.copy(image_prev);
	_image_next.copy(image_next);
	if (vector->width() == _width && vector->height() == _height) {
		_vector.copy(vector);
	} else {
		// Projection of small vector field to the scaled plane which has same range of images
		_vector.reset(_width, _height);
		for (int y = 0; y < _height; y++) {
			int Y = (int)floor(y * vector->height() / _height);
			for (int x = 0; x < _width; x++) {
				int X = (int)floor(x * vector->width() / _width);
				_vector.ref(x, y) = vector->get(X, Y);
			}
		}
	}
	_image_compensated.reset(image_prev->width(), image_prev->height());
}

MotionCompensation::~MotionCompensation(void) // Destructor
{
}


MotionCompensation &
MotionCompensation::copy(const MotionCompensation &copy)
{
	motion_compensated = copy.motion_compensated;
	_width = copy._width;
	_height = copy._height;
	_image_prev.copy(copy._image_prev);
	_image_next.copy(copy._image_next);
	_image_compensated.copy(copy._image_compensated);
	_vector.copy(copy._vector);
	return *this;
}

MotionCompensation &
MotionCompensation::set(int width, int height, const double *image_prev, const double *image_next, const VECTOR_2D *vector)
{
	if (width <= 0 || height <= 0) {
		_width = 0;
		_height = 0;
		_image_prev.reset(0, 0);
		_image_next.reset(0, 0);
		_vector.reset(0, 0);
		_image_compensated.reset(0, 0);
		return *this;
	}
	if (image_prev == nullptr || image_next == nullptr || vector == nullptr) {
		throw std::invalid_argument("pointer is NULL");
	}

	motion_compensated = false;
	_width = width;
	_height = height;
	_image_prev.reset(_width, _height, image_prev);
	_image_next.reset(_width, _height, image_next);
	_vector.reset(_width, _height, vector);
	_image_compensated.reset(_width, _height);
	return *this;
}

MotionCompensation &
MotionCompensation::set(int width, int height, const double *image_prev, const double *image_next, int v_width, int v_height, const VECTOR_2D *vector)
{
	if (width <= 0 || height <= 0
	    || v_width <= 0 || v_height <= 0) {
		_width = 0;
		_height = 0;
		_image_prev.reset(0, 0);
		_image_next.reset(0, 0);
		_vector.reset(0, 0);
		_image_compensated.reset(0, 0);
		return *this;
	}
	if (image_prev == nullptr || image_next == nullptr || vector == nullptr) {
		throw std::invalid_argument("pointer is NULL");
	}

	motion_compensated = false;
	_width = width;
	_height = height;
	_image_prev.reset(_width, _height, image_prev);
	_image_next.reset(_width, _height, image_next);
	_vector.reset(_width, _height);
	// Projection of vectors to the scaled plane which has same range of images
	for (int y = 0; y < height; y++) {
		int Y = (int)floor(y * v_height / height);
		for (int x = 0; x < width; x++) {
			int X = (int)floor(x * v_width / width);
			_vector.ref(x, y) = vector[v_width * Y + X];
		}
	}
	_image_compensated.reset(_width, _height);
	return *this;
}

MotionCompensation &
MotionCompensation::set(const ImgVector<double> &image_prev, const ImgVector<double> &image_next, const ImgVector<VECTOR_2D> &vector)
{
	if (image_prev.isNULL() || image_next.isNULL() || vector.isNULL()) {
		throw std::invalid_argument("variable is empty");
	}

	motion_compensated = false;
	_width = image_prev.width();
	_height = image_prev.height();
	_image_prev.copy(image_prev);
	_image_next.copy(image_next);
	if (vector.width() == _width && vector.height() == _height) {
		_vector.copy(vector);
	} else {
		// Projection of small vector field to the scaled plane which has same range of images
		_vector.reset(_width, _height);
		for (int y = 0; y < _height; y++) {
			int Y = (int)floor(y * vector.height() / _height);
			for (int x = 0; x < _width; x++) {
				int X = (int)floor(x * vector.width() / _width);
				_vector.ref(x, y) = vector.get(X, Y);
			}
		}
	}
	_image_compensated.reset(_width, _height);
	return *this;
}

MotionCompensation &
MotionCompensation::set(const ImgVector<double> *image_prev, const ImgVector<double> *image_next, const ImgVector<VECTOR_2D> *vector)
{
	if (image_prev == nullptr || image_next == nullptr || vector == nullptr) {
		throw std::invalid_argument("pointer is NULL");
	} else if (image_prev->isNULL() || image_next->isNULL() || vector->isNULL()) {
		throw std::invalid_argument("variable is empty");
	}

	motion_compensated = false;
	_width = image_prev->width();
	_height = image_prev->height();
	_image_prev.copy(image_prev);
	_image_next.copy(image_next);
	if (vector->width() == _width && vector->height() == _height) {
		_vector.copy(vector);
	} else {
		// Projection of small vector field to the scaled plane which has same range of images
		_vector.reset(_width, _height);
		for (int y = 0; y < _height; y++) {
			int Y = (int)floor(y * vector->height() / _height);
			for (int x = 0; x < _width; x++) {
				int X = (int)floor(x * vector->width() / _width);
				_vector.ref(x, y) = vector->get(X, Y);
			}
		}
	}
	_image_compensated.reset(_width, _height);
	return *this;
}


int
MotionCompensation::width(void) const
{
	return _width;
}

int
MotionCompensation::height(void) const
{
	return _height;
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

double
MotionCompensation::get_image_compensated(int n)
{
	if (motion_compensated == false) {
		this->create_image_compensated();
	}
	assert(0 <= n && n < _width * _height);
	return _image_compensated.get(n);
}

double
MotionCompensation::get_image_compensated(int x, int y)
{
	if (motion_compensated == false) {
		this->create_image_compensated();
	}
	printf("x %d, y %d, _width %d, _height %d\n", x, y, _width, _height);
	assert(0 <= x && x < _width && 0 <= y && y < _height);
	return _image_compensated.get(x, y);
}


// return Reference
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


ImgVector<double> &
MotionCompensation::ref_image_compensated(void)
{
	if (motion_compensated == false) {
		this->create_image_compensated();
	}
	return _image_compensated;
}

double &
MotionCompensation::ref_image_compensated(int x, int y)
{
	if (motion_compensated == false) {
		this->create_image_compensated();
	}
	return _image_compensated[_width * y + x];
}

double &
MotionCompensation::operator[](int n) // Get reference to the compensated image[n]
{
	if (motion_compensated == false) {
		this->create_image_compensated();
	}
	return _image_compensated[n];
}


/*
  void MotionCompensation::create_image_compensated(void)

  Make entirely compensated image.
  All pixels are compensated by intensity of image_prev and motion vector.
*/
void
MotionCompensation::create_image_compensated(void)
{
	_image_compensated.reset(_width, _height);
	for (int y = 0; y < _height; y++) {
		for (int x = 0; x < _width; x++) {
			VECTOR_2D v = _vector.get(x, y);
			if (x + v.x < 0 || _width <= x + v.x
			    || y + v.y < 0 || _height <= y + v.y) {
				continue;
			}
			_image_compensated.set(x + v.x, y + v.y, _image_prev.get(x, y));
		}
	}
	motion_compensated = true;
}


/*
  void MotionCompensation::create_image_masked_compensated(ImgVector<bool> *mask)

  Make compensated image. It only compensate the pixel of masked.
  ImgVector<bool> *mask should mean mask of compensated image.
  If mask(x, y) == true then the pixel would be compensated
  and if mask(x, y) == false then the pixel hold the original (image_next) intensity.
*/
void
MotionCompensation::create_image_masked_compensated(ImgVector<bool> *mask)
{
	if (mask != nullptr) {
		_image_compensated.copy(_image_next); // Initialize with Original image (image_next)
		for (int y = 0; y < _height; y++) {
			for (int x = 0; x < _width; x++) {
				VECTOR_2D v = _vector.get(x, y);
				if (mask->get(x + v.x, y + v.y) == false) {
					continue;
				}
				_image_compensated.set(x + v.x, y + v.y, _image_prev.get(x, y));
			}
		}
		motion_compensated = true;
	}
}

