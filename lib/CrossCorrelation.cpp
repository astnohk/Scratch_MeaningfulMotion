#include "CrossCorrelation.h"




ImgStatistics::ImgStatistics(void)
{
	_width = 0;
	_height = 0;
	_data = nullptr;
}

ImgStatistics::ImgStatistics(const ImgStatistics &copy)
{
	_width = 0;
	_height = 0;
	_data = nullptr;
	try {
		_data = new double[copy._width * copy._height];
	}
	catch (const std::bad_alloc &bad) {
		fprintf(stderr, "ImgStatistics::ImgStatistics(const ImgStatistics &) error : memory allocation\n");
		_data = nullptr;
		return;
	}
	_width = copy._width;
	_height = copy._height;
	for (int i = 0; i < _width * _height; i++) {
		_data[i] = copy._data[i];
	}
}

ImgStatistics::ImgStatistics(int W, int H, double *Img)
{
	_width = 0;
	_height = 0;
	_data = nullptr;
	if (W > 0 && H > 0) {
		try {
			_data = new double[W * H]();
		}
		catch (const std::bad_alloc &bad) {
			fprintf(stderr, "ImgStatistics::ImgStatistics(int, int, double *) error : memory allocation\n");
			return;
		}
		_width = W;
		_height = H;
		if (Img != nullptr) {
			for (int i = 0; i < W * H; i++) {
				_data[i] = Img[i];
			}
		}
	}
}

ImgStatistics::~ImgStatistics(void)
{
	delete _data;
}

void
ImgStatistics::set(int W, int H, double *Img)
{
	_width = 0;
	_height = 0;
	delete[] _data;
	_data = nullptr;
	if (W > 0 && H > 0) {
		try {
			_data = new double[W * H]();
		}
		catch (const std::bad_alloc &bad) {
			fprintf(stderr, "ImgStatistics::set(int, int, double *) error : memory allocation\n");
			return;
		}
		_width = W;
		_height = H;
		if (Img != nullptr) {
			for (int i = 0; i < W * H; i++) {
				_data[i] = Img[i];
			}
		}
	}
}

ImgStatistics &
ImgStatistics::copy(const ImgStatistics &copy)
{
	if (this != &copy) {
		double *tmp_data = nullptr;
		try {
			tmp_data = new double[copy._width * copy._height];
		}
		catch (const std::bad_alloc &bad) {
			fprintf(stderr, "ImgStatistics::ImgStatistics(const ImgStatistics &) error : memory allocation\n");
			return *this;
		}
		_width = copy._width;
		_height = copy._height;
		delete[] _data;
		_data = tmp_data;
		for (int i = 0; i < _width * _height; i++) {
			_data[i] = copy._data[i];
		}
	}
	return *this;
}

ImgStatistics &
ImgStatistics::operator=(const ImgStatistics &copy)
{
	if (this != &copy) {
		double *tmp_data = nullptr;
		try {
			tmp_data = new double[copy._width * copy._height];
		}
		catch (const std::bad_alloc &bad) {
			fprintf(stderr, "ImgStatistics::ImgStatistics(const ImgStatistics &) error : memory allocation\n");
			return *this;
		}
		_width = copy._width;
		_height = copy._height;
		delete[] _data;
		_data = tmp_data;
		for (int i = 0; i < _width * _height; i++) {
			_data[i] = copy._data[i];
		}
	}
	return *this;
}

double &
ImgStatistics::image(int x, int y)
{
	return _data[_width * y + x];
}

int
ImgStatistics::width(void)
{
	return _width;
}

int
ImgStatistics::height(void)
{
	return _height;
}

double
ImgStatistics::mean(void)
{
	double sum = 0.0;

	for (int y = 0; y < _height; y++) {
		for (int x = 0; x < _width; x++) {
			sum += _data[_width * y + x];
		}
	}
	return sum / (_width * _height);
}

double
ImgStatistics::mean(int center_x, int center_y, int window_width, int window_height)
{
	double sum = 0.0;
	int x_tmp, y_tmp;

	for (int y = 0; y < window_height; y++) {
		y_tmp = center_y + y - (window_height - 1) / 2;
		if (y_tmp < 0 || _height <= y_tmp) {
			continue;
		}
		for (int x = 0; x < window_width; x++) {
			x_tmp = center_x + x - (window_width - 1) / 2;
			if (x_tmp < 0 || _width <= x_tmp) {
				continue;
			}
			sum += _data[_width * y_tmp + x_tmp];
		}
	}
	return sum / (window_width * window_height);
}

double
ImgStatistics::variance(void)
{
	double sum = 0.0;
	double mu;

	mu = this->mean();
	for (int n = 0; n < _width * _height; n++) {
		sum += (_data[n] - mu) * (_data[n] - mu);
	}
	return sum;
}

double
ImgStatistics::std_deviation(void)
{
	double sum = 0.0;
	double mu;

	mu = this->mean();
	for (int n = 0; n < _width * _height; n++) {
		sum += (_data[n] - mu) * (_data[n] - mu);
	}
	return sqrt(sum);
}

double
ImgStatistics::variance(int center_x, int center_y, int window_width, int window_height)
{
	double sum = 0.0;
	double mu;
	int x_tmp, y_tmp;

	mu = this->mean(center_x, center_y, window_width, window_height);
	for (int y = 0; y < window_height; y++) {
		y_tmp = center_y + y - (window_height - 1) / 2;
		if (y_tmp < 0 || _height <= y_tmp) {
			continue;
		}
		for (int x = 0; x < window_width; x++) {
			x_tmp = center_x + x - (window_width - 1) / 2;
			if (x_tmp < 0 || _width <= x_tmp) {
				continue;
			}
			sum += (_data[_width * y_tmp + x_tmp] - mu) * (_data[_width * y_tmp + x_tmp] - mu);
		}
	}
	return sum;
}

double
ImgStatistics::std_deviation(int center_x, int center_y, int window_width, int window_height)
{
	double sum = 0.0;
	double mu;
	int x_tmp, y_tmp;

	mu = this->mean(center_x, center_y, window_width, window_height);
	for (int y = 0; y < window_height; y++) {
		y_tmp = center_y + y - (window_height - 1) / 2;
		if (y_tmp < 0 || _height <= y_tmp) {
			continue;
		}
		for (int x = 0; x < window_width; x++) {
			x_tmp = center_x + x - (window_width - 1) / 2;
			if (x_tmp < 0 || _width <= x_tmp) {
				continue;
			}
			sum += (_data[_width * y_tmp + x_tmp] - mu) * (_data[_width * y_tmp + x_tmp] - mu);
		}
	}
	return sqrt(sum);
}




CrossCorrelation::CrossCorrelation(void)
{
}

CrossCorrelation::CrossCorrelation(const CrossCorrelation &copy) // copy constructor
{
	_width = copy._width;
	_height = copy._height;
	_img0.copy(copy._img0);
	_img1.copy(copy._img1);
}

CrossCorrelation::CrossCorrelation(ImgStatistics &img0, ImgStatistics &img1)
{
	if (img0.width() != img1.width()
	    || img0.height() != img1.height()) {
		return;
	} else {
		_img0.copy(img0);
		_img1.copy(img1);
		_width = img0.width();
		_height = img0.height();
	}
}

CrossCorrelation::~CrossCorrelation(void)
{
}

CrossCorrelation &
CrossCorrelation::copy(const CrossCorrelation &copy)
{
	if (this != &copy) {
		_img0.copy(copy._img0);
		_img1.copy(copy._img1);
		_width = copy._width;
		_height = copy._height;
	}
	return *this;
}

CrossCorrelation &
CrossCorrelation::operator=(const CrossCorrelation &copy)
{
	if (this != &copy) {
		_img0.copy(copy._img0);
		_img1.copy(copy._img1);
		_width = copy._width;
		_height = copy._height;
	}
	return *this;
}

int
CrossCorrelation::width(void)
{
	return _width;
}

int
CrossCorrelation::height(void)
{
	return _height;
}

double
CrossCorrelation::NCC(int x, int y, int window_width, int window_height)
{
	double std_deviation;
	double sum = 0.0;

	std_deviation = _img0.std_deviation(x, y, window_width, window_height) * _img1.std_deviation(x, y, window_width, window_height);
	for (int y = 0; y < window_height; y++) {
		for (int x = 0; x < window_width; x++) {
			sum += (_img0.image(x, y) - _img0.mean(x, y, window_width, window_height))
			    * (_img0.image(x, y) - _img0.mean(x, y, window_width, window_height));
		}
	}
	return sum / std_deviation;
}

double
CrossCorrelation::TruncatedNCC(int x, int y, int window_width, int window_height)
{
	return std::min(1.0, 1.0 - this->NCC(x, y, window_width, window_height));
}

