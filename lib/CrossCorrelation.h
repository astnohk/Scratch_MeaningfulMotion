#include <algorithm>
#include <cmath>
#include <cstdio>
#include <new>

#ifndef nullptr
#define nullptr 0
#endif




#ifndef LIB_CrossCorrelation
#define LIB_CrossCorrelation
class ImgStatistics
{
	private:
		int _width;
		int _height;
		double *_data;
	public:
		ImgStatistics(void);
		ImgStatistics(const ImgStatistics &copy);
		ImgStatistics(int W, int H, double *Img);
		~ImgStatistics(void);
		void set(int W, int H, double *Img);
		ImgStatistics& copy(const ImgStatistics &copy);
		ImgStatistics& operator=(const ImgStatistics &copy);
		double& image(int x, int y);
		int width(void);
		int height(void);
		double mean();
		double mean(int start_x, int start_y, int end_x, int end_y);
		double variance();
		double variance(int x, int y, int window_width, int window_height);
		double std_deviation(void);
		double std_deviation(int x, int y, int window_width, int window_height);
};


class CrossCorrelation
{
	private:
		int _width;
		int _height;
		ImgStatistics _img0;
		ImgStatistics _img1;
	public:
		CrossCorrelation(void);
		CrossCorrelation(const CrossCorrelation &copy);
		CrossCorrelation(ImgStatistics &img0, ImgStatistics &img1);
		~CrossCorrelation(void);
		CrossCorrelation& copy(const CrossCorrelation &copy);
		CrossCorrelation& operator=(const CrossCorrelation &copy);
		int width(void);
		int height(void);
		double& operator[](int n);
		double NCC(int x, int y, int window_width, int window_height);
		double TruncatedNCC(int x, int y, int window_width, int window_height);
};
#endif

