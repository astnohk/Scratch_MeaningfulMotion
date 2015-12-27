#include "../ImgClass/ImgStatistics.h"
#include "../lib/ImgStruct.h"



#ifndef LIB_HOG_struct
#define LIB_HOG_struct

class HOG // Histograms of Oriented Gradients
{
	private:
		bool orient_signed;
		int bins;
		int width;
		int height;
		Histogram *hist;
	public:
		HOG(void);
		HOG(const HOG &copy);
		HOG(const bool init_signed, const int init_width, const int init_height, const int init_bins);
		HOG& copy(const HOG &copy);
		bool reset(const bool init_signed, const int init_width, const int init_height, const int init_bins);
		~HOG(void);
		void free(void);
		void setSign(const bool init_signed);
		// Read
		bool Signed(void) const;
		int Bins(void) const;
		int Width(void) const;
		int Height(void) const;
		double Hist(const int x, const int y, const int bin) const;
		const Histogram *Data(void) const;
		const Histogram *Data(const int x, const int y) const;
		// control Histogram
		bool AddHist(const int x, const int y, const int bin, const double val);
};

struct HOG_PARAM
{
	int Bins;
	bool Dense;
	bool SignedOrient;
	HOG_PARAM(void);
	void set_default(const char* name);
	void set_value(const char* name, const void* value);
};

#endif

