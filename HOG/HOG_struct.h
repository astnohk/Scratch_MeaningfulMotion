#include "../lib/ImgStruct.h"



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
		HOG(bool init_signed, int init_width, int init_height, int init_bins);
		HOG& copy(const HOG &copy);
		bool reset(bool init_signed, int init_width, int init_height, int init_bins);
		~HOG(void);
		void free(void);
		void setSign(bool init_signed);
		// Read
		bool Signed(void) const;
		int Bins(void) const;
		int Width(void) const;
		int Height(void) const;
		double Hist(int x, int y, int bin) const;
		const Histogram *Data(void) const;
		const Histogram *Data(int x, int y) const;
		// control Histogram
		bool AddHist(int x, int y, int bin, double val);
};

