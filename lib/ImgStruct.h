#include <new>

#ifndef nullptr
#define nullptr 0
#endif


#ifndef LIB_ImgStruct
#define LIB_ImgStruct
struct SIZE
{
	int width;
	int height;
	SIZE(void);
	SIZE(int w, int h);
	void reset(void);
	void set_size(const SIZE &size);
	void set_size(const SIZE *size);
};

struct COORDINATE
{
	int x;
	int y;
	COORDINATE(void);
	COORDINATE(int ix, int iy);
};

struct COORDINATE_3D
{
	double x;
	double y;
	double z;
	COORDINATE_3D(void);
	COORDINATE_3D(double ix, double iy, double iz);
	void set(double sx, double sy, double sz);
};

struct VECTOR_2D
{
	double x;
	double y;
	VECTOR_2D(void);
	VECTOR_2D(double ix, double iy);
	void reset(void);
};

class Histogram
{
	private:
		int _bins;
		double *_hist;
	public:
		Histogram(void);
		Histogram(const Histogram &copy);
		explicit Histogram(int init_bins);
		Histogram& copy(const Histogram &copy);
		Histogram& reset(int init_bins);
		~Histogram(void);
		void free(void);
		// Read
		const double* data(void) const;
		int bins(void) const;
		double get(int bin) const;
		// Control
		bool add(int bin, double val);
};
#endif

