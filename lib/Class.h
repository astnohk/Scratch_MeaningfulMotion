#include <string>

/*
#ifndef nullptr
#define nullptr 0
#endif
*/


#ifndef LIB_Class
#define LIB_Class

class ERROR
{
	private:
		std::string FunctionName;
		std::string ErrorFunctionName;
		std::string ValueName;
		std::string FileName;
	public:
		// Error data set
		ERROR(void);
		explicit ERROR(const char *name);
		const char* OutputFunctionName(void);
		void Function(const char *name);
		void Value(const char *name);
		void File(const char *name);
		// Error output
		void Others(const char *error);
		void OthersWarning(const char *error);
		void Malloc(void);
		void FunctionFail(void);
		void PointerNull(void);
		void ValueIncorrect(void);
		void ImageSize(void);
		void FileRead(void);
		void FileWrite(void);
};

class ATAN2_DIV_PI
{
	private:
		int width;
		int height;
		double *table;
	public:
		ATAN2_DIV_PI(void);
		ATAN2_DIV_PI(const ATAN2_DIV_PI &copy);
		ATAN2_DIV_PI(int W, int H);
		~ATAN2_DIV_PI(void);

		int Width(void) const;
		int Height(void) const;
		const double* Data(void) const;
		void reset(int W, int H);
		double val(int y, int x) const;
};

#endif

