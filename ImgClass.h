template <typename T>
class ImgVector
{
	private:
		T *data;
		int width;
		int height;
	public:
		ImgVector(void);
		ImgVector(int W, int H, T *array);
		~ImgVector(void);
		void reset(int W, int H, T *array);
		void set(int x, int y, T value);
		T get(int x, int y) const;
		T operator[](int n) const;
};

