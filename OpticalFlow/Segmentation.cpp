#include "Segmentation.h"




template <class T>
ImgVector<T> *
Segmentation_Merging(const ImgVector<T>& img)
{
	ImgVector<T>* img_segment = nullptr;

	try {
		img_segment = new ImgVector<int>(img.width(), img.height());
	}
	catch (const std::bad_alloc& bad) {
		fprintf(stderr, "ImgVector<T>* Segmentation_Merging(const ImgVector<T>&) : Cannot allocate memory\n");
		throw;
	}
	for (int i = 0; i < img.size(); i++) {
		(*img_segment)[i] = i;
	}
	for (int y = 0; y < img.height(); y++) {
		for (int x = 0; x < img.width(); x++) {
			img_segment
		}
	}
	return img_segment;
}

