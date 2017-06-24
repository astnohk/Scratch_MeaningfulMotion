#CC=g++
CC=g++ -std=c++11
WARNING=-Wall -Wextra
LIBES=-lm -lX11
OPTION=-O2 -fopenmp
MACROS=


LIBRARY_CFILES = lib/Class.cpp lib/ImgLibrary.cpp lib/ImgStruct.cpp lib/Library.cpp lib/Struct.cpp lib/ExtVector.cpp
IMGCLASS_CFILES = ImgClass/BlockMatching.cpp ImgClass/CrossCorrelation.cpp ImgClass/ImgStatistics.cpp ImgClass/Lab.cpp ImgClass/RGB.cpp ImgClass/Segmentation.cpp
MEANINGFUL_CFILES = MeaningfulAlignments/Detection.cpp MeaningfulAlignments/Exclusive.cpp
OPTICALFLOW_CFILES = OpticalFlow/MultiResolution.cpp OpticalFlow/MEstimator.cpp OpticalFlow/Affine_MultipleMotion.cpp OpticalFlow/OpticalFlow.cpp OpticalFlow/OpticalFlow_BlockMatching.cpp OpticalFlow/Affine_BlockMatching.cpp
HOG_CFILES = HOG/HOG.cpp HOG/HOG_struct.cpp HOG/HOG_match.cpp
PLOT_CFILES = Plot/Plotting.cpp Plot/Plot_X11.cpp Plot/Plot_X11_Struct.cpp
PNM_CFILES = pnm_cpp_lib/pnm.cpp pnm_lib_cpp/pnm_double.cpp pnm_lib_cpp/pnm_library.cpp

CFILES = main.cpp Scratch_Struct.cpp Scratch_MeaningfulMotion.cpp $(LIBRARY_CFILES) $(IMGCLASS_CFILES) $(MEANINGFUL_CFILES) $(OPTICALFLOW_CFILES) $(HOG_CFILES) $(PLOT_CFILES) $(PNM_CFILES)


LIBRARY_OFILES = Class.o ImgLibrary.o ImgStruct.o Library.o Struct.o ExtVector.o
IMGCLASS_OFILES = BlockMatching.o CrossCorrelation.o ImgStatistics.o Lab.o RGB.o Segmentation.o
MEANINGFUL_OFILES = Detection.o Exclusive.o
OPTICALFLOW_OFILES = MultiResolution.o MEstimator.o Affine_MultipleMotion.o OpticalFlow.o OpticalFlow_BlockMatching.o Affine_BlockMatching.o
HOG_OFILES = HOG.o HOG_struct.o HOG_match.o
PLOT_OFILES = Plotting.o Plot_X11.o Plot_X11_Struct.o
PNM_OFILES = pnm.o pnm_double.o pnm_library.o

OFILES = main.o Scratch_Struct.o Scratch_MeaningfulMotion.o $(LIBRARY_OFILES) $(IMGCLASS_OFILES) $(MEANINGFUL_OFILES) $(OPTICALFLOW_OFILES) $(HOG_OFILES) $(PLOT_OFILES) $(PNM_OFILES)

OUTNAME = Scratch_MeaningfulMotion




Scratch_MeaningfulMotion: $(OFILES)
	$(CC) $(WARNING) $(LIBES) $(OPTION) -o $(OUTNAME) $^


main.o: main.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^

Class.o: lib/Class.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^

Struct.o: lib/Struct.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^

BlockMatching.o: ImgClass/BlockMatching.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^

CrossCorrelation.o: ImgClass/CrossCorrelation.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^

ImgLibrary.o: lib/ImgLibrary.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^

ImgStatistics.o: ImgClass/ImgStatistics.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^

Lab.o: ImgClass/Lab.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^

RGB.o: ImgClass/RGB.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^

Segmentation.o: ImgClass/Segmentation.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^

ImgStruct.o: lib/ImgStruct.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^

Library.o: lib/Library.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^

ExtVector.o: lib/ExtVector.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^

Scratch_Struct.o: Scratch_Struct.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^

Scratch_MeaningfulMotion.o: Scratch_MeaningfulMotion.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^

Detection.o: MeaningfulAlignments/Detection.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^

Exclusive.o: MeaningfulAlignments/Exclusive.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^


MultiResolution.o: OpticalFlow/MultiResolution.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^

MEstimator.o: OpticalFlow/MEstimator.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^

Affine_MultipleMotion.o: OpticalFlow/Affine_MultipleMotion.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^

Affine_BlockMatching.o: OpticalFlow/Affine_BlockMatching.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^

OpticalFlow.o: OpticalFlow/OpticalFlow.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^

OpticalFlow_BlockMatching.o: OpticalFlow/OpticalFlow_BlockMatching.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^

HOG.o: HOG/HOG.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^

HOG_match.o: HOG/HOG_match.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^

HOG_struct.o: HOG/HOG_struct.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^

Plotting.o: Plot/Plotting.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^

Plot_X11.o: Plot/Plot_X11.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^

Plot_X11_Struct.o: Plot/Plot_X11_Struct.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^

pnm.o: pnm_lib_cpp/pnm.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^

pnm_double.o: pnm_lib_cpp/pnm_double.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^

pnm_library.o: pnm_lib_cpp/pnm_library.cpp
	$(CC) $(WARNING) $(OPTION) $(MACROS) -c $^




debug: $(CFILES)
	$(CC) $(WARNING) $(LIBES) -g -O2 $(MACROS) -o $(OUTNAME) $^

debugmp: $(CFILES)
	$(CC) $(WARNING) $(LIBES) $(OPTION) $(MACROS) -g -O2 -o $(OUTNAME) $^

clean:
	rm -f $(OFILES)
	find -name "*.gch" -exec rm {} +

