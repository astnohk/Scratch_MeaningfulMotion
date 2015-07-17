CC = g++
WARNING = -Wall -Wextra
LIBES = -lm -lX11
OPTION = -fopenmp
OPTION_STDCPP11 = -std=c++11

CFILES = main.cpp Class.cpp Struct.cpp Scratch_MeaningfulMotion.cpp Detection.cpp Exclusive.cpp Plotting.cpp Library.cpp ImgLibrary.cpp Affine_MultipleMotion.cpp OpticalFlow_MultipleMotion.cpp MultiResolution.cpp MEstimator.cpp HOG.cpp HOG_class.cpp HOG_match.cpp Plot_X11.cpp pnm.cpp pnm_double.cpp pnm_library.cpp
OFILES = main.o Class.o Struct.o Scratch_MeaningfulMotion.o Detection.o Exclusive.o Plotting.o Library.o ImgLibrary.o Affine_MultipleMotion.o OpticalFlow_MultipleMotion.o MultiResolution.o MEstimator.o HOG.o HOG_class.o HOG_match.o Plot_X11.o pnm.o pnm_double.o pnm_library.o

OUTNAME = Scratch_MeaningfulMotion


Scratch_MeaningfulMotion: $(OFILES)
	$(CC) $(WARNING) $(LIBES) $(OPTION) -O2 -o $(OUTNAME) $^

main.o: main.cpp
	$(CC) $(WARNING) $(OPTION) -c $^

Class.o: Class.cpp
	$(CC) $(WARNING) $(OPTION) -c $^

Struct.o: Struct.cpp
	$(CC) $(WARNING) $(OPTION) -c $^

Scratch_MeaningfulMotion.o: Scratch_MeaningfulMotion.cpp
	$(CC) $(WARNING) $(OPTION) -c $^

Detection.o: Detection.cpp
	$(CC) $(WARNING) $(OPTION) -c $^

Exclusive.o: Exclusive.cpp
	$(CC) $(WARNING) $(OPTION) -c $^

Plotting.o: Plotting.cpp
	$(CC) $(WARNING) $(OPTION) -c $^

Library.o: Library.cpp
	$(CC) $(WARNING) $(OPTION) -c $^

ImgLibrary.o: ImgLibrary.cpp
	$(CC) $(WARNING) $(OPTION) -c $^

Affine_MultipleMotion.o: Affine_MultipleMotion.cpp
	$(CC) $(WARNING) $(OPTION) -c $^

OpticalFlow_MultipleMotion.o: OpticalFlow_MultipleMotion.cpp
	$(CC) $(WARNING) $(OPTION) -c $^

MultiResolution.o: MultiResolution.cpp
	$(CC) $(WARNING) $(OPTION) -c $^

MEstimator.o: MEstimator.cpp
	$(CC) $(WARNING) $(OPTION) -c $^

HOG.o: HOG.cpp
	$(CC) $(WARNING) $(OPTION) -c $^

HOG_class.o: HOG_class.cpp
	$(CC) $(WARNING) $(OPTION) -c $^

HOG_match.o: HOG_match.cpp
	$(CC) $(WARNING) $(OPTION) -c $^

Plot_X11.o: Plot_X11.cpp
	$(CC) $(WARNING) $(OPTION) -c $^

pnm.o: pnm.cpp
	$(CC) $(WARNING) $(OPTION) -c $^

pnm_double.o: pnm_double.cpp
	$(CC) $(WARNING) $(OPTION) -c $^

pnm_library.o: pnm_library.cpp
	$(CC) $(WARNING) $(OPTION) -c $^


debug: $(CFILES)
	$(CC) $(WARNING) $(LIBES) -g -O2 -o $(OUTNAME) $^

debugmp: $(CFILES)
	$(CC) $(WARNING) $(LIBES) $(OPTION) -g -O2 -o $(OUTNAME) $^

clean:
	rm -f $(OFILES)

