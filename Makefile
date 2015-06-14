CC = g++
OPTION = -W -Wall -fopenmp -std=c++11
LOADLIBES = -lm -lX11 -fopenmp
OPTION_DEBUG = -W -Wall -std=c++11
LOADLIBES_DEBUG = -lm -lX11

CFILES = main.cpp Class.cpp Struct.cpp Scratch_MeaningfulMotion.cpp Detection.cpp Exclusive.cpp Plotting.cpp Library.cpp ImgLibrary.cpp Affine_MultipleMotion.cpp OpticalFlow_MultipleMotion.cpp MultiResolution.cpp MEstimator.cpp Plot_X11.cpp pnm.cpp pnm_double.cpp pnm_library.cpp
OFILES = main.o Class.o Struct.o Scratch_MeaningfulMotion.o Detection.o Exclusive.o Plotting.o Library.o ImgLibrary.o Affine_MultipleMotion.o OpticalFlow_MultipleMotion.o MultiResolution.o MEstimator.o Plot_X11.o pnm.o pnm_double.o pnm_library.o

OUTNAME = Scratch_MeaningfulMotion


Scratch_MeaningfulMotion: $(OFILES)
	$(CC) $(LOADLIBES) -O2 -o $@ $^

main.o: main.cpp
	$(CC) $(OPTION) -c $^

Class.o: Class.cpp
	$(CC) $(OPTION) -c $^

Struct.o: Struct.cpp
	$(CC) $(OPTION) -c $^

Scratch_MeaningfulMotion.o: Scratch_MeaningfulMotion.cpp
	$(CC) $(OPTION) -c $^

Detection.o: Detection.cpp
	$(CC) $(OPTION) -c $^

Exclusive.o: Exclusive.cpp
	$(CC) $(OPTION) -c $^

Plotting.o: Plotting.cpp
	$(CC) $(OPTION) -c $^

Library.o: Library.cpp
	$(CC) $(OPTION) -c $^

ImgLibrary.o: ImgLibrary.cpp
	$(CC) $(OPTION) -c $^

Affine_MultipleMotion.o: Affine_MultipleMotion.cpp
	$(CC) $(OPTION) -c $^

OpticalFlow_MultipleMotion.o: OpticalFlow_MultipleMotion.cpp
	$(CC) $(OPTION) -c $^

MultiResolution.o: MultiResolution.cpp
	$(CC) $(OPTION) -c $^

MEstimator.o: MEstimator.cpp
	$(CC) $(OPTION) -c $^

Plot_X11.o: Plot_X11.cpp
	$(CC) $(OPTION) -c $^

pnm.o: pnm.cpp
	$(CC) $(OPTION) -c $^

pnm_double.o: pnm_double.cpp
	$(CC) $(OPTION) -c $^

pnm_library.o: pnm_library.cpp
	$(CC) $(OPTION) -c $^

debug: $(CFILES)
	$(CC) -g $(OPTION_DEBUG) $(LOADLIBES_DEBUG) -O2 -o $(OUTNAME) $^

debugmp: $(CFILES)
	$(CC) -g $(OPTION_DEBUG) $(LOADLIBES) -O2 -o $(OUTNAME) $^

clean:
	rm -f $(OFILES)

