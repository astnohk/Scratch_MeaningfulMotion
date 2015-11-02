CC = g++
WARNING = -Wall -Wextra
LIBES = -lm -lX11
OPTION = -fopenmp
#OPTION = -fopenmp -std=c++11


LIBRARY_CFILES = lib/Class.cpp lib/ImgLibrary.cpp lib/ImgStruct.cpp lib/Library.cpp lib/Struct.cpp lib/Vector.cpp
IMGCLASS_CFILES = ImgClass/CrossCorrelation.cpp ImgClass/ImgStatistics.cpp
MEANINGFUL_CFILES = MeaningfulAlignments/Detection.cpp MeaningfulAlignments/Exclusive.cpp
OPTICALFLOW_CFILES = OpticalFlow/MultiResolution.cpp OpticalFlow/MEstimator.cpp OpticalFlow/Affine_MultipleMotion.cpp OpticalFlow/OpticalFlow_MultipleMotion.cpp
MOTIONCOMPENSATION_CFILES = MotionCompensation/MotionCompensation.cpp
HOG_CFILES = HOG/HOG.cpp HOG/HOG_struct.cpp HOG/HOG_match.cpp
PLOT_CFILES = Plot/Plotting.cpp Plot/Plot_X11.cpp Plot/Plot_X11_Struct.cpp
PNM_CFILES = PNM/pnm.cpp PNM/pnm_double.cpp PNM/pnm_library.cpp

CFILES = main.cpp Scratch_Struct.cpp Scratch_MeaningfulMotion.cpp $(LIBRARY_CFILES) $(IMGCLASS_CFILES) $(MOTIONCOMPENSATION_CFILES) $(MEANINGFUL_CFILES) $(OPTICALFLOW_CFILES) $(HOG_CFILES) $(PLOT_CFILES) $(PNM_CFILES)


LIBRARY_OFILES = Class.o ImgLibrary.o ImgStruct.o Library.o Struct.o Vector.o
IMGCLASS_OFILES = CrossCorrelation.o ImgStatistics.o
MEANINGFUL_OFILES = Detection.o Exclusive.o
OPTICALFLOW_OFILES = MultiResolution.o MEstimator.o Affine_MultipleMotion.o OpticalFlow_MultipleMotion.o
MOTIONCOMPENSATION_OFILES = MotionCompensation.o
HOG_OFILES = HOG.o HOG_struct.o HOG_match.o
PLOT_OFILES = Plotting.o Plot_X11.o Plot_X11_Struct.o
PNM_OFILES = pnm.o pnm_double.o pnm_library.o

OFILES = main.o Scratch_Struct.o Scratch_MeaningfulMotion.o $(LIBRARY_OFILES) $(IMGCLASS_OFILES) $(MEANINGFUL_OFILES) $(OPTICALFLOW_OFILES) $(MOTIONCOMPENSATION_OFILES) $(HOG_OFILES) $(PLOT_OFILES) $(PNM_OFILES)

OUTNAME = Scratch_MeaningfulMotion


Scratch_MeaningfulMotion: $(OFILES)
	$(CC) $(WARNING) $(LIBES) $(OPTION) -O2 -o $(OUTNAME) $^

Scratch_MeaningfulMotion.h: lib/ImgClass.h lib/ImgStatistics.h lib/ImgStruct.h lib/Vector.h PNM/pnm.h

main.o: main.cpp
	$(CC) $(WARNING) $(OPTION) -c $^
main.cpp: Scratch_MeaningfulMotion.h

Class.o: lib/Class.cpp
	$(CC) $(WARNING) $(OPTION) -c $^
Class.cpp: lib/Class.h

Struct.o: lib/Struct.cpp
	$(CC) $(WARNING) $(OPTION) -c $^
Struct.cpp: lib/Struct.h lib/Class.h

CrossCorrelation.o: ImgClass/CrossCorrelation.cpp
	$(CC) $(WARNING) $(OPTION) -c $^
CrossCorrelation.cpp: ImgClass/CrossCorrelation.h

ImgLibrary.o: lib/ImgLibrary.cpp
	$(CC) $(WARNING) $(OPTION) -c $^
ImgLibrary.cpp: Scratch_MeaningfulMotion.h

ImgStatistics.o: ImgClass/ImgStatistics.cpp
	$(CC) $(WARNING) $(OPTION) -c $^
ImgStatistics.cpp: ImgClass/ImgStatistics.h

ImgStruct.o: lib/ImgStruct.cpp
	$(CC) $(WARNING) $(OPTION) -c $^
ImgStruct.cpp: lib/ImgStruct.h

Library.o: lib/Library.cpp
	$(CC) $(WARNING) $(OPTION) -c $^
Library.cpp: Scratch_MeaningfulMotion.h

Vector.o: lib/Vector.cpp
	$(CC) $(WARNING) $(OPTION) -c $^
Vector.cpp: lib/Vector.h

Scratch_Struct.o: Scratch_Struct.cpp
	$(CC) $(WARNING) $(OPTION) -c $^
Scratch_Struct.cpp: Scratch_MeaningfulMotion.h

Scratch_MeaningfulMotion.o: Scratch_MeaningfulMotion.cpp
	$(CC) $(WARNING) $(OPTION) -c $^
Scratch_MeaningfulMotion.cpp: Scratch_MeaningfulMotion.h OpticalFlow/Affine_MultipleMotion.h OpticalFlow/OpticalFlow_MultipleMotion.h HOG/HOG.h

Detection.o: MeaningfulAlignments/Detection.cpp
	$(CC) $(WARNING) $(OPTION) -c $^
Detection.cpp: Scratch_MeaningfulMotion.h

Exclusive.o: MeaningfulAlignments/Exclusive.cpp
	$(CC) $(WARNING) $(OPTION) -c $^
Exclusive.cpp: Scratch_MeaningfulMotion.h

OpticalFlow/MultiResolution.h: ImgClass/ImgClass.h

MultiResolution.o: OpticalFlow/MultiResolution.cpp
	$(CC) $(WARNING) $(OPTION) -c $^
MultiResolution.cpp: OpticalFlow/MultiResolution.h Scratch_MeaningfulMotion.h

MultiResolution.h: lib/Vector.h

MEstimator.o: OpticalFlow/MEstimator.cpp
	$(CC) $(WARNING) $(OPTION) -c $^
MEstimator.cpp: MEstimator.h

MotionCompensation.o: MotionCompensation/MotionCompensation.cpp
	$(CC) $(WARNING) $(OPTION) -c $^
MotionCompensation.cpp: MotionCompensation/MotionCompensation.h

Affine_MultipleMotion.o: OpticalFlow/Affine_MultipleMotion.cpp
	$(CC) $(WARNING) $(OPTION) -c $^
Affine_MultipleMotion.cpp: Scratch_MeaningfulMotion.h OpticalFlow/Affine_MultipleMotion.h

OpticalFlow/Affine_MultipleMotion.h: OpticalFlow/MEstimator.h OpticalFlow/MultiResolution.h lib/Struct.h

OpticalFlow_MultipleMotion.o: OpticalFlow/OpticalFlow_MultipleMotion.cpp
	$(CC) $(WARNING) $(OPTION) -c $^
OpticalFlow_MultipleMotion.cpp: Scratch_MeaningfulMotion.h OpticalFlow/OpticalFlow_MultipleMotion.h

HOG.o: HOG/HOG.cpp
	$(CC) $(WARNING) $(OPTION) -c $^
HOG.cpp: HOG/HOG.h Scratch_MeaningfulMotion.h

HOG.h: HOG/HOG_struct.h lib/Vector.h PNM/pnm.h

HOG_match.o: HOG/HOG_match.cpp
	$(CC) $(WARNING) $(OPTION) -c $^
HOG_match.cpp: HOG/HOG.h MotionCompensation/MotionCompensation.h Scratch_MeaningfulMotion.h

HOG_struct.h: ImgClass/ImgStatistics.h lib/ImgStruct.h

HOG_struct.o: HOG/HOG_struct.cpp
	$(CC) $(WARNING) $(OPTION) -c $^
HOG_struct.cpp: HOG/HOG_struct.h

Plotting.o: Plot/Plotting.cpp
	$(CC) $(WARNING) $(OPTION) -c $^
Plotting.cpp: Scratch_MeaningfulMotion.h

Plot_X11.o: Plot/Plot_X11.cpp
	$(CC) $(WARNING) $(OPTION) -c $^
Plot_X11.cpp: Plot/Plot_X11.h Scratch_MeaningfulMotion.h

Plot_X11_Struct.o: Plot/Plot_X11_Struct.cpp
	$(CC) $(WARNING) $(OPTION) -c $^
Plot_X11_Struct.cpp: Plot/Plot_X11.h Scratch_MeaningfulMotion.h

pnm.o: PNM/pnm.cpp
	$(CC) $(WARNING) $(OPTION) -c $^
pnm.cpp: PNM/pnm.h

pnm_double.o: PNM/pnm_double.cpp
	$(CC) $(WARNING) $(OPTION) -c $^
pnm_double.cpp: PNM/pnm.h

pnm_library.o: PNM/pnm_library.cpp
	$(CC) $(WARNING) $(OPTION) -c $^
pnm_library.cpp: PNM/pnm.h




debug: $(CFILES)
	$(CC) $(WARNING) $(LIBES) -g -O2 -o $(OUTNAME) $^

debugmp: $(CFILES)
	$(CC) $(WARNING) $(LIBES) $(OPTION) -g -O2 -o $(OUTNAME) $^

clean:
	rm -f $(OFILES)

