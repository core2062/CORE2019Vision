CXX=arm-raspbian10-linux-gnueabihf-g++
DEPS_CFLAGS=-Iinclude -Iinclude/opencv -Iinclude
DEPS_LIBS=-Llib -lwpilibc -lwpiHal -lcameraserver -lntcore -lcscore -lopencv_dnn -lopencv_ml -lopencv_objdetect -lopencv_shape -lopencv_stitching -lopencv_superres -lopencv_videostab -lopencv_calib3d -lopencv_features2d -lopencv_highgui -lopencv_videoio -lopencv_imgcodecs -lopencv_video -lopencv_photo -lopencv_imgproc -lopencv_flann -lopencv_core -lwpiutil
EXE=multiCameraServerExample
DESTDIR?=/home/pi/

.PHONY: clean build install

build: ${EXE}

install: build
	cp ${EXE} runCamera ${DESTDIR}

clean:
	rm ${EXE} *.o

OBJS=GripPipeline.o main.o 

GripPipeline.o: GripPipeline.cpp
	${CXX} -pthread -g -Og -c -o GripPipeline.o ${CXXFLAGS} ${DEPS_CFLAGS} GripPipeline.cpp

main.o: GripPipeline.h main.cpp
	${CXX} -pthread -g -Og -c -o main.o ${CXXFLAGS} ${DEPS_CFLAGS} main.cpp

${EXE}: ${OBJS}
	${CXX} -pthread -g -o $@ $^ ${DEPS_LIBS} -Wl,--unresolved-symbols=ignore-in-shared-libs

