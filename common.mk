SHELL=/bin/bash

headerfilenames=constant.h utils.h ConfigFileParser.h Cartesian4D.h wavefunction.h ParameterMap.h Observable.h Laser.h InputStuff.h fftw3.h mat.h matrix.h
libraryfilenames=Cartesian4D.cpp Parameter.cpp InputStuff.cpp ConfigFileParser.cpp Laser.cpp 
commonfilenames=$(headerfilenames) $(libraryfilenames)
docfilenames=README INSTALL


MATLABDIR=/opt/matlab
MATCPPOPTS=-I$(MATLABDIR)/extern/include -L$(MATLABDIR)/bin/glnxa64
MATLIBS=-lmat -lmx -lut -licuio -licuuc -licudata -licui18n -lMTwister
#MATLIBS=-lmat -lmx -lut
OPTFLAGS=-fopenmp -lfftw3_mpi  -lfftw3 -lm  -O3 -march=k8 -m64 -msse -msse2 -fomit-frame-pointer -mfpmath=sse -funroll-loops
#OPTFLAGS=-DDEBUG -march=k8 -m64 -msse -msse2 -mfpmath=sse -O

LIBSRCDIRNAME=libsrc
LIBDIRNAME=lib
INCDIRNAME=include

CFLAGS=-fno-common $(OPTFLAGS)

define lib_msg
	@echo ' &&&&&&&&&&&&&&&&&&&&&&&&&&&  '
	@echo ' We all are attosecond fans !'
	@echo '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& '
endef

