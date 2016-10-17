
# $Id: Makefile 482 2006-10-25 16:07:34Z camilo $
include ../../common.mk

#$ module load icc/15.0.0   
#$ module load impi/4.1.3.048   
#$ module load fftw/icc/impi/double/3.3.3

# make single simple rule for creating each target 
# to have all source files listed in one rule
# to make tar file out of them
#
# here specify the default target:

INCDIR=../../$(INCDIRNAME)
LIBDIR=../../$(LIBDIRNAME)
LIBSRCDIR=../../$(LIBSRCDIRNAME)

headers=$(addprefix $(INCDIR)/, $(headerfilenames))
libraryfiles=$(addprefix $(LIBSRCDIR)/, $(libraryfilenames))

# tell make where libnpsf.a is lying
libnpsf.a: $(LIBDIR)/libnpsf.a
# tell it how to build it, if it not there
$(LIBDIR)/libnpsf.a:
	$(MAKE) -C $(LIBDIR)

he:  Main.cpp libnpsf.a  
	mpicxx  $(CFLAGS) -I$(INCDIR) Main.cpp  -g -L$(LIBDIR) -lnpsf  -lfftw3  -o $@
	$(call gettar)
	$(call lib_msg, $@)
doc:
	$(DOXYGEN)

clean: 
	rm -f ../../libsrc/libnpsf.a
	rm -f ../../lib/libnpsf.a
	rm -f ../../libsrc/*.o
	rm -f ../../libsrc/*.d

.PHONY: all clean

