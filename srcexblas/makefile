## ============================================================
## INTEL COMPILERS
## ============================================================
#
#CC = mpicc
#CFLAGS = -Wall -g -openmp -I.  -I${MKLROOT}/include
#CFLAGS = -Wall -openmp -I. -I${MKLROOT}/include
#CLFLAGS = -Wall -fopenmp -I. -I${MKLROOT}/include
#
#CLINKER = mpicc
#FLINKER = mpif77
#LDFLAGS = -openmp
#LIBLIST = -L. -lhbio -lclock -lsparsenew -lvector -lm -lc
#LIBLIST = -L. -lhbio -lclock -lvector -lm -lc
#LIBLIST = -L. -lhbio -lclock -lm -lc
#LIBLIST = -L. -lsparse -lvector -lclock -lm -lc
#
##LIBMKL = -L$(MKLROOT)/lib/intel64 $(MKL_FMULTIS_INTEL)
#LIBMKL = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread

# ============================================================
# GNU COMPILERS
# ============================================================

CC = mpicxx
CFLAGS = -std=c++11 -mavx -fabi-version=0 -Wall -fopenmp -g -I. -I${MKLROOT}/include -I${HOME}/libs
CLFLAGS = -Wall -fopenmp -I. -I${MKLROOT}/include

CLINKER = mpicxx
LDFLAGS = -fopenmp
LIBLIST = -L. -lhbio -lclock -lsparsenew -lvector -lm -lc
LIBLIST = -L. -lhbio -lclock -lvector -lm -lc
LIBLIST = -L. -lhbio -lclock -lm -lc
LIBLIST = -L. -lsparse -lvector -lclock -lm -lc

LIBMKL = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread

# ============================================================

AR = ar
ARFLAGS = ru

RL = ranlib

# ============================================================

OBJS_CLOCK  = reloj.o 
OBJS_VECTOR = ScalarVectors.o
OBJS_SPARSE = hb_io.o SparseProduct.o

OBJS = $(OBJS_CLOCK) $(OBJS_VECTOR) $(OBJS_SPARSE) 

# ============================================================

default: libclock.a libvector.a libsparse.a CG_MPI

libshared.a : $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RL) $(RLFLAGS) $@

libclock.a : $(OBJS_CLOCK)
	$(AR) $(ARFLAGS) $@ $?
	$(RL) $(RLFLAGS) $@

libvector.a : $(OBJS_VECTOR)
	$(AR) $(ARFLAGS) $@ $?
	$(RL) $(RLFLAGS) $@

libsparse.a : $(OBJS_SPARSE)
	$(AR) $(ARFLAGS) $@ $?
	$(RL) $(RLFLAGS) $@

CG_MPI: CG_MPI.o ToolsMPI.o matrix.o 
	$(CLINKER) $(LDFLAGS) -o CG_MPI CG_MPI.o ToolsMPI.o matrix.o $(LIBMKL) $(LIBLIST)

# ============================================================

.c.o:
	echo compiling
	$(CC) $(CFLAGS) -c $*.c

clean:
	rm -f *.o *.a CG_MPI

# ============================================================
