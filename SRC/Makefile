#######################################################################
#
#  This makefile creates a library for distributed SuperLU.
#  The files are organized as follows:
#
#       ALLAUX  -- Auxiliary routines called from all precisions
#       DSLUSRC -- Double precision real serial SuperLU routines
#       DPLUSRC -- Double precision real parallel SuperLU routines
#       ZSLUSRC -- Double precision complex serial SuperLU routines
#       ZPLUSRC -- Double precision complex parallel SuperLU routines
#
#  The library can be set up to include routines for any combination
#  of the two precisions.  To create or add to the library, enter make
#  followed by one or more of the precisions desired.  Some examples:
#       make double
#       make double complex16
#  Alternatively, the command
#       make
#  without any arguments creates a library of all two precisions.
#  The library is called
#       superlu.a
#  and is created at the next higher directory level.
#
#  To remove the object files after the library is created, enter
#       make clean
#
#######################################################################
include ../make.inc

DSUPERLULIB = libsuperlu_dist.a

INCLUDEDIR += -I$(SUPERLU_ROOT)/SRC

#
# Precision independent routines
#
ALLAUX 	= sp_ienv.c etree.c sp_colorder.c get_perm_c.c \
	  colamd.c mmd.c comm.c memory.c util.c superlu_grid.c \
	  pxerr_dist.c superlu_timer.c symbfact.c \
	  psymbfact.c psymbfact_util.c get_perm_c_parmetis.c mc64ad_dist.c \
	  xerr_dist.c smach_dist.c dmach_dist.c \
	  superlu_dist_version.c comm_tree.c

ifeq ($(HAVE_CUDA),TRUE)
ALLAUX += gpublas_utils.c pdgstrs_lsum_cuda.c pdgstrs_lsum_cuda_dlink.c
endif

#
# Routines literally taken from SuperLU, but renamed with suffix _dist
#
DSLUSRC	= dlangs_dist.c dgsequ_dist.c dlaqgs_dist.c dutil_dist.c \
	  dmemory_dist.c dmyblas2_dist.c dsp_blas2_dist.c dsp_blas3_dist.c
ZSLUSRC	= dcomplex_dist.c zlangs_dist.c zgsequ_dist.c zlaqgs_dist.c \
	  zutil_dist.c zmemory_dist.c zmyblas2_dist.c \
	  zsp_blas2_dist.c zsp_blas3_dist.c

#
# Routines for double precision parallel SuperLU
DPLUSRC = pdgssvx.c pdgssvx_ABglobal.c \
	  dreadhb.c dreadrb.c dreadtriple.c dreadMM.c dbinary_io.c \
	  pdgsequ.c pdlaqgs.c dldperm_dist.c pdlangs.c pdutil.c \
	  pdsymbfact_distdata.c ddistribute.c pddistribute.c \
	  pdgstrf.c dstatic_schedule.c pdgstrf2.c pdGetDiagU.c \
	  pdgstrs.c pdgstrs1.c pdgstrs_lsum.c pdgstrs_Bglobal.c \
	  pdgsrfs.c pdgsmv.c pdgsrfs_ABXglobal.c pdgsmv_AXglobal.c \
	  dreadtriple_noheader.c
#
# Routines for double complex parallel SuperLU
ZPLUSRC = pzgssvx.c pzgssvx_ABglobal.c \
	  zreadhb.c zreadrb.c zreadtriple.c zreadMM.c zbinary_io.c\
	  pzgsequ.c pzlaqgs.c zldperm_dist.c pzlangs.c pzutil.c \
	  pzsymbfact_distdata.c zdistribute.c pzdistribute.c \
	  pzgstrf.c zstatic_schedule.c pzgstrf2.c pzGetDiagU.c \
	  pzgstrs.c pzgstrs1.c pzgstrs_lsum.c pzgstrs_Bglobal.c \
	  pzgsrfs.c pzgsmv.c pzgsrfs_ABXglobal.c pzgsmv_AXglobal.c \
	  zreadtriple_noheader.c

# ifeq ($(HAVE_COMBBLAS),TRUE)
# DPLUSRC += d_c2cpp_GetHWPM.cpp
# ZPLUSRC += z_c2cpp_GetHWPM.cpp
# endif

PZGSTRF_SRC = pzgstrf.c zSchCompUdt-2Ddynamic.c zlook_ahead_update.c
PDGSTRF_SRC = pdgstrf.c dlook_ahead_update.c dSchCompUdt-2Ddynamic.c
			  

ALLAUXOBJ=$(patsubst %.c, %.o, $(ALLAUX))
DSLUOBJ=$(patsubst %.c, %.o, $(DSLUSRC))
ZSLUOBJ=$(patsubst %.c, %.o, $(ZSLUSRC))
DPLUOBJ=$(patsubst %.c, %.o, $(DPLUSRC))
ZPLUOBJ=$(patsubst %.c, %.o, $(ZPLUSRC))

ALLOBJ = ${ALLAUXOBJ} ${DSLUOBJ} ${ZSLUOBJ} ${DPLUOBJ} ${ZPLUOBJ}

ifeq ($(HAVE_SW), TRUE)
	HOST_DPLUDIR = sw/host/double
	HOST_ZPLUDIR = sw/host/doublecomplex
	SLAVE_DPLUDIR = sw/slave/double
	SLAVE_ZPLUDIR = sw/slave/doublecomplex

	HOST_DPLUSRC = $(shell find $(HOST_DPLUDIR) -name "*.c")
	HOST_ZPLUSRC = $(shell find $(HOST_ZPLUDIR) -name "*.c")
	SLAVE_DPLUSRC = $(shell find $(SLAVE_DPLUDIR) -name "*.c")
	SLAVE_ZPLUSRC = $(shell find $(SLAVE_ZPLUDIR) -name "*.c")

	HOST_DPLUOBJ = $(patsubst %.c, %.o, $(HOST_DPLUSRC))
	HOST_ZPLUOBJ = $(patsubst %.c, %.o, $(HOST_ZPLUSRC))
	SLAVE_DPLUOBJ = $(patsubst %.c, %.o, $(SLAVE_DPLUSRC))
	SLAVE_ZPLUOBJ = $(patsubst %.c, %.o, $(SLAVE_ZPLUSRC))

	HOST_OBJ = ${HOST_DPLUOBJ} ${HOST_ZPLUOBJ}
	SLAVE_OBJ = ${SLAVE_DPLUOBJ} ${SLAVE_ZPLUOBJ}

	PZGSTRF_SRC += sw/host/zlook_ahead_update_sw.c sw/host/zSchCompUdt-2Ddynamic-sw.c 
	PDGSTRF_SRC += 
endif

all:  double complex16

config_h:
	rm -f superlu_dist_config.h
ifeq ($(XSDK_INDEX_SIZE),64)
		printf "#define XSDK_INDEX_SIZE 64\n" >> superlu_dist_config.h
else
		printf "/* #define XSDK_INDEX_SIZE 64 */\n" >> superlu_dist_config.h
endif
ifeq ($(SLU_HAVE_LAPACK),TRUE)
		printf "#define SLU_HAVE_LAPACK TRUE\n" >> superlu_dist_config.h
else
		printf "/* #define SLU_HAVE_LAPACK TRUE */\n" >> superlu_dist_config.h
endif
ifeq ($(HAVE_PARMETIS),TRUE)
		printf "#define HAVE_PARMETIS TRUE\n" >> superlu_dist_config.h
else
		printf "/* #define HAVE_PARMETIS TRUE */\n" >> superlu_dist_config.h
endif
ifeq ($(HAVE_COMBBLAS),TRUE)
		printf "#define HAVE_COMBBLAS TRUE\n" >> superlu_dist_config.h
else
		printf "/* #define HAVE_COMBBLAS TRUE */\n" >> superlu_dist_config.h
endif
	printf "#if (XSDK_INDEX_SIZE == 64)\n#define _LONGINT 1\n#endif\n" >> superlu_dist_config.h


#config_h:
#	cp -f superlu_dist_config_for_makefile.h superlu_dist_config.h

double: config_h $(DSLUOBJ) $(DPLUOBJ) $(ALLAUXOBJ) ${HOST_DPLUOBJ} ${SLAVE_DPLUOBJ}
	$(ARCH) $(ARCHFLAGS) $(DSUPERLULIB) \
		$(DSLUOBJ) $(DPLUOBJ) $(ALLAUXOBJ)  ${HOST_DPLUOBJ} ${SLAVE_DPLUOBJ}
	$(RANLIB) $(DSUPERLULIB)

complex16: config_h $(ZSLUOBJ) $(ZPLUOBJ) $(ALLAUXOBJ) ${HOST_ZPLUOBJ} ${SLAVE_ZPLUOBJ}
	$(ARCH) $(ARCHFLAGS) $(DSUPERLULIB) \
		$(ZSLUOBJ) $(ZPLUOBJ) $(ALLAUXOBJ) ${HOST_ZPLUOBJ} ${SLAVE_ZPLUOBJ}
	$(RANLIB) $(DSUPERLULIB)

pdgstrf.o: $(PDGSTRF_SRC)
	$(CC) $(CFLAGS) $(CDEFS) $(BLASDEF) $(INCLUDEDIR) -c pdgstrf.c $(VERBOSE)

pzgstrf.o: $(PZGSTRF_SRC)
	$(CC) $(CFLAGS) $(CDEFS) $(BLASDEF) $(INCLUDEDIR) -c pzgstrf.c $(VERBOSE)

# pdgstrs_lsum_cuda.o: pdgstrs_lsum_cuda.cu
# 	$(NVCC) $(CUDACFLAGS) -dc $(CDEFS) $(BLASDEF) $(INCLUDEDIR) -c pdgstrs_lsum_cuda.cu $(VERBOSE)

# pdgstrs_lsum_cuda_dlink.o:	
# 	$(NVCC) $(CUDACFLAGS) -dlink $(CUDALIBS) -o pdgstrs_lsum_cuda_dlink.o pdgstrs_lsum_cuda.o $(VERBOSE)	

${ALLOBJ} : %.o : %.c
	$(CC) $(CFLAGS) $(SW_HOST_FLAGS) $(CDEFS) $(BLASDEF) $(INCLUDEDIR) -c $< $(VERBOSE) -o $@

# %.o : %.cpp
# 	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(CDEFS) $(BLASDEF) $(INCLUDEDIR) -c $< $(VERBOSE)

# %.o : %.f
# 	$(FORTRAN) $(FFLAGS) -c $< $(VERBOSE)

${HOST_OBJ} : %.o : %.c
	$(CC) $(CFLAGS) $(SW_HOST_FLAGS) $(CDEFS) $(BLASDEF) $(INCLUDEDIR) -c $< $(VERBOSE) -o $@

${SLAVE_OBJ} : %.o : %.c
	$(CC) $(CFLAGS) $(SW_SLAVE_FLAGS) $(CDEFS) $(BLASDEF) $(INCLUDEDIR) -c $< $(VERBOSE) -o $@

header_file=$(shell find . -name "*.h")

install: 
	mkdir -p $(INCLUDE_DIR)
	mkdir -p $(LIB_DIR)
	cp --path $(header_file) $(INCLUDE_DIR)
	cp $(DSUPERLULIB) $(LIB_DIR)

ALLOBJ = $(shell find . -name "*.o")

clean:	
	rm -f $(ALLOBJ) $(DSUPERLULIB) superlu_dist_config.h
