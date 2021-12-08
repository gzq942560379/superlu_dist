############################################################################
#
#  Program:         SuperLU_DIST
#
#  Module:          Makefile
#
#  Purpose:         Top-level Makefile
#
#  Creation date:   September 1, 1999  version 1.0
#
#  Modified:        
#
############################################################################

include make.inc

all: lib install example

lib: superlulib

example: install
	$(MAKE) -C EXAMPLE

clean: cleanlib cleantesting

# install:
# 	( cd INSTALL; $(MAKE) )
# #	( cd INSTALL; cp lsame.c ../SRC/; \
# #	  cp dlamch.c ../SRC/; cp slamch.c ../SRC/ )

install : lib
	$(MAKE) -C SRC install

# blaslib:
# 	( cd CBLAS; $(MAKE) )

superlulib:
	$(MAKE) -C SRC

cleanlib:
	$(MAKE) -C SRC clean
	rm -rf include lib

cleantesting:
	$(MAKE) -C EXAMPLE clean
	$(MAKE) -C FORTRAN clean
	$(MAKE) -C TEST clean
