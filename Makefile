############################################################################
#
#  Program:         ScaLAPACK
#
#  Module:          Makefile
#
#  Purpose:         Top-level Makefile
#
#  Creation date:   March 20, 1995
#
#  Modified:        February 15, 2000
#
#  Send bug reports, comments or suggestions to scalapack@cs.utk.edu
#
############################################################################

include SLmake.inc

#PRECISIONS = single double complex complex16 FRC=FRC
PRECISIONS = single double complex complex16

############################################################################
#
#  The library can be set up to include routines for any combination of the
#  four PRECISIONS.  First, modify the ARCH, ARCHFLAGS, RANLIB, F77, CC,
#  F77FLAGS, CCFLAGS, F77LOADER, CCLOADER, F77LOADFLAGS, CCLOADFLAGS and
#  CDEFS definitions in SLmake.inc to match your library archiver, compiler
#  and the options to be used.
#
#  The command
#       make
#  without any arguments creates the library of precisions defined by the
#  environment variable PRECISIONS as well as the corresponding testing
#  executables,
#       make lib
#  creates only the library,
#       make exe
#  creates only the testing executables.
#       make example
#  creates only the example
#
#  The name of the library is defined in the file called SLmake.inc and
#  is created at this directory level.
#
#  To remove the object files after the library and testing executables
#  are created, enter
#       make clean
#
############################################################################

all: lib
#all: lib exe example

lib: toolslib pblaslib redistlib scalapacklib

exe: pblasexe redistexe scalapackexe

clean: cleanlib cleanexe cleanexample

pblaslib:
	( cd $(PBLASdir)/SRC; $(MAKE) $(PRECISIONS) )

redistlib:
	( cd $(REDISTdir)/SRC; $(MAKE) integer $(PRECISIONS) )

scalapacklib:
	( cd $(SRCdir); $(MAKE) $(PRECISIONS) )

toolslib:
	( cd $(TOOLSdir); $(MAKE) $(PRECISIONS) )

pblasexe:
	( cd $(PBLASdir)/TESTING; $(MAKE) $(PRECISIONS) )
	( cd $(PBLASdir)/TIMING; $(MAKE) $(PRECISIONS) )

scalapackexe:
	( cd $(TESTdir)/LIN; $(MAKE) $(PRECISIONS) )
	( cd $(TESTdir)/EIG; $(MAKE) $(PRECISIONS) )

redistexe:
	( cd $(REDISTdir)/TESTING; $(MAKE) integer $(PRECISIONS) )

example:
	( cd EXAMPLE; $(MAKE) $(PRECISIONS) )

cleanexe:
	( cd $(PBLASdir)/TESTING; $(MAKE) clean )
	( cd $(PBLASdir)/TIMING; $(MAKE) clean )
	( cd $(TESTdir)/LIN; $(MAKE) clean )
	( cd $(TESTdir)/EIG; $(MAKE) clean )
	( cd $(REDISTdir)/TESTING; $(MAKE) clean )

cleanlib:
	( cd $(PBLASdir)/SRC; $(MAKE) clean )
	( cd $(SRCdir); $(MAKE) clean )
	( cd $(TOOLSdir); $(MAKE) clean )
	( cd $(REDISTdir)/SRC; $(MAKE) clean )

cleanexample:
	( cd EXAMPLE; $(MAKE) clean )

