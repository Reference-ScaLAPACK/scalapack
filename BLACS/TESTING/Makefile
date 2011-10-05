include ../../SLmake.inc

#  ---------------------------------------------------------------------
#  The file tools.f contains some LAPACK routines that the tester calls.
#  If you have ScaLAPACK, you may point to your tools library instead
#  of compiling this file.
#  ---------------------------------------------------------------------
tools = tools.o

exe : all
ctest : xCbtest
ftest : xFbtest
all : xCbtest xFbtest

obj = blacstest.o btprim.o

xCbtest: $(obj) $(tools)
	$(CC) -c $(CCFLAGS) -DBTCINTFACE Cbt.c
	$(F77LOADER) $(F77LOADFLAGS) -o $@ $(obj) $(tools) Cbt.o ../../$(SCALAPACKLIB)

xFbtest: $(obj) $(tools)
	$(F77LOADER) $(F77LOADFLAGS) -o $@ $(obj) $(tools) ../../$(SCALAPACKLIB)

#  --------------------------------------------------------------------
#  The files tools.f and blacstest.f are compiled without optimization.
#  Tools.f contains the LAPACK routines slamch and dlamch, which only
#  operate correctly for low-levels of optimization.  Blacstest.f is
#  extremely large, and optimizing it takes a long time.  More
#  importantly, the sun's f77 compiler seems to produce errors in
#  trying to optimize such a large file.  We therefore insist that it
#  also not be optimized.
#  --------------------------------------------------------------------
tools.o : tools.f
	$(F77) $(F77NO_OPTFLAGS) -c $*.f

blacstest.o : blacstest.f
	$(F77) $(F77NO_OPTFLAGS) -c $*.f

btprim.o : btprim.f
	$(F77) -c $(F77FLAGS) $*.f

clean :
	rm -f $(obj) tools.o Cbt.o xCbtest xFbtest

.f.o: ; $(F77) -c $(F77FLAGS) $*.f
