include ../../SLmake.inc

help :
	@ echo " "
	@ echo " "
	@ echo "You need to specify which INSTALL executable to build."
	@ echo "General choices are: xsize, xintface, xsyserrors"
	@ echo "MPI specific choices are: xsyserrors, xtc_CsameF77, xtc_UseMpich,"
	@ echo "                          xcmpi_sane, xfmpi_sane"
	@ echo " "
	@ echo "Here is a brief explanation of each of these routines: "
	cat README

xsize : size.o
	$(CCLOADER) $(CCLOADFLAGS) -o $@ size.o

xintface : Fintface.o Cintface.o
	$(FCLOADER) $(FCLOADFLAGS) -o $@ Fintface.o Cintface.o

xsyserrors : syserrors.o
	$(CCLOADER) $(CCLOADFLAGS) -o $@ syserrors.o

xtc_CsameFC : tc_fCsameF77.o tc_cCsameF77.o
	$(FCLOADER) $(FCLOADFLAGS) -o $@ tc_fCsameF77.o tc_cCsameF77.o

xtc_UseMpich : tc_UseMpich.o
	$(CCLOADER) $(CCLOADFLAGS) -o $@ tc_UseMpich.o

xcmpi_sane : cmpi_sane.o
	$(CCLOADER) $(CCLOADFLAGS) -o $@ cmpi_sane.o

xfmpi_sane : mpif.h fmpi_sane.o
	$(FCLOADER) $(FCLOADFLAGS) -o $@ fmpi_sane.o

clean:
	rm -f size.o Fintface.o Cintface.o syserrors.o transcomm.o \
              mpi_sane.o fmpi_sane.o tc_UseMpich.o tc_fCsameF77.o tc_cCsameF77.o

.f.o: ; $(FC) -c $(FCFLAGS) $*.f

.c.o:
	$(CC) -c $(CCFLAGS) $(CDEFS) $<
