all : lib tester

clean:
	( cd TESTING ; make clean )
	( cd SRC ; make clean )
	( cd INSTALL; make clean )

tester :
	( cd TESTING ; make $(what) )

lib :
	( cd SRC ; make $(what) )
