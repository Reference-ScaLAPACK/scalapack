all : lib tester

clean:
	( cd TESTING ; make clean )
	( cd SRC ; make clean )
	( cd INSTALL; make clean )

tester :
	( cd TESTING ; make )

lib :
	( cd SRC ; make )
