all : lib tester

clean:
	( cd TESTING ; $(MAKE) clean )
	( cd SRC ; $(MAKE) clean )

tester :
	( cd TESTING ; $(MAKE) )

lib :
	( cd SRC ; $(MAKE) )
