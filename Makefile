FOR = gfortran -std=legacy -ffixed-line-length-none
LIB = -lmathlib

ang_dist: ang_dist.for w.o
	$(FOR) -o ang_dist w.o ang_dist.for $(LIB)

w.o: w.for
	$(FOR) -c w.for
	
clean:
	rm ang_dist *.o
