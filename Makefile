FOR = gfortran -std=legacy -ffixed-line-length-none

all:ang_dist

ang_dist: ang_dist.f w.o dwig3j64.o dwig9j64.o
	$(FOR) -o ang_dist w.o dwig3j64.o dwig9j64.o ang_dist.f

w.o: w.f
	$(FOR) -c w.f

dwig3j64.o: cernlib/dwig3j64.f
	$(FOR) -c cernlib/dwig3j64.f

dwig9j64.o: cernlib/dwig9j64.f
	$(FOR) -c cernlib/dwig9j64.f
	
clean:
	rm ang_dist *.o
