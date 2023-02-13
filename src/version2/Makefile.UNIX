CC = g++
CCFLAGS = -O3  -fopenmp
LIBS = -lm -lgsl -lgslcblas -lboost_iostreams -lz




fastenloc: main.o controller.o sigCluster.o
	$(CC) $(CCFLAGS)  $^ -o $@ $(LIBS)
static: main.o controller.o sigCluster.o
	  	$(CC) $(CCFLAGS) $^ -o fastenloc.static $(LIBS) -static
main.o: main.cc
	$(CC) $(CCFLAGS) -c  main.cc
controller.o: controller.cc controller.h
	$(CC) $(CCFLAGS) -c  controller.cc
sigCluster.o: sigCluster.h sigCluster.cc
	$(CC) $(CCFLAGS) -c sigCluster.cc
clean:
	rm *.o fastenloc
