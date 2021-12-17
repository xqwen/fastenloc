CC = g++
CFLAG = -O3  -fopenmp
LIB = -lm -lgsl -lgslcblas -lboost_iostreams -lz




fastenloc: main.o controller.o sigCluster.o
	$(CC) $(CFLAG)  $^ -o $@ $(LIB)
static: main.o controller.o sigCluster.o
	  	$(CC) $(CFLAG) $^ -o fastenloc.static $(LIB) -static
main.o: main.cc
	$(CC) $(CFLAG) -c  main.cc
controller.o: controller.cc controller.h
	$(CC) $(CFLAG) -c  controller.cc
sigCluster.o: sigCluster.h sigCluster.cc
	$(CC) $(CFLAG) -c sigCluster.cc
clean:
	rm *.o fastenloc
