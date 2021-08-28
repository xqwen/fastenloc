main: main.o controller.o sigCluster.o
	g++ -fopenmp -O3 main.o controller.o sigCluster.o  -lm -lgsl -lgslcblas -lboost_iostreams -lz  -o fastenloc
static: main.o controller.o sigCluster.o
	g++ -fopenmp -O3 main.o controller.o sigCluster.o  -lm -lgsl -lgslcblas -lboost_iostreams -lz  -static -o fastenloc.static
main.o: main.cc
	g++ -c  main.cc
controller.o: controller.cc controller.h
	g++ -fopenmp   -c  controller.cc
sigCluster.o: sigCluster.h sigCluster.cc
	g++ -c sigCluster.cc 
clean:
	rm *.o fastenloc
