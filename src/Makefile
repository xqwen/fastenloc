CXX=g++
CXXFLAGS+=-O3 -Wall
INC=-isystem /opt/homebrew/include # Mac libraries installed via Homebrew
LIBS=-lm -lgsl -lgslcblas -lboost_iostreams -lz -lomp
LDFLAGS=-L/opt/homebrew/lib

objects = main.o controller.o sigCluster.o

all: fastenloc

fastenloc: $(objects)
	$(CXX) $(CXXFLAGS) $(INC) -Xpreprocessor -fopenmp $(LDFLAGS) $(objects) $(LIBS) -o fastenloc
static: $(objects)
	$(CXX) $(CXXFLAGS) $(INC) -Xpreprocessor -fopenmp $(LDFLAGS) $(objects) $(LIBS) -static -o fastenloc.static
main.o: main.cc
	$(CXX) $(CXXFLAGS) $(INC) -c  main.cc
controller.o: controller.cc controller.h
	$(CXX) $(CXXFLAGS) $(INC) -Xpreprocessor -fopenmp -c  controller.cc
sigCluster.o: sigCluster.h sigCluster.cc
	$(CXX) $(CXXFLAGS) $(INC) -c sigCluster.cc

.PHONY: clean
clean:
	rm *.o fastenloc
