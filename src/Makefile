CXX=g++
CPPFLAGS+=-isystem /opt/homebrew/include # Mac libraries installed via Homebrew
CPPFLAGS+=-Xpreprocessor -fopenmp # Mac openmp
CXXFLAGS+=-O3 -Wall
LDLIBS=-lm -lgsl -lgslcblas -lboost_iostreams -lz -lomp #gcc build doesn't need -lomp?
LDFLAGS=-L/opt/homebrew/lib

objects = main.o controller.o sigCluster.o

all: fastenloc

fastenloc: $(objects)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $(objects) $(LDLIBS) -o fastenloc
static: $(objects)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $(objects) $(LDLIBS) -static -o fastenloc.static

controller.o: controller.cc controller.h
sigCluster.o: sigCluster.h sigCluster.cc

.PHONY: clean
clean:
	rm *.o fastenloc
