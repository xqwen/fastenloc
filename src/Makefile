CXX=g++
CXXFLAGS+=-O3 -Wall
LDLIBS+=-lm -lgsl -lgslcblas -lboost_iostreams -lz

UNAME_S := $(shell uname -s)

# MacOS / Homebrew specific configuration
ifeq ($(UNAME_S),Darwin)
	CPPFLAGS+=-isystem /opt/homebrew/include # Mac libraries installed via Homebrew
	CPPFLAGS+=-Xpreprocessor -fopenmp # Mac openmp
	LDFLAGS+=-L/opt/homebrew/lib # Homebrew library location
	LDLIBS+=-lomp #only need lomp on Mac OS
endif


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
