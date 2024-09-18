CC = g++
CCFLAGS = -O3
LIBS = -lm -lgsl -lgslcblas -lz
TARGET = fastenloc

# Control OpenMP flag with USE_OPENMP (default is to not use it)
USE_OPENMP ?= 0

ifeq ($(shell uname),Darwin) # macOS
    ifeq ($(USE_OPENMP),1)
        CCFLAGS += -Xpreprocessor -fopenmp
        LIBS += -lomp
    endif
else # Linux
    ifeq ($(USE_OPENMP),1)
        CCFLAGS += -fopenmp
    endif
endif

SRCS = main.cc controller.cc sigCluster.cc
OBJS = $(SRCS:.cc=.o)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CCFLAGS) $^ -o $@ $(LIBS)

%.o: %.cc
	$(CC) $(CCFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
