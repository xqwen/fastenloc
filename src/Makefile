CC = g++
CCFLAGS = -O3
LIBS = -lm -lgsl -lgslcblas -lz
TARGET = fastenloc

ifeq ($(shell uname),Darwin) # macOS
    CCFLAGS += -Xpreprocessor -fopenmp
    LIBS += -lomp
else # Linux
    CCFLAGS += -fopenmp
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

