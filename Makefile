CXX      := g++
IDIR     := $(CURDIR)/include
RFLAGS   := -std=c++17 -O3 -g -Wall -lm -lgsl -lgslcblas -I$(IDIR) #-Wpedantic
#CPPFLAGS += $(HDF5_INCLUDE)
#LIB      := ar cr
#FLAGS    := 
#LDFLAGS  += $(HDF5_LIBS)

#INCLUDE = -I/usr/local/include
#LIBS    = -L/usr/local/lib


COMPILE_DIR = $(CURDIR)/compile/

SRC_DIR = $(CURDIR)/src/

TESTS_DIR := $(CURDIR)/tests/

TABLES_DIR := $(CURDIR)/eos_table/

#ODIR = .
#LDIR = .

#_DEPS = constants.h fermi.h nucfrmfac.h weak_rates.h nu_elastic_scatt.h
#DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

#_OBJ = main.o
#OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

#$(ODIR)/%.o: %.c $(DEPS)
#                        $(CC) -c -o $@ $< $(CFLAGS)

SOURCES := $(shell find $(SRC_DIR) -name '*.cpp')
#$(wildcard $(SRC_DIR)*/*.cpp)
OBJECTS := $(patsubst %.cpp, %.o, $(SOURCES))

.PHONY : all
all : nu_rates


.PHONY: clean
clean:
	rm -v -f $(OBJECTS)
	rm -v -f ./main
	rm -v -f ./main.o
	rm -v -f ./make_plot.o
	rm -v -f ./tests/test_NNbrem

plot: $(OBJECTS) make_plot.o
	$(CXX) $(RFLAGS) -o ./plot ./make_plot.o $(OBJECTS)

make_plot.o: make_plot.cpp
	$(CXX) $(RFLAGS) -c make_plot.cpp

nu_rates: $(OBJECTS) main.o
	$(CXX) $(RFLAGS) -o ./main ./main.o $(OBJECTS)

main.o: main.cpp
	$(CXX) $(RFLAGS) -c main.cpp

test: $(OBJECTS)
	$(CXX) $(RFLAGS) -o tests/prova tests/prova.cpp $(SOURCES)

test_integration: $(OBJECTS) tests/test_integration.cpp
	$(CXX) $(RFLAGS) -o tests/test_integration tests/test_integration.cpp $(SOURCES)

test_NN: $(OBJECTS) tests/test_NNbrem.cpp
	$(CXX) $(RFLAGS) -o tests/test_NNbrem tests/test_NNbrem.cpp $(SOURCES)

#program.o: $(SOURCES)
#	$(CXX) $(RFLAGS) -c $(SOURCES)

#$(OBJECTS): $(SOURCES)
#	$(CXX) $(RFLAGS) -c -o $@ $^

%.o: %.cpp
	$(CXX) $(RFLAGS) -c -o $@ $<
