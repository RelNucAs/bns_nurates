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

nu_rates: $(OBJECTS) main.o
	$(CXX) $(RFLAGS) -o ./main ./main.o $(OBJECTS)

main.o: main.cpp
	$(CXX) $(RFLAGS) -c main.cpp

#program.o: $(SOURCES)
#	$(CXX) $(RFLAGS) -c $(SOURCES)

#$(OBJECTS): $(SOURCES)
#	$(CXX) $(RFLAGS) -c -o $@ $^

%.o: %.cpp
	$(CXX) $(RFLAGS) -c -o $@ $<
