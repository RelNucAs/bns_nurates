CXX      := g++
RFLAGS   := -std=c++11 -O3 -g -Wall #-Wpedantic
#CPPFLAGS += $(HDF5_INCLUDE)
#LIB      := ar cr
#FLAGS    := 
#LDFLAGS  += $(HDF5_LIBS)

#INCLUDE = -I/usr/local/include
#LIBS    = -L/usr/local/lib

#INCLUDE += -I$(CURDIR)/include

COMPILE_DIR = ./compile/

FILENAME = test_eos_performance_new

SRC = $(FILENAME).cpp
OFILE = $(COMPILE_DIR)$(FILENAME)

all: compile

compile: file.o
	$(CXX) $(RFLAGS) -o $(OFILE) $(OFILE).o

file.o: $(SRC)
	$(CXX) $(RFLAGS) -c -o $(OFILE).o $(SRC)

clean:
	rm -v -f $(COMPILE_DIR)*

