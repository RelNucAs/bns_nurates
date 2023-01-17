CXX      := g++
RFLAGS   := -std=c++17 -O3 -g -Wall #-Wpedantic
#CPPFLAGS += $(HDF5_INCLUDE)
#LIB      := ar cr
#FLAGS    := 
#LDFLAGS  += $(HDF5_LIBS)

#INCLUDE = -I/usr/local/include
#LIBS    = -L/usr/local/lib

#INCLUDE += -I$(CURDIR)/include

COMPILE_DIR = ./compile/

all: compile

compile: test_eos
	$(CXX) $(RFLAGS) -o ./main test_eos.cpp

test_performance: test_eos_performance_new.cpp
	$(CXX) $(RFLAGS) -o $(COMPILE_DIR)test_eos_performance_new test_eos_performance_new.cpp

test_accuracy: test_eos_accuracy.cpp
	$(CXX) $(RFLAGS) -o $(COMPILE_DIR)test_eos_accuracy test_eos_accuracy.cpp

test_eos: test_eos.cpp
	$(CXX) $(RFLAGS) -o $(COMPILE_DIR)test_eos test_eos.cpp

eta_table: generate_eta_table.cpp
	$(CXX) $(RFLAGS) -o $(COMPILE_DIR)generate_eta_table generate_eta_table.cpp

lep_table: generate_lep_table.cpp
	$(CXX) $(RFLAGS) -o $(COMPILE_DIR)generate_lep_table generate_lep_table.cpp

complete_table: generate_complete_table.cpp
	$(CXX) $(RFLAGS) -o $(COMPILE_DIR)generate_complete_table generate_complete_table.cpp


clean:
	rm -v -f $(COMPILE_DIR)*

