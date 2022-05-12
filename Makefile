IDIR =.
CC=g++
CFLAGS=-I$(IDIR)

ODIR = .
LDIR = .

_DEPS = constants.h corrections.h weak_rates.h nu_elastic_scatt.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = check_rates.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.c $(DEPS)
			$(CC) -c -o $@ $< $(CFLAGS)

.PHONY : all
all : nu_rates

.PHONY: clean
clean:
				rm -f nu_rates *.o

nu_rates: $(OBJ)
			$(CC) -o nu_rates $^ $(CFLAGS) $(LIBS)
