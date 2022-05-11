IDIR =.
CC=g++
CFLAGS=-I$(IDIR)

ODIR = .
LDIR = .

_DEPS = constants.h corrections.h weak_rates.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o
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
