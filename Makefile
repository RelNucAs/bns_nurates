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

nu_rates: $(OBJ)
			$(CC) -o nu_rates $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
				rm nu_rates *.o
