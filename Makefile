SRC=$(wildcard *.c) 
SLIBS=/usr/lib/libnetcdf.a
OBJ=$(SRC:.c=.o)
LIBS=-lOpenCL -pthread -lm
CFLAGS=-pipe -Werror 
BIN=qrtool


default: fast

fast:
	$(MAKE) -C src/ fast
opt:
	$(MAKE) -C src/ opt
debug:
	$(MAKE) -C src/ debug
san:
	$(MAKE) -C src/ san
valgrind:
	$(MAKE) -C src/ valgrind


install: fast
	cp src/$(BIN) .
help:
	@echo "options are: valgrind san debug opt [fast]"

clean:
	rm -f $(BIN)
	$(MAKE) -C src clean
