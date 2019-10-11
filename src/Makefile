#SRC=system.c system_cl.c device.c device_cl.c main.c#$(wildcard *.c) 
SRC=$(wildcard *.c) 
SLIBS=/usr/lib/libnetcdf.a
OBJ=$(SRC:.c=.o)
LIBS=-lOpenCL -pthread -lm
CFLAGS=-pipe -Werror 
BIN=qrtool


default: fast

help:
	@echo "options are: valgrind san [debug] opt fast"

%.o : %.c %.h defines.h
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN): $(OBJ)
	$(CC) $(CFLAGS) $(LIBS) $(OBJ) $(SLIBS) -o $(BIN)

val: valgrind
valgrind: debug
	valgrind --leak-check=full  --suppressions=suppressions.conf ./a.out 

san: CFLAGS += -fsanitize=address
san: debug

sant: CFLAGS += -fsanitize=thread -fPIE -fPIC -pie
sant: debug

debug: CFLAGS += -O0 -ggdb3 
debug: $(BIN)

prof: CFLAGS += -pg -g
prof: debug

opt: CFLAGS += -O2
opt: $(BIN)

fast: CFLAGS += -O3 -fomit-frame-pointer -march=native -mtune=native
fast: $(BIN) 

clean:
	rm -f $(OBJ) $(BIN)
