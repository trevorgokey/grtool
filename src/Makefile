#SRC=system.c system_cl.c device.c device_cl.c main.c#$(wildcard *.c) 
SRC=$(wildcard *.c) 
SLIBS=/usr/lib64/libnetcdf.a
# SLIBS=-lnetcdf
OBJ=$(SRC:.c=.o)
OCL=-L/opt/cuda/lib64 -lOpenCL  
LIBS=-L/usr/lib64 -pthread -ljpeg -lmfhdf -ldf -lhdf5_hl -lhdf5 -lm -lz -lcurl
LIBS+=$(OCL)
CFLAGS=-pipe 
BIN=qrtool
LDFLAGS=

default: fast

help:
	@echo "options are: valgrind san [debug] opt fast"

%.o : %.c %.h defines.h
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN): $(OBJ)
	$(CC) $(CFLAGS) $(LIBS) $(OBJ) $(SLIBS) $(LDFLAGS) -o $(BIN)

val: valgrind
valgrind: debug
	#valgrind --leak-check=full  --suppressions=suppressions.conf ./qrtool
	valgrind --leak-check=full ./qrtool

san: CFLAGS += -fsanitize=address
san: debug

sant: CFLAGS += -fsanitize=thread -fPIE -fPIC -pie
sant: debug

debug: CFLAGS += -O0 -ggdb3 
debug: $(BIN)

prof: CFLAGS += -pg -g
prof: debug

opt: CFLAGS += -O2  -ffat-lto-objects -flto
opt: LDFLAGS += -fuse-linker-plugin -ffat-lto-objects   -Wl,-O2 -Wl,--as-needed -flto=4
opt: $(BIN)

fast: CFLAGS += -O3 -fomit-frame-pointer -march=native -mtune=native -ffat-lto-objects -flto
fast: LDFLAGS += -fuse-linker-plugin -ffat-lto-objects   -Wl,-O2 -Wl,--as-needed -flto=4
fast: $(BIN) 

clean:
	rm -f $(OBJ) $(BIN)
