CC= g++
CFLAGS = -O3 -Wno-unused-function -std=c++11 -w -march=native -fopenmp #turn on warning
CPPFLAGS= -O3 -Wno-unused-function -std=c++11 -w -march=native -fopenmp
INCLUDES=	-I.

LIBS = -lm -fopenmp -lz -lpthread

_DEPS = creads.h cseq.h khash.h kseq.h kvec.h kmeans.h
DEPS = $(pasubst %,$(IDIR)/%,$(_DEPS))

_OBJ=main.o cseq.o preprocess.o sketch.o kmeans.o matrix.o
OBJ=$(patsubst %,%,$(_OBJ))

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

MRC: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean: rm -f *.o