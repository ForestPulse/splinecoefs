GCC = gcc
G11=g++ -std=c++11

#BINDIR=$(HOME)/bin
BINDIR=/usr/local/bin

#GDAL=-I/usr/include/gdal -L/usr/lib -Wl,-rpath=/usr/lib
GDAL=-I/usr/include/gdal -L/usr/lib/x64_64-linux-gnu -Wl,-rpath=/usr/lib
LDGDAL=-lgdal

GSL=-I/usr/local/include -L/usr/local/lib -Wl,-rpath=/usr/local/lib -DHAVE_INLINE=1 -DGSL_RANGE_CHECK=0
LDGSL=-lgsl -lgslcblas

CFLAGS=-fopenmp -O3 -Wall 
#CFLAGS=-g -Wall -fopenmp 

.PHONY: all install clean

all: splinecoefs

utils: src/utils.c
	$(GCC) $(CFLAGS) -c src/utils.c -o utils.o

#usage: src/usage.c
#	$(GCC) $(CFLAGS) $(GDAL) -c src/usage.c -o usage.o $(LDGDAL)

alloc: src/alloc.c
	$(GCC) $(CFLAGS) -c src/alloc.c -o alloc.o

date: src/date.c
	$(GCC) $(CFLAGS) -c src/date.c -o date.o

dir: src/dir.c
	$(GCC) $(CFLAGS) -c src/dir.c -o dir.o

image_io: src/image_io.c
	$(GCC) $(CFLAGS) $(GDAL) -c src/image_io.c -o image_io.o $(LDGDAL)

quality: src/quality.c
	$(GCC) $(CFLAGS) -c src/quality.c -o quality.o

table: src/table.c
	$(GCC) $(CFLAGS) -c src/table.c -o table.o

stats: src/stats.c
	$(GCC) $(CFLAGS) $(GSL) $(GDAL) -c src/stats.c -o stats.o $(LDGSL) $(LDGDAL)

#write: src/write.c
#	$(G11) $(CFLAGS) $(GDAL) -c src/write.c -o write.o $(LDGDAL)

string: src/string.c
	$(GCC) $(CFLAGS) -c src/string.c -o string.o


splinecoefs: utils alloc dir string date image_io quality table stats src/_splinecoefs.c
	$(GCC) $(CFLAGS) $(GSL) $(GDAL) -o splinecoefs src/_splinecoefs.c utils.o alloc.o dir.o string.o date.o image_io.o quality.o table.o stats.o -lm $(LDGSL) $(LDGDAL)

install:
	cp splinecoefs $(BINDIR) ; chmod 755 $(BINDIR)/splinecoefs

clean:
	rm -f splinecoefs *.o
