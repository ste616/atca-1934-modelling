CC=gcc
FF=gfortran
LIBS=-lpgplot -lcpgplot -lpng -lz -lX11 -lgsl -lgslcblas
CFLAGS=-Wall

refit_1934: refit_1934.o
	${FF} -o $@ refit_1934.o ${LIBS}
