# Richard Darst, July 2011

CFLAGS=-O3

opts=-Wall -lm -shared -fPIC

_cmodels.so: cmodels.o SFMT.o
	$(CC) ${opts} ${CFLAGS} -lm cmodels.o SFMT.o -o _cmodels.so

cmodels.o: cmodels.c cmodels.h
	$(CC) ${opts} ${CFLAGS} -c cmodels.c

SFMT.o: SFMT.c SFMT.h
	$(CC) ${opts} ${CFLAGS} -DMEXP=19937 -include SFMT-params.h -c SFMT.c

#test:
#        python tests/unittests_run.py
