# Richard Darst, July 2011

CFLAGS=-O3

opts=-Wall -shared -fPIC

_cmodels.so: cmodels.o SFMT.o
	gcc ${opts} ${CFLAGS} cmodels.o SFMT.o -o _cmodels.so

cmodels.o: cmodels.c
	gcc ${opts} ${CFLAGS} -c cmodels.c

SFMT.o: SFMT.c SFMT.h
	gcc ${opts} ${CFLAGS} -DMEXP=19937 -include SFMT-params.h -c SFMT.c

#test:
#        python tests/unittests_run.py
