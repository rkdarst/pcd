# Richard Darst, July 2011

CFLAGS=-O3
IMATRIX_T=int

opts=-Wall -Wextra -lm -shared -fPIC

_default: _cmodels.so
.PHONY: _settypes

# The following line defines the type to be used for the interaction
# matrix.  To use it, make with IMATRIX_T={int,float,double}, probably
# int or float are the most useful types.
#
# To do this definition, we use two files: imatrix_t.h and
# imatrix_t.py.  imatrix_t.h has a #define IMATRIX_T <type> (and this
# is typedef'ed in cmodels.h).  imatrix_t.py has just "c_<type>" in
# it, and this is read and eval'ed in Python
#
# Both of these are very low-tech, but makes things statically defined
# (important if I ever make Python smart enough to read cmodels.h and
# automatically detect things.
_settypes:
#	Define the imatrix type for C use
	echo "#define IMATRIX_T ${IMATRIX_T}" > imatrix_t.h
	@if [ ${IMATRIX_T} = int ] ; then \
	  echo "#define IMATRIX_T_INT" >> imatrix_t.h ; \
	fi
#	Define the imatrix type for Python use
	echo "c_${IMATRIX_T}" > imatrix_t.py

_cmodels.so: cmodels.o SFMT.o
	$(CC) ${opts} ${CFLAGS} -lm cmodels.o SFMT.o -o _cmodels.so

cmodels.o: _settypes cmodels.c cmodels.h
	$(CC) ${opts} ${CFLAGS} -c cmodels.c

SFMT.o: SFMT.c SFMT.h
	$(CC) ${opts} ${CFLAGS} -DMEXP=19937 -include SFMT-params.h -c SFMT.c

#test:
#        python tests/unittests_run.py
