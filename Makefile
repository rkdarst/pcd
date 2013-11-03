# Richard Darst, July 2011

CFLAGS=-O3
#IMATRIX_T=int
IMATRIX_T=float

opts=-Wall -Wextra -lm -shared -fPIC

_default: _cmodels.so
.PHONY: _settypes
# "always" as a dep makes a rule always run.
always:

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
TYPEDEFS=imatrix_t.h imatrix_t.py
imatrix_t.h: always
#	Define the imatrix type for C use
	@if ! grep "${IMATRIX_T}" imatrix_t.h > /dev/null ; then \
	  echo "#define IMATRIX_T ${IMATRIX_T}" > imatrix_t.h ; \
	  if [ ${IMATRIX_T} = int ] ; then \
	    echo "#define IMATRIX_T_INT" >> imatrix_t.h ; \
	  fi ; \
	fi
imatrix_t.py: always
#	Define the imatrix type for Python use
	@if ! grep "c_${IMATRIX_T}" imatrix_t.py > /dev/null ; then \
	  echo "c_${IMATRIX_T}" > imatrix_t.py ; \
	fi

_cmodels.so: cmodels.o SFMT.o
	$(CC) ${opts} ${CFLAGS} cmodels.o SFMT.o -lm `pkg-config --libs glib-2.0` `pkg-config --libs gthread-2.0` -o _cmodels.so

cmodels.o: ${TYPEDEFS} cmodels.c cmodels.h
	$(CC) ${opts} ${CFLAGS} `pkg-config --cflags glib-2.0` `pkg-config --cflags gthread-2.0` -c cmodels.c

SFMT.o: SFMT.c SFMT.h
	$(CC) ${opts} ${CFLAGS} -DMEXP=19937 -include SFMT-params.h -c SFMT.c

#test:
#        python tests/unittests_run.py

support/cavity.so: support/cavity.pyx
	~/bin/cython support/cavity.pyx -a --embed=main
	gcc -I/usr/include/python2.6/ -rdynamic -fPIC -pie -lpython2.6 -g -O3 support/cavity.c -o support/cavity.so

