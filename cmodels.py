# Richard Darst, July 2011

import ctypes
import numpy
from numpy.ctypeslib import ndpointer
import os
import random

c_double = ctypes.c_double
c_void_p = ctypes.c_void_p
c_int = ctypes.c_int
c_int_p = ctypes.POINTER(c_int)

imatrix_t = c_int_p

class cGraph(ctypes.Structure):
    pass
cGraph._fields_ = [
    ("N",            c_int),
    ("Ncmty",        c_int),
    ("oneToOne",     c_int),

    ("cmty",         c_int_p),
    ("imatrix",      imatrix_t),

    ("cmtyl",        ctypes.POINTER(c_int_p)),
    ("cmtyll",       c_int_p),
    ("cmtyN",        c_int_p),
    ("randomOrder",  c_int_p),
    ("randomOrder2", c_int_p),

    #("callback", Callback),   # callback to let us get python shell from C
    #("S", ctypes.py_object),  # Pointer for function above
    ]
cGraph_p = ctypes.POINTER(cGraph)
cGraph_fields = { }
for name, t in cGraph._fields_:
    cGraph_fields[name] = True


cfuncs = (
    ("test",          c_int,    (cGraph_p, )),

    ("cmtyListAdd",      None,      (cGraph_p, c_int, c_int)),
    ("cmtyListRemove",   None,      (cGraph_p, c_int, c_int)),
    ("cmtyListInit",     None,      (cGraph_p, )),
    ("cmtyListCheck",    c_int,     (cGraph_p, )),

    ("q",                c_int,     (cGraph_p, )),
    ("entropy",          c_double,  (cGraph_p, )),
    ("mutual_information",
                         c_double, (cGraph_p, cGraph_p )),

    ("energy_naive",  c_double, (cGraph_p, c_double)),
    ("energy",        c_double, (cGraph_p, c_double)),
    ("energy_cmty",   c_double, (cGraph_p, c_double, c_int)),
    ("energy_cmty_n", c_double, (cGraph_p, c_double, c_int, c_int)),
    ("energy_cmty_cmty", c_double, (cGraph_p, c_double, c_int, c_int)),
    ("energy_n",      c_double, (cGraph_p, c_double, c_int)),

    ("minimize_naive",c_int,    (cGraph_p, c_double)),
    ("minimize",      c_int,    (cGraph_p, c_double)),
    ("combine_cmtys", c_int,    (cGraph_p, c_double)),
    ("remap_cmtys",   c_int,    (cGraph_p, )),

    )

filename = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_cmodels.so')
C = ctypes.cdll[os.path.join(os.path.dirname(__file__), filename)]

for name, restype, argtypes in cfuncs:
    getattr(C, name).restype  = restype
    getattr(C, name).argtypes = argtypes
    globals()[name] = getattr(C, name)

C.init_gen_rand.restype = None
C.init_gen_rand(random.randrange(2**32-1))



