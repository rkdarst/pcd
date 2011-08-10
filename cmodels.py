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

class cGraph(ctypes.Structure):
    pass
cGraph._fields_ = [
    ("N",            c_int),
    ("Ncmty",        c_int),
    ("oneToOne",     c_int),

    ("cmty",         c_int_p),
    ("interactions", c_int_p),

    ("cmtyl",        ctypes.POINTER(c_int_p)),
    ("cmtyll",       c_int_p),
    ("cmtyN",        c_int_p),
    ("nodeOrder",    c_int_p),
    
    #("callback", Callback),   # callback to let us get python shell from C
    #("S", ctypes.py_object),  # Pointer for function above
    ]
cGraph_p = ctypes.POINTER(cGraph)
cGraph_fields = { }
for name, t in cGraph._fields_:
    cGraph_fields[name] = True


cfuncs = (
    ("test",          c_int,    (cGraph_p, )),
    ("minimize0",     c_int,    (cGraph_p, c_double)),
    ("minimize",      c_int,    (cGraph_p, c_double)),
    ("energy",        c_double, (cGraph_p, c_double)),
    ("energy_cmty",   c_double, (cGraph_p, c_double, c_int)),
    ("combine_cmtys", c_int,    (cGraph_p, c_double)),
    ("remap_cmtys",   c_int,    (cGraph_p, )),

    ("cmtyListAdd",      None,      (cGraph_p, c_int, c_int)),
    ("cmtyListRemove",   None,      (cGraph_p, c_int, c_int)),
    ("cmtyListInit",     None,      (cGraph_p, )),
    ("cmtyListCheck",    c_int,     (cGraph_p, )),
    )
print __file__
filename = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_cmodels.so')
C = ctypes.cdll[os.path.join(os.path.dirname(__file__), filename)]
                        
for name, restype, argtypes in cfuncs:
    getattr(C, name).restype  = restype
    getattr(C, name).argtypes = argtypes
    globals()[name] = getattr(C, name)

C.init_gen_rand.restype = None
C.init_gen_rand(random.randrange(2**32-1))



