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

class _cobj(object):
    """Generic code to interface an object with ctypes using the darst method.

    Note: work in progress.
    """
    def _allocStruct(self, cModel):
        self._struct = cModel()
        self._struct_p = ctypes.pointer(self._struct)

    def _allocArray(self, name, array=None, **args):
        """Allocate an array in a way usable by both python and C.

        `name` is the name of the array to allocate: it will be
        self.`name` and self.SD.`name`.  `**args` are arguments passed
        to `numpy.zeros` to allocate the arrays.

        If you want to initilize the array to something, you have to
        do it afterwards, like this:
        self.lattsite[:] = S12_EMPTYSITE
        """
        # This code can be hard to understand, since we can't use the
        # name, but it is equivalent to this for the name `lattsite`:
        ##  self.__dict__["lattsite"] = numpy.zeros(**args)
        ##  self.SD.lattsite = self.lattsite.ctypes.data
        # Get dtype if not already a argument
        if 'dtype' not in args:
            args['dtype'] = getattr(self._struct, name)._type_
        if array is None:
            array = numpy.zeros(**args)
        self.__dict__[name] = array
        setattr(self._struct, name,
                array.ctypes.data_as(ctypes.POINTER(
                    getattr(self._struct, name)._type_)))
        #setattr(self._struct, name, array.ctypes.data)
    def _allocArrayPointers(self, pointerarray, array):
        for i, row in enumerate(array):
            pointerarray[i] = row.ctypes.data

    def __getattr__(self, attrname):
        """Wrapper to proxy attribute gets to the C SimData struct
        """
        if attrname in cGraph_fields:
            return getattr(self._struct, attrname)
        raise AttributeError("No Such Attribute: %s"%attrname)
    def __setattr__(self, attrname, value):
        """Wrapper to proxy attribute sets to the C SimData struct
        """
        if attrname in cGraph_fields:
            setattr(self._struct, attrname, value)
        else:
            self.__dict__[attrname] = value

# This defines the C structure we use
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
    ("overlapMinimize_add",
                      c_int,    (cGraph_p, c_double)),
    ("overlapMinimize_remove",
                      c_int,    (cGraph_p, c_double)),
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



