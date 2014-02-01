"""Lazy evaluation for Python.
"""


# From http://stackoverflow.com/questions/9057669/how-can-i-intercept-calls-to-pythons-magic-methods-in-new-style-classes
class Lazy(object):
    """Wrapper class providing __x__ magic methods.

    This class provides lazy access to an instance of some internal
    instance.  It is created with an argument of a function which
    returns an instance of cls.  The function is not called until some
    property of the object is called.  Then, we load the function and
    proxy all attribute lookups to that instance.  Most complexity
    comes from having to proxy the __x__ methods, which do not use
    __getattr__.

    You must make an explicit subclass and set __wraps__ to the type
    of object wrapped.  This is needed so that __metaclass__ can make
    explicit proxy functions for the magic methods.

    Usage:

    gen = lambda: networkx.complete_graph(100)
    p = NxGraphInstance(gen)    # does not call `gen` yet
    len(p)                      # calls `gen` now and returns len

    self._lazy_func: None or callable
        callable used to create the object.  None otherwise.

    self._lazy_obj: None or object.
        when object is created, it is stored here.  None otherwise.

    self._lazy_make(): method
        method to create the internal object

    self._lazy_clear(): method
        resets self._lazy_obj to None, freeing memory."""

    __wraps__  = None
    __ignore__ = set(("class", "mro", "new", "init",
                      "setattr", "getattr", "getattribute",
                      "dict",
                      "str", "repr", # printing doesn't create object.
                      "reduce_ex", "reduce",  # allow pickle support
                      "getstate", "slots",    # further pickle support.
                      ))
    #__ignore__ = "class mro new init setattr getattr getattribute dict obj func"
    _lazy_obj      = None
    _lazy_func     = None

    def __init__(self, func):
        """Initialization

        func: callable
            Function called to initialize object."""
        if self.__wraps__ is None:
            raise TypeError("base class Wrapper may not be instantiated")
        if not hasattr(func, '__call__'):
            raise TypeError('func %s is not callable'%func)
        self._lazy_func = func
    def _lazy_make(self):
        """Make the internal object unconditionally."""
        self._lazy_obj = self._lazy_func()
        if not isinstance(self._lazy_obj, self.__wraps__):
            raise ValueError("obj must be instance of %s, not %s"%(
                self.__wraps__.__name__, self._lazy_obj.__class__.__name__))
    def _lazy_clear(self):
        """Clear the saved object"""
        self._lazy_obj = None
    def _lazy_is_loaded(self):
        """True if internal object is created."""
        return self._lazy_obj is not None

    # provide lazy access to regular attributes of wrapped object
    def __getattr__(self, name):
        if self._lazy_func is None:
            # This should only be used during unpickling.  During
            # unpickling, '__func' may not yet be set on the
            # dictionary.
            return
        if self._lazy_obj is None:
            self._lazy_make()
        #print "getattr: %s"%name
        return getattr(self._lazy_obj, name)

    # create proxies for wrapped object's double-underscore
    # attributes, except those in __ignore__
    class __metaclass__(type):
        def __init__(cls, name, bases, dct):
            def make_lazy(name):
                def lazy(self, *args):
                    #print dir(self)
                    #print 'lazy: %s'%name
                    if self._lazy_obj is None:
                        self._lazy_make()
                    return getattr(self._lazy_obj, name)
                return lazy

            type.__init__(cls, name, bases, dct)
            if cls.__wraps__:
                ignore = set("__%s__" % n for n in cls.__ignore__)
                for name in dir(cls.__wraps__):
                    if name.startswith("__"):
                        if name not in ignore and name not in dct:
                            #print name
                            setattr(cls, name, property(make_lazy(name)))

    # Needed for pickle support.  Otherwise tries for __slots__.
    def __getstate__(self):
        """Pickle support, must be defined."""
        return self.__dict__
    def __str__(self):
        return '<%s object at %#x (%s)>'%(self.__class__.__name__,
                          id(self), 'loaded' if self._lazy_is_loaded() else 'not loaded')


# These proxies must be explicitly defined in a globally importable
# way.  So any class you might want to proxy, put it here.

import networkx
class NxGraphLazy(Lazy):
    __wraps__ = networkx.Graph



def _test_creator():
    """Helper function for _test.  Must be in globals() for pickling."""
    print "making graph"
    return networkx.complete_graph(10)

def _test():
    from functools import partial
    p = NxGraphLazy(_test_creator)
    print str(p)
    print repr(p)

    assert p._lazy_obj is None
    assert len(p) == 10
    assert p._lazy_obj is not None
    p._lazy_clear()
    assert p._lazy_obj is None

    import cPickle as pickle
    import pickle

    # Test pickling and unpickling.
    print type(p)
    #print p.__dict__
    #from fitz import interactnow
    print "pickling..."
    pkl = pickle.dumps(p)
    print "unpickling..."
    p2 = pickle.loads(pkl)
    print type(p2)
    assert isinstance(p2, NxGraphLazy)
    #print len(p2)
    assert p2._lazy_obj is None
    print [x for x in p2]
    assert p2._lazy_obj is not None

    p2._lazy_clear()
    assert p2._lazy_obj is None

    p = NxGraphLazy(lambda: 1)
    #print len(p)  # should fail.

if __name__ == '__main__':
    _test()
