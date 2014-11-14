

def cache(storage, key, regenerate=False, use_cache=False, clear_cache=False,
          verbosity=0):
    """Caching decorator"""
    def cacher(func):
        if key is None:
            key_ = (func.func_name, ) + \
                  tuple(zip(func.func_code.co_varnames, func.func_defaults))
        else:
            key_ = key
        if key_ in storage:
            # Already in cache
            val = storage[key_]
            if verbosity > 0:
                print "Cached"
            return val
        if verbosity > 0:
            print "Not cached, generating"
        val = func()
        storage[key_] = val

        return val
    return cacher

