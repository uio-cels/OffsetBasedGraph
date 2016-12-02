def takes(*targs):
    def dec(func):
        def new_func(*args, **kwargs):
            i = 0
            for arg, targ in zip(args[1:], targs):
                assert isinstance(arg, targ), "%s-%s: %s is not %s" % (func.__name__, i, arg, targ)
                i += 1
            return func(*args, **kwargs)
        return new_func
    return dec
