def takes(*targs):
    def dec(func):
        def new_func(*args, **kwargs):
            for arg, targ in zip(args, targs):
                assert isinstance(arg, targ), "%s is not %s" % (arg, targ)
            return func(*args, **kwargs)
        return new_func
    return dec
