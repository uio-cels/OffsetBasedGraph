def takes(*targs):
    def dec(func):
        def new_func(*args, **kwargs):
            for arg in zip(args, targs):
                assert isinstance(arg, targs)
            return func(*args, **kwargs)
        return new_func
    return dec
