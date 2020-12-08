# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.2
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.





from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('Angpow_swig_py', [dirname(__file__)])
        except ImportError:
            import Angpow_swig_py
            return Angpow_swig_py
        if fp is not None:
            try:
                _mod = imp.load_module('Angpow_swig_py', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    Angpow_swig_py = swig_import_helper()
    del swig_import_helper
else:
    import Angpow_swig_py
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


try:
    import weakref
    weakref_proxy = weakref.proxy
except:
    weakref_proxy = lambda x: x


class wallTimer(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, wallTimer, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, wallTimer, name)
    __repr__ = _swig_repr
    def __init__(self): 
        this = Angpow_swig_py.new_wallTimer()
        try: self.this.append(this)
        except: self.this = this
    def start(self): return Angpow_swig_py.wallTimer_start(self)
    def stop(self): return Angpow_swig_py.wallTimer_stop(self)
    __swig_destroy__ = Angpow_swig_py.delete_wallTimer
    __del__ = lambda self : None;
wallTimer_swigregister = Angpow_swig_py.wallTimer_swigregister
wallTimer_swigregister(wallTimer)

class ClassFunc1D(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, ClassFunc1D, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, ClassFunc1D, name)
    __repr__ = _swig_repr
    def __init__(self): 
        this = Angpow_swig_py.new_ClassFunc1D()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = Angpow_swig_py.delete_ClassFunc1D
    __del__ = lambda self : None;
ClassFunc1D_swigregister = Angpow_swig_py.ClassFunc1D_swigregister
ClassFunc1D_swigregister(ClassFunc1D)

class Py_ClassFunc1D(ClassFunc1D):
    __swig_setmethods__ = {}
    for _s in [ClassFunc1D]: __swig_setmethods__.update(getattr(_s,'__swig_setmethods__',{}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, Py_ClassFunc1D, name, value)
    __swig_getmethods__ = {}
    for _s in [ClassFunc1D]: __swig_getmethods__.update(getattr(_s,'__swig_getmethods__',{}))
    __getattr__ = lambda self, name: _swig_getattr(self, Py_ClassFunc1D, name)
    __repr__ = _swig_repr
    def __init__(self): 
        if self.__class__ == Py_ClassFunc1D:
            _self = None
        else:
            _self = self
        this = Angpow_swig_py.new_Py_ClassFunc1D(_self, )
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = Angpow_swig_py.delete_Py_ClassFunc1D
    __del__ = lambda self : None;
    def get_value(self, *args): return Angpow_swig_py.Py_ClassFunc1D_get_value(self, *args)
    def __disown__(self):
        self.this.disown()
        Angpow_swig_py.disown_Py_ClassFunc1D(self)
        return weakref_proxy(self)
Py_ClassFunc1D_swigregister = Angpow_swig_py.Py_ClassFunc1D_swigregister
Py_ClassFunc1D_swigregister(Py_ClassFunc1D)

class CheFunc(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, CheFunc, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, CheFunc, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = Angpow_swig_py.new_CheFunc(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = Angpow_swig_py.delete_CheFunc
    __del__ = lambda self : None;
    def orderFunc(self): return Angpow_swig_py.CheFunc_orderFunc(self)
    def ChebyshevCoeffFFT(self): return Angpow_swig_py.CheFunc_ChebyshevCoeffFFT(self)
    def ChebyshevTransform(self, *args): return Angpow_swig_py.CheFunc_ChebyshevTransform(self, *args)
CheFunc_swigregister = Angpow_swig_py.CheFunc_swigregister
CheFunc_swigregister(CheFunc)

class CheAlgo(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, CheAlgo, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, CheAlgo, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = Angpow_swig_py.new_CheAlgo(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = Angpow_swig_py.delete_CheAlgo
    __del__ = lambda self : None;
    def ClenshawCurtisWeightsFast(self): return Angpow_swig_py.CheAlgo_ClenshawCurtisWeightsFast(self)
    def InverseChebyshevCoeffFFT(self): return Angpow_swig_py.CheAlgo_InverseChebyshevCoeffFFT(self)
    def InverseChebyshevTransform(self): return Angpow_swig_py.CheAlgo_InverseChebyshevTransform(self)
    def ComputeIntegralUnscaled(self): return Angpow_swig_py.CheAlgo_ComputeIntegralUnscaled(self)
CheAlgo_swigregister = Angpow_swig_py.CheAlgo_swigregister
CheAlgo_swigregister(CheAlgo)

class SwigPyIterator(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, SwigPyIterator, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, SwigPyIterator, name)
    def __init__(self, *args, **kwargs): raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = Angpow_swig_py.delete_SwigPyIterator
    __del__ = lambda self : None;
    def value(self): return Angpow_swig_py.SwigPyIterator_value(self)
    def incr(self, n=1): return Angpow_swig_py.SwigPyIterator_incr(self, n)
    def decr(self, n=1): return Angpow_swig_py.SwigPyIterator_decr(self, n)
    def distance(self, *args): return Angpow_swig_py.SwigPyIterator_distance(self, *args)
    def equal(self, *args): return Angpow_swig_py.SwigPyIterator_equal(self, *args)
    def copy(self): return Angpow_swig_py.SwigPyIterator_copy(self)
    def next(self): return Angpow_swig_py.SwigPyIterator_next(self)
    def __next__(self): return Angpow_swig_py.SwigPyIterator___next__(self)
    def previous(self): return Angpow_swig_py.SwigPyIterator_previous(self)
    def advance(self, *args): return Angpow_swig_py.SwigPyIterator_advance(self, *args)
    def __eq__(self, *args): return Angpow_swig_py.SwigPyIterator___eq__(self, *args)
    def __ne__(self, *args): return Angpow_swig_py.SwigPyIterator___ne__(self, *args)
    def __iadd__(self, *args): return Angpow_swig_py.SwigPyIterator___iadd__(self, *args)
    def __isub__(self, *args): return Angpow_swig_py.SwigPyIterator___isub__(self, *args)
    def __add__(self, *args): return Angpow_swig_py.SwigPyIterator___add__(self, *args)
    def __sub__(self, *args): return Angpow_swig_py.SwigPyIterator___sub__(self, *args)
    def __iter__(self): return self
SwigPyIterator_swigregister = Angpow_swig_py.SwigPyIterator_swigregister
SwigPyIterator_swigregister(SwigPyIterator)

class std_vector_CheFunc(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, std_vector_CheFunc, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, std_vector_CheFunc, name)
    __repr__ = _swig_repr
    def iterator(self): return Angpow_swig_py.std_vector_CheFunc_iterator(self)
    def __iter__(self): return self.iterator()
    def __nonzero__(self): return Angpow_swig_py.std_vector_CheFunc___nonzero__(self)
    def __bool__(self): return Angpow_swig_py.std_vector_CheFunc___bool__(self)
    def __len__(self): return Angpow_swig_py.std_vector_CheFunc___len__(self)
    def pop(self): return Angpow_swig_py.std_vector_CheFunc_pop(self)
    def __getslice__(self, *args): return Angpow_swig_py.std_vector_CheFunc___getslice__(self, *args)
    def __setslice__(self, *args): return Angpow_swig_py.std_vector_CheFunc___setslice__(self, *args)
    def __delslice__(self, *args): return Angpow_swig_py.std_vector_CheFunc___delslice__(self, *args)
    def __delitem__(self, *args): return Angpow_swig_py.std_vector_CheFunc___delitem__(self, *args)
    def __getitem__(self, *args): return Angpow_swig_py.std_vector_CheFunc___getitem__(self, *args)
    def __setitem__(self, *args): return Angpow_swig_py.std_vector_CheFunc___setitem__(self, *args)
    def append(self, *args): return Angpow_swig_py.std_vector_CheFunc_append(self, *args)
    def empty(self): return Angpow_swig_py.std_vector_CheFunc_empty(self)
    def size(self): return Angpow_swig_py.std_vector_CheFunc_size(self)
    def clear(self): return Angpow_swig_py.std_vector_CheFunc_clear(self)
    def swap(self, *args): return Angpow_swig_py.std_vector_CheFunc_swap(self, *args)
    def get_allocator(self): return Angpow_swig_py.std_vector_CheFunc_get_allocator(self)
    def begin(self): return Angpow_swig_py.std_vector_CheFunc_begin(self)
    def end(self): return Angpow_swig_py.std_vector_CheFunc_end(self)
    def rbegin(self): return Angpow_swig_py.std_vector_CheFunc_rbegin(self)
    def rend(self): return Angpow_swig_py.std_vector_CheFunc_rend(self)
    def pop_back(self): return Angpow_swig_py.std_vector_CheFunc_pop_back(self)
    def erase(self, *args): return Angpow_swig_py.std_vector_CheFunc_erase(self, *args)
    def __init__(self, *args): 
        this = Angpow_swig_py.new_std_vector_CheFunc(*args)
        try: self.this.append(this)
        except: self.this = this
    def push_back(self, *args): return Angpow_swig_py.std_vector_CheFunc_push_back(self, *args)
    def front(self): return Angpow_swig_py.std_vector_CheFunc_front(self)
    def back(self): return Angpow_swig_py.std_vector_CheFunc_back(self)
    def assign(self, *args): return Angpow_swig_py.std_vector_CheFunc_assign(self, *args)
    def resize(self, *args): return Angpow_swig_py.std_vector_CheFunc_resize(self, *args)
    def insert(self, *args): return Angpow_swig_py.std_vector_CheFunc_insert(self, *args)
    def reserve(self, *args): return Angpow_swig_py.std_vector_CheFunc_reserve(self, *args)
    def capacity(self): return Angpow_swig_py.std_vector_CheFunc_capacity(self)
    __swig_destroy__ = Angpow_swig_py.delete_std_vector_CheFunc
    __del__ = lambda self : None;
std_vector_CheFunc_swigregister = Angpow_swig_py.std_vector_CheFunc_swigregister
std_vector_CheFunc_swigregister(std_vector_CheFunc)

class std_vector_double(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, std_vector_double, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, std_vector_double, name)
    __repr__ = _swig_repr
    def iterator(self): return Angpow_swig_py.std_vector_double_iterator(self)
    def __iter__(self): return self.iterator()
    def __nonzero__(self): return Angpow_swig_py.std_vector_double___nonzero__(self)
    def __bool__(self): return Angpow_swig_py.std_vector_double___bool__(self)
    def __len__(self): return Angpow_swig_py.std_vector_double___len__(self)
    def pop(self): return Angpow_swig_py.std_vector_double_pop(self)
    def __getslice__(self, *args): return Angpow_swig_py.std_vector_double___getslice__(self, *args)
    def __setslice__(self, *args): return Angpow_swig_py.std_vector_double___setslice__(self, *args)
    def __delslice__(self, *args): return Angpow_swig_py.std_vector_double___delslice__(self, *args)
    def __delitem__(self, *args): return Angpow_swig_py.std_vector_double___delitem__(self, *args)
    def __getitem__(self, *args): return Angpow_swig_py.std_vector_double___getitem__(self, *args)
    def __setitem__(self, *args): return Angpow_swig_py.std_vector_double___setitem__(self, *args)
    def append(self, *args): return Angpow_swig_py.std_vector_double_append(self, *args)
    def empty(self): return Angpow_swig_py.std_vector_double_empty(self)
    def size(self): return Angpow_swig_py.std_vector_double_size(self)
    def clear(self): return Angpow_swig_py.std_vector_double_clear(self)
    def swap(self, *args): return Angpow_swig_py.std_vector_double_swap(self, *args)
    def get_allocator(self): return Angpow_swig_py.std_vector_double_get_allocator(self)
    def begin(self): return Angpow_swig_py.std_vector_double_begin(self)
    def end(self): return Angpow_swig_py.std_vector_double_end(self)
    def rbegin(self): return Angpow_swig_py.std_vector_double_rbegin(self)
    def rend(self): return Angpow_swig_py.std_vector_double_rend(self)
    def pop_back(self): return Angpow_swig_py.std_vector_double_pop_back(self)
    def erase(self, *args): return Angpow_swig_py.std_vector_double_erase(self, *args)
    def __init__(self, *args): 
        this = Angpow_swig_py.new_std_vector_double(*args)
        try: self.this.append(this)
        except: self.this = this
    def push_back(self, *args): return Angpow_swig_py.std_vector_double_push_back(self, *args)
    def front(self): return Angpow_swig_py.std_vector_double_front(self)
    def back(self): return Angpow_swig_py.std_vector_double_back(self)
    def assign(self, *args): return Angpow_swig_py.std_vector_double_assign(self, *args)
    def resize(self, *args): return Angpow_swig_py.std_vector_double_resize(self, *args)
    def insert(self, *args): return Angpow_swig_py.std_vector_double_insert(self, *args)
    def reserve(self, *args): return Angpow_swig_py.std_vector_double_reserve(self, *args)
    def capacity(self): return Angpow_swig_py.std_vector_double_capacity(self)
    __swig_destroy__ = Angpow_swig_py.delete_std_vector_double
    __del__ = lambda self : None;
std_vector_double_swigregister = Angpow_swig_py.std_vector_double_swigregister
std_vector_double_swigregister(std_vector_double)

# This file is compatible with both classic and new-style classes.


