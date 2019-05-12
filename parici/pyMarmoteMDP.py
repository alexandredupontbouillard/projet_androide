# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.12
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
if _swig_python_version_info >= (2, 7, 0):
    def swig_import_helper():
        import importlib
        pkg = __name__.rpartition('.')[0]
        mname = '.'.join((pkg, '_pyMarmoteMDP')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_pyMarmoteMDP')
    _pyMarmoteMDP = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_pyMarmoteMDP', [dirname(__file__)])
        except ImportError:
            import _pyMarmoteMDP
            return _pyMarmoteMDP
        try:
            _mod = imp.load_module('_pyMarmoteMDP', fp, pathname, description)
        finally:
            if fp is not None:
                fp.close()
        return _mod
    _pyMarmoteMDP = swig_import_helper()
    del swig_import_helper
else:
    import _pyMarmoteMDP
del _swig_python_version_info

try:
    _swig_property = property
except NameError:
    pass  # Python < 2.2 doesn't have 'property'.

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_setattr_nondynamic(self, class_type, name, value, static=1):
    if (name == "thisown"):
        return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name, None)
    if method:
        return method(self, value)
    if (not static):
        if _newclass:
            object.__setattr__(self, name, value)
        else:
            self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)


def _swig_setattr(self, class_type, name, value):
    return _swig_setattr_nondynamic(self, class_type, name, value, 0)


def _swig_getattr(self, class_type, name):
    if (name == "thisown"):
        return self.this.own()
    method = class_type.__swig_getmethods__.get(name, None)
    if method:
        return method(self)
    raise AttributeError("'%s' object has no attribute '%s'" % (class_type.__name__, name))


def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except __builtin__.Exception:
    class _object:
        pass
    _newclass = 0

class SwigPyIterator(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, SwigPyIterator, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, SwigPyIterator, name)

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _pyMarmoteMDP.delete_SwigPyIterator
    __del__ = lambda self: None

    def value(self):
        return _pyMarmoteMDP.SwigPyIterator_value(self)

    def incr(self, n=1):
        return _pyMarmoteMDP.SwigPyIterator_incr(self, n)

    def decr(self, n=1):
        return _pyMarmoteMDP.SwigPyIterator_decr(self, n)

    def distance(self, x):
        return _pyMarmoteMDP.SwigPyIterator_distance(self, x)

    def equal(self, x):
        return _pyMarmoteMDP.SwigPyIterator_equal(self, x)

    def copy(self):
        return _pyMarmoteMDP.SwigPyIterator_copy(self)

    def next(self):
        return _pyMarmoteMDP.SwigPyIterator_next(self)

    def __next__(self):
        return _pyMarmoteMDP.SwigPyIterator___next__(self)

    def previous(self):
        return _pyMarmoteMDP.SwigPyIterator_previous(self)

    def advance(self, n):
        return _pyMarmoteMDP.SwigPyIterator_advance(self, n)

    def __eq__(self, x):
        return _pyMarmoteMDP.SwigPyIterator___eq__(self, x)

    def __ne__(self, x):
        return _pyMarmoteMDP.SwigPyIterator___ne__(self, x)

    def __iadd__(self, n):
        return _pyMarmoteMDP.SwigPyIterator___iadd__(self, n)

    def __isub__(self, n):
        return _pyMarmoteMDP.SwigPyIterator___isub__(self, n)

    def __add__(self, n):
        return _pyMarmoteMDP.SwigPyIterator___add__(self, n)

    def __sub__(self, *args):
        return _pyMarmoteMDP.SwigPyIterator___sub__(self, *args)
    def __iter__(self):
        return self
SwigPyIterator_swigregister = _pyMarmoteMDP.SwigPyIterator_swigregister
SwigPyIterator_swigregister(SwigPyIterator)

class solutionMDP(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, solutionMDP, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, solutionMDP, name)
    __repr__ = _swig_repr

    def __init__(self):
        this = _pyMarmoteMDP.new_solutionMDP()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _pyMarmoteMDP.delete_solutionMDP
    __del__ = lambda self: None

    def writeSolution(self):
        return _pyMarmoteMDP.solutionMDP_writeSolution(self)

    def setSize(self, s):
        return _pyMarmoteMDP.solutionMDP_setSize(self, s)
solutionMDP_swigregister = _pyMarmoteMDP.solutionMDP_swigregister
solutionMDP_swigregister(solutionMDP)

class feedbackSolutionMDP(solutionMDP):
    __swig_setmethods__ = {}
    for _s in [solutionMDP]:
        __swig_setmethods__.update(getattr(_s, '__swig_setmethods__', {}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, feedbackSolutionMDP, name, value)
    __swig_getmethods__ = {}
    for _s in [solutionMDP]:
        __swig_getmethods__.update(getattr(_s, '__swig_getmethods__', {}))
    __getattr__ = lambda self, name: _swig_getattr(self, feedbackSolutionMDP, name)
    __repr__ = _swig_repr

    def __init__(self):
        this = _pyMarmoteMDP.new_feedbackSolutionMDP()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _pyMarmoteMDP.delete_feedbackSolutionMDP
    __del__ = lambda self: None

    def setAction(self, a):
        return _pyMarmoteMDP.feedbackSolutionMDP_setAction(self, a)

    def setValue(self, t):
        return _pyMarmoteMDP.feedbackSolutionMDP_setValue(self, t)

    def setActionIndex(self, indice, value):
        return _pyMarmoteMDP.feedbackSolutionMDP_setActionIndex(self, indice, value)

    def getActionIndex(self, indice):
        return _pyMarmoteMDP.feedbackSolutionMDP_getActionIndex(self, indice)

    def getValueIndex(self, indice):
        return _pyMarmoteMDP.feedbackSolutionMDP_getValueIndex(self, indice)

    def writeSolution(self):
        return _pyMarmoteMDP.feedbackSolutionMDP_writeSolution(self)
feedbackSolutionMDP_swigregister = _pyMarmoteMDP.feedbackSolutionMDP_swigregister
feedbackSolutionMDP_swigregister(feedbackSolutionMDP)

class marmoteSet(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, marmoteSet, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, marmoteSet, name)
    __repr__ = _swig_repr
    UNION = _pyMarmoteMDP.marmoteSet_UNION
    PRODUCT = _pyMarmoteMDP.marmoteSet_PRODUCT
    SIMPLE = _pyMarmoteMDP.marmoteSet_SIMPLE

    def __init__(self, *args):
        this = _pyMarmoteMDP.new_marmoteSet(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _pyMarmoteMDP.delete_marmoteSet
    __del__ = lambda self: None

    def cardinal(self):
        return _pyMarmoteMDP.marmoteSet_cardinal(self)

    def isFinite(self):
        return _pyMarmoteMDP.marmoteSet_isFinite(self)

    def isSimple(self):
        return _pyMarmoteMDP.marmoteSet_isSimple(self)

    def isUnion(self):
        return _pyMarmoteMDP.marmoteSet_isUnion(self)

    def isProduct(self):
        return _pyMarmoteMDP.marmoteSet_isProduct(self)

    def totNbDims(self):
        return _pyMarmoteMDP.marmoteSet_totNbDims(self)

    def enumerate(self):
        return _pyMarmoteMDP.marmoteSet_enumerate(self)

    def test_index_decode(self):
        return _pyMarmoteMDP.marmoteSet_test_index_decode(self)

    def firstState(self, buffer):
        return _pyMarmoteMDP.marmoteSet_firstState(self, buffer)

    def nextState(self, buffer):
        return _pyMarmoteMDP.marmoteSet_nextState(self, buffer)

    def decodeState(self, index, buffer):
        return _pyMarmoteMDP.marmoteSet_decodeState(self, index, buffer)

    def index(self, buffer):
        return _pyMarmoteMDP.marmoteSet_index(self, buffer)

    def isZero(self, buffer):
        return _pyMarmoteMDP.marmoteSet_isZero(self, buffer)

    def printState(self, *args):
        return _pyMarmoteMDP.marmoteSet_printState(self, *args)
marmoteSet_swigregister = _pyMarmoteMDP.marmoteSet_swigregister
marmoteSet_swigregister(marmoteSet)

class marmoteInterval(marmoteSet):
    __swig_setmethods__ = {}
    for _s in [marmoteSet]:
        __swig_setmethods__.update(getattr(_s, '__swig_setmethods__', {}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, marmoteInterval, name, value)
    __swig_getmethods__ = {}
    for _s in [marmoteSet]:
        __swig_getmethods__.update(getattr(_s, '__swig_getmethods__', {}))
    __getattr__ = lambda self, name: _swig_getattr(self, marmoteInterval, name)
    __repr__ = _swig_repr

    def __init__(self, min, max):
        this = _pyMarmoteMDP.new_marmoteInterval(min, max)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _pyMarmoteMDP.delete_marmoteInterval
    __del__ = lambda self: None

    def isFinite(self):
        return _pyMarmoteMDP.marmoteInterval_isFinite(self)

    def isZero(self, buffer):
        return _pyMarmoteMDP.marmoteInterval_isZero(self, buffer)

    def firstState(self, buffer):
        return _pyMarmoteMDP.marmoteInterval_firstState(self, buffer)

    def nextState(self, buffer):
        return _pyMarmoteMDP.marmoteInterval_nextState(self, buffer)

    def decodeState(self, index, buffer):
        return _pyMarmoteMDP.marmoteInterval_decodeState(self, index, buffer)

    def index(self, buffer):
        return _pyMarmoteMDP.marmoteInterval_index(self, buffer)

    def printState(self, out, buffer):
        return _pyMarmoteMDP.marmoteInterval_printState(self, out, buffer)

    def enumerate(self):
        return _pyMarmoteMDP.marmoteInterval_enumerate(self)
marmoteInterval_swigregister = _pyMarmoteMDP.marmoteInterval_swigregister
marmoteInterval_swigregister(marmoteInterval)

class SCC(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, SCC, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, SCC, name)
    __repr__ = _swig_repr
    __swig_setmethods__["id"] = _pyMarmoteMDP.SCC_id_set
    __swig_getmethods__["id"] = _pyMarmoteMDP.SCC_id_get
    if _newclass:
        id = _swig_property(_pyMarmoteMDP.SCC_id_get, _pyMarmoteMDP.SCC_id_set)
    __swig_setmethods__["period"] = _pyMarmoteMDP.SCC_period_set
    __swig_getmethods__["period"] = _pyMarmoteMDP.SCC_period_get
    if _newclass:
        period = _swig_property(_pyMarmoteMDP.SCC_period_get, _pyMarmoteMDP.SCC_period_set)
    __swig_setmethods__["states"] = _pyMarmoteMDP.SCC_states_set
    __swig_getmethods__["states"] = _pyMarmoteMDP.SCC_states_get
    if _newclass:
        states = _swig_property(_pyMarmoteMDP.SCC_states_get, _pyMarmoteMDP.SCC_states_set)

    def __lt__(self, a):
        return _pyMarmoteMDP.SCC___lt__(self, a)

    def __init__(self):
        this = _pyMarmoteMDP.new_SCC()
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _pyMarmoteMDP.delete_SCC
    __del__ = lambda self: None
SCC_swigregister = _pyMarmoteMDP.SCC_swigregister
SCC_swigregister(SCC)

class sparseMatrix(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, sparseMatrix, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, sparseMatrix, name)
    __repr__ = _swig_repr

    def __init__(self, *args):
        this = _pyMarmoteMDP.new_sparseMatrix(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _pyMarmoteMDP.delete_sparseMatrix
    __del__ = lambda self: None

    def setEntry(self, row, col, val):
        return _pyMarmoteMDP.sparseMatrix_setEntry(self, row, col, val)

    def getEntry(self, arg2, arg3):
        return _pyMarmoteMDP.sparseMatrix_getEntry(self, arg2, arg3)

    def getNbElts(self, row):
        return _pyMarmoteMDP.sparseMatrix_getNbElts(self, row)

    def getCol(self, row, numCol):
        return _pyMarmoteMDP.sparseMatrix_getCol(self, row, numCol)

    def getEntryByCol(self, row, numCol):
        return _pyMarmoteMDP.sparseMatrix_getEntryByCol(self, row, numCol)

    def getTransDistrib(self, row):
        return _pyMarmoteMDP.sparseMatrix_getTransDistrib(self, row)

    def rowSum(self, row):
        return _pyMarmoteMDP.sparseMatrix_rowSum(self, row)

    def evaluateMeasure(self, *args):
        return _pyMarmoteMDP.sparseMatrix_evaluateMeasure(self, *args)

    def evaluateValue(self, v, res):
        return _pyMarmoteMDP.sparseMatrix_evaluateValue(self, v, res)

    def evaluateValueState(self, v, stateIndex):
        return _pyMarmoteMDP.sparseMatrix_evaluateValueState(self, v, stateIndex)

    def copy(self):
        return _pyMarmoteMDP.sparseMatrix_copy(self)

    def uniformize(self):
        return _pyMarmoteMDP.sparseMatrix_uniformize(self)

    def embed(self):
        return _pyMarmoteMDP.sparseMatrix_embed(self)

    def diagnose(self, out):
        return _pyMarmoteMDP.sparseMatrix_diagnose(self, out)

    def write(self, out, format):
        return _pyMarmoteMDP.sparseMatrix_write(self, out, format)

    def addToEntry(self, row, col, val):
        return _pyMarmoteMDP.sparseMatrix_addToEntry(self, row, col, val)

    def normalize(self):
        return _pyMarmoteMDP.sparseMatrix_normalize(self)
sparseMatrix_swigregister = _pyMarmoteMDP.sparseMatrix_swigregister
sparseMatrix_swigregister(sparseMatrix)

class totalRewardMDP(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, totalRewardMDP, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, totalRewardMDP, name)
    __repr__ = _swig_repr

    def __init__(self, *args):
        this = _pyMarmoteMDP.new_totalRewardMDP(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _pyMarmoteMDP.delete_totalRewardMDP
    __del__ = lambda self: None

    def writeMDP(self):
        return _pyMarmoteMDP.totalRewardMDP_writeMDP(self)

    def valueIteration(self, epsilon, maxIter):
        return _pyMarmoteMDP.totalRewardMDP_valueIteration(self, epsilon, maxIter)

    def valueIterationGS(self, epsilon, maxIter):
        return _pyMarmoteMDP.totalRewardMDP_valueIterationGS(self, epsilon, maxIter)

    def policyIteration(self, maxIter):
        return _pyMarmoteMDP.totalRewardMDP_policyIteration(self, maxIter)

    def policyIterationModified(self, epsilon, maxIter, delta, maxInIter):
        return _pyMarmoteMDP.totalRewardMDP_policyIterationModified(self, epsilon, maxIter, delta, maxInIter)

    def policyCost(self, policy, epsilon, maxIter):
        return _pyMarmoteMDP.totalRewardMDP_policyCost(self, policy, epsilon, maxIter)
totalRewardMDP_swigregister = _pyMarmoteMDP.totalRewardMDP_swigregister
totalRewardMDP_swigregister(totalRewardMDP)

# This file is compatible with both classic and new-style classes.

