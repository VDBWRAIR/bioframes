from functools import partial
import itertools as it
import string
import sys
from collections import namedtuple

PY3 = sys.version[0] == '3'
imap, ifilter, izip = (map, filter, zip) if PY3 else (it.imap, it.ifilter, it.izip)

#notin = compose(_not, operator.methodcaller('__contains__'))
#notin = compose(_not, attr('__contains__'))
#mismatches = pfilter(notin('M='))
def merge_dicts(*dict_args):
    '''
    from http://stackoverflow.com/a/26853961
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    '''
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result

def _not(x):
    return not x

def partial2(method, param):
      def t(x):
              return method(x, param)
      return t

def _id(x): return x

def apply_key_func(k, v, funcdict):
    return funcdict.get(k, _id)(v)
def compose_all(*funcs):
    return reduce(compose, funcs)

#def k_compose(outer, **okwargs):
#    ''' compose(f, g)(x) == f(g(x)) '''
#    def newfunc(*args, **ikwargs):
#        _kwargs = dict( (k, apply_key_func(k, v, okwargs)) for k, v in ikwargs.items())
#        return outer(*args, **_kwargs)
#    return newfunc

def compose(outer, inner):
    ''' compose(f, g)(x) == f(g(x)) '''
    def newfunc(*args, **kwargs):
        return outer(inner(*args, **kwargs))
    return newfunc


starcompose2 = lambda f, g: lambda x: f(*g(x))

#starcompose = partial(reduce, starcompose2) #need splat
def starcompose(*funcs):
    return reduce(starcompose2, funcs)

def compose(outer, inner):
    ''' compose(f, g)(x) == f(g(x)) '''
    def newfunc(*args, **kwargs):
        return outer(inner(*args, **kwargs))
    return newfunc
def fzip(funcs, args):
    for func, arg in izip(funcs, args):
        yield func(arg)

def dictmap(func, _dict):
    return dict( (key, func(val)) for key, val in _dict.items())


def ilen(iterable):
    return sum(1 for _ in iterable)

def reverse(collection): return collection[::-1]
pifilter = partial(partial, ifilter)
compose_list = partial(reduce, compose)
compose_all = compose(compose_list, lambda *a: a)
pmap = partial(partial, map)
pfilter = partial(partial, filter)
#TODO: could use partial2 instead
pstrip = lambda x: partial(string.split, chars=x)
psplit = lambda x: partial(string.split, sep=x)
boolint = lambda x: 1 if x else 0
dictzip = compose(dict, zip)
#ilen = compose(sum, pmap(boolint))
#Given a list of functons and names, return the result of those functions dictzipped witht the names.

#TODO:
''' dictfilter '''
def apply_each(funcs, arg):
    return fzip(funcs, it.repeat(arg))

import inspect
import types
def is_local(object):
   return isinstance(object, types.FunctionType) and object.__module__ == __name__
    #use inspect.isfunction
def get_funcs():
    return inspect.getmembers(sys.modules[__name__], \
                             predicate = lambda f: inspect.isfunction(f) and f.__module__ == __name__)
    #return inspect.getmembers(sys.modules[__name__], predicate=is_local)
    #return dict( ((name, func)) for name, func in locals().items() if is_local(name))
#      for key, value in locals().items():
#          if callable(value) and value.__module__ == __name__:
#              l.append(key)


'''
compose columns + object + getters => dict
- unzip

have fqframe return a dict of functions, excluding get_row, openframe; instead passing it to a function which
arranges the getters, applies them to a get_object function, and creates an intermediate dictionary.
This function allows for optional extras, like samframe & fastq (rather than fasta
'''


unzip = starcompose(zip, _id)

def nameddict(Name, _dict):
    ''' dict to named tuple '''
    names, values = unzip(_dict.items())
    return namedtuple(Name, names)(*values)

ppartial = partial(partial)
apply_to_object = compose(apply, ppartial)

kstarcompose2 = lambda f, g: lambda x: f(**g(x))
def kstarcompose(*funcs):
    return reduce(kstarcompose2, funcs)

#kstarcompose = partial(reduce, kstarcompose2)
