from __future__ import print_function
from functools import partial
from Bio import SeqIO
import numpy as np
from operator import attrgetter as attr
from operator import add, div
from func import partial2,  pmap, pifilter,  \
       ilen, merge_dicts,  dictzip, \
    apply_to_object, dictmap
from toolz import compose
import pandas as pd
import vcf
#TODO:

def check_np_type(_type_str):
    def inner(array):
        return array.dtype == np.dtype(_type_str)
    return inner

'''
Join fastq and SAM (merge on QNAME [and SEQ])
Join VCF and SAM (merge on POS)
Pileup
Join VCF and Pileup
'''

get_fastq = partial(SeqIO.parse, format='fastq')
get_fasta = partial(SeqIO.parse, format='fasta')
to_np_int = partial(np.array, dtype=int)
gccontent = compose(ilen, pifilter('GC'.__contains__))

minus33 = partial(add, -33)
qual_int_sanger = compose(minus33, ord)

''' Error = 10^-(Phred/10) '''
qual_to_phreds = compose(to_np_int, pmap(qual_int_sanger))
error = compose(partial(pow, 10), partial2(div, -10.0))
#don't need to map because numpy vectorizes it automatically
#TODO: handle non-sanger version
sanger_qual_str_to_error = compose(error, qual_to_phreds)
#SANGER_OFFSET = 33

'''
assert len(quality) == len(error) == len(phred_scores)
'''


#validate = scheme.validate
#TODO: could make these validations match samtools spec
#TODO: Could treat options/cigar string as their own class with their own parsing and validation.

def flatten_vcf(record):
    '''
    :param VCFRecord record: VCFRecord object
    :return dict: all fields as a flat dictionary, including those fields in rec.INFO
    '''
    _ids = ['CHROM', 'REF', 'QUAL', 'POS', 'ALT', 'FILTER', 'FORMAT', 'ID'] #'INFO'
    fields = [getattr(record, _id) for _id in _ids]
    d = dict( (_id, value) for  _id, value in zip(_ids, fields) )
    d.update(dict((_id, flatten_list(field)) for _id, field in record.INFO.items()))
    return d

def collection_as_df(lambdas, columns, collection):
    '''
    Create a pandas DataFrame by applying a series of functions to a collection.
    :param list lambdas: a list of functions which take exactly one argument (an objects in collection)
    :param list columns: (str) the column names, in order with lambdas
    :param list collection: a list of objects which are valid arguments for the lambdas.
    :return pandas.DataFrame the lambda results on the collection objects as a Matrix.
    '''
    assert len(lambdas) == len(columns), "lambdas must have same length as columns"
    '''use list here to force the evaluation of the functions. otherwise the lambda grabs the last obj evaluated from collection, as in a closure.'''
    values = (list( func(obj) for func in lambdas) for obj in collection)
    return pd.DataFrame(values, columns=columns)

def load_vcf_list(vcf_records):
    '''
    Convert a list of vcf Records to a pandas DataFrame.
    '''
    return pd.DataFrame(flatten_vcf(rec) for rec in vcf_records)

flatten_list = lambda a: a if type(a) != list else a[0]
load_vcf  = compose(load_vcf_list, vcf.Reader, open)

#TODO: properly handle [None], '-', etc.

#have one return the true/false series and the other do getitem i guess
def col_compare(df, col, value, comp):
    half = partial(comp, value)
    boolean = compose(half, df.__getitem__)
    return boolean(col)

#return_frame = partial(compose, df.__getitem__)
#simple_compare = return_frame(col_compare)

def ambiguous(df):
    return df[~df.CB.isin(list('AGCT'))]

def vcalls(df):
    return df[df.REF != df.CB]

def df_by_attrs(columns, collection):
    attr_getters = map(attr, columns)
    return collection_as_df(attr_getters, columns, collection)

#optionally dict(zip(columns, apply_each))
#TODO: fix this, why doesn't repeat work?
def makeframe(biodata):
    dicts = []
    while True:
        try:
            dicts.append(_create_dict(**biodata))
        except StopIteration:
            break
    return pd.DataFrame(dicts)

#make_dicts = pmap(biodata_to_row)
#make_fastq_frame = kstarcompose(makeframe, sequenceframes.fqframe)

def obj_to_dict(obj, names_getters_map):
    apply_to_obj = partial2(apply_to_object, obj)
    return dictmap(apply_to_obj, names_getters_map)

def _create_dict(obj_func, columns, getters, validator, dictgetters):
    obj = apply(obj_func)
    pre_dict = obj_to_dict(obj, dictzip(columns, getters))
    #optional intermediate schema here
    if dictgetters:
        extra_dict = obj_to_dict(pre_dict, dictgetters)
        full_dict = merge_dicts(pre_dict, extra_dict)
    else:
        full_dict = pre_dict
    return validator.validate(full_dict)
#    compose_all(merge_dicts, obj_to_dict, partial2(repeat, 2),
    #do = starcompose(obj_to_dict, partial2(repeat, 2), obj_to_dict)
    #return do(obj, dictzip(columns, getters)


#rewrite
#unpack_biodata = itemgetter('obj_func', 'columns', 'getters',  'validator','dictgetters')
#biodata_to_row = starcompose(_create_row, unpack_biodata)
