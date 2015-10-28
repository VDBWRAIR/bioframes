import re
import pandas as pd
import numpy as np
from itertools import groupby
from func import pmap, psplit, pstrip, compose, compose_all, merge_dicts, fzip, partial2, dictmap, starcompose
from operator import itemgetter
from functools import partial
import operator as op
from operator import add, div
from schema import Schema, Use
from itertools import ifilter
# Parse options
#from pyparsing import Regex

to_np_int = partial(np.array, dtype=int)
parse_array = compose_all(to_np_int, psplit(','), pstrip('[]'))
tabsplit = psplit('\t')


minus33 = partial(add, -33)
qual_int_sanger = compose(minus33, ord)
qual_to_phreds = compose(to_np_int, pmap(qual_int_sanger))
error = compose(partial(pow, 10), partial2(div, -10.0))
#don't need to map because numpy vectorizes it automatically
#TODO: handle non-sanger version
sanger_qual_str_to_error = compose(error, qual_to_phreds)

basic_scheme={
    'QNAME' : str,
    'FLAG' : int,
    'RNAME' : str,
    'POS' : int,
    'MAPQ' : int,
    'CIGAR' : str,
    'RNEXT' : str,
    'PNEXT' : int,
    'TLEN' : int,
    #'MRNM' : str,
    #'MRNM' : '*='.__contains__,
    #'MPOS' : int,
    #'ISIZE' : int,
    'SEQ' : str,
    'QUAL' : str,
    #'OPTIONS' : str
}
#['OPT', 'MPOS', 'MRNM', 'ISIZE']
basic_schema = Schema(dictmap(Use, basic_scheme))

options_scheme = {
    'A' : chr,
    'i' : int,
    'f' : float,
    'Z' : str,
    'H' : int, # hex
    'B' : parse_array
}

def parse_option(option_str):
    ''' alternatively:   tag, type, val = option_str.split('\t')  '''
#_join = partial(reduce, lambda a, b: a+':'+b)
#    tag_type_val = op.itemgetter(0, 2, 4)
#    _tag = Regex(r'[A-Za-z][A-Za-z0-9]')
#    _type = Regex(r'[AifZHB]')
#    _value = Regex('[^\t]')
#    full = _tag + ':' + _type + ':' + _value
    #reduce(operator.add, [tag, _type, value], ':')
#    parsed_list = full.parseString(option_str)
#    return tag_type_val(parsed_list)

    tag, _type, raw_val = psplit(':')(option_str)
    val = options_scheme[_type](raw_val)
    return tag, val

    #full = _join( [tag, _type, value ] )
    #parse_array = re.compile(r'[cCsSiIf](,[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)+').match
    ''' NOTE: samfiles use ASCII of Phred-scaled base QUALity+33 '''


#TODO: handle empty cases (unmapped reads, *)

def parse_cigar(cigar_str):
    #? makes the regex not be too greedy
    cigar_regex = r'(?:([0-9]+)([MIDNSHPX=]))+?'
    reg = re.compile(cigar_regex)
    tups = reg.findall(cigar_str)
    key, value = itemgetter(1), itemgetter(0)
    groups = groupby(sorted(tups, key=key), key)
    get_counts = pmap(compose(int, itemgetter(0)))
    sum_counts = compose(sum, get_counts)
    s = "cigar_{0}".format
    cigar_dict = dict( (s(name), sum_counts(nums)) for name, nums in groups)
    #print cigar_dict
    mismatches = sum(num for k, num in cigar_dict.items() if k not in ['cigar_M', 'cigar_='])
    return merge_dicts(cigar_dict, {'cigar_score': mismatches})

#dictmap(compose(sum, get_counts), dict(groups))

#dict(map(reverse, tups))
''' assert sum(itemgetter('M', 'I', 'S', '=', 'X')) == len(seq) == len(quality), \
    "cigar string M/I/S/=/X should sum to the length of the query sequence." '''

index = ['QNAME', 'POS', 'RNAME']
#TODO:
''' POS starts at 1 . . . but if POS is set as 0 for an unmapped read without coordinate. If POS is 0, no assumptions can be made about RNAME and CIGAR
Bit 0x4 is the only reliable place to tell whether the read is unmapped. If 0x4 is set, no
assumptions can be made about RNAME , POS , CIGAR , MAPQ , and bits 0x2, 0x100, and 0x800 .'''


flag_meanings = {
0x1  : "template having multiple segments in sequencing",
0x2  : "each segment properly aligned according to the aligner",
0x4  : "segment unmapped",
0x8  : "next segment in the template unmapped",
0x10 : "SEQ being reverse complemented",
0x20 : "SEQ of the next segment in the template being reversed",
0x40 : "the rst segment in the template",
0x80 : "the last segment in the template",
0x100: "secondary alignment",
0x200: "not passing quality controls",
0x400: "PCR or optical duplicate",
0x800: "supplementary alignment"
}


eval_flag = compose(bool, op.and_)

def flag_dict(flag):
    return dict((meaning, eval_flag(bit, flag)) for bit, meaning in flag_meanings.items())
def split_list(A, idx):
    return A[:idx], A[idx:]

sam_columns = ("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL") #optiosn


#TODO: get_record function takes a filehandle and returns a single record via SeqIO, etc.
#So functions expect a dictionary I guess
#pass
parse_options = compose(dict, pmap(parse_option)) #, tabsplit)
#readfields = compose(tabsplit, next)
line_to_dict = compose_all(dict, partial(zip, sam_columns)) #, tabsplit)
validated_dict = compose(basic_schema.validate, line_to_dict)
fields_and_options = compose(partial2(split_list, len(sam_columns)), tabsplit)
parsers = partial(fzip, [validated_dict, parse_options])
parse_fields_and_options = compose(parsers, fields_and_options)
all_but_cigar_dict = starcompose(merge_dicts, parse_fields_and_options)
get_cigar_dict = compose(parse_cigar, itemgetter('CIGAR'))
get_flag_dict = compose(flag_dict, itemgetter('FLAG'))
get_error = compose(sanger_qual_str_to_error, itemgetter('QUAL'))

def load_sam(fh):
    dicts = map(get_row, ifilter(bool, fh.read().split('\n')))
    return pd.DataFrame(dicts)
#TODO: do we really need indices? it complicates querying i tlooks like maybe where plays better with them
# .set_index(index) #, index=index, columns=columns)


def get_row(row):
    result = all_but_cigar_dict(row)
    return merge_dicts(result, get_cigar_dict(result), get_flag_dict(result), {'error' :  get_error(result)})


