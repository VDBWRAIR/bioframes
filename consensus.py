from __future__ import print_function
from itertools import imap, izip
from bioframes import bioframes
from Bio import SeqIO
import sys
from toolz.functoolz import compose
from toolz.itertoolz import map, zip, second, drop, nth, iterate, first
from functools import partial
import itertools
import sh
''' samtools mpileup -cf ref.fasta hu.bam -g | bcftools view  -'''

'''
use this to create consensus, and bioframes.py to create the VCF comparison
#NOTE: freebayes requires ALL reads to be tagged with an RG, which requires a slight change to
# tagreads.py:  https://github.com/VDBWRAIR/ngs_mapper/blob/9523d32effd268543611b60758991a99373a65f5/ngs_mapper/tagreads.py#L56-L59
'''
''' how to handle file globbing and multiple dirs at once? And: a compiled report as final task target?'''
AMBIGUITY_TABLE = { 'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', 'N': 'N', 'AC': 'M', 'AG': 'R', 'AT': 'W', 'CG': 'S', 'CT': 'Y', 'GT': 'K', 'ACG': 'V', 'ACT': 'H', 'AGT': 'D', 'CGT': 'B', 'ACGT': 'N' }
get_degen = compose(AMBIGUITY_TABLE.__getitem__, ''.join, sorted)
insert_gap = lambda s, x: s[:x]+ '-' + s[x+1:]
from operator import methodcaller as call, attrgetter
def make_dict(classes):
    return dict(zip(map(call('__name__'), classes, classes)))
#TODO: fix ambiguous base definition
#TODO: this doesn't work for flu sequences which have multiple sections, \
#    and the  VCFs have multiple references. need to group by reference
def fix_fb_df(df):
    #Freebayes only ever reports one ALT?
    df.ALT = df.ALT.apply(lambda x: x[0])
    # the vcf library reports alts as _Substitution/whatever objects. extract the string.
    df.REF, df.ALT = df.REF.apply(str), df.ALT.apply(str)
    '''#TODO: this re-definition of ambiguous bases translates mult-base
    sections (e.g. AC) into single base alts or something'''
    ambiguous = ((df.AO / df.DP.apply(float)) < 0.8)
    #have to use .loc for assignment or else get shallow copy warning
    #NOTE: temporarily removed the amiguous base business.
    #df.loc[ambiguous, 'ALT'] = list(map(get_degen, zip(df.loc[ambiguous].REF, df.loc[ambiguous].ALT)))
    df['OFF'] = df.ALT.apply(len) - df.REF.apply(len)
    return df

pluck_attr = lambda a, A: map(attrgetter(a), A)

#string_to_fasta = '>FreebayesConseunsus\n'.__add__
string_to_fasta = '>Freebayes {0}\n{1}'.format
def make_consensus(bam_file, ref_file, freebayes_vcf):
    ''':retrurn str'''
    refs = list(SeqIO.parse(ref_file, 'fasta'))
    df = bioframes.load_vcf(freebayes_vcf)
    ids, raw_segment_dfs = zip(*df.groupby('CHROM'))
    segment_dfs = map(fix_fb_df, raw_segment_dfs)
    refs = sorted(refs, key=lambda x: x.id)
    assert list(pluck_attr('id', refs)) == sorted(ids)
    ref_seqs = (str(s.seq) for s in refs) #pluck_attr('seq', refs)
    segment_dfs = sorted(segment_dfs, key=lambda x: x.CHROM.iloc[0])
    def process(original_ref, df):
         new_ref, off = reduce(swap_base, zip(df.POS, df.REF, df.ALT), (original_ref, 0))
         #currently only handle inserts, > 0
         offs, off_pos = df[df.OFF > 0].OFF, df[df.OFF > 0].POS
         pileup_positions = zero_coverage_positions(bam_file, ref_file)
         return gap_fill_ref(original_ref, new_ref, pileup_positions, offs, off_pos, off)
    #return map(process, ref_seqs, segment_dfs)
    return zip(sorted(ids), map(process, ref_seqs, segment_dfs))

def gap_fill_ref(original_ref, new_ref, pileup_positions, offs, off_pos, off):
    raw_missings =  set(xrange(1, len(original_ref)+1)) - set(pileup_positions)
    missings = reduce(drop_gaps_because_insertion_offsets, izip(offs, off_pos), raw_missings)
    gap_filled = reduce(insert_gap, imap(lambda x: x - 1, missings),  new_ref)
    return nth(off, iterate(lambda s: s[:-1] if s[-1] == '-' else s, gap_filled))

def zero_coverage_positions(bam_file, ref_file):
    pileup = sh.samtools('mpileup', bam_file, f=ref_file, _iter=True)
    return map(compose(int, second, unicode.split), pileup)

def swap_base(seq_and_offset, info):
    '''Have to keep track of the offset--the number of new bases--
    because insertions push the string back each time.'''
    (seq, offset), (POS, REF, ALT) = seq_and_offset, info
    result = seq[:POS-1+offset] + ALT + seq[POS-1+len(REF)+offset:]
    assert len(ALT) >= len(REF), "length of alt and ref not the same at position %s %s %s" % (POS, ALT, REF)
    #should add N's or - for deletion for consistency with ngs_mapper
    new_offset = offset + (len(ALT) - len(REF))
    return result, new_offset

def drop_gaps_because_insertion_offsets(A, off_pos):
    '''remove those N values that are the minimum of those greater than the offset
     position, where N is offset'''
    off, pos = off_pos
    before, after = partition(lambda x: x >= pos, A)
    return itertools.chain.from_iterable([before, drop(off,  imap(lambda x: x+off, after))])

''' set(ngs.POS.unique()) - set(fb.POS.unique())
 the above doesn't quite work because freebayes will report multi-base REFs & alts, like REF: ACG  ALT: ATA or the like.'''

def partition(pred, seq):
    t1, t2 = itertools.tee(seq)
    return itertools.ifilterfalse(pred, t1), itertools.ifilter(pred, t2)


#consensus_str = compose(string_to_fasta, make_consensus)
consensus_str = compose('\n'.join, partial(itertools.starmap, string_to_fasta), make_consensus)
main = compose(print,consensus_str)
#expects bam_file, ref_file, freebayes_vcf):

if __name__ == '__main__':
    extension = lambda x: x.split('.')[-1]
    bam, ref, vcf = sorted(sys.argv[1:], key=extension)
    main(bam, ref, vcf)
