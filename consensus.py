from bioframes import bioframes
from Bio import SeqIO
from toolz.functoolz import compose
from toolz.itertoolz import map, zip, second, drop, nth, iterate
import itertools
import sh
''' samtools mpileup -cf ref.fasta hu.bam -g | bcftools view  -'''


'''
#NOTE: freebayes requires ALL reads to be tagged with an RG, which requires a slight change to
# tagreads.py:  https://github.com/VDBWRAIR/ngs_mapper/blob/9523d32effd268543611b60758991a99373a65f5/ngs_mapper/tagreads.py#L56-L59
# (only two degenerate bases in this alignment)
plumbum?
# use luigi.ExternalTask

bamfile -> tagged bam

taggedbam, reffile
-> freebayevcf

freebayesvcf
-> freebayes-consensus

freebayes-consesnsus, ngs-consensus
-> comparison report

freeabyesvcf, ngs-vcf
-> vcf comparison report
'''
AMBIGUITY_TABLE = { 'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', 'N': 'N', 'AC': 'M', 'AG': 'R', 'AT': 'W', 'CG': 'S', 'CT': 'Y', 'GT': 'K', 'ACG': 'V', 'ACT': 'H', 'AGT': 'D', 'CGT': 'B', 'ACGT': 'N' }
get_degen = compose(AMBIGUITY_TABLE.__getitem__, ''.join, sorted)
insert_gap = lambda s, x: s[:x]+ '-' + s[x+1:]

def fix_fb_df(df):
    #Freebayes only ever reports one ALT?
    df.ALT = df.ALT.apply(lambda x: x[0])
    # the vcf library reports alts as _Substitution/whatever objects. extract the string.
    df.REF, df.ALT = df.REF.apply(str), df.ALT.apply(str)
    ambiguous = ((df.AO / df.DP.apply(float)) < 0.8)
    #have to use .loc for assignment or else get shallow copy warning
    df.loc[ambiguous, 'ALT'] = list(map(get_degen, zip(df.loc[ambiguous].REF, df.loc[ambiguous].ALT)))
    df['OFF'] = df.ALT.apply(len) - df.REF.apply(len)
    return df

def make_consensus(bam_file, ref_file, freebayes_vcf):
    ''':retrurn str'''
    fa = SeqIO.parse(ref_file, 'fasta')
    original_ref = str(next(fa).seq)
    df = fix_fb_df(bioframes.load_vcf(freebayes_vcf))
    new_ref, off = reduce(swap_base, zip(df.POS, df.REF, df.ALT), (original_ref, 0))
    #currently only handle inserts, > 0
    offs, off_pos = df[df.OFF > 0].OFF, df[df.OFF > 0].POS
    pileup_positions = zero_coverage_positions(bam_file, ref_file)
    return gap_fill_ref(original_ref, new_ref, pileup_positions, offs, off_pos, off)

def gap_fill_ref(original_ref, new_ref, pileup_positions, offs, off_pos, off):
    raw_missings =  set(xrange(1, len(original_ref)+1)) - set(pileup_positions)
    missings = reduce(drop_gaps_because_insertion_offsets, zip(offs, off_pos), raw_missings)
    gap_filled = reduce(insert_gap, map(lambda x: x - 1, missings),  new_ref)
    return nth(off, iterate(lambda s: s[:-1] if s[-1] == '-' else s, gap_filled))

def zero_coverage_positions(bam_file, ref_file):
    pileup = sh.samtools('mpileup', bam_file, f=ref_file, _iter=True)
    return map(compose(int, second, unicode.split), pileup)

def swap_base(seq_and_offset, info):
    '''Have to keep track of the offset--the number of new bases--
    because insertions push the string back each time.'''
    (seq, offset), (POS, REF, ALT) = seq_and_offset, info
    result = seq[:POS-1+offset] + ALT + seq[POS-1+len(REF)+offset:]
    assert len(ALT) >= len(REF)
    #should add N's or - for deletion for consistency with ngs_mapper
    new_offset = offset + (len(ALT) - len(REF))
    return result, new_offset

def drop_gaps_because_insertion_offsets(A, off_pos):
    '''remove those N values that are the minimum of those greater than the offset
     position, where N is offset'''
    off, pos = off_pos
    before, after = partition(lambda x: x >= pos, A)
    return itertools.chain.from_iterable([before, drop(off,  map(lambda x: x+off, after))])

''' set(ngs.POS.unique()) - set(fb.POS.unique())
 the above doesn't quite work because freebayes will report multi-base REFs & alts, like REF: ACG  ALT: ATA or the like.'''

def partition(pred, seq):
    t1, t2 = itertools.tee(seq)
    return itertools.ifilterfalse(pred, t1), itertools.ifilter(pred, t2)
