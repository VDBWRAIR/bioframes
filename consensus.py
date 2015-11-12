from bioframes import bioframes
from Bio import SeqIO
from toolz.functoolz import compose
from toolz.itertoolz import map, zip, second, drop, nth, iterate
import itertools
import sh
''' samtools mpileup -cf ref.fasta hu.bam -g | bcftools view  -'''

#NOTE: freebayes requires ALL reads to be tagged with an RG, which requires a slight change to
# tagreads.py:  https://github.com/VDBWRAIR/ngs_mapper/blob/9523d32effd268543611b60758991a99373a65f5/ngs_mapper/tagreads.py#L56-L59
# (only two degenerate bases in this alignment)
get_degen = {('A', 'C') : 'M', ('C', 'A') : 'M', \
             ('A', 'G') : 'R', ('G', 'A') : 'R'}.__getitem__


AMBIGUITY_TABLE = { 'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', 'N': 'N', 'AC': 'M', 'AG': 'R', 'AT': 'W', 'CG': 'S', 'CT': 'Y', 'GT': 'K', 'ACG': 'V', 'ACT': 'H', 'AGT': 'D', 'CGT': 'B', 'ACGT': 'N' }
get_degen = compose(AMBIGUITY_TABLE.__getitem__, ''.join, sorted)


df = bioframes.load_vcf('freebayes.vcf')
ref_file = 'A12x2520_Houston_Mesoni/KC807175_Mesonivirus.fasta'
fa = SeqIO.parse(ref_file, 'fasta')
bam_file = 'A12x2520_Houston_Mesoni/A12x2520_Houston_Mesoni.bam'
#This times freebayes is only ever reporting one ALT
df.ALT = df.ALT.apply(lambda x: x[0])

# the vcf library reports alts as _Substitution/whatever objects. extract the string.
df.REF, df.ALT = df.REF.apply(str), df.ALT.apply(str)
ambiguous = ((df.AO / df.DP.apply(float)) < 0.8)

#have to use .loc for assignment or else get shallow copy warning
df.loc[ambiguous, 'ALT'] = list(map(get_degen, zip(df.loc[ambiguous].REF, df.loc[ambiguous].ALT)))


original_ref = str(next(fa).seq)

def swap_base(seq_and_offset, info):
    '''Have to keep track of the offset--the number of new bases--
    because insertions push the string back each time.'''
    (seq, offset), (POS, REF, ALT) = seq_and_offset, info
    result = seq[:POS-1+offset] + ALT + seq[POS-1+len(REF)+offset:]
    assert len(ALT) >= len(REF)
    #should add N's or - for deletion for consistency with ngs_mapper
    new_offset = offset + (len(ALT) - len(REF))
    return result, new_offset

df['OFF'] = df.ALT.apply(len) - df.REF.apply(len)
offs, off_pos = df[df.OFF > 0].OFF, df[df.OFF > 0].POS
total_offset = sum(offs)
def partition(pred, seq):
    t1, t2 = itertools.tee(seq)
    return itertools.ifilterfalse(pred, t1), itertools.ifilter(pred, t2)
# remove those N values that are the minimum of those greater than the offset
# position, where N is offset
def drop_gaps_because_insertion_offsets(A, off_pos):
    off, pos = off_pos
    before, after = partition(lambda x: x >= pos, A)
    return itertools.chain.from_iterable([before, drop(off,  map(lambda x: x+off, after))])
    #return drop(off, lazysort(after_off)
    # numpy array A > off_pos

new_ref, off = reduce(swap_base, zip(df.POS, df.REF, df.ALT), (original_ref, 0))
#
#total number of new bases
assert off == 4, off
pileup = sh.samtools('mpileup', bam_file, f=ref_file, _iter=True)
pileup_positions = map(compose(int, second, unicode.split), pileup)
raw_missings =  set(xrange(1, len(original_ref)+1)) - set(pileup_positions)
missings = set(reduce(drop_gaps_because_insertion_offsets, zip(offs, off_pos), raw_missings))
insert_gap = lambda s, x: s[:x]+ '-' + s[x+1:]
_reduced, _ = reduce(swap_base, zip(df.POS, df.REF, df.ALT), (original_ref, 0))
reduced =  reduce(insert_gap, map(lambda x: x - 1, missings),  new_ref)
reduced = nth(off, iterate(lambda s: s[:-1] if s[-1] == '-' else s, reduced))
#using the gaps defined in the original ngs_mapper consensus output
# this is not correct because it fills in where there are now inserts
start = '-------------------'
end = '----------------------'
res = start + new_ref[len(start):-len(end)] + end
assert len(res) == len(original_ref) + 4
#assert res == reduced
#print res
'''
ngs = bioframes.load_vcf('ngs_vcall.vcf')
set(ngs.POS.unique()) - set(fb.POS.unique())
 the above doesn't quite work because freebayes will report multi-base REFs & alts, like

REF: ACG  ALT: ATA

or the like.'''


