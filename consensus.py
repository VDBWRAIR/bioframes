from bioframes import bioframes
from Bio import SeqIO

# (only two degenerate bases in this alignment)
get_degen = {('A', 'C') : 'M', ('C', 'A') : 'M', \
             ('A', 'G') : 'R', ('G', 'A') : 'R'}.__getitem__

df = bioframes.load_vcf('freebayes.vcf')
fa = SeqIO.parse('A12x2520_Houston_Mesoni/KC807175_Mesonivirus.fasta', 'fasta')

#This times freebayes is only ever reporting one ALT
df.ALT = df.ALT.apply(lambda x: x[0])

# the vcf library reports alts as _Substitution/whatever objects. extract the string.
df.REF, df.ALT = df.REF.apply(str), df.ALT.apply(str)
ambiguous = ((df.AO / df.DP.apply(float)) < 0.8)

#have to use .loc for assignment or else get shallow copy warning
df.loc[ambiguous, 'ALT'] = map(get_degen, zip(df.loc[ambiguous].REF, df.loc[ambiguous].ALT))


original_ref = str(next(fa).seq)

def swap_base(seq_and_offset, tup):
    '''Have to keep track of the offset--the number of new bases--
    because insertions push the string back each time.'''
    seq, offset = seq_and_offset
    POS, REF, ALT = tup
    result = seq[:POS-1+offset] + ALT + seq[POS-1+len(REF)+offset:]
    assert len(ALT) >= len(REF)
    new_offset = offset + (len(ALT) - len(REF))
    return result, new_offset

new_ref, off = reduce(swap_base, zip(df.POS, df.REF, df.ALT), (original_ref, 0))
#total number of new bases
assert off == 4, off

#using the gaps defined in the original ngs_mapper consensus output
start = '-------------------'
end = '----------------------'
res = start + new_ref[len(start):-len(end)] + end
assert len(res) == len(original_ref) + 4
print res
'''
ngs = bioframes.load_vcf('ngs_vcall.vcf')
set(ngs.POS.unique()) - set(fb.POS.unique())
 the above doesn't quite work because freebayes will report multi-base REFs & alts, like

REF: ACG  ALT: ATA

or the like.'''


