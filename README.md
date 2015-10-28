
Create a histogram of map positions:
```python 
from bioframes import samframe
samframe.load_sam(open('780.sam')) 
df.RNAME.unique() 
dfr = df[df.RNAME == r] 
%pylab inline
dfr.hist(column='POS')
```

```python
df.MAPQ.apply(np.mean)

# find the Reference with the highest average mapping quality
# MAPQ is stored internally as a numpy array
max(ifilterfalse(lambda x: '*' == x[0], df.groupby(['RNAME'])), key=lambda x: x[1].MAPQ.apply(np.mean).mean())[0]
#reference with greatest number of mapped reads
ref_groups = ifilterfalse(lambda x: '*' == x[0], df.groupby(['RNAME']))
max(ref_groups, key=lambda x: len(x[1]))[0]
```

available fields:
```python
In [58]: df.columns
Out[58]: 
Index([u'AS', u'CIGAR', u'FLAG', u'MAPQ', u'NM', u'PCR or optical duplicate',
       u'PNEXT', u'POS', u'QNAME', u'QUAL', u'RG', u'RNAME', u'RNEXT', u'SA',
       u'SEQ', u'SEQ being reverse complemented',
       u'SEQ of the next segment in the template being reversed', u'TLEN',
       u'XS', u'cigar_D', u'cigar_H', u'cigar_I', u'cigar_M', u'cigar_S',
       u'cigar_score',
       #Flags:
       u'each segment properly aligned according to the aligner', u'error',
       u'next segment in the template unmapped',
       u'not passing quality controls', u'secondary alignment',
       u'segment unmapped', u'supplementary alignment',
       u'template having multiple segments in sequencing',
       u'the rst segment in the template',
       u'the last segment in the template'],
      dtype='object')
```

``` python
len(df.groupby(['QNAME']))
Out[62]: 15417

len(df)
Out[63]: 35289

sum(df.QNAME.duplicated())
Out[64]: 19872

len(df.groupby(['QNAME'])) - len(df)
Out[65]: -19872
```
