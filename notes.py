
from bioframes import samframe
samframe.load_sam(open('780.sam')) 
df.RNAME.unique() 
dfr = df[df.RNAME == r] 
%pylab inline
dfr.hist(column='POS')


from bioframes import samframe
samframe.load_sam(open('780.sam')) 
df.RNAME.unique() 
dfr = df[df.RNAME == r] 
%pylab inline
dfr.hist(column='POS')


from bioframes import samframe
samframe.load_sam(open('780.sam')) 
df = df.RNAME.unique() 
dfr = df[df.RNAME == r] 
%pylab inline
dfr.hist(column='POS')
df = samframe.load_sam(open('780.sam')) 
df.RNAME.unique() 
dfr = df[df.RNAME == r] 
%pylab inline
dfr.hist(column='POS')
ax = dfr.hist?
%pylab inline
df.hist(column='POS')
df.hist(column='POS', grid=False)
df.hist(column='POS', grid=True)
r = df.RNAME.unique()[0]
dfr = df[df.RNAME == r]
dfr.hist(column='POS', grid=True)
dfr.POS.unique()
len(dfr.POS.unique())
dfr.POS.cumsum()
map(lambda x: len(x[1]), drf.POS)
map(lambda x: len(x[1]), dfr.POS)
map(lambda x: len(x[1]), dfr.POS.groupby())
map(lambda x: len(x[1]), dfr.groupby(['POS']))
map(lambda x: len(x[0], x[1]), dfr.groupby(['POS']))
map(lambda x: (x[0], len(x[1])), dfr.groupby(['POS']))
nums = map(lambda x: (x[0], len(x[1])), dfr.groupby(['POS']))
pd.DataFrame(nums)
import pandas as pd
pd.DataFrame(nums)
from matplotlib import pyplot as plt
dfn = pd.DataFrame(nums)
dfn.plot(x=0, y=1)
nums = max(lambda x: (x[0], len(x[1])), dfr.groupby(['POS']), key=lambda x: x[1])
nums = max(map(lambda x: (x[0], len(x[1])), dfr.groupby(['POS'])), key=lambda x: x[1])
numbs
nums
dfr.hist(column='POS', grid=True)
for ref, group in df.groupby(['REF']):
    group.hist(column='POS', label=ref)
for ref, group in df.groupby(['RNAME']):
    group.hist(column='POS', label=ref)
for ref, group in df.groupby(['RNAME']):
    print ref
for ref, group in df.groupby(['RNAME']):
    group.hist(column='POS', label=ref)
for ref, group in df.groupby(['RNAME']):
    group.hist(column='POS', xlabel=ref)
for ref, group in df.groupby(['RNAME']):
    group.hist(column='POS', label="ASFA")
for ref, group in df.groupby(['RNAME']):
    group.hist(label="ASF", column='POS')
for ref, group in df.groupby(['RNAME']):
    group.hist(label=None, column='POS')
for ref, group in df.groupby(['RNAME']):
    ax = group.hist(label=None, column='POS')
type(ax)
ax.show()
df.hist(by=['RNAME'])
df.POS.hist(by=['RNAME'])
df['POS'].hist(by=['RNAME'])
df.hist(by=['RNAME'])
df.hist(by=['RNAME'], x='POS')
df.hist(by=['POS'])
df.hist(column='POS', by=['RNAME'])
ax = df.hist(column='POS', by=['RNAME'])
%pyalab inline off
%pylab inline off
%pylab inline False
%pylab inline
%pylab inline
%pylab qt
ax = df.hist(column='POS', by=['RNAME'])
df.hist?
f.hist(column='POS', by=['RNAME'])
df.hist(column='POS', by=['RNAME'])
plt.savefig('foo.png')
import pylab
df.hist(column='POS', by=['RNAME'])
pylab.savefig('foo.png')
df.hist(column='POS', by=['RNAME'])
history
