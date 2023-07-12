import pandas as pd

df = pd.read_csv('xpehh_all.out.norm', delimiter='\t')

top = df[df['NORM'].le(df['NORM'].quantile(0.01))]

data = [top['CHR'], top['BP'], top['NORM']]

headers =['CHR', 'BP', 'NORM']

bed = pd.concat(data, axis=1, keys=headers)

bed['pos'] = bed['BP']

cols = bed.columns.tolist()
cols = ['CHR', 'BP', 'pos', 'NORM']
bed = bed[cols]

bed.to_csv('top_snp_XPEHH_lowlanders.bed', header=None, index=None, sep='\t', mode='a')
