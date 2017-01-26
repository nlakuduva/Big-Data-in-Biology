#!/usr/bin/env python

import Bio
from Bio import SeqIO
import numpy as np
import os
import pandas as pd
import re

bp = os.environ['WORK'] + '/brain_promoter'

group = pd.Series({
    'SRR1647854' : 'excitatory',
    'SRR1647855' : 'excitatory',
    'SRR1647856' : 'PV',
    'SRR1647857' : 'PV',
    'SRR1647858' : 'VIP',
    'SRR1647859' : 'VIP'
})

x = pd.read_csv(
    bp + "/GSE63137/kallisto_results/tpm.tsv.gz",
    sep = "\t",
    index_col = 0,
    header = 0
)
transcriptAnnot = x.ix[:, x.dtypes == "object"]
x = np.log2(x.ix[:, group.index] + 1)
x0 = x.copy()

xRowMeans = x.mean(axis=1)
xRowSds = x.std(axis=1)

pd.crosstab(xRowMeans > xRowMeans.median(),
            xRowSds > xRowSds.median())
# col_0  False  True 
# row_0              
# False  37644   9396
# True    9396  37644

countsFile = bp + '/GSE63137/kallisto_results/est_counts.tsv.gz'
counts = pd.read_csv(countsFile, sep="\t", index_col=0, header=0)
counts = counts.ix[:, group.index].sum(axis=1)

x = x.ix[(xRowMeans > xRowMeans.median()) &
         (xRowSds > xRowSds.median()) &
         (counts > 100)]

# srrConversion = pd.read_csv(
#     bp + "/GSE63137/SrrToSrx.tsv",
#     sep = "\t",
#     index_col = 0,
#     header = 0
# ).ix[:, 0]
# if x.columns.isin(srrConversion.index).all():
#     x.columns = srrConversion.ix[x.columns]
# annot = pd.read_csv(
#     bp + "/GSE63137-GPL17021/gse63137-gpl17021_annot.tsv",
#     sep = "\t",
#     index_col = None,
#     header = 0,
#     nrows = 8
# )
# annot.index = annot["!Sample_supplementary_file_3"].str.replace(
#     r"^.*/",
#     ""
# )

# group = annot.ix[:, 0].str.replace(r"RNA-seq_(.*)_rep\d+", r"\1")
# group = group.str.replace(r"_neurons", "")


# ## -----------------------------------------------------------------
# ## keep excitatory, PV, and VIP samples
# ## -----------------------------------------------------------------
# samplesToKeep = group[
#         group.isin(["excitatory", "PV", "VIP"])].index.to_series()
# x = x.ix[:, samplesToKeep]
# annot = annot.ix[samplesToKeep]
# group = group[samplesToKeep]


## -----------------------------------------------------------------
## load transcript annotation
## -----------------------------------------------------------------
x0sym = transcriptAnnot.ix[x0.index, "gene_name"]
xsym = transcriptAnnot.ix[x.index, "gene_name"]


## -----------------------------------------------------------------
## select highest expressing transcript for each gene
## -----------------------------------------------------------------
trxMeans = pd.DataFrame({
    "trx" : x.index,
    "symbol" : xsym,
    "mean" : xRowMeans.ix[x.index]
}).ix[:, ['trx', 'symbol', 'mean']]
bestTrx = trxMeans[['symbol', 'mean']].groupby('symbol')['mean'].idxmax()

# trxMaxes = trxMeans[['symbol', 'mean']].groupby('symbol').agg({'mean' : max})
# trxMaxes.describe()
# trxMeans.ix[bestTrx].describe()

xg = x.ix[bestTrx]
xg.index = bestTrx.index

bestTrxToSymbol = pd.Series(bestTrx.index, index=bestTrx)


## -----------------------------------------------------------------
## load promoter properties and locations
## -----------------------------------------------------------------
# seqs = SeqIO.parse(
#     open(bp + "/mouse_genome/grcm38promoters_res_3000_3000.fa"),
#     "fasta"
# )
# seqs = {re.sub(r"^\d+\s+", "", seq.description) : seq
#         for seq in seqs}

gtf = pd.read_csv(
    bp + "/mouse_genome/grcm38exon1s.gtf",
    sep = "\t",
    index_col = None,
    header = None
)
gtf.columns = ['chr', '5p', '3p',
               'gene_id', 'transcript_id', 'name', 'orientation']
gtf['start'] = int(-1)
gtf.ix[gtf['orientation'] == '+', 'start'] = \
        gtf.ix[gtf['orientation']=='+', '5p']
gtf.ix[gtf['orientation'] == '-', 'start'] = \
        gtf.ix[gtf['orientation'] == '-', '3p']
