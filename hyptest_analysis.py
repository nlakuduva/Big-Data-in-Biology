#!/Usr/bin/env python3 
# ^^ implements portability because it tells to look in all paths 

import pandas as pd
from Bio import SeqIO
#output from stdout of a command line passed in as a string and returns as string  
from subprocess import check_output
#adjusts a string to be read as a file
from io import StringIO
from scipy import stats
import numpy as np

hyptest = pd.read_csv("grcm38_hyp_tests.tsv", header=0, index_col=0, sep='\t')
DEG = hyptest[hyptest["t_pv"]>0]
nonDEG = hyptest[hyptest["t_pv"]<=0]

fasta_DEG_IDs = []
fasta_nonDEG_IDs = []
with open("grcm38promoters_res_-9000_1000.fa", "rU") as fasta_handle:
    for record in SeqIO.parse(fasta_handle, 'fasta'):
        if record.id in DEG.index:
            fasta_DEG_IDs.append(record)
        if record.id in nonDEG.index:
            fasta_nonDEG_IDs.append(record)
    
with open("DEG.fa", "w") as deg_handle:
    SeqIO.write(fasta_DEG_IDs, deg_handle, "fasta")

with open("nonDEG.fa","w") as non_deg_handle:
    SeqIO.write(fasta_nonDEG_IDs, non_deg_handle, "fasta")

# Create dataframes out of the outputs of nucfreqs.R with k=1
nfDEG = lambda k: check_output("nucfreqs.R "+str(k)+" DEG.fa", shell=True).decode('UTF-8')
nfDEG_dataframe= lambda k: pd.read_csv(StringIO(nfDEG(k)), sep="\t", index_col=0, header=0)

nfnonDEG = lambda k: check_output("nucfreqs.R "+str(k)+" nonDEG.fa", shell=True).decode('UTF-8')
nfnonDEG_dataframe= lambda k: pd.read_csv(StringIO(nfnonDEG(k)), sep="\t", index_col=0, header=0)


def create_df(k):
    vals=nfDEG_dataframe(k).columns.values
    df_DEG=nfDEG_dataframe(k)
    df_nonDEG=nfnonDEG_dataframe(k)
    df=[]
    for val in vals:
        tres=stats.ttest_ind(df_DEG[val], df_nonDEG[val])
        temp=[tres[0],tres[1],np.mean(df_DEG[val]),np.std(df_DEG[val],ddof=1),np.mean(df_nonDEG[val]),np.std(df_nonDEG[val],ddof=1)]
        df.append(temp)
    
    df = pd.DataFrame(df, columns=['tscore','pval','meanDEG','stdDEG','meannonDEG','stdnonDEG'], index=vals)
    return df        
    
    
    
#tuple with t-score, p-value
#tresults = stats.ttest_ind(nf1DEG_dataframe['A'], nf1nonDEG_dataframe['A'])
#print tresults
 
bedtools_output = check_output("bedtools intersect -wa -wb -a ../GSE63137/atac-seq/GSE63137_ATAC-seq_PV_neurons_HOMER_peaks.bed -b ../mouse_genome/tss_grcm38_res.bed", shell=True).decode('UTF-8')
bedtools_dataframe = pd.read_csv(StringIO(bedtools_output), sep="\t", index_col=0, header=0)
