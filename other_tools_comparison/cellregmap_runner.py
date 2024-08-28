# to run
# /usr/bin/time -o /directflow/SCCGGroupShare/projects/anncuo/OneK1K/saige_eqtl/cellregmap_comparison/logs/RPL23A.cis.chr17.cellregmap.logs.txt -v /share/ScratchGeneral/anncuo/jupyter/conda_notebooks/envs/cellregmap_notebook/bin/python cellregmap_runner.py /directflow/SCCGGroupShare/projects/anncuo/OneK1K/saige_eqtl/cellregmap_comparison/input_files/RPL23A.cis.chr17.input.txt /directflow/SCCGGroupShare/projects/anncuo/OneK1K/saige_eqtl/cellregmap_comparison/pvals/
# runner at (local copy of): https://github.com/annacuomo/SAIGE_QTL_analyses/other_tools_comparison/run_cellregmap.qsub

import sys

import itertools
import pandas as pd
from numpy import array, split, cumsum, zeros, append

from cellregmap import run_association_fast

input_file = str(sys.argv[1])
output_path = str(sys.argv[2])

def get_groups_from_smf(smf_df):
    n_samples = smf_df.shape[0]
    donors = smf_df['individual'].unique()
    n_donors = len(donors)
    n_cells = array([],dtype=int)
    for donor in donors:
        n_cells = append(n_cells, array(smf_df[smf_df['individual']==donor].shape[0], dtype=int))
    groups = split(range(n_samples), cumsum(n_cells))[:-1]
    return groups

def get_block_hK_from_groups(groups):
    n_samples = len(list(itertools.chain.from_iterable(groups)))
    hM = zeros((n_samples, len(groups)))
    for i, idx in enumerate(groups):
        hM[idx, i] = 1.0
    return hM

# get gene from filename
gene = input_file.split('/')[-1].split('.')[0]

# combined file for each gene including gene, covs, snps
input_df = pd.read_csv(input_file, sep='\t')

# expression of relevant gene
y = input_df['gene_sct']

# sample covariates
W = input_df[['age','sex','pc1','pc2','pc3','pc4','pc5','pc6','pf1','pf2']]

# cell contexts
C = input_df[['pf1','pf2']]

# genotypes
G = input_df.values[:,17:input_df.shape[1]]

# sample mapping file (cells to donors)
input_df['cell'] = input_df.index
smf_df = input_df[['individual','cell']]

# kinship (blocks)
groups = get_groups_from_smf(smf_df)
hK = get_block_hK_from_groups(groups)

# run
pvals = run_association_fast(y=y, W=W.values, E=C.values, G=G, hK=hK)[0]

# save
variants = input_df.columns[17:(input_df.shape[1]-1)]
pvals_df = pd.DataFrame({'variants': variants, 'pvals': pvals})
output_file = f'{output_path}/{gene}_pvals.csv'
pvals_df.to_csv(output_file)
