# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 20:40:29 2021

@author: slaiad
"""

#%%data
# DB_PATH = 'Database/fasta/uniprot_homo_sapiens.fasta'
DB_PATH = 'Database/fasta/uniprot_homo_canonical_isoform.fasta'
# uniprot_homo_canonical_isoform.fasta
SPEC_PATH = 'Dataset/PXD022999/191122_MK_SIO13_P2-GM1.mzXML'
# SPEC_PATH = 'Dataset/PXD022999/01650b_BD2-TUM_first_pool_54_01_01-2xIT_2xHCD-1h-R2.mzXML'
#%%
MIN_LENGTH_OF_PEP = 6
MAX_LENGTH_OF_PEP = 60
MIN_HK_SPEC_LEN = 2
EXTRACT_AA_TOL = 0.01
NUM_CLUSTERS_1ST_LAYER = 7
MISSED_CLEAVAGES = 2
TAG_SCORE_THRES = 1.6

AA_SCORE_THRES = 0.5

# MASS_TOL = 0.05#Da
MASS_TOL = 10*1e-6#ppm

MAX_TAGS_NUM = 100