import pandas as pd
import subprocess
import os

CHROM = 18

def get_df(chrom=CHROM):
    # Load the population reference file
    df = pd.read_csv(f'tsv/gnomad_chr{chrom}_filtered.tsv', sep='\t', 
                     names=['pos', 'ref', 'alt', 'af_eu', 'an_eu', 'af_me', 'an_me'],
                     na_values='.')

    # Cleanup population data
    df = df.dropna(subset=['af_me', 'af_eu'])
    df = df[(df['ref'].str.len() == 1) & (df['alt'].str.len() == 1)]
    df = df[df['an_me'] > 50]
    df['delta'] = (df['af_eu'] - df['af_me']).abs()
    df = df[df['delta'] > 0.1]  # only sites with >10 pp diff btwn ME/EU

    # Handle AJ subject data
    subject_file = f'tsv/hg002_chr{chrom}_genotypes.tsv'
    sub = pd.read_csv(subject_file, sep='\t', names=['pos', 'ref', 'alt', 'gt'])

    # Merge ref and subject
    df = pd.merge(df, sub[['pos', 'gt']], on='pos', how='left')
    
    df['gt'] = df['gt'].fillna('0/0')
    # Dosage | 0/0 -> 0, 0/1 -> 1, 1/1 -> 2
    df['dosage'] = df['gt'].apply(lambda x: str(x).count('1'))
    
    # needed for transition probability scaling
    df['dist'] = df['pos'].diff().fillna(100)

    return df
