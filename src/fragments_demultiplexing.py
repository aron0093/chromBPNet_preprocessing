import os

import random
import string

import numpy as np
import pandas as pd

from utils import write_chunk, create_chunk_generator

from joblib import Parallel, delayed
from tqdm.auto import tqdm

# Make pseudo-replicates from fragment df
def assign_insertion_sites(fragments_df):

    '''
    Splits fragment file into two random subsets of Tn5 sites.
    '''

    start_df, end_df = fragments_df.copy(), fragments_df.copy()

    start_df['end'] = start_df['start']+1
    start_df = start_df.sample(frac=1)

    end_df['start'] = start_df['end']-1
    end_df = end_df.sample(frac=1)

    pseudo_rep_1 = pd.concat([start_df.iloc[:int(start_df.shape[0]/2)], 
                              end_df.iloc[int(end_df.shape[0]/2):]])
    pseudo_rep_2 = pd.concat([start_df.iloc[int(start_df.shape[0]/2):], 
                              end_df.iloc[:int(end_df.shape[0]/2)]])

    return pseudo_rep_1, pseudo_rep_2

# Write out fragment file chunks
def _demux_fragments(chunk, 
                     barcodes, 
                     allowed_chrs,
                     output_loc,
                     make_pseudoreps):
    
    # Gen random string
    rand_string = ''.join(random.choice(string.ascii_uppercase)\
                    for i in range(10))

    chunk = chunk.loc[chunk.chr.isin(allowed_chrs)]

    if barcodes is not None:
        chunk = chunk.loc[chunk.barcodes.isin(barcodes)]
    write_chunk(chunk, os.path.join(output_loc, rand_string+'.txt'))
    
    if make_pseudoreps:
        pr1, pr2 = assign_insertion_sites(chunk)
        write_chunk(pr1, os.path.join(output_loc, rand_string+'_pseudorep1.txt'))
        write_chunk(pr2, os.path.join(output_loc, rand_string+'_pseudorep2.txt'))
    
    return rand_string

# Process fragment file - parallel over chunks
def demux_fragments(fragment_fil, 
                    barcodes=None,
                    output_loc='./',
                    chunksize=10000,
                    skiprows=None,
                    n_jobs=-1,
                    make_pseudoreps=True):

    '''
    Write out fragment file chunks. Useful when demuxing frag files.
    This function will produce multiple files that will be collated.

    ARGS
        fragment_fil: 
            tab separated fragment file with atleast 4 columns.
            chr, start, end, barcodes <- no headers are expected.
        barcodes: (default: None)
            Array-like of barcodes to count. If None then all.
        output_loc: os.path (default: './')
            Directory where processed chunks will be stored.
            Defaults to current directory.
        chunksize: int (default: 10000)
            Number of rows to read at once.
        skiprows: int (default: None)
            Number of rows to skip from top.
        n_jobs: int (default: -1)
            Number of cores to use.   
        make_pseudoreps: bool (default:True)
            Make pseudoreplicates by fragment splitting.  
    '''  
           
    # Check if output dir exists
    os.makedirs(output_loc, exist_ok=True)

    # Allowed chromsomes
    allowed_chrs = ["chr1", "chr10", "chr11", 
                    "chr12", "chr13", "chr14", 
                    "chr15", "chr16", "chr17", 
                    "chr18", "chr19", "chr2", 
                    "chr20", "chr21", "chr22", 
                    "chr3", "chr4", "chr5", 
                    "chr6", "chr7", "chr8", 
                    "chr9", "chrX", "chrY"]

    # Create chunk gen
    chunk_gen = create_chunk_generator(fragment_fil, 
                                       chunksize=chunksize, 
                                       skiprows=skiprows)

    # Process fragment file  
    track_chunks = []
    track_chunks.append(Parallel(n_jobs=n_jobs, 
                                 backend='threading')(delayed(_demux_fragments)(chunk, 
                                                                                barcodes,
                                                                                allowed_chrs,
                                                                                output_loc,
                                                                                make_pseudoreps,
                                                                                ) \
                                                      for chunk in tqdm(chunk_gen, 
                                                                        unit='chunks', 
                                                                        desc='Chunking fragments'
                                                                        )))

    return track_chunks[0]


       

