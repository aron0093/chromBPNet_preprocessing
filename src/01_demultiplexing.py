import os
import gzip
import random
import string

import numpy as np
import pandas as pd

from joblib import Parallel, delayed
from tqdm.auto import tqdm

import logging
import warnings

import argparse

# Read barcodes per cell grouping
def read_barcode_files(*barcode_fils):

    '''
    Read in barcode files per cell group w. no headers.
    '''

    barcodes = []
    for barcode_fil in tqdm(barcode_fils, 
                            unit='groups', 
                            desc='Reading barcode files'):
        barcodes_ = pd.read_csv(barcode_fil, 
                                header=None, 
                                names=['barcodes'])
        barcodes_['groups'] = os.path.basename(barcode_fil)
        barcodes.append(barcodes_)
    barcodes_df = pd.concat(barcodes)

    return barcodes_df

# Util functions
def _write_chunk(frag_df, output_loc, nam):

    frag_df.to_csv(os.path.join(output_loc, 
                                '{}.txt'.format(nam)),
                   sep='\t', header=False, index=False)

def _cat_files(fils, output_loc, output_nam):

    out_path = os.path.join(output_loc, output_nam)

    with open(out_path,'w') as f:
        for chunk in fils:
            with open(os.path.join(output_loc, chunk)) as infile:
                f.write(infile.read())
    f.close()

# Function to check fragment file formatting
def _check_fragment_file_formatting(fragment_fil):

    '''
    Return first line idx without header text.
    '''

    try: 
        gzip.open(fragment_fil,'rt').read(1)
        fil = gzip.open(fragment_fil,'rt')
    except: fil = open(fragment_fil, 'r')

    for idx, line in enumerate(fil):
        if line[0] == '#':
            continue
        else:
            logging.info('Skipping lines {} since they look like header text'.format(idx))
            return idx

# Make pseudo-replicates from fragment df
def _assign_insertion_sites(frag_df):

    '''
    Splits fragment file into two random subsets of Tn5 sites.
    '''

    start_df, end_df = frag_df.copy(), frag_df.copy()

    start_df['end'] = start_df['start']+1
    start_df = start_df.sample(frac=1)

    end_df['start'] = start_df['end']-1
    end_df = end_df.sample(frac=1)

    pseudo_rep_1 = pd.concat([start_df.iloc[:int(start_df.shape[0]/2)], 
                              end_df.iloc[int(end_df.shape[0]/2):]])
    pseudo_rep_2 = pd.concat([start_df.iloc[int(start_df.shape[0]/2):], 
                              end_df.iloc[:int(end_df.shape[0]/2)]])

    return pseudo_rep_1, pseudo_rep_2

# Function to demultiplex - create temp chunks
def _demultiplex_fragment_file(frag_df,
                               barcodes_df,
                               frag_nam, 
                               output_loc):
    
    '''
    Internal function for fragment file demultiplexing.
    '''

    # Check if output dir exists
    os.makedirs(output_loc, exist_ok=True)

    # Merge barcodes files to annotate group
    frag_df = frag_df.merge(barcodes_df, 
                            how='left', 
                            left_on='barcodes', 
                            right_on='barcodes')

    # Output group-wise files for each chunk
    for group in frag_df.groups.unique():

        # Output fragments
        frag_df_ = frag_df.loc[frag_df.groups==group, 
                               ['chr', 'start', 'end', 'barcodes']]

        # Skip empty chunks
        if frag_df_.shape[0]==0:
            continue

        # Add random string to chunk since outputs are not in order
        rand_string = ''.join(random.choice(string.ascii_uppercase)\
                      for i in range(10))
        _write_chunk(frag_df_, output_loc, '{}_{}_{}'.format(group, frag_nam, rand_string))

        # Output pseudo-reps
        pseudo_rep_1, pseudo_rep_2 = _assign_insertion_sites(frag_df_)

        _write_chunk(pseudo_rep_1, output_loc, '{}_pseudorep1_{}_{}'.format(group, frag_nam, rand_string))
        _write_chunk(pseudo_rep_1, output_loc, '{}_pseudorep2_{}_{}'.format(group, frag_nam, rand_string))
            
# Demultiplex fragment file - parallel over chunks
def demultiplex_fragment_file(fragment_fil, 
                              barcodes_df, 
                              output_loc=None,
                              chunk_size=10000,
                              n_jobs=-1):

    '''
    Demultiplex fragment file into cell-type groups.
    This function will produce multiple files
    depending on chunk_size which have to be collated
    per group.

    ARGS
        fragment_fil: 
            tab separated fragment file with atleast 4 columns.
            chr, start, end, barcodes <- no headers are expected.
        barcodes_df: pandas.DataFrame
            DataFrame with two columns.
            barcodes, groups
        output_loc: os.path (default: None)
            Directory where demultiplexed chunks will be stored.
            By default outputs are stored alongside input fragment files.
        chunk_size: int
            Number of rows to read. Adjust based on available memory.
        n_jobs: int (default: -1)
            Number of cores to use.    
    '''  
    
    # Read fragment file 
    start_row = _check_fragment_file_formatting(fragment_fil)

    # Create output params
    frag_nam = os.path.basename(fragment_fil)  
    if output_loc is None:
        output_loc = os.path.dirname(fragment_fil)
    if len(os.listdir(output_loc)) > 0:
        logging.warning('Remove outputs from any previous runs to avoid duplications.')
        
    # Make chunk generator
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=pd.errors.ParserWarning)
        chunk_gen = pd.read_csv(fragment_fil,
                                index_col=False,
                                header=None,
                                names=['chr', 'start', 'end', 'barcodes'],
                                sep='\t', 
                                compression='infer',
                                chunksize=chunk_size, 
                                skiprows=start_row)

        # Demultiplex fragment file  
        Parallel(n_jobs=n_jobs, 
                backend='threading')(delayed(_demultiplex_fragment_file)(chunk, 
                                                                        barcodes_df,
                                                                        frag_nam, 
                                                                        output_loc)\
                                    for chunk in tqdm(chunk_gen, 
                                                        unit='chunks', 
                                                        desc='Demultiplexing {}'\
                                                        .format(frag_nam)))

# Function to collate temp chunks
def collate_fragment_file(fil_loc, 
                          group_labels=None,
                          remove_chunks=True):

    '''
    Demultiplex fragment file into cell-groups.
    This function will produce multiple files
    depending on chunk_size which have to be collated
    per group.

    ARGS
        fil_loc: os.path
            Directory where demultiplexed chunks are stored.
        group_labels: list (default: None)
            Cell-groups to collate. By default all groups are collated.
        remove_chunks: bool (default: True)
            Remove chunks after collation.
    ''' 

    # Collate by group
    fils = os.listdir(fil_loc)
    if group_labels is None:
        _ = [fil.split('_') for fil in fils]
        for fil in _:
            try: assert len(fil) !=3
            except: raise ValueError('Group labels could not be inferred automatically.')
        group_labels = np.unique([fil[0] for fil in _])

    # Collate chunks per group
    for group in tqdm(group_labels, unit='groups', desc='Collating demultiplexed chunks'):

        # Collate fragments
        group_fils = [fil for fil in fils if group in fil and 'pseudorep' not in fil]

        fil_nam = group+'_fragments.txt'
        _cat_files(group_fils, fil_loc, fil_nam)

        # Collate pseudoreps
        pseudorep1_fils = [fil for fil in fils if group in fil and 'pseudorep1' in fil]

        fil_nam = group+'_pseudorep1.txt'
        _cat_files(pseudorep1_fils, fil_loc, fil_nam)

        pseudorep2_fils = [fil for fil in fils if group in fil and 'pseudorep2' in fil]

        fil_nam = group+'_pseudorep2.txt'
        _cat_files(pseudorep2_fils, fil_loc, fil_nam)

        # Remove chunks
        if remove_chunks:
            remove_fils = group_fils + pseudorep1_fils + pseudorep2_fils
            for chunk in remove_fils:
                os.system('rm {}'.format(os.path.join(fil_loc, chunk)))

if __name__=='__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('fragment_files', type=str)
    parser.add_argument('barcodes', type=str)
    parser.add_argument('-o', '--output_dir',  default='./', type=str)
    parser.add_argument('--chunksize', default=5000000, type=int)
    parser.add_argument('--num_cores', default=-1, type=int)
    parser.add_argument('--cleanup', action='store_true')

    args = parser.parse_args()

    # Read in barcodes (dir or file)
    if os.path.isfile(args.barcodes):
        barcode_fils = [args.barcodes]
    elif os.path.isdir(args.barcodes):
        barcode_fils = os.listdir(args.barcodes)
        barcode_fils = [os.path.join(args.barcodes, fil) for fil in barcode_fils]

    barcodes_df = read_barcode_files(*barcode_fils)

    # Read in fragments (dir or files)
    if os.path.isfile(args.fragment_files):
        fragment_fils = [args.fragment_files]
    elif os.path.isdir(args.fragment_files):
        fragment_fils = os.listdir(args.fragment_files)
        fragment_fils = [os.path.join(args.fragment_files, fil) for fil in fragment_fils]

    # Run demultiplexing
    for fil in fragment_fils:
        demultiplex_fragment_file(fragment_fil=fil, 
                                  barcodes_df=barcodes_df, 
                                  output_loc=args.output_dir,
                                  chunk_size=args.chunksize,
                                  n_jobs=args.num_cores)

    # Collate data
    collate_fragment_file(fil_loc=args.output_dir, 
                          group_labels=barcodes_df.groups.unique().tolist(),
                          remove_chunks=args.cleanup)