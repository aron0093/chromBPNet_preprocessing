import os
import gzip

import pandas as pd

import warnings

# Function to check fragment file formatting
def check_formatting(fragment_fil):

    '''
    Return first line idx without header text.
    '''

    zipped=False
    try: 
        gzip.open(fragment_fil,'rt').read(1)
        fil = gzip.open(fragment_fil,'rt')
        zipped=True
    except: fil = open(fragment_fil, 'r')

    for idx, line in enumerate(fil):
        if line[0] == '#':
            continue
        else:
            return idx, zipped

# Concatenate fragment files
def concatenate_fragments(*fragment_fils, 
                          output_path,
                          remove_fils=False):

    '''
    Concatenate fragment files into one.
    '''

    with open(output_path,'w') as f:
        for fil in fragment_fils:
            with open(fil) as infile:
                f.write(infile.read())
    f.close()

    if remove_fils:
        for fil in remove_fils:
            return_code=os.system('rm {}'.format(fil))

# Sort fragment file
def sort_fragments(fragment_fil, 
                   output_path,
                   remove_unsorted=False,
                   n_jobs=-1):

    '''
    Sort fragment file.

    ARGS
        fragment_fil: 
            tab separated fragment file with atleast 4 columns.
            chr, start, end, barcodes <- no headers are expected.
        output_path: os.path
            path at which sorted file will be created.
        remove_unsorted: bool (default: True)
            Remove unsorted files after sorting.
        n_jobs: int (default: -1)
            Number of cores to use.  
    ''' 
       
    if n_jobs < 0:
        num_threads = os.cpu_count() + n_jobs + 1
    elif n_jobs > 0:
        num_threads = n_jobs

    buffer_size = '90%'

    # Sort fils
    return_code = os.system('sort -k 1,1 -k 2,2n -u '+ \
                            '-S {} --parallel={} {} -o {}'.format(buffer_size,
                                                                  num_threads,
                                                                  fragment_fil,
                                                                  output_path))
    
    # Remove pre-sorted files
    if remove_unsorted and return_code==0:
        return_code=os.system('rm {}'.format(fragment_fil))

# Read in fragment file by chunks
def create_chunk_generator(fragment_fil, 
                           chunksize=5000000, 
                           skiprows=None):

    '''
    Create a chunk generator for fragment files.

    ARGS
        fragment_fil: 
            tab separated fragment file with atleast 4 columns.
            chr, start, end, barcodes <- no headers are expected.
        chunksize: int (default: 10000)
            Number of rows to read at once.
        skiprows: int (default: None)
            Number of rows to skip from top.
    '''

    # Read fragments file
    chunk_gen = pd.read_csv(fragment_fil,
                            index_col=False,
                            header=None,
                            names=['chr', 'start', 'end', 'barcodes'],
                            sep='\t', 
                            compression='infer',
                            chunksize=chunksize, 
                            skiprows=skiprows)
    return chunk_gen

# Write out fragment data in common format
def write_chunk(fragments_df, output_path):

    fragments_df.to_csv(output_path,
                         sep='\t', 
                         header=False, 
                         index=False)