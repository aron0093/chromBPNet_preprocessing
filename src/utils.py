import pandas as pd
from joblib import Parallel, delayed
from tqdm.auto import tqdm

import warnings

def _count_fragments(chunk):
    return chunk.value_counts('barcodes').reset_index()

def count_fragments(fragment_fil, 
                    barcodes=None,
                    chunksize=10000,
                    skiprows=None,
                    n_jobs=-1):
    '''
    Count number of fragments per cell barcode.

    ARGS
        fragment_fil: 
            tab separated fragment file with atleast 4 columns.
            chr, start, end, barcodes <- no headers are expected.
        barcodes: (default: None)
            Array-like of barcodes to count. If None then all.
        chunksize: int (default: 10000)
            Number of rows to read at once.
        skiprows: int (default: None)
            Number of rows to skip from top.
        n_jobs: int (default: -1)
            Number of cores to use.    
    '''  

    # Read fragments file
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=pd.errors.ParserWarning)
        chunk_gen = pd.read_csv(fragment_fil,
                                index_col=False,
                                header=None,
                                names=['chr', 'start', 'end', 'barcodes'],
                                sep='\t', 
                                compression='infer',
                                chunksize=chunksize, 
                                skiprows=skiprows)

        # Process fragment file  
        track_chunks = []
        track_chunks.append(Parallel(n_jobs=n_jobs, 
                                     backend='threading')(delayed(_count_fragments)(chunk) \
                                                          for chunk in tqdm(chunk_gen, 
                                                                            unit='chunks', 
                                                                            desc='Counting fragments'
                                                                            )))

    fragment_counts = pd.concat(track_chunks[0]).groupby('barcodes').sum()
    if barcodes is not None:
        fragment_counts = fragment_counts.loc[barcodes]

    return fragment_counts