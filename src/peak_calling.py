import os
import logging

from pathlib import Path

# Call peaks with MACS3 and prep ChromBPNet peak set
def call_peaks(fragment_fil, 
               q_threshold=0.01, 
               genome_size='hs',
               shift=-75,
               extsize=150,
               top_n_peaks=300000,
               output_loc='./',
               n_jobs=-1):

    '''
    Calls peaks with macs3.

    ARGS
        fragment_fil: 
            tab separated fragment file with atleast 4 columns.
            chr, start, end, barcodes <- no headers are expected.
        q_threshold: float (default: 0.01)
            q-valure threshold for peak selection.
        genome_size: {'hs', 'mm', 'ce', 'dm'} (default: 'hs')
            Approximate genome size. See MACS3 docs.
        shift: int (default: -75)

        extsize: int (default: 150)

        top_n_peaks: int (default: 300000)
            Maximum number of peaks to select.
        output_loc: os.path (default: './')
            Directory where processed data will be stored.
            Defaults to current directory.
        n_jobs: int (default: -1)
            Number of cores to use.  
    ''' 

    # Check if output dir exists
    os.makedirs(output_loc, exist_ok=True)

    # Call peaks with MACS3
    # TODO: Use MACS3 python functions instead of calling command line
    return_code = os.system('macs3 callpeak -B -f BED --nomodel --SPMR --keep-dup all --call-summits '+ \
                            '-t {}  -n {}  -q {} -g {} --shift {} --extsize {} --outdir {}'.format(fragment_fil, 
                                                                                                    Path(fragment_fil).stem, 
                                                                                                    q_threshold, 
                                                                                                    genome_size,
                                                                                                    shift,
                                                                                                    extsize,
                                                                                                    output_loc))
        
    # Get top peaks

    # Peak file names
    peak_fil_in = os.path.join(output_loc, Path(fragment_fil).stem+'_peaks.narrowPeak')
    peak_fil_out = os.path.join(output_loc, Path(fragment_fil).stem+'_peaks_sorted.narrowPeak')

    # Multi-threading
    if n_jobs < 0:
        num_threads = os.cpu_count() + n_jobs + 1
    elif n_jobs > 0:
        num_threads = n_jobs

    # Sort and select top n peaks    
    return_code = os.system('sort -k 8gr,8gr {} --parallel={} | head -n {} | '.format(peak_fil_in,
                                                                                      num_threads,
                                                                                      top_n_peaks) + \
                            'sort -k 1,1 -k2,2n --parallel={} > {}'.format(num_threads, peak_fil_out))
    
    if return_code!=0:
        logging.warning('Return code was not 0. Process was likely not completed.')