import os
import re
import gzip
import random
import string

from pathlib import Path

import numpy as np
import pandas as pd

from joblib import Parallel, delayed
from tqdm.auto import tqdm

import logging
import warnings

import argparse

# Read barcodes per cell grouping
def read_barcode_files(*barcode_fils, infer_names=True):

    '''
    Read in barcode files per cell group w. no headers.
    Group names are inferred from file name unless
    infer_names is False. Then groups are expected in second
    column of a tab separated file.
    '''

    barcodes = []
    for barcode_fil in tqdm(barcode_fils, 
                            unit='groups', 
                            desc='Reading barcode files'):
        if infer_names:
            barcodes_ = pd.read_csv(barcode_fil, 
                                    header=None, 
                                    names=['barcodes'])
            
            # Cell group name is inferred from barcode file name
            barcodes_['groups'] = Path(barcode_fil).stem
        else:
            # Cell group labels present in barcodes file
            barcodes_ = pd.read_csv(barcode_fil, 
                                    header=None, 
                                    sep='\t',
                                    names=['barcodes', 'groups'])
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
            with open(chunk) as infile:
                f.write(infile.read())
    f.close()

def _infer_labels(fil_list):

    _ = [os.path.basename(fil).split('_') for fil in fil_list]
    for fil in _:
        try: assert len(fil) !=3
        except: raise ValueError('Group labels could not be inferred automatically.')
    group_labels = np.unique([fil[0] for fil in _])
    
    return group_labels

# Function to check fragment file formatting
def _check_fragment_file_formatting(fragment_fil):

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
            logging.info('Skipping lines {} since they look like header text'.format(idx))
            return idx, zipped

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

# Function to process - create temp chunks
def _process_fragment_file(frag_df,
                           barcodes_df,
                           frag_nam, 
                           output_loc,
                           demultiplex=True,
                           make_pseudoreps=True):
    
    '''
    Internal function for fragment file processing.
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

    # Read barcodes and filter frag files
    if barcodes_df is not None:
        # Merge barcodes files to annotate group
        frag_df = frag_df.merge(barcodes_df, 
                                how='left', 
                                left_on='barcodes', 
                                right_on='barcodes')

    # If no barcodes or no demux then set groups to dummy value
    if barcodes_df is None or not demultiplex:
        frag_df['groups'] = ''
        group_labels=['']
    else:
        group_labels=frag_df.groups.unique()
        
    # Output group-wise files for each chunk <- effect of no demux if dummy value
    track_chunks = []
    for group in group_labels:

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
        if make_pseudoreps:
            frag_df_ = frag_df_.loc[frag_df_.chr.isin(allowed_chrs)]
            pseudo_rep_1, pseudo_rep_2 = _assign_insertion_sites(frag_df_)

            _write_chunk(pseudo_rep_1, output_loc, 
                         '{}_pseudorep1_{}_{}'.format(group, frag_nam, rand_string))
            _write_chunk(pseudo_rep_1, output_loc, 
                         '{}_pseudorep2_{}_{}'.format(group, frag_nam, rand_string))
        
        track_chunks.append(rand_string)
    
    return track_chunks   
            
# Process fragment file - parallel over chunks
def process_fragment_file(fragment_fil, 
                          barcodes_df, 
                          output_loc='./',
                          chunk_size=5000000,
                          demultiplex=True,
                          make_pseudoreps=True,
                          n_jobs=-1):

    '''
    Process fragment file into cell-type groups.
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
        output_loc: os.path (default: './')
            Directory where processed chunks will be stored.
            Defaults to current directory.
        chunk_size: int (default: 5000000)
            Number of rows to read. Adjust based on available memory.
        demultiplex: bool (default: True)
            Demultiplex fragment files by cell groups
        make_pseudoreps: bool (default:True)
            Make pseudoreplicates by fragment splitting.
        n_jobs: int (default: -1)
            Number of cores to use.    
    '''  
    
    # Read fragment file 
    start_row, zipped = _check_fragment_file_formatting(fragment_fil)

    # Create output params
    frag_nam = Path(fragment_fil).stem
    if zipped:
        frag_nam = Path(frag_nam).stem
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

        # Process fragment file  
        track_chunks = []
        track_chunks.append(Parallel(n_jobs=n_jobs, 
                            backend='threading')(delayed(_process_fragment_file)(chunk, 
                                                                                 barcodes_df,
                                                                                 frag_nam, 
                                                                                 output_loc,
                                                                                 demultiplex,
                                                                                 make_pseudoreps)\
                                                for chunk in tqdm(chunk_gen, 
                                                                  unit='chunks', 
                                                                  desc='Processing {}'\
                                                                  .format(frag_nam))))
    return track_chunks

# Function to collate temp chunks
def collate_fragment_files(*frag_chunks,
                           group_labels=None,
                           collate_pseudoreps=True,
                           output_loc='./',
                           remove_chunks=False):

    '''
    Collate processed chunks into
    group wise fragment files. Also creates
    pseudo-replicates for ChromBPNet.

    ARGS
        *frag_chunks: 
            Fragment file chunks to concatenate
        group_labels: list (default: None)
            Cell-groups to collate. By default all files are collated.
        collat_pseudoreps: bool (default:True)
            If pseudorep chunks are to be collated aswell.
        output_loc: os.path (default: './')
            Directory where processed data will be stored.
            Defaults to current directory.
        remove_chunks: bool (default: True)
            Remove chunks after collation.
    ''' 

    # Collate by group
    if group_labels is None:
        # #TODO: Kinda shady - "_" delimited could cause issues.
        # logging.warning('Group labels not supplied. Trying to infer group labels')
        # group_labels = _infer_labels(frag_chunks)
        
        # Collate all if no group labels
        group_labels = ['']

    # Collate chunks per group
    track_cats = []
    for group in tqdm(group_labels, unit='groups', desc='Collating processed chunks'):

        # All groups fils
        group_fils_all = [fil for fil in frag_chunks if group in fil]

        # Collate fragments
        group_fils = [fil for fil in group_fils_all if 'pseudorep' not in fil]

        fil_nam = group+'_fragments.txt'
        _cat_files(group_fils, output_loc, fil_nam)
        track_cats.append(fil_nam)

        # Collate pseudoreps
        if collate_pseudoreps:

            if len([fil for fil in group_fils_all if 'pseudorep' in fil])!=2*len(group_fils):
                logging.warning('Number of pseudoreplicate files is less than fragment files.')

            pseudorep_fils = {'pseudorep1': [fil for fil in group_fils_all if 'pseudorep1' in fil],
                              'pseudorep2': [fil for fil in group_fils_all if 'pseudorep2' in fil]
                             }
            if len(pseudorep_fils['pseudorep1'])!=len(pseudorep_fils['pseudorep2']):
                logging.warning('Number of pseudoreplicates is not equal to each other (1 and 2)')

            for key, value in pseudorep_fils.items():
                fil_nam = group+'_{}.txt'.format(key)
                _cat_files(value, output_loc, fil_nam)
                track_cats.append(fil_nam)

            fil_nam = group+'_pseudorepT.txt'
            _cat_files([os.path.join(output_loc, group+'_pseudorep1.txt'),
                        os.path.join(output_loc, group+'_pseudorep2.txt')],
                        output_loc, fil_nam)
            track_cats.append(fil_nam)

        # Remove chunks
        if remove_chunks:
            remove_fils = group_fils 
            if collate_pseudoreps:
                remove_fils+= pseudorep_fils['pseudorep1'] + \
                              pseudorep_fils['pseudorep2']
            for chunk in remove_fils:
                return_code=os.system('rm {}'.format(chunk))
    
    return track_cats

# Sort fragment files
def sort_fragment_files(*fragment_files, 
                        output_loc='./',
                        chunk_size=5000000,
                        n_jobs=-1,
                        remove_unsorted=False):

    '''
    Sort collated group wise fragment files.

    ARGS
        fragment_files: 
            List of fragment_files to sort.
        output_loc: os.path (default: './')
            Directory where processed data will be stored.
            Defaults to current directory.
        chunk_size: int (default: 5000000)
            Number of rows to read. Adjust based on available memory.
        n_jobs: int (default: -1)
            Number of cores to use.  
        remove_unsorted: bool (default: True)
            Remove unsorted files after sorting.
    ''' 
       
    if n_jobs < 0:
        num_threads = os.cpu_count() + n_jobs + 1
    elif n_jobs > 0:
        num_threads = n_jobs

    buffer_size = '90%'

    # Sort fils
    track_sorted = []
    for fil in tqdm(fragment_files, 
                    unit='files', 
                    desc='Sorting fragment files.'):
        output_fil = os.path.join(output_loc, Path(fil).stem+'_sorted.txt')
        return_code = os.system('sort -k 1,1 -k 2,2n -u '+ \
                                '-S {} --parallel={} {} -o {}'.format(buffer_size,
                                                                      num_threads,
                                                                      fil,
                                                                      output_fil))
        track_sorted.append(output_fil)
    
        # Remove pre-sorted files
        if remove_unsorted and return_code==0:
            return_code=os.system('rm {}'.format(fil))

    return track_sorted

# Call peaks with MACS3 and prep ChromBPNet peak set
def call_peaks(*fragment_files, 
               chr_order=None,
               blacklist=None,
               q_threshold=0.01, 
               genome_size='hs',
               shift=-75,
               extsize=150,
               top_n_peaks=300000,
               min_overlap=0.5,
               n_jobs=-1,
               output_loc='./'):

    '''
    Calls peaks with macs3. Intersect peak sets from
    pseudoreplicates and remove blacklist regions. Outputs
    peak set for training ChromBPNet.

    ARGS
        fragment_files: 
            List of fragment_files to call peaks on.
        chr_order: path
            Order of chromosome sizes. See ChromBPNet wiki.
        blacklist: path
            Blacklist peak regions. See ChromBPNet wiki.
        q_threshold: float (default: 0.01)
            q-valure threshold for peak selection.
        genome_size: {'hs', 'mm', 'ce', 'dm'} (default: 'hs')
            Approximate genome size. See MACS3 docs.
        shift: int (default: -75)
        
        extsize: int (default: 150)

        top_n_peaks: int (default: 300000)
            Maximum number of peaks to select.
        min_overlap: float (0.5)
            Minimum overlap to expect between peak sets.
        output_loc: os.path (default: './')
            Directory where processed data will be stored.
            Defaults to current directory.
        n_jobs: int (default: -1)
            Number of cores to use.  
    ''' 

    # Check inputs
    if chr_order is None:
        raise ValueError('Chromsome sizes not provided.')

    # Call peaks with MACS3
    # TODO: Make programatic if possible than using command line
    peak_fils = []
    for fil in fragment_files:
        if 'pseudorep' in fil:
            return_code = os.system('macs3 callpeak -B -f BED --nomodel --SPMR --keep-dup all --call-summits '+ \
                                    '-t {}  -n {}  -q {} -g {} --shift {} --extsize {} --outdir {}'.format(fil, 
                                                                                                           Path(fil).stem, 
                                                                                                           q_threshold, 
                                                                                                           genome_size,
                                                                                                           shift,
                                                                                                           extsize,
                                                                                                           output_loc))
        
            # Get top peaks
            peak_fil_in = os.path.join(output_loc, Path(fil).stem+'_peaks.narrowPeak')
            peak_fil_out = os.path.join(output_loc, Path(fil).stem+'_peaks_sorted.narrowPeak')

            if n_jobs < 0:
                num_threads = os.cpu_count() + n_jobs + 1
            elif n_jobs > 0:
                num_threads = n_jobs

            return_code = os.system('sort -k 8gr,8gr {} --parallel={} | head -n {} | '.format(peak_fil_in,
                                                                                              num_threads,
                                                                                              top_n_peaks) + \
                                    'sort -k 1,1 -k2,2n --parallel={} > {}'.format(num_threads, peak_fil_out))
            
            if return_code!=0:
                logging.warning('Return code was not 0. Process was likely not completed.')
            
            peak_fils.append(peak_fil_out)

    # Get matching peak files
    fil_nams = np.unique([fil.split('pseudorep')[0] for fil in peak_fils])

    # Intersect matching peak sets
    # TODO Replace with pybedtools instead of commandline tool
    for nam in tqdm(fil_nams, unit='pseudoreplicates', desc='Intersecting pseudoreps'):
        peak_fils_ = [fil for fil in peak_fils if nam in fil]
        if len(peak_fils_)!=3:
            raise ValueError('Matching pseudoreps were not found.')

        peak_fils_T = [fil for fil in peak_fils_ if 'pseudorepT' in fil][0]
        peak_fils_1 = [fil for fil in peak_fils_ if 'pseudorep1' in fil][0]
        peak_fils_2 = [fil for fil in peak_fils_ if 'pseudorep2' in fil][0]

        intersect_peak_fil = peak_fils_T.replace('pseudorepT', 
                                                 'pseudorep_overlap')
        return_code = os.system('bedtools intersect -u -a {} -b {} -g '.format(peak_fils_T, 
                                                                               peak_fils_1)+ \
                                '{} -f {} -F {} '.format(chr_order, 
                                                         min_overlap, 
                                                         min_overlap) + \
                                '-e -sorted | bedtools intersect -u -a stdin -b {} -g {} '.format(peak_fils_2, 
                                                                                                  chr_order)+ \
                                '-f {} -F {} -e -sorted > {}'.format(min_overlap, 
                                                                     min_overlap, 
                                                                     intersect_peak_fil))
        if return_code!=0:
            logging.warning('Return code was not 0. Process was likely not completed.')

        # Remove blacklist peaks
        if blacklist is not None:
            filtered_peak_file=intersect_peak_fil.replace('pseudorep_overlap', 
                                                          'pseudorep_overlap_filtered')
            return_code = os.system('bedtools intersect -v -a '+ \
                                    '{} -b {} > {}'.format(intersect_peak_fil,
                                                           blacklist, 
                                                           filtered_peak_file))
            if return_code!=0:
                logging.warning('Return code was not 0. Process was likely not completed.')
        else:
            filtered_peak_file = intersect_peak_fil

        # Make ChromBPNet input peak file
        bpnet_peak_file = filtered_peak_file.replace('pseudorep_overlap', 
                                                     'pseudorep_bpnet')
        with open(filtered_peak_file, "r") as f_in, open(bpnet_peak_file, "w") as f_out:
            for line in f_in:
                line_split = line.strip().split("\t")
                chro, start, peak = line_split[0], line_split[1], line_split[9]

                midpoint = int(start) + int(peak)
                start_new = midpoint - 500; end_new = midpoint + 500

                try: assert(end_new - start_new == 1000)
                except: raise ValueError('Could not make peak size 1000.')

                f_out.write(f"{chro}\t{start_new}\t{end_new}\t.\t.\t.\t.\t.\t.\t500\n")

# Make bigwig files from fragments
def make_bigwig(*fragment_files, output_loc='./', chr_ord=None, genome_fi=None):

    raise NotImplementedError()

    # script_loc="/users/salil512/chrombpnet/chrombpnet/helpers/preprocessing/reads_to_bigwig/reads_to_bigwig.py"
    # fasta_file="/srv/scratch/salil512/fastas/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
    # chro_sizes_file="/mnt/lab_data3/salil512/adult_heart/0_raw_data/other/GRCh38_EBV.chrom.sizes.tsv"

    # echo "${dataset} starting"
    # python3.8 ${script_loc} -g ${fasta_file} -ifrag ${frag_file} -c ${chro_sizes_file} -o ${out_file} -d ATAC
    # echo "${dataset} done"

if __name__=='__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('fragment_files', type=str, 
                        help='Tab separated fragment file or directory containing multiple files.')
    parser.add_argument('-b','--barcodes', default=None, type=str, 
                        help='Barcodes file or directory containing multiple files. One file per cell group.')
    parser.add_argument('-c','--chromosomes', default=None, type=str, 
                        help='File containing chromosome sizes.')
    parser.add_argument('-n','--blacklist', default=None, type=str, 
                        help='Files containing blacklist regions')                    
    parser.add_argument('-o', '--output_dir',  default='./outputs/', type=str)
    parser.add_argument('--chunksize', default=5000000, type=int)
    parser.add_argument('--threads', default=-1, type=int)
    parser.add_argument('--cleanup', default=False, action='store_true')
    parser.add_argument('--nodemux', default=False, action='store_true')
    parser.add_argument('--nocat', default=False, action='store_true')
    parser.add_argument('--nosplit', default=False, action='store_true')    
    parser.add_argument('--nosort', default=False, action='store_true')
    parser.add_argument('--nopeaks', default=False, action='store_true')

    args = parser.parse_args()

    # Make output directory if doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Read in fragments (dir or files)
    if os.path.isfile(args.fragment_files):
        fragment_fils = [args.fragment_files]
        args.fragment_files = os.path.dirname(args.fragment_files)
    elif os.path.isdir(args.fragment_files):
        fragment_fils = os.listdir(args.fragment_files)
        fragment_fils = [os.path.join(args.fragment_files, fil) for fil in fragment_fils]

    # Read in barcodes (dir or file)
    if args.barcodes is not None:
        if os.path.isfile(args.barcodes):
            barcode_fils = [args.barcodes]
        elif os.path.isdir(args.barcodes):
            barcode_fils = os.listdir(args.barcodes)
            barcode_fils = [os.path.join(args.barcodes, fil) for fil in barcode_fils]

        barcodes_df = read_barcode_files(*barcode_fils)
        group_labels = barcodes_df.groups.unique()
    else: 
        # If no barcodes then assume each frag file is already demuxed.
        # Make dummy "group labels" from the file itself.
        barcodes_df = None
        group_labels = []
        args.nodemux = True
        for fil in fragment_fils:
            idx, zipped = _check_fragment_file_formatting(fil)
            if zipped:
                group_labels.append(Path(frag_fil).stem.replace('.gz','').replace('.gzip',''))
            else:
                group_labels.append(Path(frag_fil).stem)

    # Run processing
    track_chunks = []
    for fil in fragment_fils:
        demultiplex = not args.nodemux
        make_pseudoreps = not args.nosplit

        # If no demux and no make_pseudoreps then dont chunk file
        if demultiplex or make_pseudoreps:
            track_chunks.extend(process_fragment_file(fragment_fil=fil, 
                                                      barcodes_df=barcodes_df, 
                                                      output_loc=args.output_dir,
                                                      chunk_size=args.chunksize,
                                                      demultiplex=demultiplex,
                                                      make_pseudoreps=make_pseudoreps,
                                                      n_jobs=args.threads))
        else:
            # If no demuxing or splitting then cat the frag files
            track_chunks = group_labels

    # Compile file list for collation
    if demultiplex or make_pseudoreps:
        # If these chunks are produced at output
        chunk_dir = args.output_dir
    else:
        # If this then chunks are in the input dir
        chunk_dir = args.fragment_files

    # Collate data
    if not args.nocat:

        frag_chunks = []
        for fil in os.listdir(chunk_dir):
            for chunk in track_chunks:
                if chunk in fil:
                    frag_chunks.append(os.path.join(chunk_dir,fil))
                    
        track_cats = collate_fragment_files(*frag_chunks, 
                                            group_labels=group_labels,
                                            output_loc=args.output_dir,
                                            remove_chunks=args.cleanup)
    
    # Raise warning if running sort with nocat
    if args.nocat and not args.nosort:
        # If no catting then either sort chunks or fragment_fils
        if len(track_chunks) != 0:
            track_cats = track_chunks
        else:
            track_cats = group_labels
        logging.warning('Running sorting w/o concatenation. Resorting will be required after concat.')

    # Compile list for sorting
    if demultiplex or make_pseudoreps:
        cat_dir = args.output_dir
    else:
        cat_dir = args.fragment_files

    if not args.nocat:
        frag_files = []
        for fil in os.listdir(cat_dir):
            # If skip all previous steps then filtered by group label here
            for nam in track_cats:
                if nam in fil:
                    frag_files.append(os.path.join(cat_dir, fil))

    # Sort data
    if not args.nosort:
        track_sorted = sort_fragment_files(*frag_files,
                                           output_loc=args.output_dir,
                                           n_jobs=args.threads,
                                           remove_unsorted=args.cleanup)
        print(track_sorted)
    
    # Call peaks and make peakset
    if not args.nopeaks:
        try: track_sorted
        except: track_sorted=fragment_fils
        call_peaks(*track_sorted, 
                   chr_order=args.chromosomes, 
                   blacklist=args.blacklist,
                   n_jobs=args.threads, 
                   output_loc=args.output_dir)

    