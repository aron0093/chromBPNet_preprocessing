import os
import logging

# Filter blacklist regions
def filter_blacklist(intersect_peak_fil, blacklist):

    '''    
    Remove blacklist regions. Outputs
    peak set for training ChromBPNet.

    ARGS
        peak_fil: 
            Intersected pseudorep file.
        blacklist: path
            Blacklist peak regions. See ChromBPNet wiki.
    '''

    # Remove blacklist peaks
    # TODO: Use pybedtools instead
    filtered_peak_file=intersect_peak_fil.replace('pseudorep_overlap', 
                                                  'pseudorep_overlap_filtered')
    return_code = os.system('bedtools intersect -v -a '+ \
                            '{} -b {} > {}'.format(intersect_peak_fil,
                                                    blacklist, 
                                                    filtered_peak_file))
    if return_code!=0:
        logging.warning('Return code was not 0. Process was likely not completed.')

# Make ChromBPNet input
def make_chrombpnet_input(peak_fil):

    # Make ChromBPNet input peak file
    bpnet_peak_file = peak_fil.replace('pseudorep_overlap', 
                                       'pseudorep_bpnet')
    with open(peak_fil, "r") as f_in, open(bpnet_peak_file, "w") as f_out:
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