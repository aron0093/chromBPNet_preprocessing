import os
import logging

def intersect_peaks(*peak_files, chr_order, min_overlap=0.5):

    '''
    Intersect peak sets from pseudoreplicates.

    ARGS
        peak_files: 
            List of pseudorep peak_files to intersect.
            Ordered as ps1, ps2 and psT
        chr_order: path
            Order of chromosome sizes. See ChromBPNet wiki.
        min_overlap: float (0.5)
            Minimum overlap to expect between peak sets.
    ''' 

    # Check inputs
    if chr_order is None:
        raise ValueError('Chromsome sizes not provided.')


    # Intersect matching peak sets
    # TODO Replace with pybedtools instead of commandline tool
    intersect_peak_fil = peak_files[-1].replace('pseudorepT', 
                                                'pseudorep_overlap')

    return_code = os.system('bedtools intersect -u -a {} -b {} -g '.format(peak_files[-1], 
                                                                           peak_files[0])+ \
                            '{} -f {} -F {} '.format(chr_order, 
                                                     min_overlap, 
                                                     min_overlap) + \
                            '-e -sorted | bedtools intersect -u -a stdin -b {} -g {} '.format(peak_files[1], 
                                                                                              chr_order)+ \
                            '-f {} -F {} -e -sorted > {}'.format(min_overlap, 
                                                                 min_overlap, 
                                                                 intersect_peak_fil))
    if return_code!=0:
        logging.warning('Return code was not 0. Process was likely not completed.')

