import os
import chrombpnet

import argparse

# Read arguments
parser = argparse.ArgumentParser()

parser.add_argument('fragments_file', type=str)
parser.add_argument('peak_file', type=str)
parser.add_argument('-b','--bias_model', default=0.5, type=str)
parser.add_argument('-g', '--genome_fasta', type=str)
parser.add_argument('-c','--chromosomes', default=None, type=str, 
                    help='File containing chromosome sizes.')
parser.add_argument('-n','--blacklist', default=None, type=str, 
                    help='Files containing blacklist regions')
parser.add_argument('-f', '--fold_jsons', type=str)                    
parser.add_argument('-o', '--output_dir',  default='./outputs', type=str)
parser.add_argument('--label', default='', type=str)
parser.add_argument('--nobias', default=False, action='store_true')

args = parser.parse_args()

# Get GC matched negative peak set
return_code = os.system('chrombpnet prep nonpeaks -g {} '.format(args.genome_fasta)+ \
                        '-p {} -c {} -fl {} '.format(args.peak_file, 
                                                     args.chromosomes, 
                                                     args.fold_jsons)+ \
                        '-br {} -o {}'.format(args.blacklist, 
                                              os.path.join(args.output_dir, 'data')))

# Train bias model
if not args.nobias:
    bias_model_path = os.path.join(args.output_dir, 'bias_model/')
    return_code = os.system('chrombpnet bias pipeline '+ \
                            '-ifrag {} '.format(args.fragments_file) + \
                            '-d "ATAC" '+ \
                            '-g {} '.format(args.genome_fasta) + \
                            '-c {} '.format(args.chromosomes) + \
                            '-p {} '.format(args.peak_file) + \
                            '-n {} '.format(os.path.join(args.output_dir, 'data')+'_negatives.bed') + \
                            '-fl {} '.format(args.fold_jsons) + \
                            '-b 0.5 '+ \
                            '-o {} '.format(bias_model_path)+ \
                            '-fp {}'.format(args.label))
else:
    bias_model_path = args.bias_model

# Train ChromBPNet
return_code = os.system('chrombpnet pipeline '+ \
                        '-ifrag {} '.format(args.fragments_file) + \
                        '-d "ATAC" '+ \
                        '-g {} '.format(args.genome_fasta) + \
                        '-c {} '.format(args.chromosomes) + \
                        '-p {} '.format(args.peak_file) + \
                        '-n {} '.format(os.path.join(args.output_dir, 'data')+'_negatives.bed') + \
                        '-fl {} '.format(args.fold_jsons) + \
                        '-b {} '.format(bias_model_path)+ \
                        '-o {} '.format(os.path.join(args.output_dir, 'chrombpnet_model/'))+ \
                        '-fp {}'.format(args.label))
