# Wrapper around the nano-CUT&Tag data analysis pipeline
# 27/06/2022
# The pipeline takes fastq files with multiplexed nano-CUT&Tag data as input
# It performs:
#   1. Demultiplexing of the fastq files
#   2. Run cellranger-atac
#   3. Perform custom cell picking (cellranger fails in this step)
#   4. Construct seurat object with cell x 5kb_bin matrix
#   5. Construct seurat object with cell x peak matrix

import argparse
import glob
import os, sys
import re
import yaml

parser = argparse.ArgumentParser(description='Preprocessing pipeline for single-cell nano-CUT&Tag data analysis')
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--fastq',  type=str, required=True, help='Path to folder with fastq files')
parser.add_argument('--barcodes',  type=str, required=True,nargs="+", help='list of space separated barcodes to be used for demultiplexing [e.g. ATAGAGGC TATAGCCT CCTATCCT]')
parser.add_argument('--modalities',  type=str, required=True,nargs="+", help='list of space separated modalities corresponding to the barcodes [e.g. ATAC H3K27ac H3K27me3]')
parser.add_argument('--cellranger_ref', type=str,default = '/data/ref/cellranger-atac/refdata-cellranger-atac-mm10-2020-A-2.0.0/', help = 'path to cellranger reference folder')
parser.add_argument('--threads',  type=int, default = 1, help='Number of threads')
parser.add_argument('--genome',  type=str, default = 'mm10', help='Genome to use for Gene activity scores [only mm10 supported for now]')
parser.add_argument('--tempdir',  type=str, default = '~/temp', help='Path to temp directory for sort command')
parser.add_argument('--snakeargs', default = " ",  type=str, nargs="+", help='Optional arguments to be passed to snakemake, must be formated --snakeargs="--PLACE_ARGS_HERE --MORE_ARGS" [e.g. --snakeargs="--dryrun --printshellcmds"')
args = parser.parse_args()

def get_sample_id_from_fastq(files):
    sample_names = [os.path.basename(x) for x in files]
    sample_names = [re.split("_S[0-9]_",x)[0] for x in sample_names]
    sample_names = list(set(sample_names))
    return(sample_names)

def match_barcode_to_modality(modalities, barcodes):
    if len(modalities) != len(barcodes):
        sys.stderr.write("*** ERROR: Number of modalities does not match the numbe of barcodes\n")
        sys.stderr.write("Modalities: {}\nBarcodes:{}\n".format(modalities,barcodes))
        sys.exit()
    return {modalities[i]: barcodes[i] for i,x in enumerate(modalities)}

# Parse some arguments
files     = glob.glob(args.fastq + "/*.fastq") + glob.glob(args.fastq + "/*.fq") + glob.glob(args.fastq + "/*.fastq.gz") + glob.glob(args.fastq + "/*.fq.gz")
files     = [os.path.abspath(x) for x in files]
files     = [x for x in files if '_I1_' not in x] # Remove I1 read from input
sample_id = get_sample_id_from_fastq(files)


if(len(files) == 0):
    sys.stderr.write("*** ERROR: No Fastq files in folder: " + args.fastq + "\n")
    sys.exit()

if (sum([x.endswith('.gz') for x in files]) != len(files)):
    sys.stderr.write('*** ERROR: all files must be either gziped or not; conflicting files: \n')
    sys.stderr.write("\n".join(files) + "\n")
    sys.exit()

if(len(sample_id) != 1):
    sys.stderr.write("*** ERROR: multiple or 0 sequencing runs found in folder: " + args.fastq + "\n")
    sys.stderr.write("Make sure the folder contains data from one sequencing run\n")
    sys.stderr.write("Run names: " + ",".join(samples) + "\n")
    sys.exit()



# Create config
SNAKEFILE = os.path.dirname(os.path.realpath(__file__)) + "/workflow/Snakefile"
modality_barcode_dict = match_barcode_to_modality(args.modalities, args.barcodes)

CONFIG = {}
CONFIG['cellranger_ref'] = args.cellranger_ref
CONFIG['input_folder']   = args.fastq
CONFIG['fastq']          = files
CONFIG['sample']         = sample_id
CONFIG['modalities']     = modality_barcode_dict
CONFIG['genome']         = args.genome
CONFIG['tempdir']        = args.tempdir
with open('config.yaml', 'w') as outfile:
    yaml.dump(CONFIG, outfile, default_flow_style=False)

os.system("snakemake --snakefile {snakefile} --cores {threads} --configfile {configfile} -p {pass_args}".format(
    snakefile=SNAKEFILE,
    threads = args.threads,
    configfile = 'config.yaml',
    pass_args = " ".join(args.snakeargs)))