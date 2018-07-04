# module to trim Illumina reads as part of kmer profiling
# Matthew J. Neave 04.07.2018

import subprocess
import os

def trimmer(args, stem):

    """
    uses trimmomatic to trim fastq Illunina reads (paired end)
    adapters will be detected in the default CSIRO location
    :param args: argparse arguments used to grab the F and R reads
    :param stem: stem of the filename used for output file naming
    :return: returns the names of the F and R trimmed files
    :raise: errors will be raised if the adapters can't be found or \
    if trimmomatic is not available on the path
    """
    # first check that the adapters are accessible
    try:
        adapter_fl = open("/apps/trimmomatic/0.36/adapters/TruSeq3-PE-2.fa")
        adapter_path = "/apps/trimmomatic/0.36/adapters/TruSeq3-PE-2.fa"
        print("~~~ found adapter files ~~~")
    except:
        print("\n~~~ could not find adapter file for trimmomatic! ~~~\n")
        raise

    # check that modules have been loaded
    try:
        subprocess.call(["trimmomatic", "-version"])
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            print("\ntrimmomatic could not be found: try 'module load trimmomatic'\n")
            raise

    ## TRIM READS ##
    print("~~~ beginning trimming with trimmomatic ~~~")

    subprocess.check_output(["trimmomatic", "PE", "-threads", args.threads, args.forward_reads[0], args.reverse_reads[0], "-baseout", stem+".fastq.gz", "ILLUMINACLIP:" + adapter_path + ":2:30:10", "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:20", "MINLEN:50"])

    forward_trimmed = args.forward_reads[0].rstrip("R1.fastq.gz") + "1P.fastq.gz"
    reverse_trimmed = args.reverse_reads[0].rstrip("R2.fastq.gz") + "2P.fastq.gz"

    return(forward_trimmed, reverse_trimmed)
