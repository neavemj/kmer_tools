#!/usr/bin/env python

# take raw NGS reads, trim, calc kmer peak, and draw graph
# Matthew J. Neave 21.06.2018

import argparse
import subprocess # use to pass command line calls
import os, sys
import pandas as pd
import matplotlib
matplotlib.use("Agg") # avoids 'no display' error
import matplotlib.pyplot as plt

# argparse to collect command line arguments

parser = argparse.ArgumentParser("take raw NGS reads, trim, calculate 31-mer peaks, and produce a graph\nnote: requires trimmomatic and bbmap to be present on the path\nnote: requires pandas and matplotlib to be installed in python\nnote: the bbmap step can use a lot of memory (around 70 Gb for a HiSeq lane)\nnote: the stem of input files should be separated by an underscore\n")

parser.add_argument('-1', '--forward_reads', type = str,
        nargs=1, help = "raw fastq forward reads")
parser.add_argument('-2', '--reverse_reads', type = str,
        nargs=1, help = "raw fastq reverse reads")
parser.add_argument('-t', '--threads', type = str,
        nargs="?", default="16", help = "threads for trimming and khist [default 16]")

if len(sys.argv) == 1:  # if no args are given
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check that trimmomatic adapters are present (CSIRO modules location)

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

try:
    tmp = subprocess.call(["khist.sh", "-version"])
except OSError as e:
    if e.errno == os.errno.ENOENT:
        print("\nkhist could not be found: try 'module load bbmap'\n")
        raise

# grab stem name of files for later

stem = args.forward_reads[0].split("_")[0]

## TRIM READS ##
print("~~~ beginning trimming with trimmomatic ~~~")

#subprocess.check_output(["trimmomatic", "PE", "-threads", args.threads, args.forward_reads[0], args.reverse_reads[0], "-baseout", stem+".fastq.gz", "ILLUMINACLIP:" + adapter_path + ":2:30:10", "LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:20", "MINLEN:50"])

forward_trimmed = args.forward_reads[0].rstrip("R1.fastq.gz") + "1P.fastq.gz"
reverse_trimmed = args.reverse_reads[0].rstrip("R2.fastq.gz") + "2P.fastq.gz"

# now use bbmap's khist to draw kmer profile

print("~~~ beginning kmer profile with bbmap ~~~")

#subprocess.call(["khist.sh", "in=" + forward_trimmed, "in2=" + reverse_trimmed, "khist=" + stem + ".khist.txt", "threads=16", "k=31"])

# create a quick graph from the khist file
# zoom is set to usually work (but might not always)

print("~~~ drawing kmer figure  ~~~")

khist = pd.read_csv(stem + ".khist.txt", sep="\t")

khist.plot(x="#Depth", y="Unique_Kmers", xlim=[0, 10000], ylim=[0, 3000])

plt.savefig(stem + ".khist.png", dpi=300)


