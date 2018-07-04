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

from quality_check_trim import trimmer

# argparse to collect command line arguments

parser = argparse.ArgumentParser("take raw NGS reads, trim, calculate 31-mer peaks, and produce a graph\nnote: requires trimmomatic and bbmap to be present on the path\nnote: requires pandas and matplotlib to be installed in python\nnote: the bbmap step can use a lot of memory (around 70 Gb for a HiSeq lane)\nnote: the stem of input files should be separated by an underscore\n")

parser.add_argument('-1', '--forward_reads', type = str,
        nargs=1, help = "fastq forward reads")
parser.add_argument('-2', '--reverse_reads', type = str,
        nargs=1, help = "fastq reverse reads")
parser.add_argument('--trim', action = 'store_true',
        help = "perform trimming with trimmomatic?")
parser.add_argument('-t', '--threads', type = str,
        nargs="?", default="16", help = "threads for trimming and khist [default 16]")

if len(sys.argv) == 1:  # if no args are given
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

print(args)

# check bbmap has been loaded

#try:
#    tmp = subprocess.call(["khist.sh", "-version"])
#except OSError as e:
#    if e.errno == os.errno.ENOENT:
#        print("\nkhist could not be found: try 'module load bbmap'\n")
#        raise

# grab stem name of files for later

stem = args.forward_reads[0].split("_")[0]

if args.trim:
    F_read, R_read = trimmer(args, stem)
else:
    F_read = args.forward_reads
    R_read = args.reverse_reads

# now use bbmap's khist to draw kmer profile

print("~~~ beginning kmer profile with bbmap ~~~")

subprocess.call(["khist.sh", "in=" + F_read, "in2=" + R_read, "khist=" + stem + ".khist.txt", "threads=16",
 "k=31"])

# create a quick graph from the khist file
# zoom is set to usually work (but might not always)

print("~~~ drawing kmer figure  ~~~")

khist = pd.read_csv(stem + ".khist.txt", sep="\t")

khist.plot(x="#Depth", y="Unique_Kmers", xlim=[0, 10000], ylim=[0, 3000])

plt.savefig(stem + ".khist.png", dpi=300)


