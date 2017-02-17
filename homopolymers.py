#!/usr/bin/env python

# returns positions of homopolymers for each contig

import sys
import os
import argparse
import threading
import time
from collections import namedtuple

def progress(fp, fs, fin):
    progress = ["-", "\\", "|", "/"]
    prog = 0
    while not fin.isSet():
        if sys.stderr.isatty():
            sys.stderr.write("\r{:3.0f}% {:s}\b"
                             .format(fp.tell()/float(fs)*100,
                                     progress[prog]))
            sys.stderr.flush()
            prog = (prog + 1)%4
            time.sleep(0.1)
        else:
            sys.stderr.write("{:3.0f}% "
                             .format(fp.tell()/float(fs)*100))
            time.sleep(5)
    else:
        if sys.stderr.isatty():
            sys.stderr.write("\r100% *")
        sys.stderr.write("\n")
    return


Homopolymer = namedtuple("Homopolymer", "seqname base length beg end")

# ----- command line parsing -----
parser = argparse.ArgumentParser(
    description="Finds positions of homopolyomers.")
parser.add_argument("file", type=str,
                    help="Fasta file.")
parser.add_argument("min_length", type=int,
                    help="Minimum length of homopolymer.")
parser.add_argument("-s", "--case_sensitive",
                    dest="case", action="store_true",
                    help="Case sensitive bases.")
parser.add_argument("-b", "--bases", type=str,
                    help="Comma separated list of bases to consider.")
parser.set_defaults(case=False)

args = parser.parse_args()
# ----- end command line parsing -----

fasta_size = os.path.getsize(args.file)
fasta = open(args.file)

sys.stderr.write("Reading {:s}...\n".format(args.file))
fin = threading.Event()
pthread = threading.Thread(name = "progress",
                           target = progress,
                           args = (fasta, fasta_size, fin))

homopolymers = []
name = ""
pos = 0
run = 0
last_base = ""
try:
    pthread.start()
    for line in fasta:
        if line[0] == '>':
            if run >= args.min_length:
                    homopolymers.append(
                        Homopolymer(name, last_base,
                                    run, pos-run+1, pos))
            name = line[1:-1]
            pos = 0
            run = 0
            last_base = ""
            continue
        for base in line[:-1]:
            pos += 1
            if not args.case: base = base.upper()
            if last_base == "":
                last_base = base
                run = 1
            elif base == last_base:
                run += 1
            else:
                if run >= args.min_length:
                    homopolymers.append(
                        Homopolymer(name, last_base,
                                    run, pos-run, pos-1))
                run = 1
                last_base = base
    if run >= args.min_length:
        homopolymers.append(
            Homopolymer(name, last_base,
                        run, pos-run+1, pos))
    fin.set()
    pthread.join()
except KeyboardInterrupt:
    fin.set()
    pthread.join()
    isolate_file.close()
    sys.stderr.write("\n")
    sys.exit(1)
fasta.close()

if args.bases is not None:
    if args.case:
        bases = set(args.bases.split(","))
    else:
        bases = set(args.bases.upper().split(","))
else:
    bases = None

for homo in homopolymers:
    if bases is None or homo.base in bases:
        sys.stdout.write("{:s}\t{:s}\t{:d}\t{:d}\t{:d}\n"
                         .format(*homo))

