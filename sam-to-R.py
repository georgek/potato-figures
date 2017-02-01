#!/usr/bin/env python

import sys
import argparse

def cigar_to_aln_len(cigar):
    current_number = ""
    total_length = 0
    for c in cigar:
        if c in "0123456789":
            current_number += c
        elif c in "MD=X":
            total_length += int(current_number)
            current_number = ""
        else:
            current_number = ""
    return total_length


def parse_opts(opt):
    opts = {}
    for text_opt in opt.split():
        [name,typ,value] = text_opt.split(":")
        if typ == "A":
            opts[name] = value
        elif typ == "i":
            opts[name] = int(value)
        elif typ == "f":
            opts[name] = float(value)
        elif typ == "Z":
            opts[name] = value
        elif typ == "H":
            opts[name] = bytearray(value.decode("hex"))
        elif typ == "B":
            arr = value.split(",")
            arrtyp = arr[0]
            arr = arr[1:]
            if arrtyp == "f":
                arr = [float(f) for f in arr]
            else:
                arr = [int(f) for f in arr]
            opts[name] = arr
    return opts

# ----- command line parsing -----
parser = argparse.ArgumentParser(
    prog="sam-to-svg",
    description="Turns SAM file into SVG showing alignments.")
parser.add_argument("names", type=str,
                    help="Name of alignments, comma separated.")
parser.add_argument("sorted_sam_file", type=str, nargs="+",
                    help="Sorted SAM files.")
parser.add_argument("-l", "--ref_length", type=int,
                    help="Length of reference (taken from SAM file if header is present.")
parser.add_argument("-m", "--min_length", type=int,
                    help="Minimum length of alignment to draw.")
parser.add_argument("-i", "--min_identity", type=float,
                    help="Minimum percent identity of match.")

args = parser.parse_args()
# ----- end command line parsing -----

names = args.names.split(",")
length = args.ref_length
lanes = []
lanelens = []

if len(names) < len(args.sorted_sam_file):
    sys.exit("Error: fewer names than SAM files given.")

n = 0
minlane = 0
for sam_file in args.sorted_sam_file:
    samf = open(sam_file)
    reflengths = {}
    refname = ""
    for line in samf:
        if line[0] == "@":
            if line[1:3] == "SQ":
                split = line[3:].split()
                for bit in split:
                    if bit[0:2] == "SN":
                        refname = bit[3:]
                    if bit[0:2] == "LN":
                        reflengths[refname] = int(bit[3:])
            elif line[1:3] == "PG":
                split = line[3:].split()
                
        else:
            [qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,seq,qual,opt] = line.split(None,11)
            aln_len = cigar_to_aln_len(cigar)
            opts = parse_opts(opt)
            if rname == "*":
                continue
            if args.min_identity is not None and "NM" not in opts:
                sys.exit("No NM value in SAM file, can't calculate percent identity.")
            if (args.min_length is None or aln_len >= args.min_length) \
               and (args.min_identity is None
                    or (float(aln_len - opts["NM"])/aln_len)*100 >= args.min_identity):
                # sys.stderr.write("Match id: {:f}\n".format((float(aln_len - opts["NM"])/aln_len)*100))
                inserted = False
                for lane in lanes[minlane:]:
                    if len(lane) == 0:
                        lane.append((int(pos),int(pos)+aln_len,n))
                        inserted = True
                        break
                    elif lane[-1][1] < int(pos):
                        lane.append((int(pos),int(pos)+aln_len,n))
                        inserted = True
                        break
                if not inserted:
                    lanes.append([(int(pos),int(pos)+aln_len,n)])
                    lanelens.append(reflengths[rname])
    minlane = len(lanes)
    samf.close()
    n += 1

sys.stdout.write("{:s}\t{:s}\t{:s}\t{:s}\n".format("track","file","beg","end"))

y = 0
for lane,lanelen in zip(lanes,lanelens):
    sys.stdout.write("{:d}\t{:s}\t{:d}\t{:d}\n".format(y, "space", 1, lanelen))
    for aln in lane:
        # if aln[2] in reflengths:
        #     reflength = reflengths[aln[2]]
        # elif args.ref_length is not None:
        #     reflength = args.ref_length
        # else:
        #     sys.stderr.write("Error: {:s} not found in SAM header and --ref_length not given.\n"
        #                      .format(aln[2]))
        #     continue
        
        sys.stdout.write("{:d}\t{:s}\t{:d}\t{:d}\n".format(y, names[aln[2]], aln[0], aln[1]))
        # sys.stdout.write("    <rect x=\"{:f}\" y=\"{:d}\" width=\"{:f}\" height=\"10\" />\n"
        #                  .format((float(aln[0])/reflength)*100,
        #                          y,
        #                          (float(aln[1])/reflength)*100))
    y += 1

