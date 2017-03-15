#!/usr/bin/env python

import sys

current_tid = None
current_tlen = 0
last_start = 0
last_end = 0
seg_lengths = []
seg_covs    = []

for line in sys.stdin:
    (tid,tstart,tend,tlen,
     seqid,seqstart,seqend,seqlen,
     pident,length,evalue) = line.split()
    tstart = int(tstart)
    tend = int(tend)
    tlen = int(tlen)
    seqstart = int(seqstart)
    seqend = int(seqend)
    seqlen = int(seqlen)
    if tid != current_tid:
        if current_tid is not None:
            sys.stdout.write("{:s}\t{:d}\t{:d}\t{:f}\t{:f}\n"
                             .format(current_tid,current_tlen,
                                     sum(seg_lengths),
                                     float(sum(seg_lengths))/
                                     current_tlen,
                                     sum(seg_covs)/len(seg_covs)))
        current_tid = tid
        current_tlen = tlen
        last_start = tstart
        last_end = tend
        seg_lengths = [last_end - last_start + 1]
        seg_covs = [1.0]
    else:
        if tstart < last_start:
            sys.exit("Input file not sorted by transcript start!")
        if tstart <= last_end:
            if tend > last_end:
                # extend segment
                seg_lengths[-1] = tend - last_start + 1
                seg_covs[-1] = (  seg_covs[-1]     * (tstart-last_start)
                                + (seg_covs[-1]+1) * (last_end-tstart+1)
                                + 1                * (tend-last_end)) \
                                      / seg_lengths[-1]
                last_end = tend
            else:
                # embed segment
                seg_covs[-1] = (  seg_covs[-1]     * (tstart-last_start)
                                + (seg_covs[-1]+1) * (tend-tstart+1)
                                + seg_covs[-1]     * (last_end-tend)) \
                                      / seg_lengths[-1]
        else:
            # new segment
            last_start = tstart
            last_end = tlen
            seg_lengths.append(last_end - last_start + 1)
            seg_covs.append(1.0)

if current_tid is not None:
    sys.stdout.write("{:s}\t{:d}\t{:d}\t{:f}\t{:f}\n"
                     .format(current_tid,current_tlen,
                             sum(seg_lengths),
                             float(sum(seg_lengths))/
                             current_tlen,
                             sum(seg_covs)/len(seg_covs)))
