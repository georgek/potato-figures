#!/usr/bin/env python

import sys
import argparse

from collections import namedtuple

from difflib import SequenceMatcher

class Blast:
    """A blast alignment record"""
    def __init__(self,
                 tid, tbeg, tend, tlen,
                 sid, sbeg, send, slen,
                 pident, length, evalue):
        self.tid = tid
        self.tbeg = int(tbeg)
        self.tend = int(tend)
        self.tlen = int(tlen)
        self.sid = sid
        self.sbeg = int(sbeg)
        self.send = int(send)
        self.slen = int(slen)
        self.pident = float(pident)
        self.length = int(length)
        self.evalue = float(evalue)
        
    def __repr__(self):
        return "Aln of {:s} and {:s}, length {:d}".format(
            self.tid, self.sid, self.length)


def sign(x):
    if x > 0:
        return 1
    elif x < 0:
        return -1
    elif x == 0:
        return 0
    else:
        return x


def uniq(list):
    last = None
    uniqlist = []
    for item in list:
        if last is None or item != last:
            uniqlist.append(item)
            last = item
    return uniqlist


def distance(a,b):
    return 1 - SequenceMatcher(None, a, b).ratio()


def nonorderliness(seq, key):
    """Measures how much the sequence is out of order according to a key"""
    return min(distance(seq, sorted(seq,None,key,False)),
               distance(seq, sorted(seq,None,key,True)))


def nonlinearness(seq):
    """Measure how nonlinear the sequence is,
linear means the items appear in any order but are not repeated"""
    return distance(sorted(uniq(seq)), uniq(sorted(seq)))


def inclusivetoexclusive(chunks):
    """Converts [beg,end] chunks to [beg,end)"""
    return [(beg,end+1) if beg<end else (beg,end-1) for beg,end in chunks]


def nonorientedness(chunks):
    """Calculates nonorientedness of list of chunks where chunks are [beg,end)"""
    total = sum([abs(end-beg) for beg,end in chunks])
    inorder = float(total-abs(sum([end-beg for beg,end in chunks])))/2
    if total > 0:
        return inorder/total
    else:
        return 0.0


Score = namedtuple("Score", ["nonlinearliness", "nonorientedness", "nonorderliness"])


def score_transcript(alns):
    """Calculates badness of transcript alignment."""
    nl = nonlinearness([b.sid for b in alns])
    # sys.stderr.write("nonlinearness {:s}: {:f}\n".format(
    #     alns[0].tid, nl))
    alns_per_contig = split_list_on_key(alns, lambda a: a.sid)
    no = []
    nd = []
    for contig_alns in alns_per_contig:
        no.append(nonorientedness(zip(
            [b.sbeg for b in contig_alns],
            [b.send for b in contig_alns])))
        # sys.stderr.write("nonorientedness {:s}: {:f}\n".format(
        #     contig_alns[0].sid, no))
        nd.append(nonorderliness(contig_alns, lambda a: a.sbeg))
        # sys.stderr.write("nonorderliness {:s}: {:f}\n".format(
        #     contig_alns[0].sid, nd))
    return Score(nl, sum(no)/len(no), sum(nd)/len(nd))


def split_list_on_key(l, key):
    splits = [[]]
    last_key = None
    for item in sorted(l, None, key):
        if last_key is None or key(item) == last_key:
            splits[-1].append(item)
        else:
            splits.append([item])
        last_key = key(item)
    return splits


def read_alns(file):
    alns = []
    for line in file:
        alns.append(Blast(*line.split()))
    return alns


def alns_from_file(filename):
    with open(filename) as f:
        alns = read_alns(f)
    return alns


def main():
    last_tid = None
    alns = []
    for line in sys.stdin:
        aln = Blast(*line.split())
        if last_tid is None or aln.tid == last_tid:
            alns.append(aln)
        else:
            sys.stdout.write("{:s}\t{:f}\t{:f}\t{:f}\n".format(
                last_tid, *score_transcript(alns)))
            alns = [aln]
        last_tid = aln.tid
        
if __name__ == '__main__':
    main()
