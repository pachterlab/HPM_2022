#!/usr/bin/env python3

import random
import os
import sys

def subsample_reads(ctl_path, aff_path, outputs):

    outputs = outputs.split()
    print(outputs)

    # Read pns for control and affected cohorts
    pns_ctl = [pn.rstrip('\n') for pn in open(ctl_path)]
    pns_aff = [pn.rstrip('\n') for pn in open(aff_path)]

    # Read full set of reads
    reads = {}
    p_1 = 'reads/expected_overlaps/200'
    for pn in pns_ctl:
        with open(os.path.join(p_1, f'{pn}_wildtype.R1.fq')) as fh:
            lines = fh.readlines()
            reads[pn] = [lines[i].rstrip('\n') for i in range(1, len(lines), 4)]

    for pn in pns_aff:
        with open(os.path.join(p_1, f'{pn}_affected.R1.fq')) as fh:
            lines = fh.readlines()
            reads[pn] = [lines[i].rstrip('\n') for i in range(1, len(lines), 4)]

    for output in sorted(list(outputs), key=lambda x: -int(x.split(':')[-1])):
        with open(output, 'w') as outfile:
            N = int(output.split(':')[-1])
            # For every N we subsample the previous list of reads
            for pn, r in reads.items():
                reads[pn] = random.sample(r, N)
                outfile.write('\n'.join(f'>{pn}.{i}\n{read}' for i, read in enumerate(r)))
                outfile.write('\n')

if __name__ == '__main__':
    subsample_reads(sys.argv[1], sys.argv[2], sys.argv[3])
