#!/usr/bin/env python3

# Precursor table output, from parquet to tsv, could possibly run directly
# inside DIA-NN process

import sys
import os

from pyarrow import compute as pc
from pyarrow import parquet as pq
from pyarrow import csv as pcsv
from pyarrow import Table, array

precursors = pq.read_table(sys.argv[1])

if sys.argv[2] == 'trypsin':
    miscl_re = '[KR][^P]'
elif sys.argv[2] == 'trypsinp':
    miscl_re = '[KR][A-Z]'

inputfn = sys.argv[3]
normalizing = sys.argv[4]

precursors = precursors.append_column('missed_cleavages', pc.count_substring_regex(precursors['Stripped.Sequence'], miscl_re))

samples, files = [], []
with open(inputfn) as fp:
    header = next(fp).strip().split('\t')
    for line in fp:
        vals = line.strip().split('\t')
        files.append(os.path.splitext(os.path.basename(vals[0]))[0])
        samples.append(vals[1])
fnsamples = Table.from_arrays([array(files), array(samples)], names=['Run', 'sample'])
sorted_p = precursors.join(fnsamples, 'Run', join_type='inner').sort_by('Run').sort_by('sample')


# Find which rows contain "next sample/file" in sorted table
new_file_row_breaks, new_sample_row_breaks = [], []
old_sample = False
for row, sample in enumerate(sorted_p['sample']):
    if sample != old_sample:
        new_sample_row_breaks.append(row)
        old_sample = sample
new_sample_row_breaks.append(row + 1)
old_fn = False
for row, fn in enumerate(sorted_p['Run']):
    if fn != old_fn:
        new_file_row_breaks.append(row)
        old_fn = fn
new_file_row_breaks.append(row + 1)

# scenarios:
# - multiple files (>maxslicesize), samples
# - multiple files (>maxslicesize),  one sample (so all need to be in same slice)
# - multiple files (<maxslicesize), all in one slice
# - single file

write_opt = pcsv.WriteOptions(delimiter='\t', quoting_style='none', quoting_header='none')
rowcount = row
# 100 samples is reasonable plot size
maxslicesize = 100
start_ix, start, end, chunk = 0, 0, 0, 0
max_ix = len(new_file_row_breaks)
while end < rowcount:
    # max e.g 100 files in a slice of the table, get last table row start + max
    end_ix = start_ix + maxslicesize
    if end_ix < max_ix:
        max_fn_row = new_file_row_breaks[end_ix]
    else:
        # last row reached
        max_fn_row = new_file_row_breaks[max_ix - 1]
    # make start the previous end
    start_ix = end_ix
    start = end-start
    if max_fn_row in new_sample_row_breaks:
        # lucky, sample break is right at file break (e.g. 1 file per sample)
        end = max_fn_row
    elif sample_sublist := [x for x in new_sample_row_breaks if start < x <= max_fn_row]:
        # there is no sample break at file break, so take a previous sample break
        end = sample_sublist[-1]
        pass
    elif end := next(x + 1 for x in new_sample_row_breaks if x > max_fn_row):
        # there are more files than maxslice for a given sample, so there is
        # still more files with that sample even you have taken a slice of 100 files.
        # In that case we select the next sample break
        pass
    pcsv.write_csv(sorted_p.slice(start, end-start), f'precursors_{normalizing}_{chunk}.txt', write_opt)
    chunk += 1

# Also write full table for user output
pcsv.write_csv(precursors, 'precursors.txt', write_opt)

