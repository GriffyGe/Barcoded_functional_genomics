#!/usr/bin/python

import re
import argparse
import gzip
import pandas as pd

parser = argparse.ArgumentParser(description='Identify barcode and genomic DNA.')
parser.add_argument('-i', '--input', metavar='', type=str, required=True, help='Filename, BarSeq sequencing data processed by fastp.')
parser.add_argument('-o_bc', '--output_barcode', metavar='', type=str, required=True, help='Filename, an output table of readID and barcode.')
parser.add_argument('-o_count', '--output_count', metavar='', type=str, required=True, help='Filename, an output table of barcode and counts.')
args = parser.parse_args()

def extract_barcode(sequence):
    match_checker = re.search(r'AGTGAGCTC(.*)GAATTCGAA', sequence)
    if match_checker is not None:
        return match_checker.group(1)
    else:
        return None

if __name__ == '__main__':
    # 1. Get barcodes.
    result_file = open(args.output_barcode, "a")

    with gzip.open(args.input, 'rt') as f:
        while True:                           
            read_id = f.readline().strip()
            if read_id:
                seq = f.readline().strip()
                f.readline()
                f.readline()
                barcode = extract_barcode(seq)
                if barcode is not None and len(barcode) == 17:
                    result_file.write("%s\t%s\n" % (read_id, barcode))
            else:
                break
    result_file.close()
    
    # 2. Calculate the count of each barcode.
    df = pd.read_csv(args.output_barcode, sep='\t',
                     names=['readId', 'barcode'])

    count_bc = df['barcode'].value_counts().reset_index()
    count_bc.columns = ['barcode', 'count']

    count_bc.to_csv(args.output_count, sep='\t', index=False)