#!/usr/bin/python

import re
import argparse
import gzip

parser = argparse.ArgumentParser(description='Identify barcode and genomic DNA.')
parser.add_argument('-i', '--input', metavar='', type=str, required=True, help='Filename, Tn-Seq sequencing data processed by fastp.')
parser.add_argument('-o_bc', '--output_barcode', metavar='', type=str, required=True, help='Filename, an output table of readID and barcode.')
parser.add_argument('-o_gdna', '--output_gdna', metavar='', type=str, required=True, help='Filename, an output FASTQ file of genomic DNA flanking barcoded transposons.')
args = parser.parse_args()

def extract_barcode_gdna_ins(seq):
    get_barcode_gdna_position = re.search(r"GTTCGAATTC(\w+)GAGCTCACTTGTGTATAAGAGTCAG", seq)
    if get_barcode_gdna_position != None:
        barcode = get_barcode_gdna_position.group(1)
        insertion_site = get_barcode_gdna_position.end()
        if len(barcode) == 17 and len(seq[insertion_site:]) >= 15:
                #genomic_DNA = seq[insertion_site:]
                return (barcode, seq[insertion_site:], insertion_site)
        else:
            return None
    else:
        return None

if __name__ == '__main__':
    result_gdna = open(args.output_gdna, "a")
    result_id_bc = open(args.output_barcode, "a")

    with gzip.open(args.input, 'rt') as f:
        while True:
            read_id = f.readline().strip()
            if read_id:
                sequence = f.readline().strip()
                get_bc_gdna_ins = extract_barcode_gdna_ins(sequence)
                if get_bc_gdna_ins is not None:
                    #Generate a table of readID and barcode.
                    result_id_bc.write("%s\t%s\n" % (read_id, get_bc_gdna_ins[0]))
                    
                    #Gnerate a FASTQ file of genomic DNA flanking barcoded transposon.
                    #Write the first line: readID
                    result_gdna.write("%s\n" % read_id)
                    #Write the second line: genomic DNA
                    result_gdna.write("%s\n" % get_bc_gdna_ins[1])
                    #Write the third line
                    result_gdna.write(f.readline())
                    #Write the fourth line: sequencing quality
                    seq_quality = f.readline().strip()
                    result_gdna.write("%s\n" % seq_quality[get_bc_gdna_ins[2]:])
                else:
                    f.readline()
                    f.readline()
            else:
                break

    result_gdna.close()
    result_id_bc.close()