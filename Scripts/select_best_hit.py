#!/usr/bin/python

import pandas as pd
import re
import argparse

parser = argparse.ArgumentParser(description='Select genomic DNA with only one best hit, and match the insertion site with gene info.')
# Input files
parser.add_argument('-i', '--input', metavar='',type=str, required=True, help='Filename, blastn result file.')
parser.add_argument('-gt', '--gene_table', metavar='',type=str, required=True, help='Filename, an excel file providing gene location (i.e., scaffold, begin, end, description).')
# Output files
parser.add_argument('-o', '--output', metavar='', type=str, required=True, help='Filename, the output matching result of genomic DNA with only one best hit.')
args = parser.parse_args()

def match_gene_info(blastn_result_list):
    # Regard sstart site as the insertion site.
    insertion_site = int(blastn_result_list[8])
    scaffold = blastn_result_list[1]
    # Select searching database.
    gene_db = gene_info_dict[scaffold]
    gene_loc = gene_db[['begin','end']]
    # Search for inserted gene.
    match_gene = gene_loc[(insertion_site >= gene_loc['begin']) & (insertion_site <= gene_loc['end'])].index
    
    # The insertion is located at only one gene.
    if len(match_gene) == 1:
        line_num = match_gene[0]
        gene_begin = gene_db.loc[line_num][1]
        gene_end = gene_db.loc[line_num][2]
        gene_desc = gene_db.loc[line_num][3]
        gene_id = re.findall(r'ID=(.*?);', gene_desc)[0]
        gene_name = re.findall(r'Name=(.*?);', gene_desc)[0]
        match_result = '%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n' % (blastn_result_list[0], int(scaffold), 
                                                                     insertion_site, float(blastn_result_list[2]), 
                                                                     int(blastn_result_list[3]), gene_begin, gene_end, 
                                                                     gene_desc, gene_id, gene_name)
        return (match_result,scaffold)
    
    # The insertion is located at intergenic spacer(igs).
    elif len(match_gene) == 0:
        return ('igs', scaffold)
    
    # The insertion is located at overlapping region of two genes.
    elif len(match_gene) >= 1:
        return ('olp', scaffold)


if __name__ == '__main__':
    # 1. Select the best hit of each genomic DNA.
    # Read blastn mapping results as a dataframe.
    df = pd.read_csv(args.input, sep='\t',
                 names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])

    # Filter for the best hits for each qseqid(readID)
    find_all_best_hits = df[df.groupby('qseqid')['pident'].transform(max) == df['pident']]
    # Count the number of the best hits of each qseqid.
    count_best_hits = find_all_best_hits.groupby('qseqid')['pident'].value_counts()
    # Get qseqids with >= two best hits and return results as a list of qseqids.
    multi_best_hit_qseqids = count_best_hits[count_best_hits >= 2].index.get_level_values('qseqid').tolist()
    # Get qseqids with only one best hit.
    unique_best_hit = find_all_best_hits[~find_all_best_hits['qseqid'].isin(multi_best_hit_qseqids)]
    not_unique_best_hit = find_all_best_hits[find_all_best_hits['qseqid'].isin(multi_best_hit_qseqids)]

    # Save the results.
    unique_best_hit.to_csv('.\\unique_best_hit.txt', sep='\t', index=False, header=False)
    not_unique_best_hit.to_csv('.\\not_unique_best_hit.txt', sep='\t', index=False, header=False)
    
    # 2. Match insertion site with gene info.
    # Generate gene_info_dict serving as gene_info_database to search for inserted genes.
    gene_table = pd.read_excel(args.gene_table)
    # NZ_CP096227.1 = '1', NZ_CP096228.1 = '2', NZ_CP096229.1 = '3'
    gene_info_chr = gene_table[gene_table['scaffoldId'] == 'NZ_CP096227.1']
    gene_info_p1 = gene_table[gene_table['scaffoldId'] == 'NZ_CP096228.1']
    gene_info_p2 = gene_table[gene_table['scaffoldId'] == 'NZ_CP096229.1']
    gene_info_dict = {'1':gene_info_chr, '2':gene_info_p1, '3':gene_info_p2}
    
    # Create a new file to write matching result.
    result_file = open(args.output, "a")

    # Build a dict to calculate the number of insertions located at 'only one gene', 'intergenic spacer' and 'overlapping region'.
    loc_dict = {'1':[0, 0, 0], '2':[0, 0, 0], '3':[0, 0, 0]}

    with open('.\\unique_best_hit.txt', 'r') as f:
        while True:
            blastn_result = f.readline().strip()
            if blastn_result:
                res_list = blastn_result.split("\t")
                match_res = match_gene_info(res_list)
                if match_res[0] == 'igs':
                    loc_dict[match_res[1]][1] = loc_dict[match_res[1]][1] + 1
                elif match_res[0] == 'olp':
                    loc_dict[match_res[1]][2] = loc_dict[match_res[1]][2] + 1
                else:
                    loc_dict[match_res[1]][0] = loc_dict[match_res[1]][0] + 1
                    result_file.write(match_res[0])
            else:
                break

    result_file.close()

    print('For chromosome, the number of insertions located at:\n'+
          '  only one gene: %d (%f).\n' % (loc_dict['1'][0], loc_dict['1'][0]/sum(loc_dict['1']))+
          '  intergenic region: %d (%f).\n' % (loc_dict['1'][1], loc_dict['1'][1]/sum(loc_dict['1']))+
          '  overlapping region of genes: %d (%f).\n' % (loc_dict['1'][2], loc_dict['1'][2]/sum(loc_dict['1']))+
          '-----------------------------------------------------------------------\n'+
          'For pXAC33, the number of insertions located at:\n'+
          '  only one gene: %d (%f).\n' % (loc_dict['2'][0], loc_dict['2'][0]/sum(loc_dict['2']))+
          '  intergenic region: %d (%f).\n' % (loc_dict['2'][1], loc_dict['2'][1]/sum(loc_dict['2']))+
          '  overlapping region of genes: %d (%f).\n' % (loc_dict['2'][2], loc_dict['2'][2]/sum(loc_dict['2']))+
          '-----------------------------------------------------------------------\n'+
          'For pXAC64, the number of insertions located at:\n'+
          '  only one gene: %d (%f).\n' % (loc_dict['3'][0], loc_dict['3'][0]/sum(loc_dict['3']))+
          '  intergenic region: %d (%f).\n' % (loc_dict['3'][1], loc_dict['3'][1]/sum(loc_dict['3']))+
          '  overlapping region of genes: %d (%f).\n' % (loc_dict['3'][2], loc_dict['3'][2]/sum(loc_dict['3']))
          )