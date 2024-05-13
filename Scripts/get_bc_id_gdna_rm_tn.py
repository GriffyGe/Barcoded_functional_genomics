#!/usr/bin/python

import re
import sys
import gzip

def extract_barcode_gdna_position(seq):
    get_barcode_gdna_position = re.search(r"GTTCGAATTC(\w+)GAGCTCACTTGTGTATAAGAGTCAG", seq)
    if get_barcode_gdna_position != None:
        barcode = get_barcode_gdna_position.group(1)
        position = get_barcode_gdna_position.end()
        if len(barcode) == 17 and len(seq[position:]) >= 15:
                return (barcode,position, seq[position:])
        else:
            return None
    else:
        return None

#建立result文件，后面得到的结果直接写入该文件
result_rm_tn = open(sys.argv[2], "a")
result_id_bc_gdna = open(sys.argv[3], "a")

#以1个read的4行为单位进行处理
with gzip.open(sys.argv[1], 'rt') as f:
    while True:
        read_id = f.readline().strip() #读取1条read的第1行
        if read_id:
            sequence = f.readline().strip() #读取1条read的第2行
            get_bc_pos = extract_barcode_gdna_position(sequence)
            if get_bc_pos is not None:
                #生成id and barcode_gdna文件
                result_id_bc_gdna.write("%s\t%s\t%s\n" % (read_id, get_bc_pos[0], get_bc_pos[2])) #qseqid, barcode, gDNA

                #生成截去transposon的fastq文件
                result_rm_tn.write("%s\n" % read_id)
                result_rm_tn.write("%s\n" % sequence[get_bc_pos[1]:])
                result_rm_tn.write(f.readline())
                seq_quality = f.readline().strip()
                result_rm_tn.write("%s\n" % seq_quality[get_bc_pos[1]:])  #截短测序质量
            else:
                f.readline()   #读取1条read的第3行
                f.readline()   #读取1条read的第4行
        else:
            break

result_rm_tn.close()
result_id_bc_gdna.close()