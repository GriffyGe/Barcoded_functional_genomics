# Barcoded Functional Genome Analysis

This is an instruction for analyzing functional genome of *Xanthomonas citri* pv. *citri* CQ13 using RB-TnSeq[^1].

The workflow consists of:

- Tn-Seq analysis
    - Pre-process sequencing data.[`Fastp`]
    - Get barcodes and FASTA files with transposon-truncated reads.[`get_bc_gdna.py`&`SeqKit`]
    - Map to the reference genome.[`blastn-short`]
    - Select best hit and Match insertion site with gene info.[``]
    - Match barcodes with consistently-inserted genes. [``]

- BarSeq analysis
    - Pre-process sequencing data.[`Fastp`]
    - Get barcodes.[``]
    - Calculate the count of each barcode.[``]
    - Merge barcodes with consistently-inserted genes.[``]
    - Compute strain fitness and gene fitness.[``]

# Requirements

- Miniconda (optional, for packages installation)
- Fastp
- SeqKit
- BLAST
- Python
- R, Rstudio
- Operating system: Linux
- SSH tool: finalshell (optional, better to have if using remote server)

# Tn-Seq analysis

## Pre-process sequencing data

This step processes sequencing data through quality control, trimming of adapters using fastp.

`fastp -i Read1.fq.gz -I Read2.fq.gz -o example_read1.fq.gz -O example_read2.fq.gz -h fastp_report.html`

Then, combine pair-end sequencing data into one FASTQ file.

`cat example_read1.fq.gz example_rea21.fq.gz > example.fq.gz`

## Get barcodes and FASTA files with transposon-truncated reads(only genomic DNA remained)

This step identifies barcodes and genomic DNA flanking the barcoded transposon accroding to the following model. Only reads with 17bp barcode and genomic DNA with length 15bp or more will be considered. Finally, convert the genomic DNA FASTQ file to FASTA for mapping.

`GTTCGAATTCNNNNNNNNNNNNNNNNNGAGCTCACTTGTGTATAAGAGTCAG`

This model is designed based on the 3' end of our transposon. 

- `GTTCGAATTC` is a 10bp sequence located upstream of the barcode on the transposon. 
- `NNNNNNNNNNNNNNNNN` is the 17bp barcode.
- `GAGCTCACTTGTGTATAAGAGTCAG` is the junction between barcode and genomic DNA.

**Helper Message**

Run `python .\Scripts\get_bc_gdna.py -h` on Linux terminal

Helper massage is shown:

```
usage: get_bc_gdna.py [-h] -i  -o_bc  -o_gdna

Identify barcode and genomic DNA.

optional arguments:
  -h, --help            show this help message and exit
  -i , --input          Filename, Tn-Seq sequencing data processed by fastp.
  -o_bc , --output_barcode
                        Filename, an output table of readID and barcode.
  -o_gdna , --output_gdna
                        Filename, an output FASTQ file of genomic DNA flanking barcoded transposons.
```
Run `get_bc_gdna.py` on examlple (**Linux** terminal)

`python .\Scripts\get_bc_gdna.py -i .\Example_file\example.fq.gz -o_bc .\Example_file\table_id_bc.txt -o_gdna .\Example_file\gdna.fq`

In order to perform blastn with FASTA file, **convert FASTQ file to FASTA file with SeqKit**.

`seqkit fq2fa gdna.fq -o gdna.fa`

## Map to the reference genome

This step maps the genomic DNA FASTA file to the [<u>Xac CQ13 reference genome</u>](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_023205075.1/) downloaded from NCBI. Here, we use the genome sequences instead of genome coding sequences, considering the crucial regulatory roles that ncRNAs play in cellular processes.

1. Create database from CQ13 genome sequences (a multi-FASTA file: 1 chromosome + 2 plasmids).

`makeblastdb -in CQ13.fna -dbtype nucl -parse_seqids -out db/CQ13`

2. Perform BLASTN search against the database.

`blastn -task blastn-short -query gdna.fa -db db/CQ13 -evalue 1e-5 -outfmt 6 -out blastn_short_result.txt -num_threads 20 -mt_mode 1`

blastn application options:
```
-task blastn-short: blastn supports "blastn-short" task, optimized for sequences less than 30 nucleotides.
-query <filename>: Query FASTA filename.(gdna.fa)
-out <filename>: Output mapping result filename.(blastn_short_result.txt)
-db <database_name>: BLASTN database name.(db/CQ13)
-evalue <value>: Expect value (E) for saving hits, usually use 1e-5, default value is 10.
-outfmt <string>: Set a specified output format.(outfmt 6)
-mt_mode 1 -num_threads <int>: Enable multi-threading by setting the -mt_mode to 1 (default value is 0). Then, specify the number of threads used with the "-num_threads" option(-num_threads 20).
```

### [BLASTN Tabular Outformat 6](https://www.metagenomics.wiki/tools/blast/blastn-output-format-6) (Click for more info about outfmt 6)

**Column Headers**

`qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore`

| Col | Field | Brief description |
| ---: | :--- | :--- |
| 1 | qseqid | query sequence id |
| 2 | sseqid | subject or reference genome sequence id |
| 3 | pident | percentage of identical matches |
| 4 | length |alignment length |
| 5 | mismatch | number of mismatches |
| 6 | gapopen | number of gap openings |
| 7 | qstart | start of alignment in query |
| 8 | qend | end of alignment in query |
| 9 | sstart | start of alignment in subject |
| 10 | send | end of alignment in subject |
| 11 | evalue | expect value |
| 12 | bitscore | bit score |

## Select best hit and Match insertion site with gene info

1. Select readID(qsedid) with only one best hit. 
1. Match <u>insertion sites</u> (i.e., **sstart** of the blastn result) with the positional information of genes. 

- Input files
  - BLASTN mapping result.
  - Gene table, positional information of genes.

    Gene table is generated by extracting `scaffoldId(1st col), begin(4th col), end(5th col), attribute(9th col)` from [<u>gff file</u>](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_023205075.1/) downloaded from NCBI.

- Output files
  - Blastn result of reads with only one best hit will be generated at the working directory.
  - Blastn result of reads with multiple best hits (>=2) will be generated at the working directory.
  - Matching result of insertion sites of reads with only one best hit.

**Helper Message**

Run `python .\Scripts\select_best_hit.py -h` on Linux terminal

Helper massage is shown:

```
usage: select_best_hit.py [-h] -i  -gt  -o

Select genomic DNA with only one best hit, and match the insertion site with gene info.

optional arguments:
  -h, --help           show this help message and exit
  -i , --input         Filename, blastn result file.
  -gt , --gene_table   Filename, an excel file providing gene location (i.e., scaffold, begin, end, description).
  -o , --output        Filename, the output matching result of genomic DNA with only one best hit.
```

Run `select_best_hit.py` on examlple (**Linux** terminal)

`python .\Scripts\select_best_hit.py -i .\Example_file\blastn_short_result.txt -gt .\Example_file\gene_table.xlsx -o .\Example_file\best_hit_match_gene_res.txt`

The following information will be printed on the terminal.
```
For chromosome, the number of insertions located at:
  only one gene: 5 (1.000000).
  intergenic region: 0 (0.000000).
  overlapping region of genes: 0 (0.000000).
-----------------------------------------------------------------------
For pXAC33, the number of insertions located at:
  only one gene: 1 (0.500000).
  intergenic region: 1 (0.500000).
  overlapping region of genes: 0 (0.000000).
-----------------------------------------------------------------------
For pXAC64, the number of insertions located at:
  only one gene: 6 (0.750000).
  intergenic region: 2 (0.250000).
  overlapping region of genes: 0 (0.000000).
```


# Example file

An example file folder can be found in [<u>here</u>](./Example_file/).

# Reference

[^1]: Wetmore, K. M. et al. Rapid quantification of mutant fitness in diverse bacteria by sequencing randomly bar-coded transposons. mBio 6, e00306-00315 (2015).

