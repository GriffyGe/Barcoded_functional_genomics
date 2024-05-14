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

Run `python get_bc_gdna.py -h` on Linux terminal

Helper massage is shown:

```
usage: get_bc_gdna.py [-h] -i  -o_bc  -o_gdna

Identify barcode and genomic DNA.

optional arguments:
  -h, --help  show this help message and exit
  -i          Filename, Tn-Seq sequencing data processed by fastp.
  -o_bc       Filename, an output table of readID and barcode.
  -o_gdna     Filename, an output FASTQ file of genomic DNA flanking barcoded transposons.
```
Run `get_bc_gdna.py` on examlple (**Linux** terminal)

`python .\Scripts\get_bc_gdna.py -i .\Example_file\example.fq.gz -o_bc .\Example_file\table_id_bc.txt -o_gdna .\Example_file\gdna.fq`

In order to perform blastn with FASTA file, **convert FASTQ file to FASTA file with SeqKit**.

`seqkit fq2fa gdna.fq -o gdna.fa`

## Map to the reference genome

This step maps the genomic DNA FASTA file to the [Xac CQ13 reference genome](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_023205075.1/) downloaded from NCBI.

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
-outfmt <string>: Set a specified output format.((outfmt 6)[#tabular_outformat_6])
-mt_mode 1 -num_threads <int>: Enable multi-threading by setting the -mt_mode to 1 (default value is 0). Then, specify the number of threads used with the "-num_threads" option(-num_threads 20).
```

[### Tabular outformat 6]()




# Example file

An example file folder can be found in [<u>here</u>](./Example_file/).

# Reference

[^1]: Wetmore, K. M. et al. Rapid quantification of mutant fitness in diverse bacteria by sequencing randomly bar-coded transposons. mBio 6, e00306-00315 (2015).

