# Barcoded Functional Genome Analysis

This is an instruction for analyzing functional genome of *Xanthomonas citri* using RB-TnSeq[^1].

The workflow consists of:

- Tn-Seq analysis
    - Pre-process sequencing data.[`Fastp`]
    - Get barcodes and FASTA files with transposon-truncated reads.[``]
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

This step identifies barcodes and genomic DNA flanking the barcoded transposon accroding to the following model. Only reads with 17bp barcode and genomic DNA with length 15bp or more will be considered.

`GTTCGAATTCNNNNNNNNNNNNNNNNNGAGCTCACTTGTGTATAAGAGTCAG`

This model is designed as a result of the 3' end of our transposon. 

- `GTTCGAATTC` is a 10bp sequence located upstream of the barcode on the transposon. 
- `NNNNNNNNNNNNNNNNN` is the 17bp barcode.
- `GAGCTCACTTGTGTATAAGAGTCAG` is the junction between barcode and genomic DNA.

### Input file

- FASTQ file processed by fastp

### Run `get_bc_id_gdna_rm_tn.py` on example

`python get_bc_id_gdna_rm_tn.py out_read.fq.gz results/re_tn_R2.fq results/id_bc_gdna.txt`

### Output file

- 

# Example file

An example file folder can be found in [<u>here</u>](./Example_file/).

# Reference

[^1]: Wetmore, K. M. et al. Rapid quantification of mutant fitness in diverse bacteria by sequencing randomly bar-coded transposons. mBio 6, e00306-00315 (2015).

