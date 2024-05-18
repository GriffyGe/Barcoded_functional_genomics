# 1. Merge barcodes with genomic DNA with gene info. ====
library(tidyverse)
readId_barcodes <- read_delim("Example_file/table_id_bc.txt",
                              delim = "\t",
                              col_names = c("qseqid", "barcode"))
gdna_gene_info <- read_delim("Example_file/best_hit_match_gene_res.txt", 
                              col_names = c("qseqid", "sseqid", "insertion_site",
                                            "pident", "length", "begin", "end", "desc",
                                            "gene_ID", "gene_name"))

# Make the format of column "qseqid" of the two dataframes same.
# @DP8480006254TRL1C021R05601207665/1 -> DP8480006254TRL1C021R05601207665/1
# @A00265:1138:HCFHGDSX5:4:1122:32298:10394 1:N:0:AAAGATAC+TCTGAAAC -> A00265:1138:HCFHGDSX5:4:1122:32298:10394
readId_barcodes$qseqid <- 
  if_else(str_detect(readId_barcodes$qseqid, pattern = " "),
          str_match(readId_barcodes$qseqid,pattern = "@(.*) ")[,2],
          str_match(readId_barcodes$qseqid, pattern = "@(.*)")[,2])


bc_gene_info <- merge(readId_barcodes, gdna_gene_info, by="qseqid")

# 2. Match barcode with unique inserted gene. ==== 
## 2.1 Compute the number of genes that each barcode has been matched with.====
list_gdna_bc_info <- split(bc_gene_info, bc_gene_info$barcode)
result <- data.frame()
not_unique <- c()
for (i in 1:length(list_gdna_bc_info)) {
  bc <- names(list_gdna_bc_info[i])
  df <- list_gdna_bc_info[[i]]
  # Get the number of genes that a barcode mapped to.
  gene_num <- length(split(df, df$gene_ID))
  if (gene_num == 1){
    sub <- df[1,-1]
    result <- rbind(result, sub)
  }
  else{
    not_unique <- append(not_unique, i)
  }
}
## 2.2 Generate a table of barcodes with unique gene.====
write_delim(result,"Example_file/barcode_with_unique_gene.txt")

## 2.3 Generate a table of barcodes without unique gene====
list_not_unique <- list_gdna_bc_info[not_unique]
library(plyr)
library(rio)
result_not_unique <- ldply(list_not_unique, data.frame)
export(result_not_unique,"Example_file/barcode_with_not_unique_gene.xlsx")

## 2.4 For barcodes without unique genes, use primary mapping sites as the "unique" sites.====
## Step1: Use PivotTable of Excel to compute the number of reads of each gene for each barcode.
## Step2: Calculate the mapping frequency of each gene for each barcode.
data <- import("Example_file/barcode_with_not_unique_gene.xlsx")
gene_dict <- list()

for (i in 1:nrow(data)){
  if (!grepl("XCQ", data[i,1])){
    barcode <-  data[i,1]
    total <-  data[i,2]
    gene_dict$barcode <- data.frame()   
    names(gene_dict)[length(gene_dict)] <- barcode
  }
  else{
    gene_id <- data[i,1]
    match_time <- data[i,2]
    match_frq <- data[i,2]/total
    sub <- data.frame("gene_id" = gene_id,
                      "match_time" = match_time,
                      "match_frq" = match_frq)
    gene_dict[[length(gene_dict)]] = rbind(gene_dict[[length(gene_dict)]], sub)
  }
}


not_unique_data <- import("Example_file/barcode_with_not_unique_gene.xlsx",sheet=2)

usable_result <- data.frame()

for (i in 1:length(gene_dict)){
  usable_barcode <- names(gene_dict[i])
  df <- arrange(gene_dict[[i]], desc(match_frq))
  #unique standard:1.match time>=10; 2.match_frq>=0.75; 3.2nd most match frq <= 1/8 * 1st most match frq
  if ((df[1,2]>=10) & (df[1,3]>=0.75) & (df[2,3]<=0.125*df[1,3])){
    usable_id = df[1,1]
    sub <- not_unique_data %>% 
      filter(barcode == usable_barcode & gene_ID == usable_id) %>% 
      distinct(barcode, .keep_all = TRUE) %>% 
      select(-c(.id, qseqid))
    usable_result <- rbind(usable_result, sub)
  }
}

unique_data <- read_delim("Example_file/barcode_with_unique_gene.txt")
final_unique_data <- rbind(unique_data, usable_result)
write_delim(final_unique_data,"Example_file/final_unique_barcode_and_gene.txt")

# 3. Select central insertions (10%~90%)====
result <- import("Example_file/final_unique_barcode_and_gene.txt") 
for (i in 1:nrow(result)){
  result$f[i] <- ((result$insertion_site[i] - result$begin[i] + 1)/ (result$end[i] - result$begin[i] + 1))
}

result2 <-  result[result$f >= 0.1 & result$f <= 0.9,]

export(result2,"Example_file/final_unique_in_central_barcode.xlsx")