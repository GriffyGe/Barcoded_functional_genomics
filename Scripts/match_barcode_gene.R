# Merge barcodes with genomic DNA with gene info ====
library(tidyverse)
readId_barcodes <- read_delim("Example_file/table_id_bc.txt",
                              delim = "\t",
                              col_names = c("qseqid", "barcode"))
gdna_gene_info <- read_delim("Example_file/best_hit_match_gene_res.txt", 
                              col_names = c("qseqid", "sseqid", "insertion_site",
                                            "pident", "length", "begin", "end", "desc",
                                            "gene_ID", "gene_name"))

# Make the format of column "qseqid" of the two dataframe same.
#@DP8480006254TRL1C021R05601207665/1 -> DP8480006254TRL1C021R05601207665/1
#@A00265:1138:HCFHGDSX5:4:1122:32298:10394 1:N:0:AAAGATAC+TCTGAAAC -> A00265:1138:HCFHGDSX5:4:1122:32298:10394
gdna_bc$qseqid <- str_match(gdna_bc$qseqid, pattern = "@(.*)[ /]")[,2]

bc_gene_info <- merge(gdna_bc, match_gene_info, by="qseqid")

#step9:同一barcode插入相同基因的不同reads合并====
#(由于基因存在重复序列导致比对出现重复结果的已经删除，所以剩下的只要gene_id一致即可)
list_gdna_bc_info <- split(gdna_bc_gene_info, gdna_bc_gene_info$barcode)
result <- data.frame()
not_unique <- c()
for (i in 1:length(list_gdna_bc_info)) {
  bc <- names(list_gdna_bc_info[i]) #获得barcode
  df <- list_gdna_bc_info[[i]] #获得上述barcode对应的gdna以及其他gene_info
  gene_num <- length(split(df, df$gene_ID)) #获得对应基因的个数
  if (gene_num == 1){
    sub <- df[1,-3]
    result <- rbind(result, sub)
  }
  else{
    not_unique <- append(not_unique, i)
  }
}

write_delim(result,"139/139_unique_barcode_and_gene_blastn.txt")

#生成表格查看not unique的序列的情况(excel数据透视表num)
list_not_unique <- list_gdna_bc_info[not_unique]
library(plyr)
library(rio)
result_not_unique <- ldply(list_not_unique, data.frame)
export(result_not_unique,"139/139_not_unique_barcode_and_gene_blastn.xlsx")


#step10:同一barcode对应不同gene，在excel中利用数据透视表查看筛选====
#根据num，参考文献中的unique的标准进行筛选
#第二次得到的not unique的barcode数量太多，需要写代码进行筛选
library(rio)
data <- import("139/139_not_unique_barcode_and_gene_blastn.xlsx")

gene_dict <- list()
for (i in 1:nrow(data)){
  if (!grepl("XCQ", data[i,1])){  #XAC for 306, XCQ for CQ13
    barcode <-  data[i,1]
    total <-  data[i,2]
    gene_dict$barcode <- data.frame()  #先创建元素后命名    
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


not_unique_data <- import("139/139_not_unique_barcode_and_gene_blastn.xlsx",sheet=2)

usable_result <- data.frame()

for (i in 1:length(gene_dict)){
  usable_barcode <- names(gene_dict[i])
  df <- arrange(gene_dict[[i]], desc(match_frq))
  #unique:1.match time>=10; 2.match_frq>=0.75; 3.2nd most match frq <= 1/8 * 1st most match frq
  if ((df[1,2]>=10) & (df[1,3]>=0.75) & (df[2,3]<=0.125*df[1,3])){
    usable_id = df[1,1]
    sub <- not_unique_data %>% 
      filter(barcode == usable_barcode & gene_ID == usable_id) %>% 
      distinct(barcode, .keep_all = TRUE) %>% 
      select(qseqid,barcode,sseqid,insertion_site,pident,length,
             begin,end,desc,gene_ID,gene_name)
    usable_result <- rbind(usable_result, sub)
  }
}

unique_data <- read_delim("139/139_unique_barcode_and_gene_blastn.txt")
final_unique_data <- rbind(unique_data, usable_result)
write_delim(final_unique_data,"139/139_blastn_final_unique_barcode_and_gene.txt")

#step11:找到central insertions (10%~90%)====
result <- import("139/139_blastn_final_unique_barcode_and_gene.txt") 
for (i in 1:nrow(result)){
  #都使用sstart，如果sstart > send，说明是转座子反着插入
  result$f[i] <- ((result$insertion_site[i] - result$begin[i] + 1)/ (result$end[i] - result$begin[i] + 1))
}

result2 <-  result[result$f >= 0.1 & result$f <= 0.9,]

export(result2,"139/139_blastn_unique_in_central_barcode.xlsx")