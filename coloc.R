library(TwoSampleMR)
library(data.table)
library(tidyverse)
library(ieugwasr)
library(magrittr)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(BSgenome)
library(stringr)
library(VariantAnnotation)
library(gassocplot2)
library(gwasglue)
library(coloc)

ensembl <-  useEnsembl(biomart='ensembl', dataset="hsapiens_gene_ensembl")

for (qqq in 1 : length(exp_dat_ids)) { 
  exp_dat_id <- exp_dat_ids[qqq]
  exp <- exps[qqq]
  
  d1<- fread(paste0(getwd(),"/",FileNames[qqq])) |> as.data.frame()
  
  d1 <- d1[!is.na(d1$SNP),]
  
  symbol <- toupper(sub("\\_.*", "", exp_dat_id))
  
  range <- getBM(attributes=c( 'start_position','chromosome_name'),
                 filters=c('hgnc_symbol'),
                 values=list(symbol),
                 mart=ensembl)
  if(nrow(range) > 0){
    
  data1 <- d1 %>% filter(CHROM==range$chromosome_name,GENPOS >= range$start_position-1000000,GENPOS <= range$start_position+1000000)

  data2 <- afc %>% filter(chr==range$chromosome_name,pos >= range$start_position-1000000, pos <= range$start_position+1000000)

  data1 <- data1 %>% dplyr::select(SNP, CHROM, GENPOS, ALLELE1, ALLELE0, BETA, A1FREQ,N,SE)
    
      if(nrow(data1) == 0 ){next}              

  data2 <- afc %>% dplyr::select(rsids,chr , pos, alt, ref, beta, af_alt,n,sebeta)
  
  colnames(data1) <- c("SNP", 'CHR', "BP", 'A1', 'A2', 'BETA','MAF','N', 'SE')
  colnames(data2) <- c("SNP", 'CHR', "BP", 'A1', 'A2', 'BETA','MAF','N', 'SE')
  
  data <- merge(data1,data2,by="SNP")
  data <- data[!duplicated(data$SNP),]
  
  data <- data %>% filter((A1.x==A1.y&A2.x==A2.y)|(A1.x==A2.y&A2.x==A1.y)) 
  data <- data %>% mutate(BETA.y = ifelse(A1.x==A1.y,BETA.y,-BETA.y))
  
  
  data$VAR.x <- data$SE.x^2
  data$VAR.y <- data$SE.y^2
  data <- data[data$VAR.x!=0 & data$VAR.y!=0 ,]
  
  data1 <- data[,c("BETA.x","VAR.x","SNP")]
  data2 <- data[,c("BETA.y","VAR.y","SNP","MAF.y","N.y")]
  colnames(data1) <- c("beta","varbeta","snp")
  colnames(data2) <- c("beta","varbeta","snp","MAF","N")
  data1 <- as.list(data1)
  data2 <- as.list(data2)
  
  data1$type <- "cc"
  data2$type <- "quant"
  
  res <- coloc.abf(data1,data2,p1=1e-4,p2=1e-4,p12=1e-5)
  
  write.csv(res$results,
            paste0("/Users/xuhongxuan/Documents/cy_ukb_finn_coloc/",symbol,"_coloc.csv"))
  if(max(res$results$SNP.PP.H4)>0.8){
  write.csv(max(res$results$SNP.PP.H4),
            paste0("/Users/xuhongxuan/Documents/cy_ukb_finn_coloc/max/",symbol,"_colocmax.csv"))}
  
}}
