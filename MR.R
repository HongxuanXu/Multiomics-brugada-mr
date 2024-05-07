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

setwd(dir="")  #set your wd

#calculate F statiscs and r2
get_f<-function(dat,F_value=10){
  log<-is.na(dat$eaf.exposure)
  log<-unique(log)
  if(length(log)==1)
  {if(log==TRUE){
    print("lacking eaf, cannot calculate F statistics")
    return(dat)}
  }
  if(is.null(dat$beta.exposure[1])==T || is.na(dat$beta.exposure[1])==T){print("lacking beta, cannot calculate F statistics")
    return(dat)}
  if(is.null(dat$se.exposure[1])==T || is.na(dat$se.exposure[1])==T){print("lacking se, cannot calculate F statistics")
    return(dat)}
  if(is.null(dat$samplesize.exposure[1])==T || is.na(dat$samplesize.exposure[1])==T){print("lacking sample size ，cannot calculate F statistics")
    return(dat)}
  
  
  if("FALSE"%in%log && is.null(dat$beta.exposure[1])==F && is.na(dat$beta.exposure[1])==F && is.null(dat$se.exposure[1])==F && is.na(dat$se.exposure[1])==F && is.null(dat$samplesize.exposure[1])==F && is.na(dat$samplesize.exposure[1])==F){
    R2<-(2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$beta.exposure^2))/((2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$beta.exposure^2))+(2*(1-dat$eaf.exposure)*dat$eaf.exposure*(dat$se.exposure^2)*dat$samplesize.exposure))
    F<- (dat$samplesize.exposure-2)*R2/(1-R2)
    dat$R2<-R2
    dat$F<-F
    dat<-subset(dat,F>F_value)
    return(dat)
  }
}

outcomeid<-fread("your directory") #read outcome data into R

FileNames <- list.files(paste0(getwd()),pattern=".txt") #get exposure list for loop
exp_dat_ids <- FileNames
exps <- FileNames

dir.create(path = "/Users/xuhongxuan/Documents/brugada/ukb")  #create results directory

####UKBPPP as an example######

for (qaq in 1:length(exp_dat_ids)) { 
  exp_dat_id <- exp_dat_ids[qaq]
  
  exp <- exps[qaq]
  
  d3<- fread(paste0(getwd(),"/",FileNames[qaq]))
  
  d3 <- d3[!d3$SNP == '',]

  d3<-subset(d3,d3$`Pval`<1e-6)
  
  if(nrow(d3) == 0 ){next}
  

    
    d3$PHENO <- FileNames[qaq]
    
    d3 <- d3[!is.na(d3$SNP),]
  
    if(nrow(d3) > 0){
      
      
      exp_data<-format_data(as.data.frame(d3),
                            type="exposure",
                            chr_col = "CHROM",
                            phenotype_col = "PHENO",
                            snp_col = "SNP",
                            beta_col = "BETA",
                            se_col = "SE",
                            pval_col = "Pval",
                            samplesize_col = "N",
                            eaf_col = "A1FREQ",
                            effect_allele_col = "ALLELE1",
                            other_allele_col = "ALLELE0",
                            pos_col = "GENPOS")

        
        skip_to_next <- FALSE
        
        tryCatch(        clump <-  ld_clump_local(dat = tibble(rsid = exp_data$SNP,
                                                               pval = exp_data$pval.exposure),
                                                  clump_r2 = 0.1,
                                                  clump_p = 1,
                                                  clump_kb = 1000,
                                                  plink_bin = "/Users/xuhongxuan/Documents/plink/plink",
                                                  bfile = "/Users/xuhongxuan/Documents/g1000_eur/g1000_eur"), error = function(e) { skip_to_next <<- TRUE})
        

        if(skip_to_next) { next }   
        exp_data <- exp_data |> 
          filter(SNP %in% 
                   clump$rsid)
      
      
      
      outcome  <- outcomeid |> filter(SNP %in% exp_data$SNP)
      
      if(nrow(ouecoome)>0){  
        dat <- TwoSampleMR::harmonise_data(
          exposure_dat = exp_data,
          outcome_dat = outcome)
        
        ####回文的直接去除
        dat <-subset(dat,mr_keep==TRUE)
        
        
        
        if(nrow(dat)>0){ 
          dat <- get_f(dat, F_value = 10)
          
          
          res <- TwoSampleMR::mr(dat, method_list= c("mr_ivw" ,
                                                     "mr_wald_ratio"))
          
          print(paste0(exp,"_No.SNP_",res$nsnp[1]))
          
          results <- TwoSampleMR::generate_odds_ratios(res)
          
          results$estimate <- paste0(
            format(round(results$or, 2), nsmall = 2), " (", 
            format(round(results$or_lci95, 2), nsmall = 2), "-",
            format(round(results$or_uci95, 2), nsmall = 2), ")")
          resdata <- dat
          write.csv(dat,file = paste0("/Users/xuhongxuan/Documents/brugada/ukb/dat/",exp,"-dat.csv"))
          
          names(resdata)
          
          
          write.csv(results,
                    paste0("/Users/xuhongxuan/Documents/brugada/ukb/res/",exp,"-res.csv"))
          
          
          if(nrow(dat)>2){
            res_hete <- TwoSampleMR::mr_heterogeneity(dat)
            res_plei <- TwoSampleMR::mr_pleiotropy_test(dat)
            res_leaveone <- mr_leaveoneout(dat)  # 
            
            ######steiger检验######
            dat$r.exposure <- get_r_from_bsen(b = dat$beta.exposure,
                                              dat$se.exposure,
                                              dat$samplesize.exposure)
            dat$r.outcome <- get_r_from_bsen(b = dat$beta.outcome,
                                             dat$se.outcome,
                                             dat$samplesize.outcome)
            res_steiger <- mr_steiger(
              p_exp = dat$pval.exposure,
              p_out = dat$pval.outcome,
              n_exp = dat$samplesize.exposure,
              n_out = dat$samplesize.outcome,
              r_exp = dat$r.exposure,
              r_out = dat$r.outcome
            )
            res_steiger <- directionality_test(dat)
            
            
            
            
            
            #         res_presso <- TwoSampleMR::run_mr_presso(dat,
            #                                             NbDistribution = 100)
            # [["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]
            #       sink(paste0("/Users/xuhongxuan/Documents/brugada/ukb/",exp,"_PRESSO.txt"),
            #                        append=FALSE,split = FALSE) 
            #       print(res_presso)
            #       sink()
            #       print(res_presso)
            
            
            
            
            p1 <- mr_scatter_plot(res, dat)
            p1[[1]]
            pdf(paste0("/Users/xuhongxuan/Documents/brugada/ukb/sca/",exp,"_scatter.pdf"))
            print(p1[[1]])
            dev.off()
            
            res_single <- mr_singlesnp(dat)
            p2 <- mr_forest_plot(res_single)
            pdf(paste0("/Users/xuhongxuan/Documents/brugada/ukb/forest/",exp,"_forest.pdf"))
            print(p2[[1]])
            dev.off()
            
            p3 <- mr_funnel_plot(res_single)
            pdf(paste0("/Users/xuhongxuan/Documents/brugada/ukb/funnel/",exp,"_funnel.pdf"))
            print(p3[[1]])
            dev.off()
            
            res_loo <- mr_leaveoneout(dat)
            pdf(paste0("/Users/xuhongxuan/Documents/brugada/ukb/lou/",exp,"_leave_one_out.pdf"))
            print(mr_leaveoneout_plot(res_loo))
            dev.off()
            res3 <- res[1:row_number(res),]
            
            
          

            # Main result 
            res4 <- as.data.frame(tidyr::pivot_wider(
              res3,names_from ="method",names_vary = "slowest",
              values_from = c("b","se","pval") ))
            # Heterogeneity statistics
            res_hete2 <- tidyr::pivot_wider(
              res_hete,names_from ="method",names_vary = "slowest",
              values_from = c("Q","Q_df","Q_pval") ) %>% 
              dplyr::select( -id.exposure,-id.outcome,-outcome,-exposure)
            # Horizontal pleiotropy
            res_plei2 <- dplyr::select(res_plei,
                                       egger_intercept,se,pval)
            ##steiger
            res_steiger2 <- dplyr::select(res_steiger,
                                          correct_causal_direction,steiger_pval)
            
            
            # Merge
            res_ALL <- cbind(res4, res_hete2, res_plei2,res_steiger2)
            
            write.csv(res_ALL,file = paste0("/Users/xuhongxuan/Documents/brugada/ukb/resall/",exp,".csv"))
          }}}}

fs <- list.files("/Users/xuhongxuan/Documents/brugada/ukb/resall", pattern = ".csv",full.names = TRUE) 
df <- map_dfr(fs, read.csv)
write.csv(df,"/Users/xuhongxuan/Documents/brugada/ukb/list/resall.csv")
