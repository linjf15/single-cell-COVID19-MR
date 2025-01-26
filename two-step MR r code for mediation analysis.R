library(TwoSampleMR)
library(dplyr)
library(MRInstruments)
library(MRPRESSO)
library(ieugwasr)
library(daff)
library(fdrtool)
library(coloc)
library(ggplot2)
library(tidyr)
library(locuscomparer)
library(tidyverse)
library(data.table)
library(gwasglue)
library(SCEnt)

#TWO-STEP MR for confounding factors 
#step 1 exposure upload
mvmr_iv <- read.csv("D:/covid19/mvmr_iv.csv") # all IVs from top MR results 
#smoking/obesity outcome selection: CigarettesPerDay(ieu-b-25), BMI ("ieu-b-40")
out_dat<- extract_outcome_data(snps =mvmr_iv$SNP,outcomes = "ieu-b-25") #for BMI as outcome, just replace "ieu-b-25" with "ieu-b-40"
out_dat_id<- "smoking" #for BMI as outcome, just replace "smoking" with "BMI"
#har_data
Mydata <- harmonise_data(
  exposure_dat=mvmr_iv,
  outcome_dat=out_dat,
  action= 2)

Mydata$rsq.exposure<- (get_r_from_bsen(b=Mydata$beta.exposure,
                                       Mydata$se.exposure,
                                       Mydata$samplesize.exposure))^2
Mydata$rsq.outcome<- (get_r_from_bsen(b=Mydata$beta.outcome,
                                      Mydata$se.outcome,
                                      Mydata$samplesize.outcome))^2
Mydata$F_statistics<-(Mydata$beta.exposure)^2/(Mydata$se.exposure)^2
Mydata<- subset(Mydata,Mydata$F_statistics>=10)

res <- mr(Mydata,method_list = c("mr_wald_ratio","mr_ivw"))
res<- res[order(res$pval),]
res$FDR<- p.adjust(res$pval,"BH")
res_final<- subset(res,res$pval<0.05) 

#Directionality test using MR Steiger filtering
st <- psych::r.test( 
  n = Mydata$samplesize.exposure, 
  n2 = Mydata$samplesize.outcome, 
  r12 = sqrt(Mydata$rsq.exposure), 
  r34 = sqrt(Mydata$rsq.outcome))

Mydata$steiger_dir <- Mydata$rsq.exposure > Mydata$rsq.outcome
Mydata$steiger_pval <- pnorm(-abs(st$z)) * 2 
Mydata_res<- left_join(res_final,Mydata,by=c("id.exposure","exposure"))

setwd("H:\\covid\\reply_to_ebio\\two step mr for smoking & obesity")
write.table(Mydata, paste0("Mydata ",out_dat_id,".txt"), row.names = FALSE, sep = '\t', quote = FALSE)
write.table(res_final, paste0("res_final ",out_dat_id,".txt"), row.names = FALSE, sep = '\t', quote = FALSE)
write.table(Mydata_res, paste0("Mydata_res ",out_dat_id,".txt"), row.names = FALSE, sep = '\t', quote = FALSE)
write.table(out_dat, paste0("out_dat ",out_dat_id,".txt"), row.names = FALSE, sep = '\t', quote = FALSE)
write.table(mvmr_iv, paste0("Exp_dat ",out_dat_id,".txt"), row.names = FALSE, sep = '\t', quote = FALSE)
write.table(res, paste0("res_all_method ",out_dat_id,".txt"), row.names = FALSE, sep = '\t', quote = FALSE)

#step 2 smoking/bmi as expo
Exp_dat<- extract_instruments("ieu-b-25",p1=0.001,p2=0.001,r2=0.01,kb=5000,clump = TRUE) #for BMI ivs selection, just replace "ieu-b-25" with "ieu-b-40"
expo_name<- "smoking"  #for BMI as expo, just replace "smoking" with "BMI"

w=1  
for(w in 1:4)
{
  setwd("D:\\covid19\\20231204\\OUTCOME")
  outcome_name <- list.files()  
  outcome_name<- outcome_name[w]
  outcome_id<- fread(file =paste0(outcome_name),sep=",")
  outcome_id<- as.data.frame(outcome_id)
  
  covid_Outcome_DAT<- format_data(
    outcome_id,
    type = "outcome",
    snps =Exp_dat$SNP,
    eaf_col = "effect_allele_freq",
    pval_col = "p",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect.allele",
    other_allele_col = "other.allele",
    chr_col = "chr",
    pos_col = "pos",
    samplesize_col = "n")
  
  Mydata <- harmonise_data(
    exposure_dat=Exp_dat,
    outcome_dat=covid_Outcome_DAT,
    action= 2) 
  
  Mydata$rsq.exposure<- (get_r_from_bsen(b=Mydata$beta.exposure,
                                         Mydata$se.exposure,
                                         Mydata$samplesize.exposure))^2
  Mydata$rsq.outcome<- (get_r_from_bsen(b=Mydata$beta.outcome,
                                        Mydata$se.outcome,
                                        Mydata$samplesize.outcome))^2
  Mydata$F_statistics<-(Mydata$beta.exposure)^2/(Mydata$se.exposure)^2
  Mydata<- subset(Mydata,Mydata$F_statistics>=10)
  
  res <- mr(Mydata)
  
  #Directionality test using MR Steiger filtering
  st <- psych::r.test( 
    n = Mydata$samplesize.exposure, 
    n2 = Mydata$samplesize.outcome, 
    r12 = sqrt(Mydata$rsq.exposure), 
    r34 = sqrt(Mydata$rsq.outcome))
  
  Mydata$steiger_dir <- Mydata$rsq.exposure > Mydata$rsq.outcome
  Mydata$steiger_pval <- pnorm(-abs(st$z)) * 2 
  
  setwd("D:\\covid19\\20231204\\step 2 res")
  write.table(Mydata, paste0("Mydata ",expo_name," to ",outcome_name,".txt"), row.names = FALSE, sep = '\t', quote = FALSE)
  write.table(res, paste0("res ",expo_name," to ",outcome_name,".txt"), row.names = FALSE, sep = '\t', quote = FALSE)
  write.table(Exp_dat, paste0("Exp_dat ",expo_name," to ",outcome_name,".txt"), row.names = FALSE, sep = '\t', quote = FALSE)
  write.table(covid_Outcome_DAT, paste0("Outcome_DAT ",expo_name," to ",outcome_name,".txt"), row.names = FALSE, sep = '\t', quote = FALSE)
 }

#mediation analysis for immune related trait
setwd("H:\\covid\\reply to ebio\\mediation_analysis_cell_count")
cell_count_id <- read_excel("H:/covid/reply_to_ebio/mediation_analysis_cell_count/cell count id.xlsx", 
                            sheet = "absolute")
i=1
for(i in 1:nrow(cell_count_id)){
  dir.create(as.character(cell_count_id[i,6]))
}

#step 1 CELL COUNT as outcome SELECTION
i=1
for(i in 1:nrow(cell_count_id)){
  
  id<- as.character(cell_count_id[i,5]) #trait gwas id 
  traits<- as.character(cell_count_id[i,4])  #trait name
  
  CELL_COUNT <- extract_outcome_data(
    snps = mvmr_iv$SNP,
    outcomes = id)
  
  Mydata <- harmonise_data(
    exposure_dat=mvmr_iv,
    outcome_dat=CELL_COUNT,
    action= 2)
  
  Mydata$rsq.exposure<- (get_r_from_bsen(b=Mydata$beta.exposure,
                                         Mydata$se.exposure,
                                         Mydata$samplesize.exposure))^2
  Mydata$rsq.outcome<- (get_r_from_bsen(b=Mydata$beta.outcome,
                                        Mydata$se.outcome,
                                        Mydata$samplesize.outcome))^2
  Mydata$F_statistics<-(Mydata$beta.exposure)^2/(Mydata$se.exposure)^2
  Mydata<- subset(Mydata,Mydata$F_statistics>=10)
  
  res <- mr(Mydata,method_list = c("mr_wald_ratio","mr_ivw"))
  res<- res[order(res$pval),]
  res$FDR<- p.adjust(res$pval,"BH")
  res_final<- subset(res,res$FDR<0.05) 
  
  #Directionality test using MR Steiger filtering
  st <- psych::r.test( 
    n = Mydata$samplesize.exposure, 
    n2 = Mydata$samplesize.outcome, 
    r12 = sqrt(Mydata$rsq.exposure), 
    r34 = sqrt(Mydata$rsq.outcome))
  
  Mydata$steiger_dir <- Mydata$rsq.exposure > Mydata$rsq.outcome
  Mydata$steiger_pval <- pnorm(-abs(st$z)) * 2 
  Mydata_res<- left_join(res_final,Mydata,by=c("id.exposure","exposure"))
  setwd("H:\\covid\\reply_to_ebio\\mediation_analysis_cell_count\\2025 new mediation analysis\\step1 results")
  write.table(Mydata, paste0("Mydata ","(exp_to_",traits,")",".txt"), row.names = FALSE, sep = '\t', quote = FALSE)
  write.table(res, paste0("res ","(exp_to_",traits,")",".txt"), row.names = FALSE, sep = '\t', quote = FALSE)
  write.table(res_final, paste0("res_final ","(exp_to_",traits,")",".txt"), row.names = FALSE, sep = '\t', quote = FALSE)
  write.table(Mydata_res, paste0("Mydata_res ","(exp_to_",traits,")",".txt"), row.names = FALSE, sep = '\t', quote = FALSE)
  write.table(mvmr_iv, paste0("Exp_dat ","(exp_to_",traits,")",".txt"), row.names = FALSE, sep = '\t', quote = FALSE)
  write.table(CELL_COUNT, paste0("Outcome_DAT ","(exp_to_",traits,")",".txt"), row.names = FALSE, sep = '\t', quote = FALSE)
}

#step 2 
i=1
for(i in 1:nrow(cell_count_id)){
  id<- as.character(cell_count_id[i,5])
  traits<- as.character(cell_count_id[i,4])
  
  CELL_COUNT<- extract_instruments(id,p1=0.001,p2=0.001,r2=0.01,kb=10000,clump = TRUE)
  
  w=1
  for(w in 1:4)
  {
    setwd("D:\\covid19\\20231204\\OUTCOME")
    outcome_name <- list.files()  
    outcome_name<- outcome_name[w]
    outcome_id<- fread(file =paste0(outcome_name),sep=",")
    outcome_id<- as.data.frame(outcome_id)
    
    covid_Outcome_DAT<- format_data(
      outcome_id,
      type = "outcome",
      snps =CELL_COUNT$SNP,
      eaf_col = "effect_allele_freq",
      pval_col = "p",
      beta_col = "beta",
      se_col = "se",
      effect_allele_col = "effect.allele",
      other_allele_col = "other.allele",
      chr_col = "chr",
      pos_col = "pos",
      samplesize_col = "n")
    
    Mydata <- harmonise_data(
      exposure_dat=CELL_COUNT,
      outcome_dat=covid_Outcome_DAT,
      action= 2) 
    
    Mydata$rsq.exposure<- (get_r_from_bsen(b=Mydata$beta.exposure,
                                           Mydata$se.exposure,
                                           Mydata$samplesize.exposure))^2
    Mydata$rsq.outcome<- (get_r_from_bsen(b=Mydata$beta.outcome,
                                          Mydata$se.outcome,
                                          Mydata$samplesize.outcome))^2
    Mydata$F_statistics<-(Mydata$beta.exposure)^2/(Mydata$se.exposure)^2
    Mydata<- subset(Mydata,Mydata$F_statistics>=10)
    
    res <- mr(Mydata)
    
    #Directionality test using MR Steiger filtering
    st <- psych::r.test( 
      n = Mydata$samplesize.exposure, 
      n2 = Mydata$samplesize.outcome, 
      r12 = sqrt(Mydata$rsq.exposure), 
      r34 = sqrt(Mydata$rsq.outcome))
    
    Mydata$steiger_dir <- Mydata$rsq.exposure > Mydata$rsq.outcome
    Mydata$steiger_pval <- pnorm(-abs(st$z)) * 2 
    
    setwd("H:\\covid\\reply_to_ebio\\mediation_analysis_cell_count\\2025 new mediation analysis\\step2 results")
    write.table(Mydata, paste0("Mydata ",traits," to ",outcome_name,".txt"), row.names = FALSE, sep = '\t', quote = FALSE)
    write.table(res, paste0("res ",traits," to ",outcome_name,".txt"), row.names = FALSE, sep = '\t', quote = FALSE)
    write.table(CELL_COUNT, paste0("Exp_dat ",traits," to ",outcome_name,".txt"), row.names = FALSE, sep = '\t', quote = FALSE)
    write.table(covid_Outcome_DAT, paste0("Outcome_DAT ",traits," to ",outcome_name,".txt"), row.names = FALSE, sep = '\t', quote = FALSE)
  } 
}
