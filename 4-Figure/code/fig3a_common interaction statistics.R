library(tidyverse)
library(readxl)
library(magrittr)

# Updated from original submission 08/14/24 SS

## ------- Input -------
setwd("~/Documents/NCI/anti_pdl1_project/documents/github/nature communication revision 2/4. Figure/Code") #CHANGE
source('input.R') 

## ------- Sample random interactions -------
# ("size" is the # of cci as inferred is from IRIS in each cohort)
set.seed(20892)
lirics_cci = list(gide=sample(colnames(ICB_input[[1]]$Gide_pre), size=88), 
                  liu=sample(colnames(ICB_input[[1]]$Liu_pre), size=51), 
                  riaz=sample(colnames(ICB_input[[1]]$Riaz_pre), size=125), 
                  puch=sample(colnames(ICB_input[[1]]$Puch_pre), size=128), 
                  auslander=sample(colnames(ICB_input[[1]]$Auslander_pre), size=95))
                  
lirics_cci = unlist(lirics_cci)
background <- data.frame(table(lirics_cci))

## load IRIS interactions
AR_LIRICS = feature_RDI

## ------- Make CELL LR DF -------
Feature_df = AR_LIRICS[c(2,3,5,6,8)] %>% unlist(., recursive=F, use.names=T) %>% plyr::ldply(., cbind) %>% set_colnames(c('cohort','interaction')) %>%
  mutate(cohort = gsub("\\.(Liu_pre|Riaz_pre|Gide_pre|Auslander_pre|Puch_pre)", "", .$cohort)) %>% data.frame(.,str_split_fixed(.$interaction, '\\_', 4)) %>%
  set_colnames(c('cohorts','interaction','Lcell','Rcell','L','R')) %>% mutate(GenePair=paste(L,R,sep='_')) 

## save interaction df
for_paper = Feature_df %>% group_by(interaction, Lcell, Rcell, L, R, GenePair) %>% distinct(.) %>% dplyr::summarise(count = n())
for_paper = for_paper %>% set_colnames(c('interaction','ligand cell','receptor cell','ligand gene','receptor gene','ligand-receptor pair', '# of cohort models'))

## ------- statistics for paper -------
overlap_fg = for_paper %>% subset(`# of cohort models` >= 2) %>% nrow(.)
notoverlap_fg = nrow(for_paper) - overlap_fg

overlap_bg = background %>% subset(Freq >= 2) %>% nrow(.)
notoverlap_bg = nrow(background) - overlap_bg

### 
mat = matrix(data=c(overlap_fg,notoverlap_fg,overlap_bg,notoverlap_bg),nrow=2)
mat
x = fisher.test(mat, alternative = 'g')
x$p.value
x$estimate
