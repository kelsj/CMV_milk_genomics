library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(Hmisc)
library(ggcorrplot)
library(ggrepel)
library(readxl)
theme_set(theme_bw())

setwd("~/Documents/postdoc/milk_genomics/CMV/manuscript_drafts/code_to_share/CMV_milk_genomics/R/")

# load data - CMV aligned reads per sample
# Table S1
cmv.reads = read_excel("~/Documents/postdoc/milk_genomics/CMV/manuscript_drafts/supp_tables/CMV_milk_genomics_suppTables_v2.xlsx",
                      sheet=1)

## Figure 1B: proportion samples CMV+/-
ggplot(cmv.reads %>% mutate(grp = "all")) + 
  geom_bar(aes(x=grp,fill=milk.cmv.status,y=)) +
  scale_fill_manual(values=c("orange","slateblue")) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + xlab("") + ylab("Number of samples") +
  scale_y_continuous(expand = c(0,0), breaks=c(0,50,100,150,200,250)) + 
  scale_x_discrete(expand = c(0,0)) + 
  labs(fill="")
ggsave("../figs/fig1B.png",width=2,height=3)

## Figure 1C: proportion CMV-mapped reads for CMV+ samples
ggplot(cmv.reads %>% filter(cmv.prop>0)) + 
  geom_histogram(aes(x=cmv.prop,y=..count..),fill="slateblue",bins=15) +
  scale_x_log10(breaks=c(1e-7,1e-6,1e-5,1e-4)) + 
  xlab("Proportion CMV-mapped reads") + ylab("Count of milk samples")
ggsave("../figs/fig1C.png",width=3,height=3)

## Figure 1D: density of CMV-mapped reads across CMV genome
both.align = readRDS("../data/SMS_WGS_CMValigned_readPos.rds")

ggplot(both.align %>% mutate(pos.kb = pos1/1000)) + geom_density(aes(x=pos.kb,y=..density..),bw=0.5,color="slateblue") + 
  xlab("Position on CMV genome (kb)") + ylab("Density of CMV-mapped reads")
ggsave("../figs/fig1D.png",width=4,height=3)


## Figure 1E: confusion matrix of CMV calls from qPCR vs. shotgun sequencing






