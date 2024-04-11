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
cmv.reads = read_excel("../data/CMV_milk_genomics_suppTables_v4.xlsx",sheet=1)

## Figure 1B: proportion samples CMV+/- 
ggplot(cmv.reads %>% mutate(grp = "all")) + 
  geom_bar(aes(x=grp,fill=shotgun.cmv.status,y=)) +
  scale_fill_manual(values=c("orange","slateblue")) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + xlab("") + ylab("Number of samples") +
  scale_y_continuous(expand = c(0,0), breaks=c(0,50,100,150,200,250)) + 
  scale_x_discrete(expand = c(0,0)) + 
  labs(fill="")
ggsave("../figs/fig1B.png",width=2,height=3)

## Figure 1C: proportion CMV-mapped reads for CMV+ samples
ggplot(cmv.reads %>% filter(shotgun.cmv.prop>0)) + 
  geom_histogram(aes(x=shotgun.cmv.prop,y=..count..),fill="slateblue",bins=15) +
  scale_x_log10(breaks=c(1e-7,1e-6,1e-5,1e-4)) + 
  xlab("Proportion CMV-mapped reads") + ylab("Count of milk samples")
ggsave("../figs/fig1C.png",width=3,height=3)

## Figure 1D: density of CMV-mapped reads across CMV genome
both.align = readRDS("../data/SMS_WGS_CMValigned_readPos.rds")

ggplot(both.align %>% mutate(pos.kb = pos1/1000)) + geom_density(aes(x=pos.kb,y=..density..),bw=0.5,color="slateblue") + 
  xlab("Position on CMV genome (kb)") + ylab("Density of CMV-mapped reads")
ggsave("../figs/fig1D.png",width=4,height=3)


## Figure 1E: confusion matrix of CMV calls from qPCR vs. shotgun sequencing
ggplot(cmv.reads %>% filter(qpcr.cmv.status!="NA")) + 
  geom_confmat(aes(x=qpcr.cmv.status,y=shotgun.cmv.status)) + 
  scale_fill_gradient(low="#FEFE62",high="#D35FB7") + 
  xlab("qPCR") + ylab("Shotgun sequencing (both types)")
ggsave("../figs/fig1E.png", width=4,height=3)


## Figure 1F: comparison of qPCR vs. shotgun viral load in CMV+ samples
ggplot(cmv.reads %>% mutate(qpcr.copies.ml = as.numeric(qpcr.copies.ml))) + 
  geom_point(aes(x=qpcr.copies.ml,y=shotgun.cmv.prop)) +
  scale_x_log10() + scale_y_log10() + 
  xlab("CMV copies/mL (qPCR)") +
  ylab("Proportion CMV-mapped reads\n(shotgun sequencing)")
ggsave("../figs/fig1F.png",width=3.5,height=3)





