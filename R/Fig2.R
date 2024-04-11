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

## load differential gene expression results - Table S4
table.s4 = read_excel("../data/CMV_milk_genomics_suppTables_v4.xlsx",sheet=4)

## Figure 2A
# calculated distribution expected P-values (uniform distr.)
table.s4$Pexp = (1:nrow(table.s4))/nrow(table.s4)
# add color label for sig genes
table.s4 = table.s4 %>% mutate(lab=ifelse(padj<0.05,gene_name,NA), fdr.pass = ifelse(padj<0.05,"yes","no"))
# make qqplot
ggplot(table.s4) + geom_abline() + 
  geom_point(aes(x=-log10(Pexp), y=-log10(pvalue), color=fdr.pass)) +
  geom_text_repel(aes(x=-log10(Pexp), y=-log10(pvalue),label=lab),fontface="italic",max.overlaps = 15,box.padding = 0.1) +
  scale_color_manual(values=c("darkgray","magenta")) +
  theme(legend.position="None") + 
  xlab("Expected P-value (-log10)") + ylab("Observed P-value (-log10)")
ggsave("../figs/fig2A.png",width=4,height=4)


## Figure 2B - volcano plot
ggplot(table.s4) + 
  geom_point(aes(x=log2FoldChange, y=-log10(pvalue), color=fdr.pass)) +
  geom_text_repel(aes(x=log2FoldChange, y=-log10(pvalue),label=lab),fontface="italic") +
  xlab("Log2 fold change in CMV+") + ylab("P-value (-log10)") +
  scale_color_manual(values=c("darkgray","magenta")) +
  theme(legend.position="None")
ggsave("../figs/fig2B.png",width=4,height=4)


## Figure 2C
# load data - immune/epi cell expr. from Nyquist et al. vs. CMV logFC (our data)
fig2c.dat = readRDS("../data/immune_epithelial_cellExpr_vs_CMVlogFC.rds")
# plot
ggplot(fig2c.dat, aes(x=log2FoldChange, y=diff)) + geom_point() +
  geom_text_repel(aes(x=log2FoldChange,y=diff,label=GENE),fontface="italic",max.overlaps = 5,
                  point.padding = 0.1, force=1.5) + 
  xlab("Log2 fold change in CMV+") +
  ylab("Difference in scaled gene expression\n(immune cells - luminal cells)")
ggsave("../figs/fig2C.png",width=4.5,height=4)




