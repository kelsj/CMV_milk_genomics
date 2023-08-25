library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(Hmisc)
library(ggcorrplot)
library(ggrepel)
library(readxl)
library(ggbeeswarm)
theme_set(theme_bw())

setwd("~/Documents/postdoc/milk_genomics/CMV/manuscript_drafts/code_to_share/CMV_milk_genomics/R/")

## Figure 3A: CMV status vs. milk kynurenine abundance

# load CMV status & addl metadata
metadata = read_excel("../data/CMV_milk_genomics_suppTables_v2.xlsx", sheet=1, na = "NA")

# load metabolite abundances
metabs = read.table("../data/metabolite_data.txt",header=T,sep="\t") %>% as_tibble()

#get residual kynurenine after regressing covariates
# 142 with all data points to calculate residuals
kyn.meta.resid = metabs %>% select(indID, KYNURENINE, TRYPTOPHAN) %>% mutate(KYN = scale(log(KYNURENINE+1))) %>%
  left_join(metadata %>% select(indID, milk.cmv.status, Center, parity, matage, mat_bmi, mat.white, avg_totalhei, gdm)) %>%
  mutate_at(vars(avg_totalhei,matage,mat_bmi),scale) %>% mutate(Center=as.factor(Center))
kyn.meta.resid$KYNresid = lm(KYN ~ Center + parity + matage + mat_bmi + mat.white + avg_totalhei + gdm, 
                             data=kyn.meta.resid)$residuals

ggplot(kyn.meta.resid, aes(x=milk.cmv.status,y=KYNresid,color=milk.cmv.status)) + 
  geom_boxplot(fill=NA,outlier.shape=NA,size=0.5) +
  geom_beeswarm(cex=3,size=0.5,alpha=0.5) +
  scale_color_manual(values=c("orange","slateblue")) +
  xlab("") + ylab("Milk kynurenine\n(normalized abundance)") + theme(legend.position="None")
ggsave("../figs/fig3A.png",width=2,height=2.5)


## Figure 3B: IDO1 TPM vs. cmv status (ENSG00000131203)
# load TPM data
tpm = read.table("../data/milk_GeneExpr_TPM.txt", header=T, sep="\t") %>% as_tibble()
ido.dat = tpm %>% filter(gene_id=="ENSG00000131203") %>% gather(key="indID",value="IDO1",2:222) %>%
  left_join(metadata %>% select(indID, milk.cmv.status)) %>% filter(!is.na(milk.cmv.status))
my_y_title = bquote("Milk"~italic("IDO1")~"expression (TPM)")

ggplot(ido.dat, aes(y=IDO1,x=milk.cmv.status,color=milk.cmv.status)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_beeswarm(size=0.5,alpha=0.6,cex=2) +
  scale_color_manual(values=c("orange","slateblue")) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  ylab(my_y_title) + xlab("") + theme(legend.position="None")
ggsave("../figs/fig3B.png",width=2,height=2.5)

## Figure 3D: IDO1 expression vs. kynurenine/tryptophan ratio

# calculate KTR ratio and get residuals after regression on covars
kyn.meta.resid = kyn.meta.resid %>% mutate(KTR=KYNURENINE/TRYPTOPHAN, KTR=scale(log(KTR)))
kyn.meta.resid$KTRresid = lm(KTR ~ Center + parity + matage + mat_bmi + mat.white + avg_totalhei + gdm, 
                             data=kyn.meta.resid)$residuals
# add IDO1 tpm to metabolite data
kyn.meta.resid = kyn.meta.resid %>% left_join(ido.dat %>% select(indID, IDO1))

# plot
my_x_title = bquote("Milk"~italic("IDO1")~"expression (TPM, log10 scale)")
ggplot(kyn.meta.resid, aes(x=IDO1,y=KTRresid,color=milk.cmv.status)) + geom_smooth(method="lm",se=F) + 
  geom_point(size=0.75) + scale_color_manual(values=c("orange","slateblue")) + 
  xlab("Milk IDO1 expression\n(TPM, log10 scale)") +
  ylab("Milk kynurenine/tryptophan \n(normalized ratio of abundances)") +
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10)) +
  theme(legend.pos="None", axis.title.x=element_text(size=rel(1.2)),axis.title.y=element_text(size=rel(1.2)))
ggsave("../figs/fig3D.png",width=3.2,height=3.3)





