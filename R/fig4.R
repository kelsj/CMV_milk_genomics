library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(Hmisc)
library(ggcorrplot)
library(ggrepel)
library(readxl)
library(ggbeeswarm)
library(lmerTest)
theme_set(theme_bw())

setwd("~/Documents/postdoc/milk_genomics/CMV/manuscript_drafts/code_to_share/CMV_milk_genomics/R/")


## load metadata
metadata = read_excel("../data/CMV_milk_genomics_suppTables_v4.xlsx", sheet=1, na = "NA")

## Figure 4A: 1-month infant microbiome PC3 vs. milk CMV status
## calculate PC3 resid w/ covariates
# add metadata
pc3.resid = metadata %>% select(indID, delivery_cat, parity, matage, inf.white, mat_bmi, avg_totalhei,
                                gdm, gbs, fecal_1mo_site, inf.microbiome.1mo.PC3, shotgun.cmv.status) %>% na.omit()
# calculate resid PC3
pc3.resid$pc3.resid = lm(inf.microbiome.1mo.PC3 ~ delivery_cat + parity + matage + inf.white + mat_bmi + avg_totalhei + gdm +
                           gbs + fecal_1mo_site, data=pc3.resid)$residuals

ggplot(pc3.resid, aes(x=shotgun.cmv.status, y=pc3.resid,color=shotgun.cmv.status)) + 
  geom_boxplot(fill=NA,outlier.shape=NA) +
  geom_beeswarm(cex=3,size=0.5,alpha=0.5) +
  scale_color_manual(values=c("orange","slateblue")) +
  xlab("") + ylab("Infant fecal microbiome PC3") + theme(legend.position="None")
ggsave("../figs/fig4A.png",width=2,height=2.5)


## Figure 4B: effect estimates of CMV status on 1 & 6month microbiome, modeled together in LMM
lmm.res = read_excel("../data/CMV_milk_genomics_suppTables_v4.xlsx", sheet=10, na = "NA")

## dot plot for 2tmpt results
# plot pw sig for quant 0
dotplot.dat = lmm.res %>% filter(q.both<0.05)
dotplot.dat$taxon = str_replace_all(dotplot.dat$taxon,"s__","")
dotplot.dat$taxon = str_replace_all(dotplot.dat$taxon,"_"," ")
dotplot.dat$taxon = factor(dotplot.dat$taxon, levels = rev(dotplot.dat %>% arrange(p.both) %>% pull(taxon)))

ggplot(dotplot.dat) + geom_point(aes(x=taxon,y=est.both)) + 
  geom_segment(aes(x=taxon,xend=taxon,y=est.both-1.96*se.both,yend=est.both+1.96*se.both)) +
  coord_flip() + geom_hline(aes(yintercept=0),linetype="dashed") +
  theme(axis.text.y=element_text(face="italic"),legend.pos="None") + 
  ylab("Estimated effect of CMV+") + xlab("")
ggsave("../figs/fig4B.png",width=4,height=2.5)


## Figure 4C: milk CMV status vs. infant B. infantis at 1 & 6 months
# load 1 & 6mo taxon CLR
inf.1mo.taxa = read.table("../data/infFecalMicrobiome_1month_taxonCLR.txt",header=T) %>% as_tibble()
inf.6mo.taxa = read.table("../data/infFecalMicrobiome_6month_taxonCLR.txt",header=T) %>% as_tibble()
# extract b.infantis, scale CLR, add collection site for 1/6mo samples
b.inf = bind_rows(inf.1mo.taxa %>% select(indID, s__Bifidobacterium_infantis) %>% 
                    mutate(mytax = scale(s__Bifidobacterium_infantis), tmpt="1 month") %>%
                    left_join(metadata %>% select(indID, fecal_1mo_site) %>% rename(collected_at=fecal_1mo_site)),
                  inf.6mo.taxa %>% select(indID, s__Bifidobacterium_infantis) %>% 
                    mutate(mytax = scale(s__Bifidobacterium_infantis), tmpt="6 months") %>%
                    left_join(metadata %>% select(indID, fecal_6mo_site) %>% rename(collected_at=fecal_6mo_site)))
# add metadata
b.inf = b.inf %>% left_join(metadata %>% 
                              select(indID, delivery_cat, parity, matage, inf.white, mat_bmi, avg_totalhei, gdm, gbs,
                                     shotgun.cmv.status, ebf6mo, compfoods6mo)) %>% na.omit()

# regress out covars (excl. CMV status) 
summary(lmer(mytax ~ shotgun.cmv.status + tmpt + delivery_cat + parity + inf.white + mat_bmi + 
               ebf6mo + collected_at + gbs + compfoods6mo + gdm + (1|indID), data = b.inf))

b.inf$tax.resid = residuals(lmer(mytax ~ tmpt + delivery_cat + parity + inf.white + mat_bmi + 
                                   ebf6mo + collected_at + gbs + compfoods6mo + gdm + (1|indID), data = b.inf))
# make plot
ggplot(b.inf, aes(x=shotgun.cmv.status, y=tax.resid,color=shotgun.cmv.status)) + 
  geom_boxplot(fill=NA,outlier.shape=NA) +
  geom_beeswarm(cex=3,size=0.5,alpha=0.5) +
  scale_color_manual(values=c("orange","slateblue")) +
  xlab("") + ylab("Bifidobacterium infantis") +
  theme(legend.position="None", axis.title.y=element_text(face="italic")) + facet_wrap(~tmpt)
ggsave("../figs/fig4C.png",width=3.5,height=2.5)





