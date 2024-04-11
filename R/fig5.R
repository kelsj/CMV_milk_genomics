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


## load metadata
metadata = read_excel("../data/CMV_milk_genomics_suppTables_v4.xlsx", sheet=1, na = "NA")

## Fig. 5A: infant growth traits @ 0, 1, 6 months vs. milk CMV status
growth.res = read_excel("../data/CMV_milk_genomics_suppTables_v4.xlsx", sheet=13, na = "NA")

mypal = c("#A0D8B3","#05BFDB","#0A4D68")
growth.cmv.sens.plot = growth.res %>% separate(trait, into=c("trait","tmpt"), convert=T) %>%
  mutate(ci0 = est-1.96*se, ci1 = est+1.96*se) %>% filter(tmpt!="3") %>%
  mutate(dodge.x = ifelse(trait=="laz",tmpt-0.15,ifelse(trait=="wlz",tmpt+0.15,tmpt))) %>%
  mutate(trait=toupper(trait))
ggplot(growth.cmv.sens.plot) + geom_hline(yintercept=0,linetype="dashed",color="gray") +
  geom_point(aes(x=dodge.x,y=est,color=trait)) + 
  geom_segment(aes(x=dodge.x,xend=dodge.x,y=ci0,yend=ci1,color=trait)) +
  geom_line(aes(x=dodge.x,y=est,color=trait)) + #facet_wrap(~trait,ncol=3) +
  xlab("Infant age (months)") + scale_x_continuous(breaks=c(0,1,6)) +
  ylab("Estimated effect of\nCMV+ milk on infant trait") + labs(color="Infant\ntrait") +
  scale_color_manual(values=mypal)
ggsave("../figs/fig5A.png",width=4,height=3)


## Fig. 5B: WLZ 1mo. vs. proportion CMV-mapped reads
# get residual wlz 1mo after correcting for covars
growth.resid = metadata %>% select(indID,wlz_0,wlz_1,wlz_6,laz_0,laz_1,laz_6,waz_0,waz_1,waz_6,inf.white,
                                 mat_bmi,income_cat,delivery_cat,gdm,shotgun.cmv.prop,shotgun.cmv.status) %>% na.omit()

growth.resid$wlz1.resid = lm(wlz_1 ~ wlz_0 + inf.white + mat_bmi + income_cat + delivery_cat + gdm,
                           data=growth.resid)$residuals

ggplot(growth.resid %>% filter(shotgun.cmv.prop>0) %>% mutate(logAb=log(shotgun.cmv.prop)),
       aes(x=logAb,y=wlz1.resid)) + 
  geom_smooth(method="lm",color="gray",linetype="dashed",alpha=0.2) + 
  geom_point(color="slateblue",size=0.75) + 
  xlab("Proportion CMV-mapped reads (log)") +
  ylab("Infant 1 month\nweight-for-length Z-score")
ggsave("../figs/fig5B.png",width=3,height=3)


## Fig 5C: LAZ 1mo. vs. proportion CMV-mapped reads
growth.resid$laz1.resid = lm(laz_1 ~ laz_0 + inf.white + mat_bmi + income_cat + delivery_cat + gdm,
                             data=growth.resid)$residuals

ggplot(growth.resid %>% filter(shotgun.cmv.prop>0) %>% mutate(logAb=log(shotgun.cmv.prop)),
       aes(x=logAb,y=laz1.resid)) + 
  geom_smooth(method="lm",color="gray",linetype="dashed",alpha=0.2) + 
  geom_point(color="slateblue",size=0.75) + 
  xlab("Proportion CMV-mapped reads (log)") +
  ylab("Infant 1 month\nlength-for-age Z-score")
ggsave("../figs/fig5C.png",width=3,height=3)



## Fig 5D: WLZ 1mo. vs. milk kynurenine
# calculate residual WLZ1 against covars
covar.incl.kynu.resid = metadata %>% select(indID,shotgun.cmv.status,wlz_0,wlz_1,wlz_6,laz_0,laz_1,laz_6,waz_0,waz_1,
                                            waz_6,Center,parity,inf.white,mat_bmi,income_cat,delivery_cat,gdm,
                                            shotgun.cmv.prop, logKynurenine) %>% 
  mutate(logKynurenine = scale(logKynurenine)) %>% na.omit() 

covar.incl.kynu.resid$wlz1.resid = glm(wlz_1 ~ wlz_0 + inf.white + mat_bmi + income_cat + delivery_cat + Center + 
                                         parity + gdm, data=covar.incl.kynu.resid)$residuals
# plot
ggplot(covar.incl.kynu.resid, aes(x=logKynurenine,y=wlz1.resid,color=shotgun.cmv.status)) +
  geom_smooth(method="lm",alpha=0.2,se=T,size=0.6,aes(fill=shotgun.cmv.status)) + 
  geom_point(size=0.75) +
  scale_color_manual(values=c("orange","slateblue")) +
  scale_fill_manual(values=c("orange","slateblue")) +
  xlab("Milk kynurenine\n(normalized abundance)") + ylab("Infant 1 month\nweight-for-length Z-score")+ 
  theme(legend.position="None") + facet_wrap(~shotgun.cmv.status,nrow=2)
ggsave("../figs/fig5D.png",width=3,height=3)


## Fig. 5E: SEM with CMV status
# 5E, 5F SEM plots were made in powerpoint; but code shown below
# load lavaan, set seed
library(lavaan)
set.seed(1234)

# get data, N=200
med.dat = metadata %>% dplyr::select(indID,wlz_0,wlz_1,shotgun.cmv.status,logKynurenine) %>% 
  na.omit() %>% mutate(cmv.bin = ifelse(shotgun.cmv.status=="CMV+",1,0)) %>%
  mutate(logKynurenine=as.numeric(scale(logKynurenine)))

# reverse CMV/kynu in 4b
model4d <- ' # direct effect
             wlz_1 ~ f*logKynurenine + z*wlz_0 
           # mediator
             cmv.bin ~ d*logKynurenine
             wlz_1 ~ e*cmv.bin
           # indirect effect (d*e)
             de := d*e
          # total effect
           total := f + (d*e) '

model4d.fit <- sem(model4d, data = med.dat2, bootstrap=1000)
lavaanPlot(model4d.fit, coefs = TRUE)
summary(model4d.fit)$pe



## Fig. 5F: SEM with CMV load
# get data, N = 76
cmv.pos.med = metadata %>% select(indID,wlz_0,wlz_1,shotgun.cmv.prop,logKynurenine) %>% 
  filter(shotgun.cmv.prop>0) %>%
  mutate(log.ra = as.numeric(scale(log10(shotgun.cmv.prop)))) %>% 
  mutate(logKynurenine=as.numeric(scale(logKynurenine))) %>% na.omit()

# relab>kynu>wlz1
pos.mod1 <- ' # direct effect
             wlz_1 ~ z*wlz_0 + c*log.ra
           # mediator
             logKynurenine ~ a*log.ra
             wlz_1 ~ b*logKynurenine
           # indirect effect (a*b)
             ab := a*b
           # total effect
            total := c + (a*b) '

pos.mod1.fit <- sem(pos.mod1, data = cmv.pos.med, bootstrap=1000)
lavaanPlot(pos.mod1.fit, coefs = TRUE)
summary(pos.mod1.fit)$pe


