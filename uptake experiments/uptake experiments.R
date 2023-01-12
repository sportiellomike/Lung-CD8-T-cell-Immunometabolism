#load in libraries
library(ggstatsplot)
library(viridis)
library(gridExtra)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(readxl)

dat<-read_xlsx('MFIdatfullnorm.xlsx')
# dat$MFI<-dat$MFI*dat$Multiplier
#assign levels
dat$Subset <- factor(dat$Subset, levels = c('DP','CD49a SP','CD103 SP','DN'))


# New facet label names for dose variable
experiment.labs <- c("Bodipy d14", "FA Uptake d14", "Mitotracker d14",'NBDG Uptake d14','TMRE d14','Mitotracker d60+','NBDG Uptake d60+','FA Uptake d60+','TMRE d60+','Bodipy d60+')
names(experiment.labs) <- c("lungbodipy", "lungfa", "lungmito",'lungNBDG','lungtmre','lungmito70','lungnbdg70','lungfa70','lungtmre70','lungbodipy70')



datlung<-subset(dat,organ=='lung')


names(experiment.labs) <- c("lungbodipy", "lungfa", "lungmito",'lungNBDG','lungtmre','lungmito70','lungnbdg70','lungfa70','lungtmre70','lungbodipy70')

datlung$organexperimentf = factor(datlung$organexperiment, levels=c('lungmito','lungtmre','lungNBDG','lungfa','lungbodipy',
                                                              'lungmito70','lungtmre70','lungnbdg70','lungfa70','lungbodipy70'))



cleandat<-datlung
cleandat$organexperimentsubset<-paste0(cleandat$organexperimentf,cleandat$Subset)
# stats
# first test outliers
outliers<-cleandat %>%
  group_by(organexperimentsubset) %>%
  identify_outliers(MFI)
outliersextreme<-subset(outliers,outliers$is.extreme==TRUE)


# seven outliers from five mice were identified as extreme, so we will remove those mice here
dim(cleandat)
outliersremoveddat<-subset(cleandat,cleandat$Mouse!='Mouse S' | cleandat$organexperimentf!='lungmito')
dim(outliersremoveddat)
outliersremoveddat<-subset(outliersremoveddat,outliersremoveddat$Mouse!='Mouse H' | outliersremoveddat$organexperimentf!='lungmito')
dim(outliersremoveddat)
outliersremoveddat<-subset(outliersremoveddat,outliersremoveddat$Mouse!='Mouse T' | outliersremoveddat$organexperimentf!='lungNBDG')
dim(outliersremoveddat)
outliersremoveddat<-subset(outliersremoveddat,outliersremoveddat$Mouse!='Mouse B' | outliersremoveddat$organexperimentf!='lungtmre70')
dim(outliersremoveddat)
outliersremoveddat<-subset(outliersremoveddat,outliersremoveddat$Mouse!='Mouse D' | outliersremoveddat$organexperimentf!='lungtmre70')
dim(outliersremoveddat)



# then look at qqq plots for normality
ggqqplot(outliersremoveddat, "MFI", facet.by = "organexperimentsubset") # they all look good so we can continue

# then do anovas
anovaresult<-outliersremoveddat %>%
  group_by(organexperimentf) %>%
  anova_test(dv=MFI,wid = Mouse,within =  Subset)
get_anova_table(anovaresult)
# lungmito, lungtmre70, and lungbodipy70 don't pass anova, so we must remove them when putting p values on plots later


stat.test <- outliersremoveddat %>%
  group_by(organexperimentf) %>%
  t_test(MFI ~ Subset,paired=T) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance("p.adj")

stat.test <- stat.test %>% add_xy_position(x = 'Subset')
# move the lungnbdg brackets down by .8 so less space is wasted on the plot
lungnbdg<-subset(stat.test,stat.test$organexperimentf=='lungNBDG')
lungnbdg$y.position<-lungnbdg$y.position-.8
stat.test<-subset(stat.test,stat.test$organexperimentf!='lungNBDG')
stat.test<-rbind(stat.test,lungnbdg)

stat.test<-subset(stat.test,stat.test$organexperimentf!= 'lungtmre70' & # removing the three that did not pass anova
                    stat.test$organexperimentf!= 'lungbodipy70' & 
                    stat.test$organexperimentf!= 'lungmito')

# set themes
theme_set(theme(axis.text = element_text(size=8),
                axis.title = element_text(size=8),
                strip.text = element_text(size=8),
                axis.text.x = element_text(angle=90),
                legend.position = 'none',
                strip.text.x = element_text(size = 8,margin = margin(.1, 0, .1, 0, "cm")),
                # strip.text.x = element_text(size = 8),
                legend.text = element_text(size=7),
                legend.title = element_text(size=7),
                # legend.position = 'bottom',panel.border=element_blank(),
                panel.background = element_blank(),
                panel.grid.major = element_line(color='light grey')
))


metabtimepoints<-ggplot(outliersremoveddat,aes(x = Subset, y = MFI)) +
  geom_bar(stat='summary',aes(fill=Subset)) +
  # geom_jitter(width = .2, alpha = .6) +
  stat_summary(geom="errorbar", fun.data="mean_sdl",fun.args=list(mult=1),
               color="black", width=.4)+
  scale_fill_viridis(discrete=T,begin=.8,end=0.1)+
  # theme(legend.position = 'none')+
  facet_wrap(~organexperimentf,nrow=2,labeller = as_labeller(experiment.labs))+
  stat_pvalue_manual(stat.test, hide.ns = T,label = "p.adj.signif",size = 6) +
  scale_y_continuous(name="Normalized MFI")

metabtimepoints
ggsave(file.path("./", paste0("metabtimepoints-lung-", Sys.Date(), ".png")), plot = metabtimepoints, width = 180,height = 200, units = "mm",dpi = 600)


datbal<-subset(dat,organ=='bal')


# New facet label names for dose variable
experiment.labs.bal <- c(
  "Mitotracker d14",
  'TMRE d14',
  'NBDG Uptake d14',
  "FA Uptake d14",
  "Bodipy d14")
names(experiment.labs.bal) <- c(
  "balmito",
  'baltmre',
  'balNBDG',
  "balfa",
  "balbodipy")

datbal$organexperimentf = factor(datbal$organexperiment, levels=c("balmito",
                                                               'baltmre',
                                                               'balNBDG',
                                                               "balfa",
                                                               "balbodipy"))


cleandat<-datbal
# first test outliers
outliers<-cleandat %>%
  group_by(organexperiment) %>%
  identify_outliers(MFI)
outliersextreme<-subset(outliers,outliers$is.extreme==TRUE)


# four outliers from three mice were identified as extreme, so we will remove tose mice here
dim(cleandat)
outliersremoveddat<-subset(cleandat,cleandat$Mouse!='Mouse E' | cleandat$organexperimentf!='baltmre')
dim(outliersremoveddat)
# then do anovas
anovaresult<-outliersremoveddat %>%
  group_by(organexperimentf) %>%
  anova_test(dv=MFI,wid = Mouse,within =  Subset)
get_anova_table(anovaresult)
# balmito doesn't pass anova


# t tests
stat.test.bal <- outliersremoveddat %>%
  group_by(organexperimentf) %>%
  t_test(MFI ~ Subset,paired=T) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance("p.adj")

stat.test.bal <- stat.test.bal %>% add_xy_position(x = 'Subset')
stat.test.bal<-subset(stat.test.bal,stat.test.bal$organexperimentf!= 'balmito' ) # removing the two that did not pass anova

# set themes
theme_set(theme(axis.text = element_text(size=8),
                axis.title = element_text(size=8),
                strip.text = element_text(size=8),
                axis.text.x = element_text(angle=90),
                legend.position = 'none',
                strip.text.x = element_text(size = 8,margin = margin(.1, 0, .1, 0, "cm")),
                legend.text = element_text(size=7),
                legend.title = element_text(size=7),
                # legend.position = 'bottom',panel.border=element_blank(),
                panel.background = element_blank(),
                panel.grid.major = element_line(color='light grey')
))


metabtimepointsbal<-ggplot(outliersremoveddat,aes(x = Subset, y = MFI)) +
  geom_bar(stat='summary',aes(fill=Subset)) +
  # geom_jitter(width = .2, alpha = .6) +
  stat_summary(geom="errorbar", fun.data="mean_sdl",fun.args=list(mult=1),
               color="black", width=.4)+
  scale_fill_viridis(discrete=T,begin=.8,end=0.1)+
  facet_wrap(~organexperimentf,nrow=1,labeller = as_labeller(experiment.labs.bal))+
  stat_pvalue_manual(stat.test.bal, hide.ns = T,label = "p.adj.signif",size = 6) +
  scale_y_continuous(name="Normalized MFI")

metabtimepointsbal
ggsave(file.path("./", paste0("metabtimepoints-bal-", Sys.Date(), ".png")), plot = metabtimepointsbal, width = 180,height = 100, units = "mm",dpi = 600)

sessionInfo()
### FIN ###