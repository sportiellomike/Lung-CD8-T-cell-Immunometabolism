# load libraries
library(ggplot2)
library(viridis)
library(readxl)
library(ggpubr)
library(rstatix)
library(dplyr)
dat<-read_excel('dat-7.13.22 day 6 post tgfb-both-reps.xlsx',sheet = 1)
dat$Plate<-factor(dat$Plate)
dat$cd8sperul<-dat$count.live.singlet.CD8Tcells/dat$volume
# DPcd49acd103
dat$freq.CD49apos=dat$CD49asp+dat$DPcd49acd103
dat$freq.CCR7pos=dat$CCR7sp+dat$DPccr7cd62L
dat$freq.CD62Lpos=dat$CD62Lsp+dat$DPccr7cd62L
dat$freq.CD103pos=dat$CD103sp+dat$DPcd49acd103
dat$cd49acountperul<-dat$freq.CD49apos*dat$cd8sperul

dat$cd49anegcd103negccr7negcd62lnegperul<-dat$CD49anegCD103negCCR7negCD62lneg*dat$cd8sperul
dat$CD49aCD103DPcountperul<-dat$DPcd49acd103*dat$cd8sperul
dat$CCR7countperul<-dat$freq.CCR7pos*dat$cd8sperul
dat$CD62Lcountperul<-dat$freq.CD62Lpos*dat$cd8sperul
dat$CCR7CD62LDPcountperul<-dat$DPccr7cd62L*dat$cd8sperul
dat$CD103countperul<-dat$freq.CD103pos*dat$cd8sperul

dat$lipid.gluc<-paste0(dat$Lipid.level,'.',dat$Glucose.level)

dat1<-subset(dat,dat$Lipid.level!=8 & dat$Glucose.level!=100)
dat2<-subset(dat,dat$Lipid.level!=8 & dat$Lipid.level!=4 & dat$Glucose.level!=100)

dat.lipid0<-subset(dat,dat$Lipid.level==0)
dat.lipid2<-subset(dat,dat$Lipid.level==2)
dat.lipid4<-subset(dat,dat$Lipid.level==4)
dat.lipid8<-subset(dat,dat$Lipid.level==8)

dat.glucose0<-subset(dat,dat$Glucose.level==0)
dat.glucose25<-subset(dat,dat$Glucose.level==25)
dat.glucose50<-subset(dat,dat$Glucose.level==50)
dat.glucose100<-subset(dat,dat$Glucose.level==100)



theme_set(theme(axis.text = element_text(size=8),
                axis.title = element_text(size=8),
                strip.text = element_text(size=8),
                axis.text.x = element_text(angle=90),
                legend.position = 'none',
                # strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm"))
                strip.text.x = element_text(size = 8,margin = margin(.1, 0, .1, 0, "cm")),
                legend.text = element_text(size=7),
                legend.title = element_text(size=7)
))


datforposter<-dat
datforposter$Glucose.num<-datforposter$Glucose.level
datforposter$Lipid.num<-datforposter$Lipid.level

datforposter$Lipid.level<-as.character(datforposter$Lipid.level)
datforposter["Lipid.level"][datforposter["Lipid.level"] == 0] <- '0x Lipid'
datforposter["Lipid.level"][datforposter["Lipid.level"] == 2] <- '2x Lipid'
datforposter["Lipid.level"][datforposter["Lipid.level"] == 4] <- '4x Lipid'
datforposter["Lipid.level"][datforposter["Lipid.level"] == 8] <- '8x Lipid'

datforposter$Glucose.level<-as.character(datforposter$Glucose.level)
datforposter["Glucose.level"][datforposter["Glucose.level"] == 0] <- '4.5 mg/mL'
datforposter["Glucose.level"][datforposter["Glucose.level"] == 25] <- '9.5 mg/mL'
datforposter["Glucose.level"][datforposter["Glucose.level"] == 50] <- '14.5 mg/mL'
datforposter["Glucose.level"][datforposter["Glucose.level"] == 100] <- '24.5 mg/mL'

datforposter$Glucose.level <- factor(datforposter$Glucose.level, levels=c('4.5 mg/mL','9.5 mg/mL','14.5 mg/mL','24.5 mg/mL'))

datforposter["Glucose.num"][datforposter["Glucose.num"] == 0] <- 4.5
datforposter["Glucose.num"][datforposter["Glucose.num"] == 25] <- 9.5
datforposter["Glucose.num"][datforposter["Glucose.num"] == 50] <- 14.5
datforposter["Glucose.num"][datforposter["Glucose.num"] == 100] <- 24.5


# make ids
vecvec<-c()

for (q in 1:10) {
  addvec<-rep(q, 16)
  vecvec<-c(vecvec,addvec)
}
datforposter$id <- vecvec


# unfaceting
datforposter$condition<-paste0(datforposter$Lipid.level,'-',datforposter$Glucose.level)

datforposter$condition <- factor(datforposter$condition, levels=c('0x Lipid-4.5 mg/mL',
                                                                  '0x Lipid-9.5 mg/mL',
                                                                  '0x Lipid-14.5 mg/mL',
                                                                  '0x Lipid-24.5 mg/mL',
                                                                  '2x Lipid-4.5 mg/mL',
                                                                  '2x Lipid-9.5 mg/mL',
                                                                  '2x Lipid-14.5 mg/mL',
                                                                  '2x Lipid-24.5 mg/mL',
                                                                  '4x Lipid-4.5 mg/mL',
                                                                  '4x Lipid-9.5 mg/mL',
                                                                  '4x Lipid-14.5 mg/mL',
                                                                  '4x Lipid-24.5 mg/mL',
                                                                  '8x Lipid-4.5 mg/mL',
                                                                  '8x Lipid-9.5 mg/mL',
                                                                  '8x Lipid-14.5 mg/mL',
                                                                  '8x Lipid-24.5 mg/mL'))

# set themes
theme_set(theme(axis.text = element_text(size=8),
                axis.title = element_text(size=8),
                strip.text = element_text(size=8),
                axis.text.x = element_text(angle=90),
                legend.position = 'none',
                # strip.text.x = element_text(margin = margin(.1, 0, .1, 0, "cm"))
                strip.text.x = element_text(size = 8,margin = margin(.1, 0, .1, 0, "cm")),
                legend.text = element_text(size=6),
                legend.title = element_text(size=6)
))

#### CD49a positivity ####
### CD49a cells per uL ###
#stats
stat1 <- datforposter %>%
  group_by(Glucose.level) %>%
  t_test(cd49acountperul ~ condition,paired=T) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance("p.adj")
stat1<-stat1[-1]

stat2 <- datforposter %>%
  group_by(Lipid.level) %>%
  t_test(cd49acountperul ~ condition,paired=T) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance("p.adj")
stat2<-stat2[-1]



stat1 <- stat1 %>% add_xy_position(x = "condition")
stat2 <- stat2 %>% add_xy_position(x = "condition")
maxstat2<-max(stat2$y.position)
stat1$y.position<-stat1$y.position+maxstat2

# plot
plot1<-ggplot(datforposter,aes(x=condition,y=cd49acountperul))+
  stat_summary(fun=mean, geom="col",aes(fill=Glucose.level))+ 
  geom_jitter(aes(color=Lipid.level),alpha=.8,size=.75)+
  scale_color_viridis("Lipid Level",discrete=T,option = 'mako',begin =.3,end=.8)+
  scale_fill_viridis("Glucose level",discrete = T,option='rocket',begin =.3,end=.8)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar",  width=0.2)+
  #stat_pvalue_manual(stat1, hide.ns = T,label = "p.adj.signif",step.increase = .1,vjust = .9,bracket.nudge.y = -100000) +
  #stat_pvalue_manual(stat2, hide.ns = T,label = "p.adj.signif",vjust = .9) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  ylab('CD49a positive cells per uL')+
  theme(axis.title.x = element_blank())+
  theme(legend.position = 'bottom',panel.border=element_blank(),panel.background = element_blank(),panel.grid.major = element_line(color='light grey'))

plot1
### Percent CD49a cells ###
#stats
stat3 <- datforposter %>%
  group_by(Glucose.level) %>%
  t_test(freq.CD49apos ~ condition,paired=T) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance("p.adj")
stat3<-stat3[-1]

stat4 <- datforposter %>%
  group_by(Lipid.level) %>%
  t_test(freq.CD49apos ~ condition,paired=T) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance("p.adj")
stat4<-stat4[-1]

stat3 <- stat3 %>% add_xy_position(x = "condition")
stat4 <- stat4 %>% add_xy_position(x = "condition")
maxstat4<-max(stat4$y.position)
stat3$y.position<-stat3$y.position+maxstat4

# plot
plot2<-ggplot(datforposter,aes(x=factor(condition),y=freq.CD49apos))+
  stat_summary(fun=mean, geom="col",aes(fill=Glucose.level))+ 
  geom_jitter(aes(color=Lipid.level),alpha=.8,size=.75)+
  scale_color_viridis("Lipid Level",discrete=T,option = 'mako',begin =.3,end=.8)+
  scale_fill_viridis("Glucose level",discrete = T,option='rocket',begin =.3,end=.8)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar",  width=0.2)+
  # facet_wrap(~Lipid.level,nrow=1)+
  #stat_pvalue_manual(stat3, hide.ns = T,label = "p.adj.signif",step.increase = .05,vjust = .9,bracket.nudge.y = -100) +
  #stat_pvalue_manual(stat4, hide.ns = T,label = "p.adj.signif",vjust = .9) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  ylab('Percent CD49a positive')+
  # xlab('Glucose concentration')+
  theme(axis.title.x = element_blank())+
  theme(legend.position = 'bottom',panel.border=element_blank(),panel.background = element_blank(),panel.grid.major = element_line(color='light grey'))

plotggarrange_CD49a<-ggarrange(plot1,plot2,common.legend = T,ncol=2,legend = 'bottom')
ggsave(paste0(Sys.Date(),'-','CD49aplots-nopvalues.png'),plot=plotggarrange_CD49a,units = 'in',width = 7,height=7,dpi=600)

############################################
#### CCR7 positivity ####
### CCR7 cells per uL ###
#stats
stat1 <- datforposter %>%
  group_by(Glucose.level) %>%
  t_test(CCR7countperul ~ condition,paired=T) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance("p.adj")
stat1<-stat1[-1]

stat2 <- datforposter %>%
  group_by(Lipid.level) %>%
  t_test(CCR7countperul ~ condition,paired=T) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance("p.adj")
stat2<-stat2[-1]

stat1 <- stat1 %>% add_xy_position(x = "condition")
stat2 <- stat2 %>% add_xy_position(x = "condition")
maxstat2<-max(stat2$y.position)
stat1$y.position<-stat1$y.position+maxstat2

# plot
plot1<-ggplot(datforposter,aes(x=factor(condition),y=CCR7countperul))+
  stat_summary(fun=mean, geom="col",aes(fill=Glucose.level))+ 
  geom_jitter(aes(color=Lipid.level),alpha=.8,size=.75)+
  scale_color_viridis("Lipid Level",discrete=T,option = 'mako',begin =.3,end=.8)+
  scale_fill_viridis("Glucose level",discrete = T,option='rocket',begin =.3,end=.8)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar",  width=0.2)+
  # facet_wrap(~Lipid.level,nrow=1)+
  #stat_pvalue_manual(stat1, hide.ns = T,label = "p.adj.signif",step.increase = .1,vjust = .9,bracket.nudge.y = -150000) +
  #stat_pvalue_manual(stat2, hide.ns = T,label = "p.adj.signif",vjust = .9) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  ylab('CCR7 positive cells per uL')+
  # xlab('Glucose concentration')+
  theme(axis.title.x = element_blank())+
  theme(legend.position = 'bottom',panel.border=element_blank(),panel.background = element_blank(),panel.grid.major = element_line(color='light grey'))

plot1
### Percent CCR7 cells ###
#stats
stat3 <- datforposter %>%
  group_by(Glucose.level) %>%
  t_test(freq.CCR7pos ~ condition,paired=T) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance("p.adj")
stat3<-stat3[-1]

stat4 <- datforposter %>%
  group_by(Lipid.level) %>%
  t_test(freq.CCR7pos ~ condition,paired=T) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance("p.adj")
stat4<-stat4[-1]

stat3 <- stat3 %>% add_xy_position(x = "condition")
stat4 <- stat4 %>% add_xy_position(x = "condition")
maxstat4<-max(stat4$y.position)
stat3$y.position<-stat3$y.position+maxstat4
# plot
plot2<-ggplot(datforposter,aes(x=factor(condition),y=freq.CCR7pos))+
  stat_summary(fun=mean, geom="col",aes(fill=Glucose.level))+ 
  geom_jitter(aes(color=Lipid.level),alpha=.8,size=.75)+
  scale_color_viridis("Lipid Level",discrete=T,option = 'mako',begin =.3,end=.8)+
  scale_fill_viridis("Glucose level",discrete = T,option='rocket',begin =.3,end=.8)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar",  width=0.2)+
  # facet_wrap(~Lipid.level,nrow=1)+
  #stat_pvalue_manual(stat3, hide.ns = T,label = "p.adj.signif",step.increase = .05,vjust = .9,bracket.nudge.y = -100) +
  #stat_pvalue_manual(stat4, hide.ns = T,label = "p.adj.signif",vjust = .9) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  ylab('Percent CCR7 positive')+
  # xlab('Glucose concentration')+
  theme(axis.title.x = element_blank())+
  theme(legend.position = 'bottom',panel.border=element_blank(),panel.background = element_blank(),panel.grid.major = element_line(color='light grey'))
plot2
plotggarrange_CCR7<-ggarrange(plot1,plot2,common.legend = T,ncol=2,legend = 'bottom')
ggsave(paste0(Sys.Date(),'-','CCR7plots-nopvalues.png'),plot=plotggarrange_CCR7,units = 'in',width = 7,height=7,dpi=600)

############################################
#### DP positivity ####
### DP cells per uL ###
#stats
stat1 <- datforposter %>%
  group_by(Glucose.level) %>%
  t_test(CD49aCD103DPcountperul ~ condition,paired=T) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance("p.adj")
stat1<-stat1[-1]

stat2 <- datforposter %>%
  group_by(Lipid.level) %>%
  t_test(CD49aCD103DPcountperul ~ condition,paired=T) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance("p.adj")
stat2<-stat2[-1]

stat1 <- stat1 %>% add_xy_position(x = "condition")
stat2 <- stat2 %>% add_xy_position(x = "condition")
maxstat2<-max(stat2$y.position)
stat1$y.position<-stat1$y.position+maxstat2
# plot
plot1<-ggplot(datforposter,aes(x=factor(condition),y=CD49aCD103DPcountperul))+
  stat_summary(fun=mean, geom="col",aes(fill=Glucose.level))+ 
  geom_jitter(aes(color=Lipid.level),alpha=.8,size=.75)+
  scale_color_viridis("Lipid Level",discrete=T,option = 'mako',begin =.3,end=.8)+
  scale_fill_viridis("Glucose level",discrete = T,option='rocket',begin =.3,end=.8)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar",  width=0.2)+
  # facet_wrap(~Lipid.level,nrow=1)+
  #stat_pvalue_manual(stat1, hide.ns = T,label = "p.adj.signif",step.increase = .15,vjust = .9,bracket.nudge.y = -90000) +
  #stat_pvalue_manual(stat2, hide.ns = T,label = "p.adj.signif",vjust = .9) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  ylab('CD49aCD103 double positive cells per uL')+
  # xlab('Glucose concentration')+
  theme(axis.title.x = element_blank())+
  theme(legend.position = 'bottom',panel.border=element_blank(),panel.background = element_blank(),panel.grid.major = element_line(color='light grey'))
plot1

### Percent DP cells ###
#stats
stat3 <- datforposter %>%
  group_by(Glucose.level) %>%
  t_test(DPcd49acd103 ~ condition,paired=T) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance("p.adj")
stat3<-stat3[-1]

stat4 <- datforposter %>%
  group_by(Lipid.level) %>%
  t_test(DPcd49acd103 ~ condition,paired=T) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance("p.adj")
stat4<-stat4[-1]

stat3 <- stat3 %>% add_xy_position(x = "condition")
stat4 <- stat4 %>% add_xy_position(x = "condition")
maxstat4<-max(stat4$y.position)
stat3$y.position<-stat3$y.position+maxstat4
# plot
plot2<-ggplot(datforposter,aes(x=factor(condition),y=DPcd49acd103))+
  stat_summary(fun=mean, geom="col",aes(fill=Glucose.level))+ 
  geom_jitter(aes(color=Lipid.level),alpha=.8,size=.75)+
  scale_color_viridis("Lipid Level",discrete=T,option = 'mako',begin =.3,end=.8)+
  scale_fill_viridis("Glucose level",discrete = T,option='rocket',begin =.3,end=.8)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar",  width=0.2)+
  # facet_wrap(~Lipid.level,nrow=1)+
  #stat_pvalue_manual(stat3, hide.ns = T,label = "p.adj.signif",vjust = .9,bracket.nudge.y = -120) +
  #stat_pvalue_manual(stat4, hide.ns = T,label = "p.adj.signif",vjust = .9) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  ylab('Percent CD49aCD103 double positive')+
  # xlab('Glucose concentration')+
  theme(axis.title.x = element_blank())+
  theme(legend.position = 'bottom',panel.border=element_blank(),panel.background = element_blank(),panel.grid.major = element_line(color='light grey'))
plot2
plotggarrange_DP<-ggarrange(plot1,plot2,common.legend = T,ncol=2,legend = 'bottom')
ggsave(paste0(Sys.Date(),'-','DPplots-nopvalues.png'),plot=plotggarrange_DP,units = 'in',width = 7,height=7,dpi=600)

sessionInfo()
### FIN ###