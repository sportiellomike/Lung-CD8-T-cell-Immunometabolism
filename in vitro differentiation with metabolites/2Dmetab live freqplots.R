# load libraries
library(ggplot2)
library(viridis)
library(readxl)
library(ggpubr)
library(rstatix)
library(dplyr)
dat<-read_excel('PercentAlive2D.xlsx')


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

dat$Tube <- factor(dat$Tube, levels=c('Tube 1',
                                                    'Tube 2',
                                                    'Tube 3',
                                                    'Tube 4',
                                                    'Tube 5',
                                                    'Tube 6',
                                                    'Tube 7',
                                                    'Tube 8',
                                                    'Tube 9',
                                                    'Tube 10',
                                                    'Tube 11',
                                                    'Tube 12',
                                                    'Tube 13',
                                                    'Tube 14',
                                                    'Tube 15',
                                                    'Tube 16'))

ggplot(dat,aes(x=Tube,y=`Percent alive`))+
  stat_summary(fun=mean, geom="col",aes(fill=Condition))+ 
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar",  width=0.2)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  ylab('Percent alive')+
  theme(axis.title.x = element_blank(),
        axis.text.x=element_blank())+
  theme(legend.position = 'bottom',panel.border=element_blank(),panel.background = element_blank(),panel.grid.major = element_line(color='light grey'))+
  scale_fill_viridis(discrete=T)
ggsave(filename = 'PercentAlive2Dconditions.png',dpi=600)

sessionInfo()
### FIN ###