# plotting metabolic model barplot
library(viridis)
library(ggplot2)
library(ggpubr)
library(reshape2)

dntable<-read.delim('tissuemodeldntable.txt', header = TRUE, sep = "\t", dec = ".")
dptable<-read.delim('tissuemodeldptable.txt', header = TRUE, sep = "\t", dec = ".")
cd49asptable<-read.delim('tissuemodelcd49asptable.txt', header = TRUE, sep = "\t", dec = ".")
cd103sptable<-read.delim('tissuemodelcd103sptable.txt', header = TRUE, sep = "\t", dec = ".")

colnames(dntable)<-paste0(colnames(dntable),'-dn')
colnames(dptable)<-paste0(colnames(dptable),'-dp')
colnames(cd49asptable)<-paste0(colnames(cd49asptable),'-cd49asp')
colnames(cd103sptable)<-paste0(colnames(cd103sptable),'-cd103sp')

names(dntable)[1]<-'rxns'
names(dptable)[1]<-'rxns'
names(cd49asptable)[1]<-'rxns'
names(cd103sptable)[1]<-'rxns'

merge1<-merge(dptable,dntable,by = 'rxns',all=T)
merge2<-merge(cd49asptable,cd103sptable,by = 'rxns',all=T)
merge3<-merge(merge1,merge2,by='rxns',all=T)
merge<-merge3

# # fatty acid oxidation
mergefao<-subset(merge,merge$`subsystems-dp`=='Fatty acid oxidation' |
                   merge$`subsystems-dn`=='Fatty acid oxidation' |
                   merge$`subsystems-cd49asp`=='Fatty acid oxidation' |
                   merge$`subsystems-cd103sp`=='Fatty acid oxidation')

mergefaononzero<-subset(mergefao,(abs(mergefao$`fbasolutiontissuev-dp`)+
                                 abs(mergefao$`fbasolutiontissuev-dn`)+
                                 abs(mergefao$`fbasolutiontissuev-cd49asp`)+
                                 abs(mergefao$`fbasolutiontissuev-cd103sp`)>0))

meltfaononzero<-mergefaononzero[c(1:6,11,16,21)]
colnames(meltfaononzero)<-c('rxns','rxnnames','subsystem','keggid','reactomeid','DP','DN','CD49a SP','CD103 SP')

meltfaononzero<-melt(meltfaononzero,variable.name = 'Subset',value.name = c('Flux'))
meltfaononzero$Subset<-factor(meltfaononzero$Subset,c('DP','CD49a SP','CD103 SP','DN'))

# prepare plot for Beta oxidation
selectfaoreactions<-c(
  'Beta oxidation of med/long chain fatty acid'
)
selectfaoreactionsindex<-meltfaononzero$rxnnames %in% selectfaoreactions
selectfaoreactionsmelt<-meltfaononzero[selectfaoreactionsindex,]


selectfaoreactionsmeltsub<-subset(selectfaoreactionsmelt,selectfaoreactionsmelt$Subset!='CD103 SP')

level_order <- c('DP','CD49a SP','CD103 SP','DN')
selectfaoreactionsmeltplotsub<-ggplot(selectfaoreactionsmeltsub,aes(fill=Subset,x=factor(Subset),y=Flux))+
  geom_col()+
  # facet_wrap(~rxnnames,strip.position = 'bottom',scales = 'free',ncol = 2)+
  theme(legend.position = 'none',panel.border=element_blank(),panel.background = element_blank(),panel.grid.major = element_line(color='light grey'))+
  scale_fill_manual(values = c("#7AD151FF", "#1FA188FF" , "#482576FF"))+
  ylab(expression(Net~Flux~(mmol~g^-1~DW^-1~h^-1)))+
  xlab('Subset')+
  theme(strip.text = element_blank())+
  ggtitle('Net flux of Beta Oxidation by Subset')+
  theme(plot.title =  element_text(hjust = 0.5),
        axis.title.x = element_blank())
selectfaoreactionsmeltplotsub

# save th eplot
ggsave(paste0(Sys.Date(),'plosonefluxplot.png'),plot=selectfaoreactionsmeltplotsub,width=85,height=85,units = 'mm',dpi=600)

sessionInfo()
### FIN ###