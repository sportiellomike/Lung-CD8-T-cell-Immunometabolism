# manipulating tissuemodeldptable and tissuemodeldntable
library(viridis)
library(ggplot2)
library(reshape2)

dntable<-read.delim('tissuemodeldntable.txt', header = TRUE, sep = "\t", dec = ".")
dptable<-read.delim('tissuemodeldptable.txt', header = TRUE, sep = "\t", dec = ".")

colnames(dntable)<-paste0(colnames(dntable),'-dn')
colnames(dptable)<-paste0(colnames(dptable),'-dp')

names(dntable)[1]<-'rxns'
names(dptable)[1]<-'rxns'

merge<-merge(dptable,dntable,by = 'rxns',all=T)
faodn<-subset(merge,merge$`subsystems-dn`=='Fatty acid oxidation')
faodp<-subset(merge,merge$`subsystems-dp`=='Fatty acid oxidation')

diffdn<-subset(faodn,faodn$`fbasolutiontissuev-dp`!=faodn$`fbasolutiontissuev-dn`)
diffdp<-subset(faodp,faodp$`fbasolutiontissuev-dp`!=faodp$`fbasolutiontissuev-dn`)


diffsubsetfao<-subset(diffdn,abs(diffdn$`fbasolutiontissuev-dn`+diffdn$`fbasolutiontissuev-dp`)>1)
write.csv(diffsubsetfao,'rxnswithdiffvectors-FAO.csv')


#######################

glydn<-subset(merge,merge$`subsystems-dn`=='Glycolysis/gluconeogenesis')
glydp<-subset(merge,merge$`subsystems-dp`=='Glycolysis/gluconeogenesis')

diffdngly<-subset(glydn,glydn$`fbasolutiontissuev-dp`!=glydn$`fbasolutiontissuev-dn`)
diffdpgly<-subset(glydp,glydp$`fbasolutiontissuev-dp`!=glydp$`fbasolutiontissuev-dn`)


diffsubset<-subset(diffdngly,abs(diffdngly$`fbasolutiontissuev-dn`+diffdngly$`fbasolutiontissuev-dp`)>1)
write.csv(diffsubset,'rxnswithdiffvectors-gly.csv')

#####
faogg<-diffsubsetfao[c(1:4,6,11)]
colnames(faogg)<-c('rxns','rxnnames','subsystem','keggid','DP','DN')
head(faogg)
melt<-melt(faogg,variable.name = 'Subset',value.name = c('Flux'))
meltplotrxns<-c('Beta oxidation fatty acid','Lipase, extracellular','Medium-Chain Acyl Coenzyme A Dehydrogenase')
meltplot1<-subset(melt,melt$rxnnames=='Beta oxidation  fatty acid')
meltplot2<-subset(melt,melt$rxnnames=='Lipase, extracellular')
meltplot3<-subset(melt,melt$rxnnames=='Medium-Chain Acyl Coenzyme A Dehydrogenase')
meltplot<-rbind(meltplot1,meltplot2,meltplot3)


p<-ggplot(meltplot,aes(fill=Subset,x=Subset,y=Flux))+
  geom_col()+
  facet_wrap(~rxnnames,strip.position = 'bottom',scales = 'free')+
  theme(legend.position = 'none',panel.border=element_blank(),panel.background = element_blank(),panel.grid.major = element_line(color='light grey'))+
  scale_fill_viridis(discrete=T,begin=.7,end=0)+
  ylab('Net Flux')


ggsave(file.path("./", paste0("faofluxplot", Sys.Date(), ".png")), plot = p, height = 4, width = 10, units = "in",dpi = 600)


sessionInfo()
