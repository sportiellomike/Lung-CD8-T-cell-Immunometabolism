#we included final output files in our github repository, but are including this for sake of reproducibility. You do not need to run this script.
#you'll have to correct all these different working directories for it to work on your machine.I renamed 
#'featurecounts.210319.goodjobadam' to 'counts for converttotpmR' that you downloaded from this repository.
#get files into variable
source('tpmfunction.R')
files<-list.files('./data/featurecounts.210319.goodjobadam/counts/',full.names = T)
fileshort<-list.files('./data/featurecounts.210319.goodjobadam/counts/')
files

#read files in
CD103_1<-read.table(files[[1]],header=T)
tablefile<-CD103_1
names(tablefile)[7] <- "Counts"
tablefile$tpm<-tpm(tablefile$Counts,tablefile$Length)
CD103_1<-tablefile

CD49a_1<-read.table(files[[2]],header=T)
tablefile<-CD49a_1
names(tablefile)[7] <- "Counts"
tablefile$tpm<-tpm(tablefile$Counts,tablefile$Length)
CD49a_1<-tablefile
DN_1<-read.table(files[[3]],header=T)
tablefile<-DN_1
names(tablefile)[7] <- "Counts"
tablefile$tpm<-tpm(tablefile$Counts,tablefile$Length)
DN_1<-tablefile

DP_1<-read.table(files[[4]],header=T)
tablefile<-DP_1
names(tablefile)[7] <- "Counts"
tablefile$tpm<-tpm(tablefile$Counts,tablefile$Length)
DP_1<-tablefile

CD103_2<-read.table(files[[5]],header=T)
tablefile<-CD103_2
names(tablefile)[7] <- "Counts"
tablefile$tpm<-tpm(tablefile$Counts,tablefile$Length)
CD103_2<-tablefile

CD49a_2<-read.table(files[[6]],header=T)
tablefile<-CD49a_2
names(tablefile)[7] <- "Counts"
tablefile$tpm<-tpm(tablefile$Counts,tablefile$Length)
CD49a_2<-tablefile

DN_2<-read.table(files[[7]],header=T)
tablefile<-DN_2
names(tablefile)[7] <- "Counts"
tablefile$tpm<-tpm(tablefile$Counts,tablefile$Length)
DN_2<-tablefile

DP_2<-read.table(files[[8]],header=T)
tablefile<-DP_2
names(tablefile)[7] <- "Counts"
tablefile$tpm<-tpm(tablefile$Counts,tablefile$Length)
DP_2<-tablefile

CD103_3<-read.table(files[[9]],header=T)
tablefile<-CD103_3
names(tablefile)[7] <- "Counts"
tablefile$tpm<-tpm(tablefile$Counts,tablefile$Length)
CD103_3<-tablefile

CD49a_3<-read.table(files[[10]],header=T)
tablefile<-CD49a_3
names(tablefile)[7] <- "Counts"
tablefile$tpm<-tpm(tablefile$Counts,tablefile$Length)
CD49a_3<-tablefile

DN_3<-read.table(files[[11]],header=T)
tablefile<-DN_3
names(tablefile)[7] <- "Counts"
tablefile$tpm<-tpm(tablefile$Counts,tablefile$Length)
DN_3<-tablefile

DP_3<-read.table(files[[12]],header=T)
tablefile<-DP_3
names(tablefile)[7] <- "Counts"
tablefile$tpm<-tpm(tablefile$Counts,tablefile$Length)
DP_3<-tablefile


#combine tpms into one df
CD103<-data.frame(CD103_1$tpm,CD103_2$tpm,CD103_3$tpm)
rownames(CD103)<-CD103_1$Geneid
CD49a<-data.frame(CD49a_1$tpm,CD49a_2$tpm,CD49a_3$tpm)
rownames(CD49a)<-CD49a_1$Geneid
DN<-data.frame(DN_1$tpm,DN_2$tpm,DN_3$tpm)
rownames(DN)<-DN_1$Geneid
DP<-data.frame(DP_1$tpm,DP_2$tpm,DP_3$tpm)
rownames(DP)<-DP_1$Geneid

#get average tpm
CD103$rowmean<-rowMeans(CD103)
CD49a$rowmean<-rowMeans(CD49a)
DP$rowmean<-rowMeans(DP)
DN$rowmean<-rowMeans(DN)

#encode this as binary and fix ensembl decimal
CD103$binary<-ifelse(CD103$rowmean > 1,1,0)
CD103$trans<-rownames(CD103)
split<-do.call(rbind, strsplit(as.character(CD103$trans),"\\."))
CD103$trans<-split[,1]
rownames(CD103)<-CD103$trans

CD49a$binary<-ifelse(CD49a$rowmean > 1,1,0)
CD49a$trans<-rownames(CD49a)
split<-do.call(rbind, strsplit(as.character(CD49a$trans),"\\."))
CD49a$trans<-split[,1]
rownames(CD49a)<-CD49a$trans

DP$binary<-ifelse(DP$rowmean > 1,1,0)
DP$trans<-rownames(DP)
split<-do.call(rbind, strsplit(as.character(DP$trans),"\\."))
DP$trans<-split[,1]
rownames(DP)<-DP$trans

DN$binary<-ifelse(DN$rowmean > 1,1,0)
DN$trans<-rownames(DN)
split<-do.call(rbind, strsplit(as.character(DN$trans),"\\."))
DN$trans<-split[,1]
rownames(DN)<-DN$trans

#map ensembl to genes
library(AnnotationDbi)
library(org.Mm.eg.db)
columns(org.Mm.eg.db)
keytypes(org.Mm.eg.db)

CD103$symbol<-mapIds(org.Mm.eg.db, keys = CD103$trans,column = c('SYMBOL'), keytype = 'ENSEMBL')
CD103$genename<-mapIds(org.Mm.eg.db, keys = CD103$trans,column = c('GENENAME'), keytype = 'ENSEMBL')
CD103$gene<-mapIds(org.Mm.eg.db, keys = CD103$trans,column = c('ENTREZID'), keytype = 'ENSEMBL')

CD49a$symbol<-mapIds(org.Mm.eg.db, keys = CD49a$trans,column = c('SYMBOL'), keytype = 'ENSEMBL')
CD49a$genename<-mapIds(org.Mm.eg.db, keys = CD49a$trans,column = c('GENENAME'), keytype = 'ENSEMBL')
CD49a$gene<-mapIds(org.Mm.eg.db, keys = CD49a$trans,column = c('ENTREZID'), keytype = 'ENSEMBL')

DN$symbol<-mapIds(org.Mm.eg.db, keys = DN$trans,column = c('SYMBOL'), keytype = 'ENSEMBL')
DN$genename<-mapIds(org.Mm.eg.db, keys = DN$trans,column = c('GENENAME'), keytype = 'ENSEMBL')
DN$gene<-mapIds(org.Mm.eg.db, keys = DN$trans,column = c('ENTREZID'), keytype = 'ENSEMBL')

DP$symbol<-mapIds(org.Mm.eg.db, keys = DP$trans,column = c('SYMBOL'), keytype = 'ENSEMBL')
DP$genename<-mapIds(org.Mm.eg.db, keys = DP$trans,column = c('GENENAME'), keytype = 'ENSEMBL')
DP$gene<-mapIds(org.Mm.eg.db, keys = DP$trans,column = c('ENTREZID'), keytype = 'ENSEMBL')

head(DP)
#you're going to have to change working directories to something that makes sense for you.
setwd("C:/Users/msportiello/Box/TRM_Metabolism_Modeling/metabinformatics/results/dp")
expressiondata<-data.frame(DP$gene,DP$symbol,DP$rowmean)
head(expressiondata)
expressiondata<-expressiondata[!is.na(expressiondata$DP.gene), ]
names(expressiondata)[1] <- "gene"
names(expressiondata)[2] <- "symbol"
names(expressiondata)[3] <- "tpm"
head(expressiondata)
dim(expressiondata)
write.table(expressiondata,'dp.genesymboltpm.txt',row.names = F,col.names =T,quote = F)
#these .txt files were then saved as excels for use in the matlab scripts.
sessionInfo()
citation()



###extra stuff if you want it you don't need. 
#gene<-expressiondata$gene
#gene<-t(gene)
#gene
#generow<-gene
#write.table(generow,'generow.txt')
#gene<-as.vector(gene)
#class(gene)
#head(gene)
#writeChar(gene,'genechar.txt')

#value<-expressiondata$value
#class(value)
#head(value)

#binary<-expressiondata$binary
#class(binary)
#head(binary)
#length(binary)
#write.table(gene,'genetable.txt',row.names=F,quote=F,col.names=F)
#write.table(value,'valuetable.txt',row.names=F,quote=F,col.names=F)
#write.table(binary,'binarytable.txt',row.names=F,quote=F,col.names=F)

########################
#you're going to have to change working directories to something that makes sense for you.
setwd("C:/Users/msportiello/Box/TRM_Metabolism_Modeling/metabinformatics/results/dn")
expressiondata<-data.frame(DN$gene,DN$symbol,DN$rowmean)
head(expressiondata)
expressiondata<-expressiondata[!is.na(expressiondata$DN.gene), ]
names(expressiondata)[1] <- "gene"
names(expressiondata)[2] <- "symbol"
names(expressiondata)[3] <- "tpm"
head(expressiondata)
dim(expressiondata)
write.table(expressiondata,'DN.genesymboltpm.txt',row.names = F,col.names =T,quote = F)
#these .txt files were then saved as excels for use in the matlab scripts.
sessionInfo()
citation()

###extra stuff if you want it you don't need. 
#gene<-expressiondata$gene
#gene<-t(gene)
#gene
#generow<-gene
#write.table(generow,'generowDN.txt')
#gene<-as.vector(gene)
#class(gene)
#head(gene)
#writeChar(gene,'genecharDN.txt')

#value<-expressiondata$value
#class(value)
#head(value)

#binary<-expressiondata$binary
#class(binary)
#head(binary)
#length(binary)
#write.table(gene,'genetableDN.txt',row.names=F,quote=F,col.names=F)
#write.table(value,'valuetableDN.txt',row.names=F,quote=F,col.names=F)
#write.table(binary,'binarytableDN.txt',row.names=F,quote=F,col.names=F)

