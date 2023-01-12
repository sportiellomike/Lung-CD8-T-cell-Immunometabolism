library(ggplot2)
library(viridis)
library(gplots)
library(DESeq2)
library(dplyr)
library(tibble)
library(ashr)
library(IHW)
#Produce plots for figure 1
# Import data from featureCounts
countdata <- read.table("rawcounts.5.6.20.tsv", header=TRUE, row.names=1)

# Convert to matrix
countdata <- as.matrix(countdata)
head(countdata)

# Assign conditions
(condition <- factor(c(rep("CD103", 3), rep("DP", 3), rep("DN",3), rep("CD49a",3))))


# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design= ~ condition)
dds



# Run the DESeq pipeline
dds <- DESeq(dds)


rld <- rlogTransformation(dds)
head(assay(rld))


res.CD103vsCD49a <- lfcShrink(dds,contrast =c('condition','CD103','CD49a'),type='ashr')
res.CD103vsDN <- lfcShrink(dds,contrast =c('condition','CD103','DN'),type='ashr')
res.CD103vsDP <- lfcShrink(dds,contrast =c('condition','CD103','DP'),type='ashr')
res.CD49avsDP <- lfcShrink(dds,contrast =c('condition','CD49a','DP'),type='ashr')
res.CD49avsDN <- lfcShrink(dds,contrast =c('condition','CD49a','DN'),type='ashr')
res.DPvsDN <- lfcShrink(dds,contrast =c('condition','DP','DN'),type='ashr')

# Get differential expression results
res.CD103vsCD49a <- res.CD103vsCD49a[res.CD103vsCD49a$baseMean>10,]
res.CD103vsDN <- res.CD103vsDN[res.CD103vsDN$baseMean>10,]
res.CD103vsDP <- res.CD103vsDP[res.CD103vsDP$baseMean>10,]
res.CD49avsDP <- res.CD49avsDP[res.CD49avsDP$baseMean>10,]
res.CD49avsDN <- res.CD49avsDN[res.CD49avsDN$baseMean>10,]
res.DPvsDN <- res.DPvsDN[res.DPvsDN$baseMean>10,]

res.CD103vsCD49a <- na.omit(res.CD103vsCD49a)
res.CD103vsDN  <- na.omit(res.CD103vsDN)
res.CD103vsDP <- na.omit(res.CD103vsDP)
res.CD49avsDP <- na.omit(res.CD49avsDP)
res.CD49avsDN  <- na.omit(res.CD49avsDN)
res.DPvsDN <- na.omit(res.DPvsDN)

res.CD103vsCD49a[!complete.cases(res.CD103vsCD49a),]
res.CD103vsDN[!complete.cases(res.CD103vsDN),]
res.CD103vsDP[!complete.cases(res.CD103vsDP),]
res.CD49avsDP[!complete.cases(res.CD49avsDP),]
res.CD49avsDN[!complete.cases(res.CD49avsDN),]
res.DPvsDN[!complete.cases(res.DPvsDN),]

#How many results have a p value of less than 0.05?
table(res.CD103vsCD49a$padj<0.05)
table(res.CD103vsDN$padj<0.05)
table(res.CD103vsDP$padj<0.05)
table(res.CD49avsDP$padj<0.05)
table(res.CD49avsDN$padj<0.05)
table(res.DPvsDN$padj<0.05)

## Order by adjusted p-value
res.CD103vsCD49a <- res.CD103vsCD49a[order(res.CD103vsCD49a$padj), ]
res.CD103vsDN <- res.CD103vsDN[order(res.CD103vsDN$padj), ]
res.CD103vsDP <- res.CD103vsDP[order(res.CD103vsDP$padj), ]
res.CD49avsDP <- res.CD49avsDP[order(res.CD49avsDP$padj), ]
res.CD49avsDN <- res.CD49avsDN[order(res.CD49avsDN$padj), ]
res.DPvsDN <- res.DPvsDN[order(res.DPvsDN$padj), ]

#look at the first few lines
head(res.CD103vsCD49a)
head(res.CD103vsDN)
head(res.CD103vsDP)
head(res.CD49avsDP)
head(res.CD49avsDN)
head(res.DPvsDN)

#IHW
#run IHW package (independent hypothesis testing) from this paper: 
#http://bioconductor.org/packages/release/bioc/vignettes/IHW/inst/doc/introduction_to_ihw.html#an-example-rna-seq-differential-expression
#https://www.nature.com/articles/nmeth.3885

res.CD103vsCD49aDF <- as.data.frame(res.CD103vsCD49a)
res.CD103vsDPDF <- as.data.frame(res.CD103vsDP)
res.CD103vsDNDF <- as.data.frame(res.CD103vsDN)
res.CD49avsDPDF <- as.data.frame(res.CD49avsDP)
res.CD49avsDNDF <- as.data.frame(res.CD49avsDN)
res.DPvsDNDF <- as.data.frame(res.DPvsDN)

#In particular, we have p-values and baseMean (i.e., the mean of normalized counts) for each gene. As argued in the DESeq2 paper, these two statistics are approximately independent under the null hypothesis. Thus we have all the ingredients necessary for a IHW analysis (p-values and covariates), which we will apply at a significance level 0.05.
ihw.res.CD103vsCD49a <- ihw(pvalue ~ baseMean,  data = res.CD103vsCD49aDF, alpha = 0.05)
ihw.res.CD103vsDP <- ihw(pvalue ~ baseMean,  data = res.CD103vsDPDF, alpha = 0.05)
ihw.res.CD103vsDN <- ihw(pvalue ~ baseMean,  data = res.CD103vsDNDF, alpha = 0.05)
ihw.res.CD49avsDP <- ihw(pvalue ~ baseMean,  data = res.CD49avsDPDF, alpha = 0.05)
ihw.res.CD49avsDN <- ihw(pvalue ~ baseMean,  data = res.CD49avsDNDF, alpha = 0.05)
ihw.res.DPvsDN <- ihw(pvalue ~ baseMean,  data = res.DPvsDNDF, alpha = 0.05)
#This returns an object of the class ihwResult.

#add the gene names to the ihw df so I can get the gene names to give to enrichr
#CD103vsCD49a
de.ihw.res.CD103vsCD49a <- ihw.res.CD103vsCD49a@df
dim(res.CD103vsCD49a)
dim(de.ihw.res.CD103vsCD49a)
genes.res.CD103vsCD49a <- rownames(res.CD103vsCD49a)
genes.res.CD103vsCD49a <- as.vector(genes.res.CD103vsCD49a)
de.ihw.res.CD103vsCD49a["gene"] <- genes.res.CD103vsCD49a
ihw.res.CD103vsCD49a.all <- de.ihw.res.CD103vsCD49a
de.ihw.res.CD103vsCD49a <- de.ihw.res.CD103vsCD49a[de.ihw.res.CD103vsCD49a[,"weighted_pvalue"]<0.05,]
dim(de.ihw.res.CD103vsCD49a)
ihw.genes.res.CD103vsCD49a <- de.ihw.res.CD103vsCD49a$gene
head(ihw.genes.res.CD103vsCD49a)
ihw.res.CD103vsCD49a.all <- data.frame(ihw.res.CD103vsCD49a.all, row.names = 8)
ihw.res.CD103vsCD49a.all['log2FoldChange'] <- as.vector(res.CD103vsCD49a$log2FoldChange)
#make set of genes that are up and down regulated
ihw.genes.res.CD103vsCD49a.up <- ihw.res.CD103vsCD49a.all[ihw.res.CD103vsCD49a.all[,"log2FoldChange"]>0.5,]
ihw.genes.res.CD103vsCD49a.up <- ihw.genes.res.CD103vsCD49a.up[ihw.genes.res.CD103vsCD49a.up[,"weighted_pvalue"]<0.05,]
ihw.genes.res.CD103vsCD49a.down <- ihw.res.CD103vsCD49a.all[ihw.res.CD103vsCD49a.all[,"log2FoldChange"]<(-0.5),]
ihw.genes.res.CD103vsCD49a.down <- ihw.genes.res.CD103vsCD49a.down[ihw.genes.res.CD103vsCD49a.down[,"weighted_pvalue"]<0.05,]
upgenes.ihw.res.CD103vsCD49a <- rownames(ihw.genes.res.CD103vsCD49a.up)
downgenes.ihw.res.CD103vsCD49a <- rownames(ihw.genes.res.CD103vsCD49a.down)

#CD103vsDP
de.ihw.res.CD103vsDP <- ihw.res.CD103vsDP@df
dim(res.CD103vsDP)
dim(de.ihw.res.CD103vsDP)
genes.res.CD103vsDP <- rownames(res.CD103vsDP)
genes.res.CD103vsDP <- as.vector(genes.res.CD103vsDP)
de.ihw.res.CD103vsDP["gene"] <- genes.res.CD103vsDP
ihw.res.CD103vsDP.all <- de.ihw.res.CD103vsDP
de.ihw.res.CD103vsDP <- de.ihw.res.CD103vsDP[de.ihw.res.CD103vsDP[,"weighted_pvalue"]<0.05,]
dim(de.ihw.res.CD103vsDP)
ihw.genes.res.CD103vsDP <- de.ihw.res.CD103vsDP$gene
head(ihw.genes.res.CD103vsDP)
ihw.res.CD103vsDP.all <- data.frame(ihw.res.CD103vsDP.all, row.names = 8)
ihw.res.CD103vsDP.all['log2FoldChange'] <- as.vector(res.CD103vsDP$log2FoldChange)
#make set of genes that are up and down regulated
ihw.genes.res.CD103vsDP.up <- ihw.res.CD103vsDP.all[ihw.res.CD103vsDP.all[,"log2FoldChange"]>0.5,]
ihw.genes.res.CD103vsDP.up <- ihw.genes.res.CD103vsDP.up[ihw.genes.res.CD103vsDP.up[,"weighted_pvalue"]<0.05,]
ihw.genes.res.CD103vsDP.down <- ihw.res.CD103vsDP.all[ihw.res.CD103vsDP.all[,"log2FoldChange"]<(-0.5),]
ihw.genes.res.CD103vsDP.down <- ihw.genes.res.CD103vsDP.down[ihw.genes.res.CD103vsDP.down[,"weighted_pvalue"]<0.05,]
upgenes.ihw.res.CD103vsDP <- rownames(ihw.genes.res.CD103vsDP.up)
downgenes.ihw.res.CD103vsDP <- rownames(ihw.genes.res.CD103vsDP.down)

#CD103vsDN
de.ihw.res.CD103vsDN <- ihw.res.CD103vsDN@df
dim(res.CD103vsDN)
dim(de.ihw.res.CD103vsDN)
genes.res.CD103vsDN <- rownames(res.CD103vsDN)
genes.res.CD103vsDN <- as.vector(genes.res.CD103vsDN)
de.ihw.res.CD103vsDN["gene"] <- genes.res.CD103vsDN
ihw.res.CD103vsDN.all <- de.ihw.res.CD103vsDN
de.ihw.res.CD103vsDN <- de.ihw.res.CD103vsDN[de.ihw.res.CD103vsDN[,"weighted_pvalue"]<0.05,]
dim(de.ihw.res.CD103vsDN)
ihw.genes.res.CD103vsDN <- de.ihw.res.CD103vsDN$gene
head(ihw.genes.res.CD103vsDN)
ihw.res.CD103vsDN.all <- data.frame(ihw.res.CD103vsDN.all, row.names = 8)
ihw.res.CD103vsDN.all['log2FoldChange'] <- as.vector(res.CD103vsDN$log2FoldChange)
#make set of genes that are up and down regulated
ihw.genes.res.CD103vsDN.up <- ihw.res.CD103vsDN.all[ihw.res.CD103vsDN.all[,"log2FoldChange"]>0.5,]
ihw.genes.res.CD103vsDN.up <- ihw.genes.res.CD103vsDN.up[ihw.genes.res.CD103vsDN.up[,"weighted_pvalue"]<0.05,]
ihw.genes.res.CD103vsDN.down <- ihw.res.CD103vsDN.all[ihw.res.CD103vsDN.all[,"log2FoldChange"]<(-0.5),]
ihw.genes.res.CD103vsDN.down <- ihw.genes.res.CD103vsDN.down[ihw.genes.res.CD103vsDN.down[,"weighted_pvalue"]<0.05,]
upgenes.ihw.res.CD103vsDN <- rownames(ihw.genes.res.CD103vsDN.up)
downgenes.ihw.res.CD103vsDN <- rownames(ihw.genes.res.CD103vsDN.down)

#CD49avsDP
de.ihw.res.CD49avsDP <- ihw.res.CD49avsDP@df
dim(res.CD49avsDP)
dim(de.ihw.res.CD49avsDP)
genes.res.CD49avsDP <- rownames(res.CD49avsDP)
genes.res.CD49avsDP <- as.vector(genes.res.CD49avsDP)
de.ihw.res.CD49avsDP["gene"] <- genes.res.CD49avsDP
ihw.res.CD49avsDP.all <- de.ihw.res.CD49avsDP
de.ihw.res.CD49avsDP <- de.ihw.res.CD49avsDP[de.ihw.res.CD49avsDP[,"weighted_pvalue"]<0.05,]
dim(de.ihw.res.CD49avsDP)
ihw.genes.res.CD49avsDP <- de.ihw.res.CD49avsDP$gene
head(ihw.genes.res.CD49avsDP)
ihw.res.CD49avsDP.all <- data.frame(ihw.res.CD49avsDP.all, row.names = 8)
ihw.res.CD49avsDP.all['log2FoldChange'] <- as.vector(res.CD49avsDP$log2FoldChange)
#make set of genes that are up and down regulated
ihw.genes.res.CD49avsDP.up <- ihw.res.CD49avsDP.all[ihw.res.CD49avsDP.all[,"log2FoldChange"]>0.5,]
ihw.genes.res.CD49avsDP.up <- ihw.genes.res.CD49avsDP.up[ihw.genes.res.CD49avsDP.up[,"weighted_pvalue"]<0.05,]
ihw.genes.res.CD49avsDP.down <- ihw.res.CD49avsDP.all[ihw.res.CD49avsDP.all[,"log2FoldChange"]<(-0.5),]
ihw.genes.res.CD49avsDP.down <- ihw.genes.res.CD49avsDP.down[ihw.genes.res.CD49avsDP.down[,"weighted_pvalue"]<0.05,]
upgenes.ihw.res.CD49avsDP <- rownames(ihw.genes.res.CD49avsDP.up)
downgenes.ihw.res.CD49avsDP <- rownames(ihw.genes.res.CD49avsDP.down)

#CD49avsDN
de.ihw.res.CD49avsDN <- ihw.res.CD49avsDN@df
dim(res.CD49avsDN)
dim(de.ihw.res.CD49avsDN)
genes.res.CD49avsDN <- rownames(res.CD49avsDN)
genes.res.CD49avsDN <- as.vector(genes.res.CD49avsDN)
de.ihw.res.CD49avsDN["gene"] <- genes.res.CD49avsDN
ihw.res.CD49avsDN.all <- de.ihw.res.CD49avsDN
de.ihw.res.CD49avsDN <- de.ihw.res.CD49avsDN[de.ihw.res.CD49avsDN[,"weighted_pvalue"]<0.05,]
dim(de.ihw.res.CD49avsDN)
ihw.genes.res.CD49avsDN <- de.ihw.res.CD49avsDN$gene
head(ihw.genes.res.CD49avsDN)
ihw.res.CD49avsDN.all <- data.frame(ihw.res.CD49avsDN.all, row.names = 8)
ihw.res.CD49avsDN.all['log2FoldChange'] <- as.vector(res.CD49avsDN$log2FoldChange)
#make set of genes that are up and down regulated
ihw.genes.res.CD49avsDN.up <- ihw.res.CD49avsDN.all[ihw.res.CD49avsDN.all[,"log2FoldChange"]>0.5,]
ihw.genes.res.CD49avsDN.up <- ihw.genes.res.CD49avsDN.up[ihw.genes.res.CD49avsDN.up[,"weighted_pvalue"]<0.05,]
ihw.genes.res.CD49avsDN.down <- ihw.res.CD49avsDN.all[ihw.res.CD49avsDN.all[,"log2FoldChange"]<(-0.5),]
ihw.genes.res.CD49avsDN.down <- ihw.genes.res.CD49avsDN.down[ihw.genes.res.CD49avsDN.down[,"weighted_pvalue"]<0.05,]
upgenes.ihw.res.CD49avsDN <- rownames(ihw.genes.res.CD49avsDN.up)
downgenes.ihw.res.CD49avsDN <- rownames(ihw.genes.res.CD49avsDN.down)

#DPvsDN
de.ihw.res.DPvsDN <- ihw.res.DPvsDN@df
dim(res.DPvsDN)
dim(de.ihw.res.DPvsDN)
genes.res.DPvsDN <- rownames(res.DPvsDN)
genes.res.DPvsDN <- as.vector(genes.res.DPvsDN)
de.ihw.res.DPvsDN["gene"] <- genes.res.DPvsDN
ihw.res.DPvsDN.all <- de.ihw.res.DPvsDN
de.ihw.res.DPvsDN <- de.ihw.res.DPvsDN[de.ihw.res.DPvsDN[,"weighted_pvalue"]<0.05,]
dim(de.ihw.res.DPvsDN)
ihw.genes.res.DPvsDN <- de.ihw.res.DPvsDN$gene
head(ihw.genes.res.DPvsDN)
ihw.res.DPvsDN.all <- data.frame(ihw.res.DPvsDN.all, row.names = 8)
ihw.res.DPvsDN.all['log2FoldChange'] <- as.vector(res.DPvsDN$log2FoldChange)
#make set of genes that are up and down regulated
ihw.genes.res.DPvsDN.up <- ihw.res.DPvsDN.all[ihw.res.DPvsDN.all[,"log2FoldChange"]>0.5,]
ihw.genes.res.DPvsDN.up <- ihw.genes.res.DPvsDN.up[ihw.genes.res.DPvsDN.up[,"weighted_pvalue"]<0.05,]
ihw.genes.res.DPvsDN.down <- ihw.res.DPvsDN.all[ihw.res.DPvsDN.all[,"log2FoldChange"]<(-0.5),]
ihw.genes.res.DPvsDN.down <- ihw.genes.res.DPvsDN.down[ihw.genes.res.DPvsDN.down[,"weighted_pvalue"]<0.05,]
upgenes.ihw.res.DPvsDN <- rownames(ihw.genes.res.DPvsDN.up)
downgenes.ihw.res.DPvsDN <- rownames(ihw.genes.res.DPvsDN.down)

#GSEA
#Enrichr
library(enrichR)
listEnrichrDbs()
#https://www.biostars.org/p/343196/
library(ggplot2)
library(ggbreak) 
library(patchwork)
# theme_set(theme_grey())  

# set themes
theme_set(theme(axis.text = element_text(size=8),
                axis.title = element_text(size=8),
                strip.text = element_text(size=8),
                axis.text.x = element_text(angle=90),
                # legend.position = 'none',
                # strip.text.x = element_text(size = 8,margin = margin(.1, 0, .1, 0, "cm")),
                legend.text = element_text(size=8),
                legend.title = element_text(size=8),
                # legend.position = 'bottom',panel.border=element_blank(),
                panel.background = element_blank(),
                panel.grid.major = element_line(color='light grey')
))


dbs <- listEnrichrDbs()
alldbs <- listEnrichrDbs()
dbs <- c( "KEGG_2019_Mouse", 
          "Reactome_2016")
viridispal2 <- viridis(n=2)

#DP vs DN
list_up <- c(upgenes.ihw.res.DPvsDN)
list_down <- c(downgenes.ihw.res.DPvsDN)

eup <- enrichr(list_up, dbs)
edown <- enrichr(list_down, dbs)
up <- eup$KEGG_2019_Mouse
down <- edown$KEGG_2019_Mouse
#kegg
up$type <- "up"
down$type <- "down"
up <- up[up$Adjusted.P.value<.05,]
up <- up[order(up$Combined.Score), ]  # sort
down <- down[down$Adjusted.P.value<0.05,]
down <- down[order(down$Combined.Score), ]  # sort
down$Combined.Score <- (-1) * down$Combined.Score
gos <- rbind(down,up)
gos <- na.omit(gos) # Diverging Barcharts
dpdnkegg<-ggplot(gos, aes(x=reorder(Term,Combined.Score), y=Combined.Score , label=Combined.Score)) + 
  geom_bar(stat='identity', aes(fill=Adjusted.P.value), width=.5,position="dodge")  +
  scale_fill_viridis() + 
    labs(subtitle="KEGG") + 
  ylab('Combined Score')+ xlab('')+
 # scale_y_break(c(-1202, -21)) + scale_y_break(c(50, 250))+
  labs(fill = expression(P[adj]))+
  #theme(axis.text.x=element_text(angle=30))+
  coord_flip()
dpdnkegg

gosfilterdpdnkegg<-gos
gosfilterdpdnkegg$comparison<-'DP vs DN'
gosfilterdpdnkegg$database<-'KEGG'
gosdpdnkegg<-gosfilterdpdnkegg
gosfilterdpdnkegg<-subset(gosfilterdpdnkegg,gosfilterdpdnkegg$Term=='Steroid biosynthesis')



#Reactome_2016
up <- eup$Reactome_2016
down <- edown$Reactome_2016
up$type <- "up"
down$type <- "down"
up <- up[up$Adjusted.P.value<.05,]
up <- up[order(up$Combined.Score), ]  # sort
down <- down[down$Adjusted.P.value<0.05,]
down <- down[order(down$Combined.Score), ]  # sort
down$Combined.Score <- (-1) * down$Combined.Score
gos <- rbind(down,up)
gos <- na.omit(gos) # Diverging Barcharts
dpdnreactome<-ggplot(gos, aes(x=reorder(Term,Combined.Score), y=Combined.Score , label=Combined.Score)) + 
  geom_bar(stat='identity', aes(fill=Adjusted.P.value), width=.5,position="dodge")  +
  scale_fill_viridis() + 
  labs(subtitle="Reactome") + 
  ylab('Combined Score')+ xlab('')+
  labs(fill = expression(P[adj]))+
    coord_flip()
dpdnreactome

gosfilterdpdnreact<-gos
gosfilterdpdnreact$comparison<-'DP vs DN'
gosfilterdpdnreact$database<-'Reactome'
gosdpdnreact<-gosfilterdpdnreact
gosfilterdpdnreact<-subset(gosfilterdpdnreact,gosfilterdpdnreact$Term=='Purine metabolism Homo sapiens R-HSA-73847' |
                             gosfilterdpdnreact$Term=='Fatty acid, triacylglycerol, and ketone body metabolism Homo sapiens R-HSA-535734' |
                             gosfilterdpdnreact$Term=='Metabolism of lipids and lipoproteins Homo sapiens R-HSA-556833' |
                             gosfilterdpdnreact$Term== 'Regulation of cholesterol biosynthesis by SREBP (SREBF) Homo sapiens R-HSA-1655829' |
                             gosfilterdpdnreact$Term== 'Cholesterol biosynthesis Homo sapiens R-HSA-191273' |
                             gosfilterdpdnreact$Term== 'Metabolism Homo sapiens R-HSA-1430728' |
                             gosfilterdpdnreact$Term== 'Metabolism of proteins Homo sapiens R-HSA-392499' |
                             gosfilterdpdnreact$Term== 'Metabolism of amino acids and derivatives Homo sapiens R-HSA-71291' |
                             gosfilterdpdnreact$Term== 'Selenoamino acid metabolism Homo sapiens R-HSA-2408522' |
                             gosfilterdpdnreact$Term== 'Selenocysteine synthesis Homo sapiens R-HSA-2408557'
                           )




#there are no up or down regulated pathways in DP vs cd49a sp
list_up <- c(upgenes.ihw.res.CD49avsDP)
list_down <- c(downgenes.ihw.res.CD49avsDP)

eup <- enrichr(list_up, dbs)
edown <- enrichr(list_down, dbs)
up <- eup$KEGG_2019_Mouse
down <- edown$KEGG_2019_Mouse
#kegg
up$type <- "up"
down$type <- "down"
up <- up[up$Adjusted.P.value<.05,]
up <- up[order(up$Combined.Score), ]  # sort
down <- down[down$Adjusted.P.value<0.05,]
down <- down[order(down$Combined.Score), ]  # sort
down$Combined.Score <- (-1) * down$Combined.Score
gos <- rbind(down,up)
gos <- na.omit(gos) # Diverging Barcharts
cd49aDPkegg<-ggplot(gos, aes(x=reorder(Term,Combined.Score), y=Combined.Score , label=Combined.Score)) + 
  geom_bar(stat='identity', aes(fill=Adjusted.P.value), width=.5,position="dodge")  +
  scale_fill_viridis() + 
  labs(subtitle="KEGG") + 
  ylab('Combined Score')+ xlab('')+
  # scale_y_break(c(-1202, -21)) + scale_y_break(c(50, 250))+
  labs(fill = expression(P[adj]))+
  #theme(axis.text.x=element_text(angle=30))+
  coord_flip()
cd49aDPkegg


#Reactome_2016
up <- eup$Reactome_2016
down <- edown$Reactome_2016
up$type <- "up"
down$type <- "down"
up <- up[up$Adjusted.P.value<.05,]
up <- up[order(up$Combined.Score), ]  # sort
down <- down[down$Adjusted.P.value<0.05,]
down <- down[order(down$Combined.Score), ]  # sort
down$Combined.Score <- (-1) * down$Combined.Score
gos <- rbind(down,up)
gos <- na.omit(gos) # Diverging Barcharts
cd49aDPreactome<-ggplot(gos, aes(x=reorder(Term,Combined.Score), y=Combined.Score , label=Combined.Score)) + 
  geom_bar(stat='identity', aes(fill=Adjusted.P.value), width=.5,position="dodge")  +
  scale_fill_viridis() + 
  labs(subtitle="Reactome") + 
  ylab('Combined Score')+ xlab('')+
  labs(fill = expression(P[adj]))+
  coord_flip()
cd49aDPreactome

#CD49a vs dn
list_up <- c(upgenes.ihw.res.CD49avsDN)
list_down <- c(downgenes.ihw.res.CD49avsDN)

eup <- enrichr(list_up, dbs)
edown <- enrichr(list_down, dbs)
up <- eup$KEGG_2019_Mouse
down <- edown$KEGG_2019_Mouse
up$type <- "up"
down$type <- "down"
up <- up[up$Adjusted.P.value<.05,]
up <- up[order(up$Combined.Score), ]  # sort
down <- down[down$Adjusted.P.value<0.05,]
down <- down[order(down$Combined.Score), ]  # sort
down$Combined.Score <- (-1) * down$Combined.Score
gos <- rbind(down,up)
gos <- na.omit(gos) # Diverging Barcharts
cd49adnkegg<-ggplot(gos, aes(x=reorder(Term,Combined.Score), y=Combined.Score , label=Combined.Score)) + 
  geom_bar(stat='identity', aes(fill=Adjusted.P.value), width=.5,position="dodge")  +
  scale_fill_viridis() + 
  labs(subtitle="KEGG") + 
  ylab('Combined Score')+ xlab('')+
  labs(fill = expression(P[adj]))+
  coord_flip()
cd49adnkegg

gosfiltercd49adnkegg<-gos
gosfiltercd49adnkegg$comparison<-'CD49a vs DN'
gosfiltercd49adnkegg$database<-'KEGG'
goscd49adnkegg<-gosfiltercd49adnkegg
gosfiltercd49adnkegg<-subset(gosfiltercd49adnkegg,gosfiltercd49adnkegg$Term=='Sphingolipid signaling pathway' |
                                gosfiltercd49adnkegg$Term=='Glycosphingolipid biosynthesis' |
                                gosfiltercd49adnkegg$Term=='Steroid biosynthesis'
                              )

#Reactome_2016
up <- eup$Reactome_2016
down <- edown$Reactome_2016
up$type <- "up"
down$type <- "down"
up <- up[up$Adjusted.P.value<.05,]
up <- up[order(up$Combined.Score), ]  # sort
down <- down[down$Adjusted.P.value<0.05,]
down <- down[order(down$Combined.Score), ]  # sort
down$Combined.Score <- (-1) * down$Combined.Score
gos <- rbind(down,up)
gos <- na.omit(gos) # Diverging Barcharts
cd49adnreactome<-ggplot(gos, aes(x=reorder(Term,Combined.Score), y=Combined.Score , label=Combined.Score)) + 
  geom_bar(stat='identity', aes(fill=Adjusted.P.value), width=.5,position="dodge")  +
  scale_fill_viridis() + 
  labs(subtitle="Reactome") + 
  ylab('Combined Score')+ xlab('')+
  labs(fill = expression(P[adj]))+
  coord_flip()
cd49adnreactome


gosfiltercd49adnreact<-gos
gosfiltercd49adnreact$comparison<-'CD49a vs DN'
gosfiltercd49adnreact$database<-'Reactome'
goscd49adnreact<-gosfiltercd49adnreact
gosfiltercd49adnreact<-subset(gosfiltercd49adnreact,gosfiltercd49adnreact$Term=='Metabolism of amino acids and derivatives Homo sapiens R-HSA-71291' |
                                gosfiltercd49adnreact$Term=='Selenoamino acid metabolism Homo sapiens R-HSA-2408522' |
                                gosfiltercd49adnreact$Term=='Selenocysteine synthesis Homo sapiens R-HSA-2408557'
)


###############
theme_set(theme(axis.text = element_text(size=7),
                axis.text.x = element_text(angle=90),
                axis.title = element_text(size=7),
                strip.text = element_text(size=12),
                legend.text = element_text(size=7),
                legend.title = element_text(size=7),
                legend.position = 'bottom',
                panel.background = element_blank(),
                panel.grid.major = element_line(color='light grey')
))

library(ggallin)
metabolismenrich<-rbind(gosfiltercd49adnkegg,gosfiltercd49adnreact,gosfilterdpdnkegg,gosfilterdpdnreact)
metabolismenrich$comparison_f = factor(metabolismenrich$comparison, levels=c('DP vs DN','CD49a vs DN'))

metabolismenrich$comparison <- factor(metabolismenrich$comparison, levels=c('DP vs DN','CD49a vs DN'))



metabolismenrichplot<-ggplot(metabolismenrich, aes(x=reorder(Term,Combined.Score), y=Combined.Score , label=Combined.Score)) + 
  geom_bar(stat='identity', aes(fill=Adjusted.P.value), width=.5,position="dodge")  +
  scale_fill_viridis() + 
  ylab('Combined Score')+ xlab('')+
  labs(fill = expression(P[adj]))+
  coord_flip()+scale_y_continuous(trans = pseudolog10_trans,breaks=c(0,10,50,100,500,1000,5000,10000,50000,100000,500000,-10,-50,-100,-500,-1000,-5000,-10000,-50000,-100000,-500000))+ #seq(from=-20000,to=0,by=10000),100,300,-100,-200,-1000
  facet_wrap(~comparison,scales = 'free_x',ncol=1)
metabolismenrichplot

pathwayanalysiscsv<-rbind(gosdpdnkegg,gosdpdnreact,goscd49adnkegg,goscd49adnreact)
write.csv(pathwayanalysiscsv,'pathway-analysis.csv')
ggsave(paste0(Sys.Date(),'metabolism-enrich-plot.png'),plot = metabolismenrichplot,units = 'mm',width=180,dpi=600)

sessionInfo()
### FIN ###