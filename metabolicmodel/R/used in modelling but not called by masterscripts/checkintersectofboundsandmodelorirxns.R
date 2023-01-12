#this script does not need to be run, as its final output 'exchangeboundsfromR.txt' 
#is included in the github repository. I am including it for sake of reproducibility.
#this script exists to change the reaction bounds to those published in DOI: 10.1039/c5mb00480b
#in their supplemental, they have excel of constraints. They used RECON2, 
#so I manually curated them to make sure the reaction names were correct before reading into R
#read in csvs and rename columns. The first csv is the exchange bounds that were manually curated, 
#the second is the list of reactions in the model
exchangebounds<-read.csv('exchangebounds.csv')
names(exchangebounds)[4] <- "notes"
modelorirxns<-read.csv('modelorirxns.csv',header = F)
names(modelorirxns)[1] <- "rxns"

#see which reactions are in both exchange bounds table and the modelori (iMM1865)
intersect<-intersect(exchangebounds$reactions,modelorirxns$rxns)
length(intersect)

#see which reactions are in the exchange bounds table but not in model ori
inexchangeboundsbutnotori<-setdiff(exchangebounds$reactions,modelorirxns$rxns)
inexchangeboundsbutnotori
length(inexchangeboundsbutnotori)
#delete notes column that I put there from when I manually curated it.
exchangebounds<-exchangebounds[1:3]

#look at how many rows there are
nrow(exchangebounds)
#subset the exchange bounds table to only keep those that are present in model ori 
exchangebounds<-exchangebounds[exchangebounds$reactions %in% intersect,]
#make sure this worked by making sure that the number of rows decreased by 
#the same number that were found in 'inexchangeboundsbutnorori' above
nrow(exchangebounds)

#save this and then read this back into matlab
write.table(exchangebounds,'exchangeboundsfromR.txt',row.names = F,col.names =T,quote = F)
#print session info
sessionInfo()
citation()
