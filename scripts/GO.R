setwd("../data/GO")

#Installer packages pr manip données
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler", version = "3.8")

#Load database
library(readr)
db <- read_csv("trinotate_annotation_report_modified.csv")

#GOstats package
source("https://bioconductor.org/biocLite.R")
biocLite("GOstats")
library(GOstats)

#Tabulate variable
###Make sure the Na's are listed as NA and not other characters (.,:; etc.)###
counts=table(na.omit(db$Gene_ontology))
counts2=table(na.omit(db$GO_short))
GO=data.frame(counts2)
#Avec labels qui buggent
barplot(counts2,
        main = 'Fonctional Gene Ontology found in Inonotus obliquus',
        xlab = 'GO number',
        ylab = 'Quantity',
        xpd = FALSE,
        axisnames = TRUE,
        names.arg = GO$Var1,
        las = 2)
#Sans labels
barplot(counts2,
        main = 'Fonctional Gene Ontology found in Inonotus obliquus',
        xlab = 'GO number',
        ylab = 'Quantity',
        xpd = FALSE,
        axisnames = TRUE,
        names.arg = FALSE,
        las = 2)

###COG###
#Extraire COG
setwd("../data/GO/COG")
library(readr)
categories <- read_csv("cog2cat42163_10-nov-2018.xls.csv")
db <- read_excel("GO_COG.xlsx", col_names = FALSE)
View(db)
counts=table(na.omit(db$X__2))
COG=data.frame(counts)
colnames(COG)=c('COG ID','Freq')
#Merge les tableaux
merged=merge(COG,categories)
PT=aggregate(merged$Freq~merged$Catergory,merged,sum)
PT$legend=c('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','X','Z')
colnames(PT)=c('levels','Frequency')
#COG graph
library(ggplot2)
p <- ggplot(data=PT, aes(x=levels, y=Frequency))+geom_bar(stat="identity", position=position_dodge(), colour='black')
p + coord_flip() + scale_x_discrete(limits = rev(levels(PT$legend))) + guides (fill = guide_legend(ncol = 1))+ylab("Frequency")+xlab('Factor class')+ggtitle("COG Representation by functional category")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black")) + geom_text(aes(label=Frequency), position=position_dodge(width=0.9), hjust=-0.25)
#COG graph en %
merged$Percent=(merged$Freq/sum(merged$Freq)*100)
datp <- data.frame(FunctionClass = factor(merged$`Catergory Code`), levels=merged$`Catergory Code`,legend = merged$Catergory,Percentage=merged$Percent)
library(ggplot2)
pp <- ggplot(data=datp, aes(x=FunctionClass, y=Percentage, fill=legend))+geom_bar(stat="identity", position=position_dodge(), colour="seashell")
pp + coord_flip() + guides (fill = guide_legend(ncol = 1))+xlab("Frequency")+ylab('Factor class')+ggtitle("COG Percentage by functional category")
#Tableau de %
Pt=aggregate(Percent~`Catergory Code`,merged,sum)
write.csv(Pt,file='Percent table.csv')

#COG/traitement
library(readxl)
BET_COG <- read_excel("COG/BET/BET-COG.xlsx",col_names = FALSE)
EBB_COG <- read_excel("COG/EBB/EBB-COG.xlsx",col_names = FALSE)
STD_COG <- read_excel("COG/STD/STD-COG.xlsx",col_names = FALSE)
counts_BET=table(na.omit(BET_COG$X__2))
counts_EBB=table(na.omit(EBB_COG$X__2))
counts_STD=table(na.omit(STD_COG$X__2))
COG_BET=data.frame(counts_BET)
COG_EBB=data.frame(counts_EBB)
COG_STD=data.frame(counts_STD)
colnames(COG_BET)=c('COG ID','Freq')
colnames(COG_EBB)=c('COG ID','Freq')
colnames(COG_STD)=c('COG ID','Freq')
merged_BET=merge(COG_BET,categories)
merged_EBB=merge(COG_EBB,categories)
merged_STD=merge(COG_STD,categories)
PT_BET=aggregate(merged_BET$Freq~merged_BET$Catergory,merged_BET,sum)
PT_EBB=aggregate(merged_EBB$Freq~merged_EBB$Catergory,merged_EBB,sum)
PT_STD=aggregate(merged_STD$Freq~merged_STD$Catergory,merged_STD,sum)
colnames(PT_BET)=c('levels','Frequency')
colnames(PT_EBB)=c('levels','Frequency')
colnames(PT_STD)=c('levels','Frequency')
library(ggplot2)
p_BET <- ggplot(data=PT_BET, aes(x=levels, y=Frequency))+geom_bar(stat="identity", position=position_dodge(), colour='Seashell')
p_BET + coord_flip() + scale_x_discrete(limits = rev(levels(PT_BET$legend))) + guides (fill = guide_legend(ncol = 1))+ylab("Frequency")+xlab('Factor class')+ggtitle("COG Representation by functional category for betuline treatment") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ geom_text(aes(label=Frequency), position=position_dodge(width=0.9), hjust=-0.25)
p_EBB <- ggplot(data=PT_EBB, aes(x=levels, y=Frequency))+geom_bar(stat="identity", position=position_dodge(), colour='Seashell')
p_EBB + coord_flip() + scale_x_discrete(limits = rev(levels(PT_EBB$legend))) + guides (fill = guide_legend(ncol = 1))+ylab("Frequency")+xlab('Factor class')+ggtitle("COG Representation by functional category for bark treatment") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ geom_text(aes(label=Frequency), position=position_dodge(width=0.9), hjust=-0.25)
p_STD <- ggplot(data=PT_STD, aes(x=levels, y=Frequency))+geom_bar(stat="identity", position=position_dodge(), colour='Seashell')
p_STD + coord_flip() + scale_x_discrete(limits = rev(levels(PT_STD$legend))) + guides (fill = guide_legend(ncol = 1))+ylab("Frequency")+xlab('Factor class')+ggtitle("COG Representation by functional category for standard treatment") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+ geom_text(aes(label=Frequency), position=position_dodge(width=0.9), hjust=-0.25)
Pt_BET=aggregate(Freq~`Catergory Code`,merged_BET,sum)
Pt_EBB=aggregate(Freq~`Catergory Code`,merged_EBB,sum)
Pt_STD=aggregate(Freq~`Catergory Code`,merged_STD,sum)
write.csv(Pt_BET,file='Percent table_BET.csv')
write.csv(Pt_EBB,file='Percent table_EBB.csv')
write.csv(Pt_STD,file='Percent table_STD.csv')
Pt=aggregate(Freq~`Catergory Code`,merged,sum)
write.csv(Pt,file='Percent table_total.csv')

###Graph abondance pfam###
#Load data
GO_pfam <- read_excel("Total/GO_pfam.xlsx", col_names = FALSE)
categories <- read_csv("Total/pfamlist48269_20-nov-2018.xls.csv")
#Faire disparaitre les NA's et merger
counts=as.data.frame(table(na.omit(GO_pfam$X__2)))
colnames(counts)=c('Pfam ID','Freq')
merged=merge(counts,categories)
#Ordoner pfam les plus de frequents et isoler 20
GO_sort=merged[order(-merged$Freq),]
GO_short=as.data.frame(GO_sort$Freq[1:20])
#Pfam en rownames
names <- paste(GO_sort$`Protein name`[1:20], sep=":")
rownames(GO_short)=names
GO_short$`Pfam ID`=NULL
GO_short$`Pfam Name`=NULL
colnames(GO_short)=c('Freq')
dat <- data.frame(FunctionClass = factor(row.names(GO_short)),legend=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'), levels=row.names(GO_short),Frequency=GO_short$Freq)
#Plot graph
library(ggplot2)
p <- ggplot(data=dat, aes(x=reorder(levels,-Frequency), y=Frequency))+geom_bar(stat="identity", position=position_dodge(), colour='Seashell')
p + coord_flip() + guides (fill = guide_legend(ncol = 1))+ylab("Frequency")+xlab('Protein family name')+ggtitle("Protein families representation by functional category")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black")) + geom_text(aes(label=Frequency), position=position_dodge(width=0.9), hjust=-0.25)
#Pfam multiples
setwd("../data/GO/pfam")
library(readr)
categories <- read_csv("Total/pfamlist48269_20-nov-2018.xls.csv")
library(readxl)
STD_GO_pfam <- read_excel("STD/STD-GO_pfam.xlsx", col_names = FALSE)
EBB_GO_pfam <- read_excel("EBB/EBB-GO_pfam.xlsx",col_names = FALSE)
BET_GO_pfam <- read_excel("BET/BET-GO_pfam.xlsx", col_names = FALSE)
BET_counts=as.data.frame(table(na.omit(BET_GO_pfam$X__2)))
EBB_counts=as.data.frame(table(na.omit(EBB_GO_pfam$X__2)))
STD_counts=as.data.frame(table(na.omit(STD_GO_pfam$X__2)))
colnames(BET_counts)=c('Pfam ID','Freq')
colnames(EBB_counts)=c('Pfam ID','Freq')
colnames(STD_counts)=c('Pfam ID','Freq')
BET_merged=merge(BET_counts,categories)
EBB_merged=merge(EBB_counts,categories)
STD_merged=merge(STD_counts,categories)
BET_GO_sort=BET_merged[order(-BET_merged$Freq),]
EBB_GO_sort=EBB_merged[order(-EBB_merged$Freq),]
STD_GO_sort=STD_merged[order(-STD_merged$Freq),]
BET_GO_short=as.data.frame(BET_GO_sort$Freq[1:20])
EBB_GO_short=as.data.frame(EBB_GO_sort$Freq[1:20])
STD_GO_short=as.data.frame(STD_GO_sort$Freq[1:20])
BET_names <- paste(BET_GO_sort$`Pfam ID`[1:20], BET_GO_sort$`Pfam Name`[1:20], sep=":")
rownames(BET_GO_short)=BET_names
EBB_names <- paste(EBB_GO_sort$`Pfam ID`[1:20], EBB_GO_sort$`Pfam Name`[1:20], sep=":")
rownames(EBB_GO_short)=EBB_names
STD_names <- paste(STD_GO_sort$`Pfam ID`[1:20], STD_GO_sort$`Pfam Name`[1:20], sep=":")
rownames(STD_GO_short)=STD_names
colnames(BET_GO_short)=c('Freq')
colnames(EBB_GO_short)=c('Freq')
colnames(STD_GO_short)=c('Freq')
dat_BET <- data.frame(FunctionClass = factor(row.names(BET_GO_short)),legend=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'), levels=row.names(BET_GO_short),Frequency=BET_GO_short$Freq)
dat_EBB <- data.frame(FunctionClass = factor(row.names(EBB_GO_short)),legend=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'), levels=row.names(EBB_GO_short),Frequency=EBB_GO_short$Freq)
dat_STD <- data.frame(FunctionClass = factor(row.names(STD_GO_short)),legend=c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'), levels=row.names(STD_GO_short),Frequency=STD_GO_short$Freq)
library(ggplot2)
p <- ggplot(data=dat_BET, aes(x=reorder(levels,-Frequency), y=Frequency,fill=legend))+geom_bar(stat="identity", position=position_dodge(), colour='Seashell')
p + coord_flip() + guides (fill = guide_legend(ncol = 1))+ylab("Frequency")+xlab('Factor class')+ggtitle("Protein families representation by functional category for betuline treatment")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_text(aes(label=Frequency), position=position_dodge(width=0.9), hjust=-0.25)
p <- ggplot(data=dat_EBB, aes(x=reorder(levels,-Frequency), y=Frequency,fill=legend))+geom_bar(stat="identity", position=position_dodge(), colour='Seashell')
p + coord_flip() + guides (fill = guide_legend(ncol = 1))+ylab("Frequency")+xlab('Factor class')+ggtitle("Protein families representation by functional category for bark treatment")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_text(aes(label=Frequency), position=position_dodge(width=0.9), hjust=-0.25)
p <- ggplot(data=dat_STD, aes(x=reorder(levels,-Frequency), y=Frequency,fill=legend))+geom_bar(stat="identity", position=position_dodge(), colour='Seashell')
p + coord_flip() + guides (fill = guide_legend(ncol = 1))+ylab("Frequency")+xlab('Factor class')+ggtitle("Protein families representation by functional category for standard treatment")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_text(aes(label=Frequency), position=position_dodge(width=0.9), hjust=-0.25)