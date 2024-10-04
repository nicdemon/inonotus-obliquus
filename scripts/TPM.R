#Load données
setwd("../data/")
library(readr)
db <- read_csv("Tab AA genes.csv", col_types = cols(`Reads BETU(TPM)` = col_number(), 
                                                      +     `Reads EBB (TPM)` = col_number(), `Reads STD (TPM)` = col_number()))
View(db)

#Faire barplot
par(mfrow=c(3,1))
barplot(db$`Reads STD (TPM)`,names.arg=db$Name,ylab='Mean expression (TPM)')
barplot(db$`Reads BETU(TPM)`,names.arg=db$Name,ylab='Mean expression (TPM)')
barplot(db$`Reads EBB (TPM)`,names.arg=db$Name,ylab='Mean expression (TPM)')
par(mfrow=c(1,1))

#Barplot de total
library(readr)
db <- read_csv("Tab AA genes.csv", col_types = cols(`Reads moyens (TPM)` = col_number()))
barplot(db$`Reads moyens (TPM)`,names.arg=db$Name,cex.lab=1.5,cex.axis=1.5,ylim=c(0,20),ylab = 'Mean expression of gene isoform (TPM)')