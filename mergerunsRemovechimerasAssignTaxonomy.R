#script written by M.Wutkowska 20190429

R.version #R version 3.5.2 (2018-12-20)

#clean the environment
rm(list = ls(all = TRUE))
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)


#set dir
getwd()
setwd("/Users/magdalenawutkowska/Dropbox/PhD/PROJECT_3/r-dada2_test/")


#read multiple runs
fencesA <- readRDS("./seqtab_fencesA.rds")
fencesB <- readRDS("./seqtab_fencesB.rds")
largeA <- readRDS("./seqtab_largeA.rds") 
largeB <-readRDS("./seqtab_largeB.rds")
tempA <- readRDS("./seqtab_tempA.rds")
tempB <- readRDS("./seqtab_tempB.rds")


#merge all the ASV tables
library(dada2) #packageVersion("dada2") [1] ‘1.11.1
project3_ASVtable <- mergeSequenceTables(fencesA, fencesB, largeA, largeB, tempA, tempB)
dim(project3_ASVtable) #472 11243
rownames(project3_ASVtable)
colnames(project3_ASVtable)


#remove chimeras
project3_ASVtable_nochim <- removeBimeraDenovo(project3_ASVtable, method="consensus", multithread=TRUE, verbose = TRUE) #Identified 3759 bimeras out of 11243 input sequences.
dim(project3_ASVtable_nochim) #472 7484
class(project3_ASVtable_nochim)


#inspect the data
rownames(project3_ASVtable_nochim)
boxplot(rowSums(project3_ASVtable_nochim))
hist(rowSums(project3_ASVtable_nochim))
boxplot(colSums(project3_ASVtable_nochim))
colnames(project3_ASVtable_nochim)


#now let's check seq lengths and remove these that are shorter than 200bp
table(nchar(getSequences(project3_ASVtable_nochim)))
seqlens <- nchar(getSequences(project3_ASVtable_nochim))
project3_ASVtable_rightlen <- project3_ASVtable_nochim[,seqlens >= 200] #53 sequences were removed
dim(project3_ASVtable_rightlen) #472 7431


#let's sum sequences from forward(A) and reverse (B) orientations into samples, to do that the first step would be to:
  #(1)convert row names into a first column

df_asvtable <- as.data.frame(project3_ASVtable_rightlen)
dim(df_asvtable)

library(dplyr)
library(tibble)
df_asvtable_tbl <- as_tibble(rownames_to_column(df_asvtable))
dim(df_asvtable_tbl)
remove_rownames(df_asvtable_tbl)
rownames(df_asvtable_tbl)
colnames(df_asvtable_tbl)

class(df_asvtable_tbl)
class(df_asvtable_tbl$rowname)
class(df_asvtable_tbl$ACGCAAGTTGCGCCCGAAGCCTTCGGGCCGAGGGCACGTCTGCATGGGCGTCACGCACAGCGTCGCCCCCACCCCACTCGTGGGGCGTGGGGCGGATTCTGGCCCCCCGTGTGCTCCCGCGCGCGGTCGGCCTAAAATCAGACCCCGTGGCCGCGAAATGCCGCGACGATTGGTGGTGTACGTGGCGGCCTCGAGCCTCCGAACATCGCGTCGCGCCTTCCGTGGCCCTCTGGAGTCAAGAGGACCCTCGAGAGCCCTCCGCCGGTGCGGAGGGGCCTCTCAACCGTTGCGACCCCATGTCAGGCGGGACTACCCGCTGAGTTTAA)


  #(2)remove whatever is after "_" in sample column
df_asvtable_tbl$rowname <- sapply(strsplit(basename(df_asvtable_tbl$rowname), "_"), `[`, 1)
as.factor(df_asvtable_tbl$rowname)


  #(3)sum values for samples, first convert wide format data frame into long and aggregate based on sample name and variable
library(reshape2)
df_asvtable_tbl_long <- melt(df_asvtable_tbl, id.vars = c("rowname"))
df_asvtable_long_sum <- aggregate(. ~ rowname + variable, df_asvtable_tbl_long, sum)


  #(4)go back to wide format
#library(reshape2)
df_asvtable_wide <- dcast(df_asvtable_long_sum, rowname ~ variable)
dim(df_asvtable_wide)
class(df_asvtable_wide)


#checked sample names (column 'rowname') with a table with Bistorta masterfile
#remove samples here in the df that do not correspond with any Bistorta data (3A,4A,5A,6A,7A,Blank):
df_asvtable_wide[,1]
df_asvtable_complete <- df_asvtable_wide[-c(140, 141, 142, 143, 144, 169), ]


#putting back samples names to the rownames
rownames(df_asvtable_complete) <- df_asvtable_complete$rowname
dim(df_asvtable_complete)
df_asvtable_complete <- df_asvtable_complete[,-1]
df_asvtable_complete <- df_asvtable_complete[, -colSums(df_asvtable_complete)<=1]
dim(df_asvtable_complete)

#reading in a file with environmental variables and Bistorta measurments
env_bistorta <- read.delim("/Users/magdalenawutkowska/Dropbox/PhD/PROJECT_3/data/env_bistorta.txt", header = TRUE, row.names = 3, sep = "\t")


#number of all detected AVSs
library(vegan) #This is vegan 2.5-4
env_bistorta$alldetectedASVs <- specnumber(df_asvtable_complete)
boxplot(env_bistorta$alldetectedASVs)
env_bistorta$noofseq <- rowSums(df_asvtable_complete)
cor.test(env_bistorta$alldetectedASVs, env_bistorta$noofseq) # 0.6290395 


#no of sequences per sample sorted from the smallest to the highest value
sort(env_bistorta$noofseq) #where should I make a cutoff? >20 000 seq?


#rarefy
df_asvtable_forrarefaction <- df_asvtable_complete[rowSums(df_asvtable_complete)>=21639,] #this removes 15 samples from the analysis
df_asvtable_removedinrarefaction <- df_asvtable_complete[rowSums(df_asvtable_complete)<21639,]

asvtable_rarefied <- rrarefy(df_asvtable_forrarefaction, 21639)
dim(asvtable_rarefied)

#make sure that dim of asvtable_rariefied and env_bistorta are the same
rownames(df_asvtable_removedinrarefaction)
rownames(env_bistorta)

env_bistorta_rarefied <- env_bistorta[-c(37, 46, 72, 82, 87, 89, 90, 93, 112, 186, 187, 192, 196, 210, 211), ] #these are the 15 samples that had less than 21639 reads
dim(asvtable_rarefied)
dim(env_bistorta_rarefied) #excellent! both 214 samples
colnames(env_bistorta_rarefied)
env_bistorta_rarefied <- env_bistorta_rarefied[, -c(35:42)]


#count number of ASVs in rarified samples and check how does that correlate with number of detected ASV before rarefaction
env_bistorta_rarefied$ASVs_rarefied <- specnumber(asvtable_rarefied)

boxplot(env_bistorta_rarefied$ASVs_rarefied)
cor.test(env_bistorta_rarefied$alldetectedASVs, env_bistorta_rarefied$ASVs_rarefied) # 0.9525898 #that correlates really strongly


# Assign taxonomy to non-zero abundance ASVs
asvtable_rarefied_nozerocolumns <- asvtable_rarefied[,colSums(asvtable_rarefied)!=0] 
dim(asvtable_rarefied_nozerocolumns) #214 6665, which means that 749ASVs were removed

tax_rarefied <- assignTaxonomy(asvtable_rarefied_nozerocolumns, "./sh_general_release_dynamic_02.02.2019.fasta", multithread=TRUE, verbose = TRUE)
summary(tax_rarefied)


# Save the output
saveRDS(asvtable_rarefied_nozerocolumns, "./asvtable_rarefied_nozerocolumns.rds") #ASV table
write.table(asvtable_rarefied_nozerocolumns, file="asvtable_rarefied_nozerocolumns.txt", row.names=TRUE, col.names=TRUE, sep = '\t')

saveRDS(tax_rarefied, "./taxtable_rarefied.rds") #taxonomy table
write.table(tax_rarefied, file="taxtable_rarefied.txt", row.names=TRUE, col.names=TRUE, sep = '\t')

write.table(env_bistorta_rarefied, file = "env_bistorta_rarefied.txt", row.names = TRUE, col.names = TRUE, sep = '\t')