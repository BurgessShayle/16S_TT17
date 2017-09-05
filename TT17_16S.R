# 20170815 phloseq the 16S data for TT17
# OTU table made in mothur

rm(list=ls()) #clears workspace 


biocpkgs_needed = c("ggplot2","genefilter", "impute",
                    "phyloseq")
pkgs_needed = c("PMA","dplyr","ade4","ggrepel","vegan")
letsinstall = setdiff(pkgs_needed, installed.packages()) 
letsinstallb = setdiff(biocpkgs_needed, installed.packages()) 
if (length(letsinstallb) > 0) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(letsinstallb)
}

if (length(letsinstall) > 0) {
  install.packages(letsinstall)
}
#Read in required libraries
##### Include Versions of libraries
library('phyloseq')
library("ggplot2")
library("PMA")
library("dplyr")
library("ade4")
library("genefilter")
library("ggrepel")
library("phyloseq")
library("vegan")
library("grid")
library("reshape2")

#set working directory
setwd("~/Desktop/16S_Seqs")

# Set plotting theme
theme_set(theme_bw())


###### Mike Lee's ideas 
## reading in the otu count table:
  raw_otu_tab<-read.table("TT16_OTU_Table_mothur_edited2.txt", header=T)
head(raw_otu_tab)
dim(raw_otu_tab)
 # make row names right
rownames(raw_otu_tab) <- raw_otu_tab$Group 
raw_otu_tab$Group = NULL
# transform
raw_otu_tab <- as.data.frame(t(raw_otu_tab))
head(raw_otu_tab)
class(raw_otu_tab$`1037`)


## reading in the taxonomy file:
  otu_taxon_tab<-as.matrix(read.csv("Master_fasta_QCd_merged.good.unique.good.filter.precluster.pick.pick.opti_mcc.unique_list.0.03.cons.taxonomy_GenusTab.csv", header=T, row.names=1, na.strings=c("", "NA")))
head(otu_taxon_tab)

## reading in a sample information file (that has info about each sample, examples of each of these files are attached):
  sample_info_tab<-read.table("frag_Info_16S_phyloseq.txt", header=T)
head(sample_info_tab)




## making the otu table phyloseq-ready:
  OTU<-otu_table(raw_otu_tab, taxa_are_rows=T)

## now for taxon table: 
  OTU_TAX<-tax_table(otu_taxon_tab)
  head(OTU_TAX)
## now for the sample information: 
  sample_data<-sample_data(sample_info_tab)
  
  

  
  
  
  
  
  

  
  

## now putting them together into a phyloseq object: 
  OTU_physeq<-phyloseq(OTU, OTU_TAX, sample_data) 
























# Assign variables for imported data
sharedfile = "Master_fasta_QCd_merged.good.unique.good.filter.precluster.pick.pick.opti_mcc.unique_list.shared"
taxfile = "Master_fasta_QCd_merged.good.unique.good.filter.precluster.pick.rdp.wang.pick.taxonomy"
mapfile = "frag_Info_phyloseq.csv"

# Import mothur data
mothur_data <- import_mothur(mothur_shared_file = sharedfile,
                             mothur_constaxonomy_file = taxfile)


# Import sample metadata
map <- read.csv(mapfile)

# convert dataframe to phyloseq format
map <- sample_data(map)

# Assign rownames to be Sample ID's
rownames(map) <- map$SampleID

# Merge mothurdata object with sample metadata
moth_merge <- merge_phyloseq(mothur_data, map)
moth_merge


# trying holmes
mothlist = system.file("extdata", "esophagus.fn.list.gz", package = "phyloseq")
mothgroup = system.file("extdata", "esophagus.good.groups.gz", package = "phyloseq")
mothtree = system.file("extdata", "esophagus.tree.gz", package = "phyloseq")
cutoff = "0.10"
esophman = import_mothur(mothlist, mothgroup, mothtree, cutoff)
print(esophman)
