# COMPLETION OF TAXA FILE

#clear the object from memory
rm(list=ls())

#### CHANGE THIS PATHS ####
source <- "your_wd"
out <- "your_output_folder"
##########################

# install required libraries
if (!require(ggplot2)){BiocManager::install("ggplot2")}
if (!require(pander)){install.packages("pander")}
if (!require(knitr)){install.packages("knitr")}
if (!require(myTAI)){install.packages("myTAI")}

#set wd
setwd(source)

## libraries
library(ggplot2)
library(pander)
library(knitr)
library(myTAI)

#read taxa file
silva18s_conf <- read.table( "ASV_pipeline_published_samples/4_output/fungi_ASV/fungi_taxa_conf_18S.tab", row.names=1, sep="\t", blank.lines.skip = FALSE)
silva18s_0.5 <- read.table( "ASV_pipeline_published_samples/4_output/fungi_ASV/fungi_taxa_0.5_18S.tab", row.names=1, sep="\t", blank.lines.skip = FALSE)


#remove assigning rows from SILVA
silva18s_0.5 <- silva18s_0.5[,1:(ncol(silva18s_0.5)-1)]
silva18s_conf <- silva18s_conf[,1:(ncol(silva18s_conf)-2)]

#rename
colnames(silva18s_0.5) <- c("domain", "subdomain1", "subdomain2", "kingdom", "order/family", "genus", "genus_species")
colnames(silva18s_conf) <- c("domain", "domain_conf", "subdomain1", "subdomain1_conf", "subdomain2", "subdomain2_conf",
                             "kingdom", "kingdom_conf", "order/family", "order/family_conf", "genus", "genus_conf", "genus_species", "genus_species_conf")


#for (i in c("silva18s_0.5", "silva18s_conf")) {
for (name in c("silva18s_0.5")) {
  
  db <- get(name)
  
  #add empty columns
  db$phylum <- NA
  db$class <- NA
  db$order <- NA
  db$family <- NA
  
  #add order and family for genus
  for (i in 1:nrow(db)) {
    genus <- as.character(db$genus[i])
    
    if(!is.na(genus)){
      tax <- taxonomy(genus)
      try(db$phylum[i] <- tax$name[tax$rank=="phylum"])
      try(db$class[i] <- tax$name[tax$rank=="class"])
      try(db$order[i] <- tax$name[tax$rank=="order"])
      try(db$family[i] <- tax$name[tax$rank=="family"])
    }
  }
  
  #add order and family where genus is unknown
  for (i in 1:nrow(db)) {
    of <- as.character(db$`order/family`[i])
    
    if(is.na(db$order[i])){
      if(is.na(db$family[i])){
        tax <- taxonomy(of)
        try(db$phylum[i] <- tax$name[tax$rank=="phylum"])
        try(db$class[i] <- tax$name[tax$rank=="class"])
        try(db$order[i] <- tax$name[tax$rank=="order"])
        try(db$family[i] <- tax$name[tax$rank=="family"])
        print(i)
      }
    }
  }
  
  assign(paste(name,"completed_unsort",sep="_"), db)
}

#change order
silva18s_0.5_completed <- silva18s_0.5_completed_unsort[,c(4,8:11,6,7)]
silva18s_conf_completed <- silva18s_conf_completed_unsort[,c(7,8,15:18,10:13)]

#output
write.table(silva18s_0.5_completed, paste(out, "fungi_taxa_SILVA18S_0.5_completed.tab", sep=""), sep="\t")
write.table(silva18s_conf_completed, paste(out, "fungi_taxa_SILVA18S_conf_completed.tab", sep=""), sep="\t")