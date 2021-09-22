#!/usr/bin/env Rscript

rm(list=ls())

args = commandArgs(trailingOnly = T)

#'---------------------------------------------------------
#'title: Lineage Assignment Windows
#'author: Kathryn Campbell
#'date: 09/09/2021
#'---------------------------------------------------------
data <- read.csv(system.file("extdata", paste("References/", args[2], "/reference_clusters.csv", sep = ""), package = "MADDOG"))
clusters <- read.csv(system.file("extdata", paste("References/", args[2], "/lineage_info.csv", sep = ""), package = "MADDOG"))

alignment<-seqinr::read.alignment(file = (paste(args[1], "/", args[1], "_withref.fasta", sep = "")), format = "fasta")

`%notin%` <- Negate(`%in%`)

test_seqs<-which(alignment$nam %notin% data$ID)

test_seq_assignment<-data.frame(ID=alignment$nam[test_seqs], lineage = NA)

calculate_mode <- function(x) {
  uniqx<-unique(na.omit(x))
  uniqx[which.max(tabulate(match(x, uniqx)))]
}


for (i in 1:length(test_seqs)){
  x<-1
  y<-1
  test<-1
  down<-data$cluster[which(data$ID == row.names(alignment)[(test_seqs[i]-1)])]
  repeat {
    down<-data$cluster[which(data$ID == row.names(alignment)[(test_seqs[i]-(1+x))])]
    x<-x+1
    if (length(down) != 0){
      break
    }
    if (x > length(ref_align)){
      test<-NA
      break
    }
  }
  up<-data$cluster[which(data$ID == row.names(alignment)[(test_seqs[i]+1)])]
  repeat {
    up<-data$cluster[which(data$ID == row.names(alignment)[(test_seqs[i]+(1+y))])]
    y<-y+1
    if (length(up) != 0){
      break
    }
    if (y > length(ref_align)){
      test<-NA
      break
    }
  }
  if (!is.na(test)){
    test<-c(down, up)
  }
  test_seq_assignment$lineage[i]<-calculate_mode(test)
}

test_seq_assignment$lineage_countries_seen<-NA
test_seq_assignment$lineage_first_seen<-NA
test_seq_assignment$lineage_last_seen<-NA

numbers<-1:length(test_seq_assignment$ID)
numbers<-numbers[-c(  which(is.na(test_seq_assignment$lineage)))]
for (i in numbers) {
  test_seq_assignment$lineage_countries_seen[i]<-clusters$country[which(clusters$cluster == test_seq_assignment$lineage[i])]
  test_seq_assignment$lineage_first_seen[i]<-clusters$year_first[which(clusters$cluster == test_seq_assignment$lineage[i])]
  test_seq_assignment$lineage_last_seen[i]<-clusters$year_last[which(clusters$cluster == test_seq_assignment$lineage[i])]
}


write.csv(test_seq_assignment, file = (paste(args[1], "/", args[1], "_assignment.csv", sep = "")), row.names = F)

