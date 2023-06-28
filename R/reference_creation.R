library(ggtree)
library(seqinr)
library(adephylo)
library(phangorn)

rm(list=ls())

args = ""

tree<-read.tree("Asian/Trees/Asian_combined_aligned.fasta.contree")
tree$tip.label <- gsub("\\..*", "", tree$tip.label, perl = T)
tree$node.comment<- gsub(".*=", "", tree$node.label, perl = T)

sequence_data<-read.csv("Asian/Outputs/sequence_data.csv")
lineage_info<-read.csv("Asian/Outputs/new_lineages.csv")
lineage_info<-data.frame(lineage = c("Asian SEA1b_C1.2", "Asian SEA1b_C1.2.1"))
alignment<-read.alignment("Datasets/Asian_N/Asian_N_aligned.fasta", format = "fasta")
alignment$nam <- gsub("\\..*", "", alignment$nam, perl = T)

distances<-as.matrix(distTips(tree, tips = "all", method = "patristic"))
lineages<-rep(lineage_info$lineage, times = 4)
reference_set<-data.frame(lineage = lineages, sequence = NA)
reference_set<-reference_set[order(reference_set$lineage),]

#sequence_data$lineage<-sequence_data$cluster

for (x in 1:length(lineage_info$lineage)) {
  test<-sequence_data$ID[which(sequence_data$lineage == lineage_info$lineage[x])]
  subset<-distances[which(colnames(distances) %in% test),which(rownames(distances) %in% test)]

  subset<-as.data.frame(subset)

  names<-names(which((apply(subset, 2, max)) == max(subset)))
  names<-c(names, names(which((apply(subset, 1, max)) == max(subset))))
  reference<-unique(names)

  test<-(subset[which(rownames(subset) == reference[1]),])
  test<-t(test)

  if (length(test) %% 2 == 0) {
    reference<-c(reference, rownames(test)[length(test)/2])
    reference<-c(reference, rownames(test)[(length(test)/2)+1])
  } else {
    reference<-c(reference, rownames(test)[(length(test)+1)/2])
    reference<-c(reference, NA)
  }

  reference<-unique(reference)

  for (i in 1:4) {
    reference_set$sequence[which(reference_set$lineage == lineage_info$lineage[x])[i]]<-reference[i]

  }
}

reference_set<-reference_set[-c(which(is.na(reference_set$sequence))),]

numbers<-which(alignment$nam %in% reference_set$sequence)

write.fasta(sequences = alignment$seq[numbers], names = alignment$nam[numbers], file.out =
              'Datasets/Asian_N/reference.fasta')

reference<-data.frame(ID = reference_set$sequence, cluster = reference_set$lineage)

write.csv(reference, file = "Datasets/Asian_N/reference.csv", row.names = F)
