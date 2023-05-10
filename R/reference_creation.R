library(ggtree)
library(seqinr)
library(adephylo)
library(phangorn)

rm(list=ls())

args = "Vac_WGS"

tree<-read.tree("~/Downloads/Phil_designation/Trees/Phil_designation_combined_aligned.fasta.contree")
tree$tip.label <- gsub("\\..*", "", tree$tip.label, perl = T)
tree$node.comment<- gsub(".*=", "", tree$node.label, perl = T)

sequence_data<-read.csv("~/Downloads/Phil_designation/Outputs/sequence_data.csv")
lineage_info<-read.csv("~/Downloads/Phil_designation/Outputs/new_lineages.csv")
alignment<-read.alignment("~/Downloads/Phil_designation/Alignment/Phil_designation_combined_aligned.fasta", format = "fasta")
alignment$nam <- gsub("\\..*", "", alignment$nam, perl = T)

distances<-as.matrix(distTips(tree, tips = "all", method = "patristic"))
lineages<-rep(lineage_info$lineage, times = 4)
reference_set<-data.frame(lineage = lineages, sequence = NA)
reference_set<-reference_set[order(reference_set$lineage),]

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

reference_set$sequence[1]<-"R2-081"

numbers<-which(alignment$nam %in% reference_set$sequence)

write.fasta(sequences = alignment$seq[numbers], names = alignment$nam[numbers], file.out =
              '~/Downloads/Phil_designation/reference.fasta')

reference<-data.frame(ID = reference_set$sequence, cluster = reference_set$lineage)

write.csv(reference, file = "~/Downloads/Phil_designation/reference.csv", row.names = F)
