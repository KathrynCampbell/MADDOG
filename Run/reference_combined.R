rm(list=ls())
args = commandArgs(trailingOnly = T)

folder<-args[1]
WGS<-args[2]

lineage_info<-read.csv("inst/extdata/References/RABV/lineage_info.csv")
new_lineage<-read.csv(paste(folder, "/reference/", folder, "_lineage_info.csv", sep = ""))

if(length(which(names(new_lineage) == "cluster")) != 0) {
  names(new_lineage)[which(names(new_lineage) == "cluster")]<-"lineage"
}

if (length(which(new_lineage$lineage %in% lineage_info$lineage)) != 0) {
  new_lineage<-new_lineage[-c(which(new_lineage$lineage %in% lineage_info$lineage)),]
}

lineage_info<-rbind(lineage_info, new_lineage)

lineage_info$lineage<-gsub("Cosmo ", "Cosmopolitan ", lineage_info$lineage)
lineage_info$lineage<-gsub("Cosmo_", "Cosmopolitan_", lineage_info$lineage)

clusters<-read.csv("inst/extdata/References/RABV/reference_clusters.csv")
new_clusters<-read.csv(paste(folder, "/reference/", folder, "_reference.csv", sep = ""))

if (length(which(new_clusters$cluster %in% clusters$cluster)) != 0) {
  new_clusters<-new_clusters[-c(which(new_clusters$cluster %in% clusters$cluster)),]
}

clusters<-rbind(clusters, new_clusters)

clusters$cluster<-gsub("Cosmo ", "Cosmopolitan ", clusters$cluster)
clusters$cluster<-gsub("Cosmo_", "Cosmopolitan_", clusters$cluster)

if (length(new_lineage$lineage) != 0) {
  sequences<-seqinr::read.alignment("inst/extdata/References/RABV/reference_aligned.fasta", format = "fasta")
  new_sequences<-seqinr::read.alignment(paste(folder, "/reference/", folder, "_reference.fasta", sep = ""), format = "fasta")

  if (length(which(new_sequences$nam %in% sequences$nam)) != 0) {
    new_sequences$seq<-new_sequences$seq[-c(which(new_sequences$nam %in% sequences$nam))]
    new_sequences$nam<-new_sequences$nam[-c(which(new_sequences$nam %in% sequences$nam))]
  }

  sequence_data<-c(sequences$seq, new_sequences$seq)

  names<-c(sequences$nam, new_sequences$nam)

  seqinr::write.fasta(sequence_data, names, file.out = "inst/extdata/References/RABV/seq.fasta")
  write.csv(clusters, "inst/extdata/References/RABV/reference_clusters.csv", row.names = F)
  write.csv(lineage_info, "inst/extdata/References/RABV/lineage_info.csv", row.names = F)
} else {
  print("No new lineages to add at this stage")
}

