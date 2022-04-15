rm(list=ls())

args = commandArgs(trailingOnly = T)

sequences<-seqinr::read.fasta(paste(args, "/", args, ".fasta", sep = ""))
assignments<-read.csv(paste(args, "/Assignment/assignment.csv", sep = ""))
write.csv(assignments, paste(args, "/", args, "_assignments.csv", sep = ""), row.names = F)

clades<-data.frame(clade=c("Africa", "Asian", "Arctic", "Bat", "Cosmopolitan", "Indian", "RAC"), present=NA)

assignments<-unique(assignments$lineage)

for (i in 1:length(clades$clade)) {
  if (length(grep(clades$clade[i], assignments)) != 0) {
    clades$present[i]<-"Y"
  }
}

numbers<-which(clades$present == "Y")

sequences1<-NA
alignment1<-seqinr::as.alignment(NA)
metadata1<-NA
sequences2<-NA
alignment2<-seqinr::as.alignment(NA)
metadata2<-NA
sequences3<-NA
alignment3<-seqinr::as.alignment(NA)
metadata3<-NA
sequences4<-NA
alignment4<-seqinr::as.alignment(NA)
metadata4<-NA
sequences5<-NA
alignment5<-seqinr::as.alignment(NA)
metadata5<-NA
sequences6<-NA
alignment6<-seqinr::as.alignment(NA)
metadata6<-NA
sequences7<-NA
alignment7<-seqinr::as.alignment(NA)
metadata7<-NA

if (1 %in% numbers) {
  sequences1<-read.csv("Datasets/Africa_N/Africa_N_sequence_data.csv")
  alignment1<-seqinr::read.alignment("Datasets/Africa_N/Africa_N_aligned.fasta", format = "fasta")
  metadata1<-read.csv("Datasets/Africa_N/Africa_N_metadata.csv")
}

if (2 %in% numbers) {
  sequences2<-read.csv("Datasets/Asian_N/Asian_N_sequence_data.csv")
  alignment2<-seqinr::read.alignment("Datasets/Asian_N/Asian_N_aligned.fasta", format = "fasta")
  metadata2<-read.csv("Datasets/Asian_N/Asian_N_metadata.csv")
}

if (3 %in% numbers) {
  sequences3<-read.csv("Datasets/Arctic_N/Arctic_N_sequence_data.csv")
  alignment3<-seqinr::read.alignment("Datasets/Arctic_N/Arctic_N_aligned.fasta", format = "fasta")
  metadata3<-read.csv("Datasets/Arctic_N/Arctic_N_metadata.csv")
}

if (4 %in% numbers) {
  sequences4<-read.csv("Datasets/Bat_N/Bat_N_sequence_data.csv")
  alignment4<-seqinr::read.alignment("Datasets/Bat_N/Bat_N_aligned.fasta", format = "fasta")
  metadata4<-read.csv("Datasets/Bat_N/Bat_N_metadata.csv")
}

if (5 %in% numbers) {
  sequences5<-read.csv("Datasets/Cosmo_N/Cosmo_N_sequence_data.csv")
  alignment5<-seqinr::read.alignment("Datasets/Cosmo_N/Cosmo_N_aligned.fasta", format = "fasta")
  metadata5<-read.csv("Datasets/Cosmo_N/Cosmo_N_metadata.csv")
}

if (6 %in% numbers) {
  sequences6<-read.csv("Datasets/Indian_N/Indian_N_sequence_data.csv")
  alignment6<-seqinr::read.alignment("Datasets/Indian_N/Indian_N_aligned.fasta", format = "fasta")
  metadata6<-read.csv("Datasets/Indian_N/Indian_N_metadata.csv")
}

if (7 %in% numbers) {
  sequences7<-read.csv("Datasets/RAC-SK_N/RAC-SK_N_sequence_data.csv")
  alignment7<-seqinr::read.alignment("Datasets/RAC-SK_N/RAC-SK_N_aligned.fasta", format = "fasta")
  metadata7<-read.csv("Datasets/RAC-SK_N/RAC-SK_N_metadata.csv")
}

metadata<-rbind(metadata1, metadata2, metadata3, metadata4, metadata5, metadata6, metadata7)
metadata<-metadata[-c(which(is.na(metadata$ID))),]

sequences<-rbind(sequences1, sequences2, sequences3, sequences4, sequences5, sequences6, sequences7)
sequences<-sequences[-c(which(is.na(sequences$ID))),]

align_nam<-c(alignment1$nam, alignment2$nam, alignment3$nam, alignment4$nam, alignment5$nam, alignment6$nam,
             alignment7$nam)

align_seq<-c(alignment1$seq, alignment2$seq, alignment3$seq, alignment4$seq, alignment5$seq, alignment6$seq,
             alignment7$seq)

alignment<-seqinr::as.alignment(nb=length(align_nam), nam = align_nam, seq = align_seq, com = NA)

all_lineage<-read.csv("inst/extdata/References/RABV/lineage_info.csv")

x<-1

while(x < 10){
  assignments<-c(assignments, all_lineage$lineage[which(all_lineage$parent %in% assignments)])
  x<-x+1
}

assignments<-unique(assignments)

x<-1

while(x < 10){
  assignments<-c(assignments, all_lineage$parent[which(all_lineage$lineage %in% assignments)])
  x<-x+1
}

assignments<-unique(assignments)

int_seq<-sequences$ID[which(sequences$cluster %in% assignments)]

int_seq<-int_seq[which(int_seq %in% alignment$nam)]
int_seq<-int_seq[which(int_seq %in% metadata$ID)]

numbers<-which(alignment$nam %in% int_seq)

query_alignment<-seqinr::read.alignment(paste(args, "/", args, ".fasta", sep = ""), format = "fasta")

seqinr::write.fasta(sequences = c(alignment$seq[numbers], query_alignment$seq),
                    names = c(alignment$nam[numbers], query_alignment$nam),
                    file.out = paste(args, "/", args, "_combined.fasta", sep = ""))

query_metadata<-read.csv(paste(args, "/", args, "_metadata.csv", sep = ""))

metadata<-metadata[which(metadata$ID %in% int_seq),]

metadata<-data.frame(ID=metadata$ID, year=metadata$year, country=metadata$country, assignment=metadata$assignment)

combined_metadata<-rbind(metadata, query_metadata)

write.csv(combined_metadata, paste(args, "/", args, "_combined_metadata.csv", sep = ""), row.names = F)
