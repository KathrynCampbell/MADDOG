rm(list=ls())

args = commandArgs(trailingOnly = T)

tree<-ape::read.tree(paste(args, "/Trees/", args, "_combined_aligned.fasta.contree", sep = ""))
ancestral<-seqinr::read.alignment(paste(args, "/Timetree/ancestral_sequences.fasta", sep = ""), format = "fasta")
metadata<-read.csv(paste(args, "/", args, "_combined_metadata.csv", sep = ""))
alignment<-seqinr::read.alignment(paste(args, "/Alignment/", args, "_combined_aligned.fasta", sep = ""), format = "fasta")

node_data<-MADDOG::node_info(tree, 70, alignment, metadata, ancestral)
seq_data<-MADDOG::seq_designation(tree, 70, alignment, metadata, ancestral)

names(node_data)[5]<-"number"

for (i in 1:length(node_data$lineage)) {
  seq_data$lineage[which(seq_data$lineage == node_data$lineage[i])]<-node_data$number[i]
}

lineage_info<-MADDOG::lineage_info(seq_data, metadata)

assignments<-read.csv(paste(args, "/", args, "_assignments.csv", sep = ""))

clades<-data.frame(clade=c("Africa", "Asian", "Arctic", "Bat", "Cosmopolitan", "Indian", "RAC"), present=NA)

assignment_names<-unique(assignments$lineage)

for (i in 1:length(clades$clade)) {
  if (length(grep(clades$clade[i], assignment_names)) != 0) {
    clades$present[i]<-"Y"
  }
}

assignments$lineage<-gsub("Cosmopolitan ", "", assignments$lineage)
assignments$lineage<-gsub("Cosmopolitan_", "", assignments$lineage)

numbers<-which(clades$present == "Y")

sequences1<-NA
sequences2<-NA
sequences3<-NA
sequences4<-NA
sequences5<-NA
sequences6<-NA
sequences7<-NA

if (1 %in% numbers) {
  sequences1<-read.csv("Datasets/Africa_N/Africa_N_sequence_data.csv")
}
if (2 %in% numbers) {
  sequences2<-read.csv("Datasets/Asian_N/Asian_N_sequence_data.csv")
}

if (3 %in% numbers) {
  sequences3<-read.csv("Datasets/Arctic_N/Arctic_N_sequence_data.csv")
}

if (4 %in% numbers) {
  sequences4<-read.csv("Datasets/Bat_N/Bat_N_sequence_data.csv")
}

if (5 %in% numbers) {
  sequences5<-read.csv("Datasets/Cosmo_N/Cosmo_N_sequence_data.csv")
}

if (6 %in% numbers) {
  sequences6<-read.csv("Datasets/Indian_N/Indian_N_sequence_data.csv")
}

if (7 %in% numbers) {
  sequences7<-read.csv("Datasets/RAC-SK_N/RAC-SK_N_sequence_data.csv")
}

sequences<-rbind(sequences1, sequences2, sequences3, sequences4, sequences5, sequences6, sequences7)
sequences<-sequences[-c(which(is.na(sequences$ID))),]

numbers<-which(clades$present == "Y")

current<-data.frame(lineage=unique(assignments$lineage), node=NA)

for (i in 1:length(current$lineage)) {
  current$node[i]<-ape::getMRCA(tree, tip = c(sequences$ID[which(sequences$cluster == current$lineage[i])]))
}

current<-current[-c(which(duplicated(current$node))),]

node_data<-node_data[-c(which(node_data$node %in% current$node)),]

updates<-data.frame(lineage = current$lineage, count = NA, node = NA)

for (i in 1:length(updates$lineage)) {
  updates$count[i]<-length(which(assignments$lineage == updates$lineage[i]))
}

updates<-updates[which(updates$count >= 10),]

for (i in 1:length(updates$lineage)) {
  updates$node[i]<-ape::getMRCA(tree, tip = c(assignments$ID[which(assignments$lineage == updates$lineage[i])]))
}

if (length(which(duplicated(updates$node))) != 0) {
  updates<-updates[-c(which(duplicated(updates$node))),]
}

updates$viable<-NA

for (i in 1:length(updates$lineage)) {
  updates$viable[i]<-
    length(which(node_data$node %in% phangorn::Descendants(tree, updates$node[i], type = "all"))) -
    (length(which(updates$node %in% phangorn::Descendants(tree, updates$node[i], type = "all"))) +
       length(which(current$node[-c(which(current$node %in% node_data$node))]
                    %in% phangorn::Descendants(tree, updates$node[i], type = "all")))
    )
}


if(length(which(updates$viable <= 0)) != 0) {
  updates<-updates[-c(which(updates$viable <= 0)),]
}

if (length(updates$lineage) != 0) {
  existing<-updates

  updates<-updates[order(updates$viable),]

  test<-data.frame(lineage = NA, count = NA, node = NA, viable = NA)

  numbers<-which(updates$viable > 0)

  x<-1

  while (length(numbers) != 0 && x < 100) {
    test<-data.frame(lineage = NA, count = NA, node = NA, viable = NA)
    for (i in 1:length(numbers)) {
      test$node<-node_data$node[which(node_data$node %in% phangorn::Descendants(tree, updates$node[numbers[i]], type = "all"))][1]
      test$lineage<-paste(updates$lineage[numbers[i]], ".1", sep = "")
      updates<-rbind(updates, test)
    }

    for (i in 1:length(updates$lineage)) {
      updates$viable[i]<-
        length(which(node_data$node %in% phangorn::Descendants(tree, updates$node[i], type = "all"))) -
        (length(which(updates$node %in% phangorn::Descendants(tree, updates$node[i], type = "all"))) +
           length(which(current$node[-c(which(current$node %in% node_data$node))]
                        %in% phangorn::Descendants(tree, updates$node[i], type = "all")))
        )
    }

    if (length(which(duplicated(updates$node))) != 0) {
      updates<-updates[-c(which(duplicated(updates$node))),]
    }
    int<-updates[-c(which(updates$lineage %in% existing$lineage)),]
    int<-rbind(int, updates[which(duplicated(updates$lineage)),])

    problem_names<-data.frame(letters = c("A1", "B1", "C1", "D1", "E1", "F1", "G1", "H1", "I1", "J1", "K1", "L1", "M1", "N1",
                                          "O1", "P1", "Q1", "R1", "S1", "T1", "U1", "V1", "W1", "X1", "Y1", "Z1"))

    int$update<-int$lineage
    int$update<-stringr::str_replace(int$update, "A1\\..\\..\\..", "B1")
    int$update<-stringr::str_replace(int$update, "B1\\..\\..\\..", "C1")
    int$update<-stringr::str_replace(int$update, "C1\\..\\..\\..", "D1")
    int$update<-stringr::str_replace(int$update, "D1\\..\\..\\..", "E1")
    int$update<-stringr::str_replace(int$update, "E1\\..\\..\\..", "F1")
    int$update<-stringr::str_replace(int$update, "F1\\..\\..\\..", "G1")
    int$update<-stringr::str_replace(int$update, "G1\\..\\..\\..", "H1")
    int$update<-stringr::str_replace(int$update, "H1\\..\\..\\..", "I1")
    int$update<-stringr::str_replace(int$update, "I1\\..\\..\\..", "J1")
    int$update<-stringr::str_replace(int$update, "J1\\..\\..\\..", "K1")
    int$update<-stringr::str_replace(int$update, "K1\\..\\..\\..", "L1")
    int$update<-stringr::str_replace(int$update, "L1\\..\\..\\..", "M1")
    int$update<-stringr::str_replace(int$update, "M1\\..\\..\\..", "N1")


    duplicates<-int$update[which(duplicated(int$update))]
    duplicates<-c(duplicates, int$update[which(int$update %in% c(sequences$cluster, existing$lineage))])
    problems<-duplicates[which(stringr::str_count(duplicates, pattern = "\\.") == 0)]
    duplicates<-duplicates[which(stringr::str_count(duplicates, pattern = "\\.") != 0)]

    while (length(duplicates != 0)) {
      for (i in 1:length(duplicates)) {
        test<-which(int$update == duplicates[i])
        x<-1
        for (j in 1:length(test)) {
          name<-unlist(stringr::str_split(int$update[test[j]], "\\."))
          name[length(name)]<-x+as.integer(name[length(name)])
          x<-(x+1)
          int$update[test[j]]<-paste(c(name), collapse='.' )
        }
        duplicates<-int$update[which(duplicated(int$update))]
        duplicates<-c(duplicates, int$update[which(int$update %in% existing$lineage)])
        duplicates<-duplicates[which(stringr::str_count(duplicates, pattern = "\\.") != 0)]
      }
    }

    while (length(problems != 0)) {
      for (i in 1:length(problems)) {
        test<-which(int$update == problems[i])
        if(length(strsplit(problems[i], "_")[[1]]) != 1){
          subclade<-strsplit(problems[i], "_")[[1]][1]
          lineage<-strsplit(problems[i], "_")[[1]][2]
          lineage<-problem_names$letters[(which(problem_names$letters == lineage))+1]
          int$update[test[1]]<-paste(subclade, lineage, sep = "_")
        } else {
          int$update[test[1]]<-problem_names$letters[(which(problem_names$letters == problems[i]))+1]
        }
      }
      duplicates<-int$update[which(duplicated(int$update))]
      duplicates<-c(duplicates, int$update[which(int$update %in% sequences$cluster)])
      problems<-duplicates[which(stringr::str_count(duplicates, pattern = "\\.") == 0)]
    }

    for (i in 1:length(int$lineage)) {
      updates$lineage[which(updates$node == int$node[i])]<-int$update[i]
    }
    for (i in 1:length(updates$lineage)) {
      updates$viable[i]<-
        length(which(node_data$node %in% phangorn::Descendants(tree, updates$node[i], type = "all"))) -
        (length(which(updates$node %in% phangorn::Descendants(tree, updates$node[i], type = "all"))) +
           length(which(current$node[-c(which(current$node %in% node_data$node))]
                        %in% phangorn::Descendants(tree, updates$node[i], type = "all")))
        )
    }
    numbers<-which(updates$viable > 0)
    x<-x+1
  }

  updates<-updates[-c(which(updates$lineage %in% sequences$cluster)),]

  node_updates<-data.frame()

  for (i in 1:length(updates$lineage)) {
    node_updates<-rbind(node_updates, node_data[which(node_data$node == updates$node[i]),])
    node_updates$lineage[i]<-updates$lineage[i]
  }
  for (i in 1:length(node_updates$node)) {
    seq_data$lineage[which(seq_data$lineage == node_updates$number[i])]<-node_updates$lineage[i]
  }
} else {
  print("no lineages to update")
}

for (i in 1:length(node_updates$node)) {
  lineage_info$lineage[which(lineage_info$lineage == node_updates$number[i])]<-node_updates$lineage[i]
}

lineage_info<-lineage_info[-c(which(lineage_info$lineage %in% 1:1000)),]

numbers<-which(assignments$ID %in% seq_data$ID[which(seq_data$lineage %in% 1:1000)])

for (i in 1:length(numbers)) {
  seq_data$lineage[which(seq_data$ID == assignments$ID[numbers[i]])]<-assignments$lineage[numbers[i]]
}

seq_data<-seq_data[-c(which(is.na(seq_data$lineage))),]

seq_data<-seq_data[-c(which(seq_data$lineage %in% 1:1000)),]

new_seq<-seq_data[which(seq_data$ID %in% assignments$ID),]

all_lineage<-read.csv("inst/extdata/References/RABV/lineage_info.csv")

all_lineage$lineage<-gsub("Cosmopolitan ", "", all_lineage$lineage)
all_lineage$lineage<-gsub("Cosmopolitan_", "", all_lineage$lineage)

all_lineage<-all_lineage[which(all_lineage$lineage %in% assignments$lineage),]

lineage_info$parent<-NA

numbers<-grep("\\.", lineage_info$lineage)

for (i in 1:length(numbers)) {
  parent<-strsplit(lineage_info$lineage[grep("\\.", lineage_info$lineage)], "\\.")[[i]]
  parent<-parent[1:length(parent)-1]

  if (length(parent) == 1) {
    lineage_info$parent[numbers[i]]<-parent
  } else {
    parent<-paste(parent[1], parent[2], sep = ".")
    lineage_info$parent[numbers[i]]<-parent
  }
}

numbers<-which(is.na(lineage_info$parent))

if (length(numbers)!=0){
  for (i in 1:length(numbers)) {
    lineage_info$parent[numbers[i]]<-
      unique(assignments$lineage[
        which(assignments$ID %in% seq_data$ID[which(seq_data$lineage == lineage_info$lineage[numbers[i]])])])[1]
  }
}


write.csv(lineage_info, paste(args, "/Outputs/new_lineages.csv", sep = ""), row.names = F)
write.csv(all_lineage, paste(args, "/Outputs/relevant_lineages.csv", sep = ""), row.names = F)
write.csv(new_seq, paste(args, "/Outputs/sequence_data.csv", sep = ""), row.names = F)

