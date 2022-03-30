rm(list=ls())

library(stringr)
library(seqinr)
library(plotly)
library(dplyr)

source("R/designation_updates.r")
args = commandArgs(trailingOnly = T)

folder<-args[1]
WGS<-args[2]

N<-seqinr::read.alignment(paste(folder, "/Alignment/", folder, "_aligned.fasta", sep = ""), format = "fasta")
N_meta<-read.csv(paste(folder, "/", folder, "_metadata.csv", sep = ""))
WGS_lineages<-read.csv("inst/extdata/References/RABV/lineage_info.csv")
WGS_meta<-read.csv(paste("inst/extdata/WGS/", WGS, "/", WGS, "_WGS_metadata.csv", sep = ""))

lineages<-data.frame(original = unique(WGS_lineages$lineage), subclade = NA)

for(i in 1:length(lineages$original)) {
  lineages$subclade[i]<-str_split(lineages$original[i], "_")[[1]][1]
}

N_subset<-N_meta[which(N_meta$assignment %in% unique(unlist(lineages$subclade))),]

N_subset<-N_subset[-c(which(N_subset$ID %in% WGS_meta$ID)),]

useful_numbers<-which(N$nam %in% N_subset$ID)
N_subset_seq<-seqinr::as.alignment(nb = length(useful_numbers),
                     nam = N$nam[useful_numbers],
                     seq = N$seq[useful_numbers],
                     com = NA)

write.fasta(sequences = N_subset_seq$seq, names = N_subset_seq$nam, file.out =
              paste(folder, "/assignment/", folder, "_WGS_assignment.fasta", sep = ""))

#############################################
#            IMPORT THE DATA                #
#############################################
#'
#'**TREE**
#'========================================================================================================
#' The tree must contain the element 'node.comment' which contains the bootstrap support/posterior support
#' And the element 'tip.label' which lists all the sequence ID's
#' These sequence ID's must match the sequence ID's in the metadata and alignment
#'=========================================================================================================
tree <- ape::read.tree(file = paste(folder, "/Trees/", folder, "_aligned.fasta.contree", sep = ""))
# Sequence names got messed up in MAFFT, need to fix these so they match metadata and alignment
# Also node comment is sometimes weird, fix it
# #KB- can replace above 2 lines with this:
tree$tip.label <- gsub("\\..*", "", tree$tip.label, perl = T)
tree$node.comment<- gsub(".*=", "", tree$node.label, perl = T)

#'**METADATA**
#'========================================================================================================
#' The metadata must contain the element 'year' which lists the collection year for each sequence
#' And the element 'ID' which lists all the sequence ID's
#' These sequence ID's must match the sequence ID's in the tree and alignment
#'=========================================================================================================
metadata <- read.csv(file = paste(folder, "/", folder, "_metadata.csv", sep = ""))

#'**ALIGNMENT**
#'========================================================================================================
#' The alignment must contain the element 'seq' which contains the sequences
#' And the element 'nam' which lists all the sequence ID's
#' These sequence ID's must match the sequence ID's in the metadata and tree
#'=========================================================================================================
alignment <- read.alignment(file = (paste(folder, "/Alignment/", folder, "_aligned.fasta", sep = "")), format = "fasta")

# Sequence names got messed up in MAFFT, need to fix these so they match metadata and alignment
# #KB- can replace above 2 lines with this:
alignment$nam <- gsub("\\..*", "", alignment$nam, perl = T)


#'**TIMETREE**
#'========================================================================================================
#'
#'=========================================================================================================
ancestral <- read.alignment(file = (paste(folder, "/Timetree/ancestral_sequences.fasta", sep = "")), format = "fasta")
ancestral$nam <- gsub("\\..*", "", ancestral$nam, perl = T)


sequence_data <- lineage_assignment(tree, min.support = 70, max.support = 100, alignment, metadata, ancestral)[[2]]
node_data <- lineage_assignment(tree, min.support = 70, max.support = 100, alignment, metadata, ancestral)[[1]]

node_data_copy<-node_data

sequence_data$previous <- NA
for (i in 1:length(sequence_data$ID)) {
  sequence_data$previous[i]<-
    metadata$assignment[which(metadata$ID == sequence_data$ID[i])]
}

previous_assignments<-data.frame(assignment = unique(sequence_data$previous), node = NA)

node_data$previous<-NA

for (i in 1:length(node_data$Node)) {
  clades<-unique(sequence_data$previous[
    which(sequence_data$ID %in% tree$tip.label[c(unlist(
      phangorn::Descendants(tree, node_data$Node[i], type = "tips")))])])

  node_data$previous[i]<-
    paste(c(clades), collapse = ", ")

}

for (i in 1:length(previous_assignments$assignment)) {
  previous_assignments$node[i]<-which(node_data$previous == previous_assignments$assignment[i])[1]
  previous_assignments$assignment[i]<-previous_assignments$assignment[i]
}

possible_names<-data.frame(names = rep(previous_assignments$assignment, 26))
previous_assignments$assignment<-paste(previous_assignments$assignment, "_A1", sep = "")

for (i in 1:length(previous_assignments$assignment)) {
  node_data$cluster[previous_assignments$node[i]]<-previous_assignments$assignment[i]
}

lineages_2<-data.frame(original=unique(metadata$assignment), subclade=NA)

for(i in 1:length(lineages_2$original)) {
  lineages_2$subclade[i]<-str_split(lineages_2$original[i], "_")[[1]][1]
}

clades<-unique(lineages_2$subclade)[-c(which((unique(lineages_2$subclade)) %in% unique(lineages$subclade)))]

previous_assignments<-previous_assignments[which(previous_assignments$assignment %in% clades),]

if (length(which(node_data$previous %in% clades)) != 0) {
  node_data<-node_data[which(node_data$previous %in% clades),]
  node_data$test <- NA
  problem_names<-data.frame(letters = c("A1", "B1", "C1", "D1", "E1", "F1", "G1", "H1", "I1", "J1", "K1", "L1", "M1", "N1",
                                        "O1", "P1", "Q1", "R1", "S1", "T1", "U1", "V1", "W1", "X1", "Y1", "Z1"))
  possible_names<-possible_names[order(possible_names$names),]
  possible_names<-paste(possible_names, problem_names$letters, sep = "_")

  for (i in 1:length(node_data$Node)) {
    test<-which(node_data$Node %in% ips::descendants(tree, node_data$Node[i], type = "all", ignore.tip = T))
    node_data$test[c(test)] <- paste(node_data$cluster[i], ".1", sep = "")
    node_data$test<-str_replace(node_data$test, "A1\\..\\..\\..", "B1")
    node_data$test<-str_replace(node_data$test, "B1\\..\\..\\..", "C1")
    node_data$test<-str_replace(node_data$test, "C1\\..\\..\\..", "D1")
    node_data$test<-str_replace(node_data$test, "D1\\..\\..\\..", "E1")
    node_data$test<-str_replace(node_data$test, "E1\\..\\..\\..", "F1")
    node_data$test<-str_replace(node_data$test, "F1\\..\\..\\..", "G1")
    node_data$test<-str_replace(node_data$test, "G1\\..\\..\\..", "H1")

    majors<-which(grepl("_", node_data$test))
    node_data$cluster[c(majors)] <- node_data$test[c(majors)]

    for (k in 1:length(possible_names)) {
      if (length(which(node_data$cluster == possible_names[k]))>1) {
        problems<-which(node_data$cluster == possible_names[k])
        problems<-problems[-c(1)]
        y=1
        for (a in 1:length(problems)) {
          letter<-which(problem_names$letters == (str_split(node_data$cluster[problems[a]], "_")[[1]][2]))
          node_data$cluster[problems[a]]<-paste((str_split(node_data$cluster[problems[a]], "_")[[1]][1]), problem_names$letters[(letter+y)], sep = "_")
          y = y+1
        }
      }
    }
    duplicates<-unique(node_data$cluster[duplicated(node_data$cluster)])
    problems<-duplicates[which(str_count(duplicates, pattern = "\\.") == 0)]
    duplicates<-duplicates[which(str_count(duplicates, pattern = "\\.") != 0)]

    for (i in 1:length(duplicates)) {
      test<-which(node_data$cluster == duplicates[i])
      test<-test[-c(1)]
      x<-1
      for (j in 1:length(test)) {
        name<-unlist(str_split(node_data$cluster[test[j]], "\\."))
        name[length(name)]<-x+as.integer(name[length(name)])
        x<-(x+1)
        node_data$cluster[test[j]]<-paste(c(name), collapse='.' )
      }
    }
  }

  unclassified<-which(!grepl("_", node_data$cluster))
  unclassified<-unclassified[c(-1)]
  for (i in 1:length(node_data$Node)) {
    test<-which(node_data$Node %in% ips::descendants(tree, node_data$Node[i], type = "all", ignore.tip = T))
    node_data$test[c(test)] <- paste(node_data$cluster[i], ".1", sep = "")
    node_data$test<-str_replace(node_data$test, "A1\\..\\..\\..", "B1")
    node_data$test<-str_replace(node_data$test, "B1\\..\\..\\..", "C1")
    node_data$test<-str_replace(node_data$test, "C1\\..\\..\\..", "D1")
    node_data$test<-str_replace(node_data$test, "D1\\..\\..\\..", "E1")
    node_data$test<-str_replace(node_data$test, "E1\\..\\..\\..", "F1")
    node_data$test<-str_replace(node_data$test, "F1\\..\\..\\..", "G1")
    node_data$test<-str_replace(node_data$test, "G1\\..\\..\\..", "H1")

    node_data$cluster[unclassified]<-node_data$test[unclassified]

    for (v in 1:length(problem_names$letters)) {
      if (length(which(node_data$cluster == problem_names$letters[v]))>1) {
        problems<-which(node_data$cluster == problem_names$letters[v])
        problems<-problems[-c(1)]
        y=1
        for (f in 1:length(problems)) {
          letter<-which(problem_names$letters == (node_data$cluster[problems[f]]))
          node_data$cluster[problems[f]]<-problem_names$letters[(letter+y)]
          y = y+1
        }
      }
    }
    duplicates<-unique(node_data$cluster[duplicated(node_data$cluster)])
    problems<-duplicates[which(str_count(duplicates, pattern = "\\.") == 0)]
    duplicates<-duplicates[which(str_count(duplicates, pattern = "\\.") != 0)]

    for (i in 1:length(duplicates)) {
      test<-which(node_data$cluster == duplicates[i])
      test<-test[-c(1)]
      x<-1
      for (j in 1:length(test)) {
        name<-unlist(str_split(node_data$cluster[test[j]], "\\."))
        name[length(name)]<-x+as.integer(name[length(name)])
        x<-(x+1)
        node_data$cluster[test[j]]<-paste(c(name), collapse='.' )
      }
    }
  }
} else {
  print("No new subclades to add")
}

for (i in 1:length(node_data$Node)) {
  node_data_copy$cluster[which(node_data_copy$Node == node_data$Node[i])]<-node_data$cluster[i]
}

node_data<-node_data_copy
node_data$number<-1:length(node_data$Node)

previous_assignments<-data.frame(assignment = unique(sequence_data$previous), node = NA)

rename<-(1:length(node_data$cluster))[-c(which(node_data$cluster %in% 1:1000))]

for (i in 1:length(node_data$cluster)) {
  sequence_data$cluster[which(sequence_data$cluster == i)] <- node_data$cluster[i]
}


write.csv(sequence_data, file = (paste(folder, "/temp_designation/sequence_data.csv", sep = "")), row.names=F)
write.csv(node_data, file = (paste(folder, "/temp_designation/node_data.csv", sep = "")), row.names=F)

#############################################
#         LINEAGE INFORMATION TABLE         #
#############################################

# Extract the place information from the metadata
for (i in 1:length(sequence_data$ID)) {
  sequence_data$Country[i] <- metadata$country[(which(metadata$ID == sequence_data$ID[i]))]
}

# Create a data frame ready to fill in information about each cluster
clusters <- sequence_data %>%
  group_by(cluster) %>%
  summarise()
clusters$country <- NA
clusters$year_first <- NA
clusters$year_last <- NA

# Add another column listing the number of sequences assigned to each cluster
clusters$n_seqs<-(sequence_data %>%
                    group_by(cluster)%>%
                    summarise(n=n()))$n

if(length(which(is.na(clusters$cluster)))!= 0){
  clusters<-clusters[-c(which(is.na(clusters$cluster))),]
}
# For each cluster, find and list the earliest collection year, the latest collection year and all the places
# that cluster has been found
for (i in 1:length(clusters$cluster)) {
  clusters$year_first[i] <- sequence_data %>%
    filter(cluster == clusters$cluster[i])%>%
    group_by(Year)%>%
    summarise()%>%
    min()
  clusters$year_last[i] <- sequence_data %>%
    filter(cluster == clusters$cluster[i])%>%
    group_by(Year)%>%
    summarise()%>%
    max()
  clusters$country[i]<-paste((sequence_data %>%
                                filter(cluster == clusters$cluster[i] & Country !="-") %>%
                                group_by(Country) %>%
                                summarise()))
}

write.csv(clusters, file = (paste(folder, "/temp_designation/_all_lineage_info.csv", sep = "")), row.names=F)

# For each lineage, calculate the pairwise distance for all the sequences allocated to each lineage
# Extract the mean and max distance for each lineage

distances<-as.matrix(adephylo::distTips(tree, tips = "all", method = "patristic"))
node_data<-node_data[-c(which(node_data$cluster %in% 1:1000)),]
lineages<-rep(node_data$cluster, times = 4)
reference_set<-data.frame(lineage = lineages, sequence = NA)
reference_set<-reference_set[order(reference_set$lineage),]

clusters<-clusters[-c(which(clusters$cluster %in% 1:1000)),]

lineage_info<-clusters
lineage_info$lineage<-lineage_info$cluster
sequence_data$lineage<-sequence_data$cluster

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

if(length(which(is.na(reference_set$sequence)))!= 0){
  reference_set<-reference_set[-c(which(is.na(reference_set$sequence))),]
}

numbers<-which(alignment$nam %in% reference_set$sequence)

write.fasta(sequences = alignment$seq[numbers], names = alignment$nam[numbers],
            file.out = (paste(folder, "/reference/", folder, "_reference.fasta", sep="")))

reference<-data.frame(ID = reference_set$sequence, cluster = reference_set$lineage)

write.csv(reference, file = (paste(folder, "/reference/", folder, "_reference.csv", sep="")), row.names = F)

write.csv(clusters, file = (paste(folder, "/reference/",  folder, "_lineage_info.csv", sep = "")), row.names=F)

