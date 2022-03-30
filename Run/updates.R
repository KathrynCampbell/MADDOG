rm(list=ls())

library(phangorn)
library(dplyr)

args = commandArgs(trailingOnly = T)

folder<-args[1]
WGS_arg<-args[2]

WGS<-read.csv("inst/extdata/References/RABV/lineage_info.csv")
assignments<-read.csv(paste(folder, "/assignment/", folder, "_assignment.csv", sep = ""))
node_data<-read.csv(paste(folder, "/temp_designation/node_data.csv", sep = ""))
tree<-ape::read.tree(paste(folder, "/Trees/", folder, "_aligned.fasta.contree", sep = ""))
WGS_seq<-read.csv(paste("inst/extdata/WGS/", WGS_arg, "/", WGS_arg, "_WGS_sequence_data.csv", sep = ""))

WGS_seq$lineage<-gsub("Cosmo_", "Cosmopolitan_", WGS_seq$lineage)

node_data<-node_data[(which(node_data$cluster %in% 1:1000)),]

WGS$cluster<-WGS$lineage

current<-data.frame(lineage = unique(WGS_seq$lineage), node = NA)
current<-current[-c(which(is.na(current$lineage))),]

for (i in 1:length(current$lineage)) {
  current$node[i]<-ape::getMRCA(tree, tip = c(WGS_seq$ID[which(WGS_seq$lineage == current$lineage[i])]))
}

updates<-data.frame(lineage = WGS$cluster, count = NA, node = NA)

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
    length(which(node_data$Node %in% Descendants(tree, updates$node[i], type = "all"))) -
    (length(which(updates$node %in% Descendants(tree, updates$node[i], type = "all"))) +
       length(which(current$node[-c(which(current$node %in% node_data$Node))]
                    %in% Descendants(tree, updates$node[i], type = "all")))
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
      test$node<-node_data$Node[which(node_data$Node %in% Descendants(tree, updates$node[numbers[i]], type = "all"))][1]
      test$lineage<-paste(updates$lineage[numbers[i]], ".1", sep = "")
      updates<-rbind(updates, test)
    }

    for (i in 1:length(updates$lineage)) {
      updates$viable[i]<-
        length(which(node_data$Node %in% Descendants(tree, updates$node[i], type = "all"))) -
        (length(which(updates$node %in% Descendants(tree, updates$node[i], type = "all"))) +
           length(which(current$node[-c(which(current$node %in% node_data$Node))]
                        %in% Descendants(tree, updates$node[i], type = "all")))
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

    duplicates<-int$update[which(duplicated(int$update))]
    duplicates<-c(duplicates, int$update[which(int$update %in% c(WGS$cluster, existing$lineage))])
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
        duplicates<-c(duplicates, int$update[which(int$update %in% WGS$cluster)])
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
      duplicates<-c(duplicates, int$update[which(int$update %in% WGS$cluster)])
      problems<-duplicates[which(stringr::str_count(duplicates, pattern = "\\.") == 0)]
    }

    for (i in 1:length(int$lineage)) {
      updates$lineage[which(updates$node == int$node[i])]<-int$update[i]
    }
    for (i in 1:length(updates$lineage)) {
      updates$viable[i]<-
        length(which(node_data$Node %in% Descendants(tree, updates$node[i], type = "all"))) -
        (length(which(updates$node %in% Descendants(tree, updates$node[i], type = "all"))) +
           length(which(current$node[-c(which(current$node %in% node_data$Node))]
                        %in% Descendants(tree, updates$node[i], type = "all")))
        )
    }
    numbers<-which(updates$viable > 0)
    x<-x+1
  }

  updates<-updates[-c(which(updates$lineage %in% WGS$cluster)),]

  node_updates<-data.frame()

  for (i in 1:length(updates$lineage)) {
    node_updates<-rbind(node_updates, node_data[which(node_data$Node == updates$node[i]),])
    node_updates$numbers[i]<-updates$lineage[i]
  }
  sequences<-read.csv(paste(folder, "/temp_designation/sequence_data.csv", sep = ""))
  for (i in 1:length(node_updates$Node)) {
    sequences$cluster[which(sequences$cluster == node_updates$cluster[i])]<-node_updates$numbers[i]
  }
} else {
  print("no lineages to update")
  sequences<-read.csv(paste(folder, "/temp_designation/sequence_data.csv", sep = ""))
}


lineage_info<-read.csv(paste(folder, "/temp_designation/_all_lineage_info.csv", sep = ""))

sequences<-sequences[which(sequences$ID %in% assignments$ID),]

numbers<-which(sequences$cluster %in% 1:1000)


for (i in 1:length(numbers)) {
  sequences$cluster[numbers[i]]<-assignments$lineage[which(assignments$ID == sequences$ID[numbers[i]])]
}

metadata<-read.csv(paste(folder, "/", folder, "_metadata.csv", sep = ""))

# Extract the place information from the metadata
for (i in 1:length(sequences$ID)) {
  sequences$Country[i] <- metadata$country[(which(metadata$ID == sequences$ID[i]))]
}

# Create a data frame ready to fill in information about each cluster
clusters <- sequences %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise()
clusters$country <- NA
clusters$year_first <- NA
clusters$year_last <- NA

# Add another column listing the number of sequences assigned to each cluster
clusters$n_seqs<-(sequences %>%
                    dplyr::group_by(cluster)%>%
                    dplyr::summarise(n=dplyr::n()))$n

if(length(which(is.na(clusters$cluster)))!= 0){
  clusters<-clusters[-c(which(is.na(clusters$cluster))),]
}

# For each cluster, find and list the earliest collection year, the latest collection year and all the places
# that cluster has been found
for (i in 1:length(clusters$cluster)) {
  clusters$year_first[i] <- sequences %>%
    filter(cluster == clusters$cluster[i])%>%
    dplyr::group_by(Year)%>%
    dplyr::summarise()%>%
    min()
  clusters$year_last[i] <- sequences %>%
    filter(cluster == clusters$cluster[i])%>%
    dplyr::group_by(Year)%>%
    dplyr::summarise()%>%
    max()
  clusters$country[i]<-paste((sequences %>%
                                filter(cluster == clusters$cluster[i] & Country !="-") %>%
                                dplyr::group_by(Country) %>%
                                dplyr::summarise()))
}

combine<-data.frame(lineage = unique(WGS_seq$lineage), country_WGS = NA,
                    country_N = NA, year_first_WGS = NA, year_first_N = NA,
                    year_last_WGS = NA, year_last_N = NA, n_WGS_seqs = NA, n_N_seqs = NA)

if(length(which((is.na(combine$lineage)))) != 0) {
  combine<-combine[-c(which(is.na(combine$lineage))),]
}


for (i in 1:length(combine$lineage)) {
  combine$country_WGS[i]<-WGS$country[which(WGS$lineage == combine$lineage[i])]
  combine$year_first_WGS[i]<-WGS$year_first[which(WGS$lineage == combine$lineage[i])]
  combine$year_last_WGS[i]<-WGS$year_last[which(WGS$lineage == combine$lineage[i])]
  combine$n_WGS_seqs[i]<-WGS$n_seqs[which(WGS$lineage == combine$lineage[i])]
}

N_add<-which(combine$lineage %in% clusters$cluster)

for (i in 1:length(N_add)) {
  combine$n_N_seqs[N_add[i]]<-clusters$n_seqs[which(clusters$cluster == combine$lineage[N_add[i]])]
  combine$country_N[N_add[i]]<-clusters$country[which(clusters$cluster == combine$lineage[N_add[i]])]
  combine$year_first_N[N_add[i]]<-clusters$year_first[which(clusters$cluster == combine$lineage[N_add[i]])]
  combine$year_last_N[N_add[i]]<-clusters$year_last[which(clusters$cluster == combine$lineage[N_add[i]])]
}

`%notin%` = Negate(`%in%`)

if(length(which(clusters$cluster[grep(WGS_arg, clusters$cluster)] %notin% combine$lineage)) != 0) {
  numbers<-which(clusters$cluster[grep(WGS_arg, clusters$cluster)] %notin% combine$lineage)
  for (i in 1:length(numbers)) {
    y<-length(combine$lineage)
    combine[(y+1),1]<-clusters$cluster[numbers[i]]
    combine[(y+1),3]<-clusters$country[numbers[i]]
    combine[(y+1),5]<-clusters$year_first[numbers[i]]
    combine[(y+1),7]<-clusters$year_last[numbers[i]]
    combine[(y+1),9]<-clusters$n_seqs[numbers[i]]
  }
}

strange<-assignments[c(which(assignments$lineage %in% clusters$cluster[-c(grep(WGS_arg, clusters$cluster))])),]

if(length(which(WGS$lineage[grep(WGS_arg, WGS$lineage)] %notin% combine$lineage)) != 0) {
  numbers<-grep(WGS_arg, WGS$lineage)[which(WGS$lineage[grep(WGS_arg, WGS$lineage)] %notin% combine$lineage)]
  for (i in 1:length(numbers)) {
    y<-length(combine$lineage)
    combine[(y+1),1]<-WGS$lineage[numbers[i]]
    combine[(y+1),3]<-WGS$country[numbers[i]]
    combine[(y+1),5]<-WGS$year_first[numbers[i]]
    combine[(y+1),7]<-WGS$year_last[numbers[i]]
    combine[(y+1),9]<-WGS$n_seqs[numbers[i]]
  }
}

write.csv(combine, file = paste(folder, "/Outputs/", folder, "_updated_lineage_info.csv", sep = ""), row.names = F)
write.csv(strange, file = paste(folder, "/Outputs/", folder, "_strange_assignments.csv", sep = ""), row.names = F)
write.csv(sequences, file = paste(folder, "/Outputs/", folder, "_sequence_data.csv", sep = ""), row.names = F)

new<-combine$lineage[which(combine$lineage %notin% WGS$lineage)]

lineages<-rep(new, times = 4)
reference_set<-data.frame(lineage = lineages, sequence = NA)
reference_set<-reference_set[order(reference_set$lineage),]

alignment<-seqinr::read.alignment(paste(folder, "/Alignment/", folder, "_aligned.fasta", sep = ""), format = "fasta")

distances<-as.matrix(adephylo::distTips(tree, tips = "all", method = "patristic"))
for (x in 1:length(new)) {
  test<-sequences$ID[which(sequences$cluster == new[x])]
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
    reference_set$sequence[which(reference_set$lineage == new[x])[i]]<-reference[i]

  }
}

if(length(which(is.na(reference_set$sequence)))!= 0){
  reference_set<-reference_set[-c(which(is.na(reference_set$sequence))),]
}

if (length(grep("subset", reference_set$sequence)) != 0) {
  reference_set<-reference_set[-c(grep("subset", reference_set$sequence)),]
}


numbers<-which(alignment$nam %in% reference_set$sequence)

seqinr::write.fasta(sequences = alignment$seq[numbers], names = alignment$nam[numbers],
            file.out = (paste(folder, "/reference/", folder, "_reference.fasta", sep="")))

reference<-data.frame(ID = reference_set$sequence, cluster = reference_set$lineage)

write.csv(reference, file = (paste(folder, "/reference/", folder, "_reference.csv", sep="")), row.names = F)

write.csv(clusters, file = (paste(folder, "/reference/",  folder, "_lineage_info.csv", sep = "")), row.names=F)

