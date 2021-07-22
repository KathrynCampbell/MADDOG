#' Lineage Assignment
#' 
#' This function tests query sequences against a reference set of sequences and provides a lineage
#' assignment for the query sequences. It also returns information about each of the lineages assigned.
#' 
#' @param sequences The query sequences to provide a lineage assignment for in fasta format
#' @param reference The reference set to compare the query sequences to (Cosmo_WGS, Cosmo_N)
#' @return A lineage assignment for each query sequence, and information about that lineage
#' @export
assign_lineages<-function(sequences, reference) {
  alignment<-seqinr::read.alignment(system.file("extdata", paste("References/", reference, "/reference_aligned.fasta", sep = ""), package = "MADDOG"), format = "fasta")
  data <- read.csv(system.file("extdata", paste("References/", reference, "/reference_clusters.csv", sep = ""), package = "MADDOG"))
  clusters <- read.csv(system.file("extdata", paste("References/", reference, "/lineage_info.csv", sep = ""), package = "MADDOG"))
  
  #SOMETHING HERE TO ADD SEQUENCES TO ALIGNMENT AND REALIGN
  
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
    down<-data$cluster[which(data$ID == alignment$nam[(test_seqs[i]-1)])]
    repeat {
      down<-data$cluster[which(data$ID == alignment$nam[(test_seqs[i]-(1+x))])]
      x<-x+1
      if (length(down) != 0){
        break
      }
    }
    up<-data$cluster[which(data$ID == alignment$nam[(test_seqs[i]+1)])]
    repeat {
      up<-data$cluster[which(data$ID == alignment$nam[(test_seqs[i]+(1+y))])]
      y<-y+1
      if (length(up) != 0){
        break
      }
    }
    test<-c(down, up)
    test_seq_assignment$lineage[i]<-calculate_mode(test)
  }

  test_seq_assignment$lineage_countries_seen<-NA
  test_seq_assignment$lineage_first_seen<-NA
  test_seq_assignment$lineage_last_seen<-NA

  for (i in 1:length(test_seq_assignment$ID)) {
    test_seq_assignment$lineage_countries_seen[i]<-clusters$country[which(clusters$cluster == test_seq_assignment$lineage[i])]
    test_seq_assignment$lineage_first_seen[i]<-clusters$year_first[which(clusters$cluster == test_seq_assignment$lineage[i])]
    test_seq_assignment$lineage_last_seen[i]<-clusters$year_last[which(clusters$cluster == test_seq_assignment$lineage[i])]
  }


  return(test_seq_assignment)
}

