#' Lineage Information
#'
#' This function tests query sequences against a reference set of sequences and provides a lineage
#' assignment for the query sequences. It also returns information about each of the lineages assigned.
#'
#' @param sequence_data The output of the seq_designation function
#' @param metadata The metadata corresponding to the sequences in the data table, including "ID" and "country"
#' @return A table of information about each of the designated lineages
#' @export
lineage_info<-function(sequence_data, metadata) {
  for (i in 1:length(sequence_data$ID)) {
    sequence_data$country[i] <- metadata$country[(which(metadata$ID == sequence_data$ID[i]))]
  }

  # Create a data frame ready to fill in information about each cluster
  clusters <- sequence_data %>%
    dplyr::group_by(lineage) %>%
    dplyr::summarise()
  clusters$country <- NA
  clusters$year_first <- NA
  clusters$year_last <- NA
  # Add another column listing the number of sequences assigned to each cluster
  clusters$n_seqs<-(sequence_data %>%
                      dplyr::group_by(lineage)%>%
                      dplyr::summarise(n=dplyr::n()))$n

  clusters<-clusters[-c(which(is.na(clusters$lineage))),]
  # For each cluster, find and list the earliest collection year, the latest collection year and all the places
  # that cluster has been found
  for (i in 1:length(clusters$lineage)) {
    clusters$year_first[i] <- sequence_data %>%
      dplyr::filter(sequence_data$lineage == clusters$lineage[i])%>%
      dplyr::group_by(year)%>%
      dplyr::summarise()%>%
      min()
    clusters$year_last[i] <- sequence_data %>%
      dplyr::filter(sequence_data$lineage == clusters$lineage[i])%>%
      dplyr::group_by(year)%>%
      dplyr::summarise()%>%
      max()
    clusters$country[i]<-paste((sequence_data %>%
                                  dplyr::filter(sequence_data$lineage == clusters$lineage[i] & country !="-") %>%
                                  dplyr::group_by(country) %>%
                                  dplyr::summarise()))
  }

  return(clusters)
}
