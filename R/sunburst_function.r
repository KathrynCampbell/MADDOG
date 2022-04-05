#' Lineage Figures - Sunburst
#'
#' This group of Lineage Figure functions produce figures to assist in interpretation of lineage data.
#'
#' This function produces a sunburst plot to visualise the hierarchical relationship of the lineages.
#' The outputs of the seq_designation, node_info and lineage_info functions are required, along with the
#' phylogenetic tree and corresponding metadata file used as input for the sequence data, node info and lineage
#' information.
#'
#' @param lineage_info The output of the lineage_info function
#' @param node_data The output of the node_info function
#' @param tree A phylogenetic tree
#' @param metadata The metadata corresponding to the sequences in the tree, including "ID" "assignment" "country" and "year"
#' @param sequence_data The output of the seq_designation function
#' @return A sunburst plot showing the hierarchal relationship of the lineages
#' @export
sunburst <- function(lineage_info, node_data, tree, metadata, sequence_data) {
  tree$tip.label <- gsub("\\..*", "", tree$tip.label, perl = T)
  tree$node.comment<- gsub(".*=", "", tree$node.label, perl = T)
  previous<-data.frame(assignment = unique(metadata$assignment), parent = "", n_seqs = NA)
  previous$parent[1]<-""
  for (i in 1:length(previous$assignment)) {
    previous$n_seqs[i]<-length(which(metadata$alignment.name == previous$assignment[i]))
  }
  node_data<-node_data[order(node_data$lineage),]

  node_data$parent<-NA
  node_data$parent[1]<-""


  for (i in 2:length(node_data$node)) {
    if (length(which(node_data$node %in% treeio::ancestor(tree, node_data$node[i]))) == 0) {
      node_data$parent[i]<-""
    } else {
      parent<-node_data$lineage[which(node_data$node %in% treeio::ancestor(tree, node_data$node[i]))]
      node_data$parent[i]<-parent[length(parent)]
    }
  }

  lineage_info$parent<-NA

  for (i in 1:length(lineage_info$lineage)) {
    lineage_info$parent[i]<-node_data$parent[which(node_data$lineage == lineage_info$lineage[i])]

  }

  lineage_info$colour<-NA

  Colours<-c("Reds","Purples","YlOrBr","PuBuGn","YlOrRd","OrRd","PuBu","Pastel1","Greens","Greys",
             "GnBu","BuGn","RdPu","Oranges","BuPu","YlGn","PuRd","YlGnBu")

  lineages<-data.frame(lineage = lineage_info$lineage, subclade = NA)

  for (i in 1:length(lineages$lineage)) {
    lineages$subclade[i]<-strsplit(lineages$lineage[i], "_")[[1]][1]
  }

  letters <- c("A1", "B1", "C1", "D1", "E1", "F1", "G1", "H1", "I1", "J1", "K1", "L1", "M1", "N1",
               "O1", "P1", "Q1", "R1", "S1", "T1", "U1", "V1", "W1", "X1", "Y1", "Z1")

  if(length(grep("_", lineage_info$lineage)) != 0) {
    if (length(which(lineages$subclade %in% letters)) != 0) {
      lineages<-lineages[-c(which(lineages$subclade %in% letters)),]
    }
  }

  clades<-unique(lineages$subclade)

  if(length(grep("\\.", clades)) != 0 ) {
    clades<-clades[-c(grep("\\.", clades))]
  }

  lineage<-lineage_info$lineage[-c(grep("_", lineage_info$lineage))]
  cols<-RColorBrewer::brewer.pal(9, "Blues")
  pal<-colorRampPalette(c(cols))
  pal<-rev(pal(length(lineage)))
  lineage_info$colour[-c(grep("_", lineage_info$lineage))]<-pal

  for (i in 1:length(clades)) {
    lineage<-grep(clades[i], lineage_info$lineage)
    cols<-RColorBrewer::brewer.pal(3, Colours[i])
    pal<-colorRampPalette(c(cols))
    pal<-rev(pal(length(lineage)))
    lineage_info$colour[(grep(clades[i], lineage_info$lineage))]<-pal
  }

  new<-plotly::plot_ly(
    labels = c(lineage_info$lineage),
    parents = c(lineage_info$parent),
    values = c(lineage_info$n_seqs),
    type = "sunburst",
    marker = list(colors = (lineage_info$colour))
  )

  return(new)
}


