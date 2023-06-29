#!/usr/bin/env Rscript

rm(list=ls())

args = commandArgs(trailingOnly = T)

#'---------------------------------------------------------
#'title: Lineage Designation
#'author: Kathryn Campbell
#'date: 05/09/2021
#'---------------------------------------------------------

source("R/seq_designation.R")
source("R/node_info.R")
source("R/lineage_info.R")
source("R/sunburst_function.r")
source("R/lineage_tree.r")
library(dplyr)
library(ggtree)

devtools::install_version('rvcheck',repos = "http://cran.us.r-project.org", version='0.1.8')


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
tree <- ape::read.tree(file = paste(args, "/Trees/", args, "_aligned.fasta.contree", sep = ""))
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
metadata <- read.csv(file = paste(args, "/", args, "_metadata.csv", sep = ""))

#'**ALIGNMENT**
#'========================================================================================================
#' The alignment must contain the element 'seq' which contains the sequences
#' And the element 'nam' which lists all the sequence ID's
#' These sequence ID's must match the sequence ID's in the metadata and tree
#'=========================================================================================================
alignment <- seqinr::read.alignment(file = (paste(args, "/Alignment/", args, "_aligned.fasta", sep = "")), format = "fasta")

# Sequence names got messed up in MAFFT, need to fix these so they match metadata and alignment
# #KB- can replace above 2 lines with this:
alignment$nam <- gsub("\\..*", "", alignment$nam, perl = T)


#'**TIMETREE**
#'========================================================================================================
#'
#'=========================================================================================================
ancestral <- seqinr::read.alignment(file = (paste(args, "/Timetree/ancestral_sequences.fasta", sep = "")), format = "fasta")
ancestral$nam <- gsub("\\..*", "", ancestral$nam, perl = T)

#############################################
#           RUN DESIGNATION                #
#############################################
sequence_designation<-seq_designation(tree, 90, alignment, metadata, ancestral)
defining_node_information<-node_info(tree, 90, alignment, metadata, ancestral)
lineage_info<-lineage_info(sequence_designation, metadata)

write.csv(sequence_designation, file = (paste(args, "/Outputs/", args, "_sequence_data.csv", sep = "")), row.names=F)
write.csv(defining_node_information, file = (paste(args, "/Outputs/", args, "_node_data.csv", sep = "")), row.names=F)

write.csv(lineage_info, file = (paste(args, "/Outputs/", args, "_lineage_info.csv", sep = "")), row.names=F)

#############################################
#               FIGURES                     #
#############################################

new<-sunburst(lineage_info, defining_node_information, tree, metadata, sequence_designation)

htmlwidgets::saveWidget(plotly::as_widget(new), (paste(args, "/Figures/", args, "_sunburst.html", sep = "")))

plot_tree<-lineage_tree(lineage_info, defining_node_information, tree, metadata, sequence_designation)

ggplot2::ggsave(paste(args, "/Figures/", args, "_lineage_tree.png", sep = ""),
       plot = plot_tree)


