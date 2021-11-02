## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----data_import, include=FALSE-----------------------------------------------
library(MADDOG)
 tree<-ape::read.tree(system.file("extdata", "Examples/Lineage_designation/Example_designation_aligned.fasta.contree", package = "MADDOG"))
 alignment<-seqinr::read.fasta(system.file("extdata", "Examples/Lineage_designation/Example_designation_aligned.fasta", package = "MADDOG"))
 alignment<-ape::as.alignment(alignment)
 metadata<-read.csv(system.file("extdata", "Examples/Lineage_designation/Example_designation_metadata.csv", package = "MADDOG"))
 ancestral<-seqinr::read.fasta(system.file("extdata", "Examples/Lineage_designation/ancestral_sequences.fasta", package = "MADDOG"))
 ancestral<-ape::as.alignment(ancestral)

## ----designation, results='hide'----------------------------------------------
sequence_designation<-seq_designation(tree, 90, alignment, metadata, ancestral)
defining_node_information<-node_info(tree, 90, alignment, metadata, ancestral)


## ----designation_outputs------------------------------------------------------
head(sequence_designation, 20)
defining_node_information

## ----lineage_info-------------------------------------------------------------
lineage_info<-lineage_info(sequence_designation, metadata) ; lineage_info

## ----data_import2, include = FALSE--------------------------------------------
sequences<-seqinr::read.fasta(system.file("extdata", "Examples/Lineage_assignment/example.fasta", package = "MADDOG"))

## ----assignment---------------------------------------------------------------
assignments<-assign_lineages(sequences, "Cosmo_N"); assignments

