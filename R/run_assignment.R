#!/usr/bin/env Rscript

rm(list=ls())

args = commandArgs(trailingOnly = T)

#'---------------------------------------------------------
#'title: Lineage Designation
#'author: Kathryn Campbell
#'date: 05/09/2021
#'---------------------------------------------------------

devtools::install_github("KathrynCampbell/MADDOG", dependencies = F)
library(MADDOG)

sequences <- seqinr::read.fasta(file = (paste(args[1], "/", args[1], ".fasta", sep = "")))

assignments<-MADDOG::assign_lineages(sequences, args[2])

write.csv(assignments, file = (paste(args[1], "/", args[1], "_assignment.csv", sep = "")))
