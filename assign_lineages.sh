#!/usr/bin/env bash

set -e

#command line inputs
echo "What is the name of the folder containing your sequences?"
read runname
echo "What reference would you like to use? Please type Cosmo_N or Cosmo_WGS"
read reference

mv $runname/*.fasta $runname/$runname".fasta"

#Lineage assignment
Rscript R/run_assignment.R $runname $reference
