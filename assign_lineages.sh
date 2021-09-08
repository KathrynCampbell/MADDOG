#!/bin/bash

set -e

#command line inputs
echo "What is the name of the folder containing your sequences?"
read runname
echo "What reference would you like to use? Please type Cosmo_N or Cosmo_WGS"
read reference
echo "What operating system are you using? Write unix or windows"
read OS

mv -vn $runname/*.fasta $runname/$runname".fasta"

#Lineage assignment
Rscript Run/run_assignment.R $runname $reference $OS
