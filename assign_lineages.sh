#!/bin/bash

set -e

#command line inputs
echo "What is the name of the folder containing your sequences?"
read runname

reference=RABV

mv -vn $runname/*.fasta $runname/$runname".fasta"

#Lineage assignment
Rscript Run/run_assignment.R $runname $reference $OS
