#!/bin/bash

set -e

#command line inputs
echo "What is the name of the folder containing your sequences?"
read runname

mv -vn $runname/*.fasta $runname/$runname".fasta"

mafft --add $runname/$runname".fasta" --reorder inst/extdata/References/RABV/reference_aligned.fasta > $runname/$runname"_withref.fasta"

#Lineage assignment
Rscript Run/windows_run_assignment.R $runname

rm $runname/$runname"_withref.fasta"
