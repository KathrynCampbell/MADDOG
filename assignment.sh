#!/bin/bash

set -e

#Initial checks
echo "It is essential to pull the MADDOG repository each time before starting a run to ensure any new and updated lineages are included. Have you pulled the repository? Type Y or N"
read answer1

if [[ "$answer1" == "Y" ]]
then

	#command line inputs
	echo "What is the name of the folder containing your sequences?"
	read runname

	mv -vn $runname/*.fasta $runname/$runname".fasta"

	mafft --add $runname/$runname".fasta" --reorder inst/extdata/References/RABV/reference_aligned.fasta > $runname/$runname"_withref.fasta"

	#Lineage assignment
	Rscript Run/windows_run_assignment.R $runname

	rm $runname/$runname"_withref.fasta"

else
	echo "Please update MADDOG before starting a run!"

fi

