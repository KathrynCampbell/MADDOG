#!/bin/bash

set -e

#command line inputs
echo "What is the name of the folder containing your sequences?"
read runname
echo "What whole genome set are you adding to? Type one of the following options: Africa Asian Arctic Bat Cosmopolitan RAC-SK Vac"
read reference

#make sure the metadata and sequences are named to work in future script
mv -vn $runname/*.fasta $runname/$runname".fasta"
mv -vn $runname/*.csv $runname/$runname"_metadata.csv"

#alignment
mkdir $runname/Alignment
mafft $runname/*.fasta > $runname/$runname"_aligned.fasta"
mv $runname/$runname"_aligned.fasta" $runname/Alignment
#mafft sometimes messes names up so they don't match the tree, in case this happens have this line (though will need to fix later):
awk '{sub(/,/,"_")}1' $runname/Alignment/$runname"_aligned.fasta" > temp.fasta && mv temp.fasta $runname/Alignment/$runname"_aligned.fasta"

#IQTree (ultrafast 1000)
mkdir $runname/Trees
iqtree -s $runname/Alignment/$runname"_aligned.fasta" -B 1000
mv $runname/Alignment/$runname"_aligned.fasta."* $runname/Trees

#Timetree
mkdir $runname/Timetree
treetime ancestral --aln $runname/Alignment/$runname"_aligned.fasta" --tree $runname/Trees/$runname"_aligned.fasta.contree" --outdir $runname/Timetree

#Make directories for R
mkdir $runname/assignment
mkdir $runname/temp_designation
mkdir $runname/reference
mkdir $runname/Outputs

#Export sequences from WGS subclades for windows_run_assignment
#Designation for non-WGS subclades
Rscript Run/N_WGS_comparable.R $runname $reference
Rscript Run/reference_combined.R $runname $reference

#Realign reference
mafft inst/extdata/References/RABV/seq.fasta > inst/extdata/References/RABV/reference_aligned.fasta

#Assign WGS lineage sequences
Rscript Run/run_assignment_partial.R $runname $reference

#Updates from assignments
Rscript Run/updates.R $runname $reference

#Add new designations to reference and realign
Rscript Run/reference_combined.R $runname $reference
mafft inst/extdata/References/RABV/seq.fasta > inst/extdata/References/RABV/reference_aligned.fasta

rm -r $runname/assignment
rm -r $runname/reference
rm -r $runname/temp_designation
