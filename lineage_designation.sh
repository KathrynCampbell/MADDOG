#! /bin/bash
set -e

#command line inputs
echo "What is the name of the folder containing your sequences and metadata?"
read runname

#make sure the metadata is named to work in future script
mv $runname/*.csv $runname/$runname"_metadata.csv"

#alignment
mkdir $runname/Alignment
mafft $runname/*.fasta > $runname/$runname"_aligned.fasta"
mv $runname/$runname"_aligned.fasta" $runname/Alignment
#mafft sometimes messes names up so they don't match the tree, in case this happens have this line (though will need to fix later):
awk '{sub(/,/,"_")}1' $runname/Alignment/$runname"_aligned.fasta" > temp.fasta && mv temp.fasta $runname/Alignment/$runname"_aligned.fasta"

#IQTree
mkdir $runname/Trees
iqtree -s $runname/Alignment/$runname"_aligned.fasta" -B 1000
mv $runname/Alignment/$runname"_aligned.fasta."* $runname/Trees

#Timetree
mkdir $runname/Timetree
treetime ancestral --aln $runname/Alignment/$runname"_aligned.fasta" --tree $runname/Trees/$runname"_aligned.fasta.contree" --outdir $runname/Timetree

#Lineage assignment
mkdir $runname/Outputs
Rscript R/run_designation.R $runname
