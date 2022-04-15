#! /bin/bash
set -e

#command line inputs
echo "What is the name of the folder containing your sequences and metadata?"
read runname

#make sure the metadata is named to work in future script
mv -vn $runname/*.csv $runname/$runname"_metadata.csv"

mv -vn $runname/*.fasta $runname/$runname".fasta"

mafft --add $runname/$runname".fasta" --reorder inst/extdata/References/RABV/reference_aligned.fasta > $runname/$runname"_withref.fasta"

#Lineage assignment
Rscript Run/windows_run_assignment.R $runname

rm $runname/$runname"_withref.fasta"

mkdir $runname/Assignment

cp $runname/$runname"_assignment.csv" $runname/Assignment/assignment.csv

Rscript Run/additions2.R $runname

#alignment
mkdir $runname/Alignment
mafft $runname/$runname"_combined.fasta" > $runname/$runname"_combined_aligned.fasta"
mv $runname/$runname"_combined_aligned.fasta" $runname/Alignment
#mafft sometimes messes names up so they don't match the tree, in case this happens have this line (though will need to fix later):
awk '{sub(/,/,"_")}1' $runname/Alignment/$runname"_cmbined_aligned.fasta" > temp.fasta && mv temp.fasta $runname/Alignment/$runname"_combined_aligned.fasta"

#IQTree
mkdir $runname/Trees
iqtree -s $runname/Alignment/$runname"_combined_aligned.fasta" -B 1000
mv $runname/Alignment/$runname"_combined_aligned.fasta."* $runname/Trees

#Timetree
mkdir $runname/Timetree
treetime ancestral --aln $runname/Alignment/$runname"_combined_aligned.fasta" --tree $runname/Trees/$runname"_combined_aligned.fasta.contree" --outdir $runname/Timetree

mkdir $runname/Outputs
mkdir $runname/Figures
Rscript Run/additions3.R $runname

rm $runname/$runname"_assignment.csv"
rm $runname/$runname"_combined.fasta"
rm $runname/$runname"_assignments.csv"
rm -r $runname/Figures/$runname"_sunburst_files"
rm temp.fasta
