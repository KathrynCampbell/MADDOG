#! /bin/bash
set -e

#command line inputs
echo "What is the name of the folder containing your sequences and metadata?"
read runname

#make sure the metadata is named to work in future script
mv -vn $runname/*.csv $runname/$runname"_metadata.csv"


#Timetree
mkdir $runname/Timetree
treetime ancestral --aln $runname/Alignment/$runname"_aligned.fasta" --tree $runname/Trees/$runname"_aligned.fasta.contree" --outdir $runname/Timetree

#Lineage assignment
mkdir $runname/Outputs
Rscript Run/run_designation_slow.R $runname
