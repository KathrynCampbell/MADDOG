#! /bin/bash
set -e

#Initial checks
echo "It is essential to pull the MADDOG repository each time before starting a run to ensure any new and updated lineages are included. Have you pulled the repository? Type Y or N"
read answer1

if [[ "$answer1" == "Y" ]]
then

	#command line inputs
	echo "What is the name of the folder containing your sequences and metadata?"
	read runname

	#make sure the metadata  is named to work in future script
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

	#IQTree
	mkdir $runname/Trees
	iqtree -s $runname/Alignment/$runname"_combined_aligned.fasta" -B 1000 -m GTR+F+R6 -nt AUTO
	mv $runname/Alignment/$runname"_combined_aligned.fasta."* $runname/Trees

	#Timetree
	mkdir $runname/Timetree
	treetime ancestral --aln $runname/Alignment/$runname"_combined_aligned.fasta" --tree $runname/Trees/$runname"_combined_aligned.fasta.contree" --outdir $runname/Timetree

	mkdir $runname/Outputs
	mkdir $runname/Figures
	Rscript Run/additions3.1.R $runname


	# Reporting
	
	Rscript Run/run_markdown.R $runname 
	mkdir $runname/Report
	mv "report.html" $runname/Report/"report.html"
	rm "report.md"
	mv "report_files" $runname/Report/"report_files"
	cp troupin.png $runname/Report/

	if [ -f $runname/"Outputs/new_lineages.csv" ]; then
	    echo "New lineages have been found within your dataset! In order for these lineages to be confirmed they need to be added to the MADDOG reference.
	    Would you like to do this? NOTE: This does not add the lineages right now but will generate a file called NEXT_STEPS.eml with further instructions.
	    Please type Y or N"
	    read answer2
	    if [[ "$answer2" == "Y" ]]
	    then 
	    	rm -r $runname/Figures/$runname"_sunburst_files"
	    	zip -r $runname/$runname".zip" $runname
	    	cp inst/extdata/updates.eml $runname/"NEXT_STEPS.eml"
	    fi
	fi

	rm $runname/$runname"_assignment.csv"
	rm $runname/$runname"_combined.fasta"
	rm $runname/$runname"_assignments.csv"
	rm -r $runname/Figures/$runname"_sunburst_files"
	rm Rplots.pdf


else
	echo "Please update MADDOG before starting a run!"

fi
