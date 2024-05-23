	#! /bin/bash
set -e

#Initial checks
echo "It is essential to pull the MADDOG repository each time before starting a run to ensure any new and updated lineages are included. Have you pulled the repository? Type Y or N"
read answer1

echo "What is the name of the folder containing your sequences and metadata?"
read runname

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