# MAD DOG: Method for Assignment, Definition and Designation Of Global Lineages

[![DOI](https://zenodo.org/badge/388405372.svg)](https://zenodo.org/badge/latestdoi/388405372)

1. Overview
1. Command Line Tool
    1. How to use
        1. Setting up the environment - Unix/MAC
        1. Setting up the environment - Windows
        1. Input data requirements
            1. Lineage designation
            1. Lineage assignment
        1. Running the code
            1. Lineage designation
            1. Lineage assignment
        1. Example
    1. Outputs
        1. Lineage designation
        1. Lineage assignment
    1. Software used
1. R Package

## Overview
The ‘Zero by 30’ global strategy, led by the WHO and partners, aims to eliminate human deaths from dog-mediated rabies by 2030. This goal requires improved and coordinated surveillance to understand the current epidemiological situation and monitor the impact of control measures. The existing system for classifying rabies virus diversity defines six global clades that circulate in domestic dog populations, broadly according to their geographic distribution. While useful for coarse taxonomic classification, this system lacks the definition needed for monitoring rabies viruses circulating on a more local scale. Here, we present an updated lineage designation and assignment tool for rabies virus based on the dynamic nomenclature used for SARS-CoV-2 by Rambaut et al. (2020). Application of this tool can be used to generate detailed information to inform control efforts and monitor progress towards the elimination of rabies virus.

## Command Line Tool

### How to use
Make sure you have this GitHub Repository cloned and saved somewhere easy to navigate to on your computer.

You'll also need to have miniconda installed. You can download this here: https://docs.conda.io/en/latest/miniconda.html 

#### Setting up the environment - Unix/MAC
Open your terminal/ command line and navigate to where this repository is saved on your computer. 

If you're not familiar with Linux commands, you can navigate to the MADDOG repository in your file/finder window. Type `cd` then a space into the terminal, then drag the folder into your terminal window which will show its filepath after the cd command. If you press enter, you should enter the MADDOG folder, which you can see by the `MADDOG` name appearing in place of the `~` in your command line.

Once you're in the repository in your terminal window, type `conda env create -f environment.yml` and press enter to create the environment. This may take a few minutes. 

Once the environment is created, to activate it type `conda activate MADDOG` into the command line and press enter any time you want to use it. `MADDOG` should then appear at the start of your command line.

#### Setting up the environment - Windows
Turn on Windows Subsystem for Linux (https://mafft.cbrc.jp/alignment/software/windowsfeatures.html)

Download Ubuntu (https://www.microsoft.com/en-gb/p/ubuntu/9nblggh4msv6?rtc=1&activetab=pivot:overviewtab)

Within Ubuntu run:
curl -sL \ "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" > "Miniconda3.sh"

Then:

`bash Miniconda3.sh`

You can then remove the installer by:

`rm Miniconda3.sh`

Finally, update conda using:

`conda update conda`

If this doesn't work, and gives a http error, you may need to install Win64 OpenSSL v1.1.1L (https://slproweb.com/products/Win32OpenSSL.html)

Then navigate to where the repository is stored on your computer. This is more difficult using this system, and may take some trial and error! To start, use `cd ~/../../mnt/c`

You can then use `ls` to view the contents of the directory, and `cd` followed by the next directory name to change directory. 

Once you've navigated to the MADDOG repo, type `conda env create -f environment.yml` and press enter to create the environment. This may take a few minutes. 

Once the environment is created, to activate it type `conda activate MADDOG` into the command line and press enter any time you want to use it. `MADDOG` should then appear at the start of your command line.

You'll also need to convert some file types so windows can use them.

Write `conda install -c conda-forge dos2unix` 

Then `dos2unix assignment.sh` and next `dos2unix designation.sh` followed by `dos2unix new_designation.sh`

#### Input data requirements
The scripts require only a few things to run, but they need to be formatted in a certain way.

You'll need to have a folder with the name of your dataset (E.g. Example_designation). It's best not to put any spaces or special characters in this. In this folder, you'll need to have a fasta file containing all the sequences. 

##### Lineage designation
For the lineage designation you'll also need the metadata for these sequences in a csv file. The metadata file must have a column called 'ID' which contains the ID's for all the sequences. These must match the sequence ID's in the fasta file. It also must have a column called 'year' that contains the collection year of each sequence. Additionally, it must have a column called 'country'. This should contain the country of origin of each sequence. If you don't have this information, you can leave the column blank but it needs to exist. If you have more detailed information about the origin of the sequence (e.g. a state) you can put this in the country column instead of the country. It also needs to have an 'assignment' column which contains any known clade assignments. If you don't have this information, you can leave this column blank but it needs to exist.

Your metadata file must have only the 'ID', 'year', 'country' and 'assignment' columns. Look at the metadata file in the example folder to see how this is set out. 

Please note - you can call your fasta and csv files anything, as long as they're in the dataset folder. However, the metadata file will be renamed when running the script to (dataset)_metadata.csv

##### Lineage assignment
For the lineage assignment you only need the folder with the fasta file. 

#### Running the code
To run the code, in your terminal/command line make sure you are in the MADDOG repository and the MADDOG environment (see Setting up the environment if you're unsure how to do either of these). If you've done this, you should see `MADDOG` appear at either end of your command line. 

##### Lineage designation
Type `sh designation.sh` and then press enter. You'll then be prompted to write the name of the folder that contains your data and press enter.

A lot of numbers and different processes should flash past on the screen - this means it's working! If you have a very large dataset this may take some time.

Lineages are designated using a reference set of designations as a starting point. Your sequences will first be assigned a lineage based on the reference set, and then checked to see if the inclusion of your sequences results in new lineages. If you wish to run an entirely new designation, not including the reference sequences, you can run `sh new_designation.sh` which will produce designation relevant only to your study and not comparable to the global set.

##### Lineage assignment
Type `sh assignment.sh` and press enter. You'll then be prompted to write the name of the folder that contains your sequences and press enter. 

The assignment will then begin. This shouldn't take too long. 

#### Example
The Example_designation and Example_assignment folders show how you need to set out your dataset folder and format your input data. 

You can test out the lineage designation code by writing `sh designation.sh` When you're in the MADDOG repository and environment and then writing `Example_designation` when prompted.  

You can test out the lineage assignment code by writing `sh assignment.sh` when you're in the MADDOG repository and environment and writing `Example_assignment` when prompted. 

The outputs will appear in the Example folders.

### Outputs

#### Lineage designation
Lineages are designated using a reference set of designations as a starting point. Your sequences will first be assigned a lineage based on the reference set, and then checked to see if the inclusion of your sequences results in new lineages. If so, these new lineages will be designated. 
**PLEASE NOTE** If new lineages are designated, these are not added to the reference set and therefore the resultant designations apply to your study only.

The lineage designation code will produce a number of outputs in folders called 'Alignment', 'Tree', 'Treetime', 'Figures' and 'Outputs'. These will contain:

**Alignment:** The assignments of all your sequences against the reference set. 

**Alignment:** The alignment of all input sequences plus sequences from all relevant existing lineages done by the MAFFT FFT-NS-2 algorithm in fasta format.

**Tree:** All outputs of running IQTREE2 with model selection and 1000 ultrafast bootstraps on the alignment produced. The .contree file is the finalised tree. The tree contains your input sequences plus all sequences from relevant lineages.

**Treetime:** This contains the reconstructed ancestral sequences for each node of the tree in ancestral_sequences.fasta.

**Outputs:** This contains 1-5 files.

The relevant_lineages file will always be included. This gives details about all the lineages the sequences are assigned to and the number of sequences assigned to each lineage. 

If your sequences result in new lineages being designated, a new_lineages file will be present. This gives details such as the countries the lineage has been seen in, the first and most recent years the lineage has been sequenced, and the number of sequences present in the lineage. 

If you have a new_lineages file, you will also have a sequence_data file which details the lineage each of your sequences has been designated. If there are no new lineages, this file won't exist and the file in the Assignments folder should be used instead.

The lineage designation process also checks fo lineages that are undersampled or emerging (all lineage defining requirements met, but only 5-9 sequences within a 5 year period). If any of these are detected within your sequence set and relevant lineages, a undersampled_emerging file will be present detailing the potentially emerging/undersampled lineages.

The lineage designation process also identifies singletons of interest that may indicate sequencing errors, very undersampled lineages or very newly emerging lineages. If any of these are present, a singletons_of_interest file will be produced, detailing the number of singletons of interest for each lineage, and where and when these singleton sequences are from.

**Figures:**This contains two files to help interpret the lineage information. 

First is a html file of a sunburst plot. This shows the hierarchal relationships of all relevant lineages (new and existing). The html file is interactive, and sections can be collapsed/ zoomed in as required. This also provides the colour scheme for all the figures. 

Second is a lineage tree with 2 panels. The left panel is a phylogenetic tree of the input sequences plus all sequences from relevant lineages, with the tips coloured by lineage (colour scheme from sunburst plot) with a lineage bar along side to easily show the positions of the lineages within the tree. The right hand panel indicates which of the tips represent your input sequences; showing how these fit in to the existing phylogeny.

A csv file detailing the metadata of all the sequences from relevant lineages included during the designation process is also produced.

#### Lineage assignment
The lineage assignment code produces a csv file detailing the assignment for each sequence, and details about the lineage assigned including the countries it has been seen it and when that lineage was first and last sequenced.

### Software used
These are all installed as part of the environment and don't need to be installed separately.

* MAFFT (https://mafft.cbrc.jp/alignment/software/)
* IQTREE (http://www.iqtree.org)
* Treetime (https://github.com/neherlab/treetime)
* R (https://www.r-project.org)
    * seqinr
    * phangorn,
    * stringr,
    * caper,
    * ips,
    * adegenet,
    * ape,
    * rgdal,
    * dplyr,
    * adephylo
    * gridExtra
    * phytools
    * castor
    * ggtree
    * treeio
    
## R Package 
This tool is also available as an R package which is capable of Lineage Designation and Assignment.
The lineage designation in the package requires the alignment, tree, ancestral reconstruction and metadata as input. 
**PLEASE NOTE** The lineage designation done in the package currently is equivalent to `new_designation.sh` - it produces a set of designation but these are not based on the reference set. It is instead used to give an idea of the diversity within your set of sequences.  

Full details of how to use the package are available in the vignette.

To install the package, run `devtools::install_github("kathryncampbell/MADDOG", build_vignettes = TRUE, force = TRUE)` within R. 

To view the vignette use `browseVignettes("MADDOG")`

**Be aware; the assign_lineages function currently is unavailable for Windows, but the command line tool can be used instead.**

 
