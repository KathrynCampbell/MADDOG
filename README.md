# Lineage Assignment

1. Overview
1. Command Line Tool
    1. How to use
        1. Setting up the environment
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

#### Setting up the environment
Open your terminal/ command line and navigate to where this repository is saved on your computer. 

If you're not familiar with Linux commands, you can navigate to the MADDOG repository in your file/finder window. Type `cd` then a space into the terminal, then drag the folder into your terminal window which will show its filepath after the cd command. If you press enter, you should enter the MADDOG folder, which you can see by the `MADDOG` name appearing in place of the `~` in your command line.

Once you're in the repository in your terminal window, type `conda env create -f environment.yml` and press enter to create the environment. This may take a few minutes. 

Once the environment is created, to activate it type `conda activate MADDOG` into the command line and press enter any time you want to use it. `MADDOG` should then appear at the start of your command line.

#### Input data requirements
The scripts require only a few things to run, but they need to be formatted in a certain way.

You'll need to have a folder with the name of your dataset (E.g. Example). It's best not to put any spaces or special characters in this. In this folder, you'll need to have a fasta file containing all the sequences. 

##### Lineage designation
For the lineage designation you'll also need the metadata for these sequences in a csv file. The metadata file must have a column called 'ID' which contains the ID's for all the sequences. These must match the sequence ID's in the fasta file. It also must have a column called 'year' that contains the collection year of each sequence. Additionally, it must have a column called 'country'. This should contain the country of origin of each sequence. If you don't have this information, you can leave the column blank but it needs to exist. If you have more detailed information about the origin of the sequence (e.g. a state) you can put this in the country column instead of the country. It also needs to have an 'assignment' column which contains any known clade assignments. If you don't have this information, you can leave this column blank but it needs to exist.

Your metadata file can have more columns and data as well, but it must contain at least the 'ID', 'year', 'country' and 'assignment' columns. Look at the metadata file in the example folder to see how this is set out. 

Please note - you can call your fasta and csv files anything, as long as they're in the dataset folder. However, the metadata file will be renamed when running the script to (dataset)_metadata.csv

##### Lineage assignment
For the lineage assignment you only need the folder with the fasta file. 

#### Running the code
To run the code, in your terminal/command line make sure you are in the MADDOG repository and the MADDOG environment (see Setting up the environment if you're unsure how to do either of these). If you've done this, you should see `MADDOG` appear at either end of your command line. 

##### Lineage designation
Type `sh lineage_designation.sh` and then press enter. You'll then be prompted to write the name of the folder that contains your data and press enter.

A lot of numbers and different processes should flash past on the screen - this means it's working! If you have a very large dataset this may take some time.

##### Lineage assignment
Type `sh assign_lineages.sh` and press enter. You'll then be prompted to write the name of the folder that contains your sequences and press enter. You'll then be prompted to write which reference set you want to use.

The assignment will then begin. This shouldn't take too long. 

#### Example
The Example_designation and Example_assignment folders show how you need to set out your dataset folder and format your input data. 

You can test out the lineage designation code by writing `sh lineage_designation.sh` When you're in the MADDOG repository and environment and then writing `Example_designation` when prompted.  

You can test out the lineage assignment code by writing `sh assign_lineages.sh` when you're in the MADDOG repository and environment and writing `Example_assignment` when prompted. You'll then need to write `Cosmo_N` when prompted.

The outputs will appear in the Example folders.

### Outputs

#### Lineage designation
The lineage designation code will produce a number of outputs in folders called 'Alignment', 'Tree', 'Treetime', 'Figures' and 'Outputs'. These will contain:

**Alignment:** The alignment of all input sequences done by the MAFFT FFT-NS-2 algorithm in fasta format.

**Tree:** All outputs of running IQTREE2 with model selection and 1000 ultrafast bootstraps on the alignment produced. The .contree file is the finalised tree.

**Treetime:** This contains the reconstructed ancestral sequences for each node of the tree in ancestral_sequences.fasta.

**Outputs:** This contains three files. the lineage_info file contains information about each lineage, including all the countries each lineage is seen in, the first and last years the lineage was sequenced, the genetic distance of each lineage and the number of sequences assigned to it from the given fasta file. It also contains a file giving information about each sequence, and information about each lineage defining node.

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
    * dplyr,
    * adephylo
    
## R Package 
This tool is also available as an R package which is capable of Lineage Designation and Assignment.
The lineage designation in the package requires the alignment, tree, ancestral reconstruction and metadata as input. 
Full details of how to use the package are available in the vignette.

To install the package, run `devtools::install_github("kathryncampbell/MADDOG", build_vignettes = TRUE, force = TRUE)` within R. 

To view the vignette use `browseVignettes("MADDOG")`

 
