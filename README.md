# **Am**plicon **Pre**processing **Pi**peline (AmPrePi)

## Introduction

The 16s rRNA gene amplicon is important to understand the influence of microbes in an environment.
There are many tools, scripts and manuals for processing and analysing these raw data.
In many cases, less experience users might be intimidated by these tools and scripts.
Here, we present a 16s rRNA pre-process pipeline, which uses the inference method to individual sample.
The script processes your raw files through multiple steps from sample inference, taxonomy annotation using dada2 and contamination removal using decontam.
This pipeline is based of well documented workflow.
Unlike other programs/scripts, our final product is a reproducible html report, which contains information, figures and phyloseq file that can be used for further analysis.
Some intermediate files are also available in the report.
The script is simple and designed for less experience users and at the same time is customizable for expert users.
This pipeline is ideal for people who would like to try different settings for optimization or have multiple projects


## Organization

AmPrePi is composed of 4 steps that can be run as a set or individually.

1.	Raw preprocessing
2.	Dada2 processing
3.	Taxonomic typing
4.	Contamination and final report

The scripts simplify the preprocessing step for analysing demultiplex 16S amplicon sequencings.
To run the pipeline, R (ver 4.1.0 and above) and correct packages need to be installed.
It is possible to run either R or R-studio, which will be covered down below.
The pipeline requires raw fastq and metafile in csv format.

## Requirements

### File format
Input fastq file could be zipped (gz format) or unzipped.
For reverse and forward reads, either the word R1 or R2 needs to be in the raw file.
The word fastq must also present.

### Metafile example

| Samples     | Sample_Type | Pool       | DNA_con    |
| ----------- | ----------- | ---------- | ---------- |
| sample-1    | Sample      | 1          | 8000       | 
| sample-2    | Sample      | 1          | 7000       | 
| Mock        | Mock        | 1          | 10         |
| Nc-1        | Negative    | 1          | 5000       |
| sample-3    | Sample      | 2          | 3456       |

**Samples** column contains the names of all the sample pairs in a comma or tab separated file.
This column is required.
The pipeline will use text before the first underscore as the name of the sample.
It is very important that the name matches with the raw sample so the pipeline can detect these files.
All names must be unique.
If this is not the case please rename the raw file and update your metafile.

Eg: The pair sample-1_R1.fastq.gz and sample-1_R2.fastq.gz should be named sample-1.


These columns are optional.
Please fill out completely.


**Sample_Type** descripes if the analysis file is either Sample, Mock or Negative.
Please only use one of the three words provided per sample.


**DNA_con** is the DNA concentation in ug/ml.
For the sample-free negative control (Negative), a concentration of 10 must be filled out. 
This column is required for contamination detection.  


Different contamination detection methods are scanned when both Sample_type and DNA_con are fully provided.
If either columns is available, only the appropriate contamination detection method is triggered.


**Pool** is used to differentiate between runs.
When fully filled out, the script will process samples based on the corresponding pool.

**Your_meta_data** can be added to the metafile.
Please watch out for comma as that might change the file format/structure.
Also avoid spaces, brackets and symbols.  


## Preparing to run

There are multiple ways to run the script.  
The pipeline uses R version 4.1.0.

### Version control with renv

The simples way to run is to use renv.

```
# Install renv
install.packages("renv")
# Load in renv library
library("renv")
# Restore renv image
renv::restore()
```

All the packages required will be install accordingly. 

Note: This renv.lock file was generated from a Window machine.
It might not work with Linux/Mac due to packages conflict.

### Conda Window

The environment yaml file called "envs_win.yaml" contains core programs with version control.
Navigate to the pipeline directory using conda terminal:

```
# Create environment for the script
conda env create -f envs_win.yml

# Activate conda environment
conda activate AMPREPIwin
```

Same thing can be achieved using the GUI.
Opens Anaconda Navigator > Environments > Import  
Select envs_win.yaml and wait for conda to create the environment.  
Once the environment is activated > Open Terminal > R


![alt text](https://github.com/MUMC-MEDMIC/AMPREPI/blob/a1610346773c30156a237d85c01ae872422b56e5/tutorial/import.png?raw=true)


You just created the base programs to install the rest of the packages.
Note: It is also possible to continue further using version control with renv.
Opens R and simply follow the guide above.  

You can use **BiocManager** to install the rest of the packages.

```
library("BiocManager")
BiocManager::install(c("dada2","decontam","DECIPHER","microbiome"))
```

### Conda Linux/Mac

```
# Create environment for the script
conda env create -f envs.yml
```

### Docker

A Docker image with all the github files can be downloaded.

```
# Pull the docker image
docker pull ngocgiangle/amprepi:version0

# Run image with input directory
docker run -it -v /your/directory/:/tmp ngocgiangle/amprepi:version0
```

Your files will be placed in the tmp folder.  
The pipeline and example files are located in the home directory.


## Running the pipeline

To start the pipeline, specify the location of the raw files and the metafile as well as the report location need to be specified.
The script will check the raw folder and compare the raw files found in the folder with the metafile.
Any miss match will not be processed.

From there you can edit the settings as you prefer.

### Windows

It is possible to customise the input using parameter.
Do as follow:

![alt text](https://github.com/MUMC-MEDMIC/AMPREPI/blob/564897be18e366a3aef5d7e45327c2b2661a7700/tutorial/knirt.png?raw=true)


### Linux/Mac

Run the command below from linux to generate the report
```
# Basic running command
Rscript -e "rmarkdown::render('AMPREPI.Rmd')"

# Change the output name 
Rscript -e "rmarkdown::render('AMPREPI.Rmd', output_file = 'new_file_name.html')"

# Keeps track with log file
Rscript -e "rmarkdown::render('AMPREPI.Rmd', output_file = 'new_file_name.html')" &> log.txt

```



## Citation:

If you have any questions I can be reached via:
Giang.le@mumc.nl
