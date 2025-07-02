# Mycobacteria Amplicon Sequencing Tool (MAST)

The Mycobacteria Amplicon Sequencing Tool (MAST) is a worklow made with nextflow used for Mycobacteria AMR prediction for amplicon sequences. The results are formatted into a customizable patient report, where both patient information and drug resistances are indicated. 

#Installation

This workflow was created and tested on macOS 14.5 (Sonoma). 

The environment for the workflow is installed using [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). Please make sure that Conda is installed. The environment will be set up, and located in the /work/conda folder when the pipeline is ran. 

The environment set up has been tested on Conda 23.7.4. 

To clone the repository, please run the code below in your root folder. At this point, the workflow assumes you are in your home directory. 

```
git clone https://github.com/guthrielab/MAST
```

##Running MAST

To run it, please specify the fastq file that is to be analyzed, and the output directory for the results.

```
nextflow run main.nf --data data.fq --outdir outdir
```
