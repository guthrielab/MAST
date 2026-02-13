[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17460753.svg)](https://doi.org/10.5281/zenodo.17460753)
 # Mycobacteria Amplicon Sequencing Tool (MAST)

The Mycobacteria Amplicon Sequencing Tool (MAST) is a worklow made with nextflow used for Mycobacteria AMR prediction for amplicon sequences. The results are formatted into a customizable patient report, where both patient information and drug resistances are indicated. 

## Installation

This workflow was created and tested on macOS 14.5 (Sonoma) and Linux Ubuntu 22.04.4 LTS. 

The environment set up has been tested on Conda 23.7.4. 

Please ensure you have [Nextflow](https://www.nextflow.io/docs/latest/install.html) on your machine before cloning the repository. Nextflow requires Java version 17 or later (up to 25) to run. Steps to install Nextflow and the required Java version is provided in the Installation page for Nextflow.

To clone the repository, please run the code below. 

```
git clone https://github.com/guthrielab/MAST
```

## Dependencies

This project uses a Conda environment to manage all dependencies.

**Packages:**
  - `python v3.10`
  - `pandas v2.3.3`
  - `biopython v1.86`
  - `python-docx v1.2.0`
  - `jinja2 v3.1.6`
  - `bcftools v1.17`
  - `samtools v1.21`
  - `ivar v1.4.3`
  - `bwa v0.7.18`
  - `freebayes v1.3.6`
  - `cutadapt v5.2`
  - `filtlong v0.2.1`
  - `bedtools v2.31.1`
  - `gsl v2.7`
  - `seqkit v2.12.0`
  - `pysam v0.22.1`

The environment for the workflow is installed using [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). Please make sure that Conda is installed. The environment is set up, and located in the /work/conda folder when the pipeline is ran. 


## Running MAST

Once you clone the MAST repository, you can launch MAST from the parent directory where it was cloned by using the commands below.

To display help:
```nextflow run MAST --help```

## MAST Tutorial

1. Download the test data from the [Zenodo link](https://doi.org/10.5281/zenodo.17460753) into the parent directory where MAST is
2. Extract the zip file and you should have a directory named `17460753` with **5 fastq.gz** files and a **patient_info.csv** file
3. Replace the patient_info.csv file in the MAST/Data/ directory with the one in the `17460753` folder that you have just downloaded
4. From the parent directory of MAST, run: `nextflow run MAST --input 17460753 --outdir tutorial_output`
5. You should have the following prompt on your screen as the pipeline analyzes your test data
   ```
   N E X T F L O W   ~  version 24.04.3

   Launching `MAST/main.nf` [angry_mahavira] DSL2 - revision: 5aa6077172

   executor >  local (40)
   [04/df42cf] runQualityTrimming (1) [100%] 5 of 5 ✔
   [8f/a0ff96] runAlignment (4)       [100%] 5 of 5 ✔
   [44/f9aa5b] runSortAndIndex (5)    [100%] 5 of 5 ✔
   [e4/d72d0a] runPrimerTrimming (5)  [100%] 5 of 5 ✔
   [12/1cab06] runVariantCalling (5)  [100%] 5 of 5 ✔
   [73/50e189] runFilterVariants (5)  [100%] 5 of 5 ✔
   [65/87e21f] runConvertToTSV (5)    [100%] 5 of 5 ✔
   [15/8a08d2] compareMutations (5)   [100%] 5 of 5 ✔
   Completed at: 11-Feb-2026 16:47:37
   Duration    : 8m 45s
   CPU hours   : 1.2
   Succeeded   : 40

Your tutorial_output directory should have 10 files in it, 5 DOCX files and 5 TSV files for each input file


## MAST Usage

To run MAST on your own dataset, please specify the folder with fastq files that is to be analyzed, and the output directory for the results.
> **Important:** Ensure that the name of the FASTQ file matches the `barcode` column entry in the `patient_info.csv` file.  
> If the names do not match, the pipeline will raise an error.

The contents of the final report can be customized using the `patient.csv` file, located in the `MAST/Data` folder.

```
nextflow run MAST \
  --input {fastq_folder} \
  --outdir {results_directory}
```

MAST accepts single-end FASTQ files as input. The pipeline is optimized for amplicon sequencing data, applying appropriate depth thresholds for filtering. Once the reads are cleaned and aligned, MAST trims primers and detects variants. These variants are then cross-referenced with the [WHO](https://www.who.int/publications/i/item/9789240082410) catalogue of mutations, and compiled into a report. 

The `--outdir` specifies the location of the final report. The report will be titled after the barcode of the fastq file. 

## Workflow
![MAST_workflow](https://github.com/guthrielab/MAST/blob/main/Data/MAST_workflow_v3.png)


## Caching and Resuming

Each process run by MAST is cached in the `/work` directory.

It is important to note that this folder contains intermediary files and can become large with time, so it is advisable to clean it up or delete the subdirectories when you have completed your analysis to avoid storage problems.

If an error occurs due to issues with the input, you do not need to re-run the entire pipeline from scratch. Instead, use the `-resume` flag when re-running the pipeline:

```
nextflow run MAST \
  --input {fastq_folder} \
  --outdir {results_directory} \
  -resume
```
