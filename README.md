 # Mycobacteria Amplicon Sequencing Tool (MAST)

The Mycobacteria Amplicon Sequencing Tool (MAST) is a worklow made with nextflow used for Mycobacteria AMR prediction for amplicon sequences. The results are formatted into a customizable patient report, where both patient information and drug resistances are indicated. 

## Installation

This workflow was created and tested on macOS 14.5 (Sonoma). 

The environment set up has been tested on Conda 23.7.4. 

To clone the repository, please run the code below in your root folder. At this point, the workflow assumes you are in your home directory. 

```
git clone https://github.com/guthrielab/MAST
```

## Dependencies

This project uses a Conda environment to manage all dependencies.

**Packages:**
- `python`
- `pandas`
- `biopython`
- `python-docx`
- `jinja2`
- `bcftools`
- `samtools`
- `ivar`
- `bwa`
- `freebayes`
- `cutadapt`
- `filtlong`
- `bedtools`
- `gsl`
- `seqkit`

The environment for the workflow is installed using [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). Please make sure that Conda is installed. The environment is set up, and located in the /work/conda folder when the pipeline is ran. 


## Running MAST

To run it, please specify the folder with fastq files that is to be analyzed, and the output directory for the results.

```
nextflow run -c nextflow.config main.nf \
  --data fastq_folder \
  --outdir results \
  --reference reference_H37RV.fasta \
  --primers tb-amplicon-primers.bed \
  --compare_mutations compare_mutations.py 
```

MAST accepts single-end FASTQ files as input. The pipeline is optimized for amplicon sequencing data, applying appropriate depth thresholds for filtering. Once the reads are cleaned and aligned, MAST trims primers and detects variants. These variants are then cross-referenced with the [WHO](https://www.who.int/publications/i/item/9789240082410) catalogue of mutations, and compiled into a report. 

The `--outdir` specifies the location of the final report. The report will be titled after the barcode of the fastq file. 

## Customizing the Report

The contents of the final report can be customized using the `patient.csv` file, located in the `/Data` folder.

> **Important:** Ensure that the name of the FASTQ file matches the `barcode` column entry in the `patient.csv` file.  
> If the names do not match, the pipeline will raise an error.

## Caching and Resuming

Each process run by MAST is cached in the `/work` directory.

If an error occurs due to issues with the input, you do not need to re-run the entire pipeline from scratch. Instead, use the `-resume` flag when re-running the pipeline:

```bash
nextflow run -c nextflow.config main.nf \
  --data fastq_folder \
  --outdir results \
  --reference reference_H37RV.fasta \
  --primers tb-amplicon-primers.bed \
  --compare_mutations compare_mutations.py \
  -resume
