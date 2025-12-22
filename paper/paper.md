---
title: 'MAST: Mycobacteria Amplicon Sequencing Tool'
tags:
- amplicon sequencing
- pathogen genomics
- bioinformatics
- antimicrobial resistance
- public health
date: "23 November 2025"
output:
  html_document:
    df_print: paged
  pdf_document: default
  word_document: default
authors:
- name: Idowu B. Olawoye
  orcid: "0000-0002-6658-9917"
  equal-contrib: true
  affiliation: 1
- name: Maxim Federov
  equal-contrib: true
  affiliation: 1
- name: Robert A. Petit III
  orcid: "0000-0002-1350-9426"
  affiliation: 2
- name: Jennifer L. Guthrie
  orcid: "0000-0001-8565-203X"
  corresponding: true
  affiliation: 1, 3
bibliography: paper.bib
affiliations:
- name: Department of Microbiology and Immunology, University of Western Ontario,
    London, ON, Canada
  index: 1
- name: Wyoming Public Health Laboratory, Wyoming Department of Health, Cheyenne,
    Wyoming, USA
  index: 2
- name: Public Health Ontario, Toronto, ON, Canada
  index: 3
---

# Summary

The Mycobacteria Amplicon Sequencing Tool (MAST) <https://github.com/guthrielab/MAST> is a modular Nextflow pipeline for antimicrobial resistance prediction and lineage classification of *Mycobacterium tuberculosis* from amplicon or whole-genome sequencing datasets. The workflow automates read processing, variant calling, and annotation to produce standardized, human-readable reports in .docx format, as well as machine-readable summary files. The pipeline is fully reproducible, scalable across computing environments, and available as open-source software on GitHub.

# Statement of need

*Mycobacterium tuberculosis*, the causative agent of tuberculosis, is responsible for an estimated 1.2 million deaths globally each year [@who2024]. Given its global burden and rising drug resistance, *M. tuberculosis* is routinely monitored through public health surveillance, with genome sequencing increasingly used to predict antimicrobial resistance and characterize circulating lineages. These data provide insight into phylogeographic distribution, transmissibility, and resistance patterns [@Meehan:2019; @Napier:2020; @Walker2013]. In addition, targeted amplicon sequencing has recently been used to rapidly detect resistance-associated variants [@Murphy2023].

To support public health laboratories and researchers conducting targeted amplicon sequencing of *M. tuberculosis*, we developed MAST (Mycobacteria Amplicon Sequencing Tool), a Nextflow-based pipeline that parallelizes the analysis of whole-genome and targeted amplicon sequencing data. The pipeline detects sequence variants and compares them against a catalogue of 11,390 high-confidence variants and their associated drug-resistance phenotypes (associated with resistance confidence grading) curated by the World Health Organization, 2nd edition [@who2023; @Walker:2022]. It also applies a lineage barcode schema to assign isolates to 24 major lineages and sub-lineages [@Napier:2020], enabling standardized lineage identification alongside antimicrobial resistance prediction.

# Implementation

MAST is an open-source software that is available on GitHub. It uses conda dependencies such as Filtlong for quality trimming, BWA-MEM, Pysam, and SAMtools for BAM manipulation, Freebayes for variant calling and a custom Python script to compare mutations against the drug-resistance and lineage databases [@Wick; @Li:2009a; @pysam; @Li:2009b; @Garrison:2012] as seen in \autoref{fig:Figure1}.
The workflow produces two output files per input sequence:

1. A user-friendly MS Word (.docx) report summarizing the predicted drug-resistance profile, lineage assignment, and select patient information provided by the user.

2. A tab-delimited summary file (.tsv) suitable for downstream parsing, integration with other tools, or storage in a database.

# Validation

We evaluated MAST against the widely used tool TB-Profiler (v6.6.5) [@Phelan:2019] using 168 publicly available *M. tuberculosis* genomes selected from NCBI, restricting to assemblies with an average sequencing depth of at least 30× and including isolates spanning major global lineages and both drug-susceptible and drug-resistant profiles. Lineage predictions and antimicrobial-resistance outputs produced by MAST were compared with TB-Profiler run using default parameters.

Lineage assignments were identical for all genomes (100% concordance). Antimicrobial-resistance predictions showed 99.7% concordance, with 350 of 351 drug-specific resistance calls matching those reported by TB-Profiler. The single discordant call corresponded to a variant responsible for drug resistance that fell below MAST’s reporting threshold due to filtering criteria, reflecting a difference in variant-handling logic rather than an incomplete mutation list (\autoref{tbl:table1}).


|**Category**                 |**MAST**    |**TB-Profiler**|**Concordance %**|
|:---------------------------:|:----------:|:-------------:|:---------------:|
|**Lineage assignment**       | 168        | 168           | 100%            |
|**AMR prediction**           | 350        | 351           | 99.7%           |

Table: Comparison of TB-Profiler and MAST demonstrating concordance in lineage assignment and antimicrobial resistance (AMR) prediction. \label{tbl:table1}

A key advantage of MAST is its ability to process both targeted amplicon data and whole-genome sequencing data, while remaining accessible to users with diverse levels of bioinformatics expertise. To support reproducibility and practical utility, we engaged domain experts to review the workflow and its applicability in routine settings. The pipeline will be maintained, supported, and openly available to the community for at least 10 years.

# Figures

![Flowchart of the MAST workflow. Green ovals represent quality control steps.\label{fig:Figure1}](MAST_workflow_v3.png)

# Funding

This work was supported by the Western University’s Frugal Biomedical Innovations Catalyst Grant (to I.B.O. and J.L.G.) and by a Natural Sciences and Engineering Research Council of Canada (NSERC) Undergraduate Student Research Award to M.F.

# References
