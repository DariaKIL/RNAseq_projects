# RNA-Seq Projects Repository

Welcome to the RNA-Seq projects repository! This repository contains a collection of pet projects focused on various RNA sequencing (RNA-Seq) data analysis approaches. Each folder within this repository corresponds to a different project, where you'll find comprehensive explanations, plots, and detailed notes that I've made for myself and for others who might find them useful.

## Repository Structure

- **One folder = One project:** Each folder contains the scripts, data, and results specific to an individual project. Explore the folders to dive into each RNA-Seq analysis project.

## Project 1: Wound Signaling in Arabidopsis

### Overview

This project, titled **Wound_Signaling_Arabidopsis**, explores the transcriptomic responses in Arabidopsis thaliana, particularly focusing on the role of the gene *IMPAIRED SUCROSE INDUCTION 1 (ISI1)* in long-distance wound signaling.

### Study Summary

- **Objective:** To assess the impact of *ISI1* on the wound response by comparing the transcriptomes of wounded *ISI1* mutants with those of wounded wild-type plants.
- **Data Source:** The dataset used in this analysis is publicly available and can be accessed through the [NCBI GEO database](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157938).
- **Experimental Design:**
  - **Conditions:** 
    1. WT - control / *ISI1* mutant plants
    2. W - wounded plants / UW unwounded plants
  - **Replicates:** Each group included three biological replicates.
  - **Sampling:** Distal leaves 13 were collected one hour after wounding leaf 8 in wounded plants, or without any treatment in unwounded plants. 
  - **Sequencing:** The study employed a single-read sequencing design.

### Results and Analysis

The detailed results, including plots and interpretations, are available at the following URL: [View Results](https://dariakil.github.io/RNAseq_projects/Wound_Signaling_Arabidopsis/Wound_Signaling_Arabidopsis.nb.html).

## Project 2: Micro-RNA Prediction

This project focuses on predicting important microRNAs (miRNAs) involved in the development of disease, specifically Acute Respiratory Distress Syndrome (ARDS), by analyzing miRNA and mRNA data from patients with ARDS. The data was sourced from the article [“T cell dysfunction in elderly ARDS patients based on miRNA and mRNA integration analysis”](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1368446/full). The objective is to identify miRNAs that might play a critical role in the pathogenesis and severity of ARDS.

**Background**

Acute Respiratory Distress Syndrome (ARDS) is a life-threatening condition that results in respiratory failure, often occurring in critically ill patients. Despite its significant impact, the molecular mechanisms underlying ARDS are not fully understood. In this study, both mRNA and miRNA were sequenced from patients with ARDS and healthy controls, with the goal of better understanding the pathogenesis of ARDS.

**Project Workflow**

This project processes the sequencing data from ARDS patients and healthy donors, using the following approach:
1. Data Processing:
- Raw sequencing data (FASTQ format) from 49 files were processed using a Snakemake pipeline.
- The pipeline uses STAR for aligning RNA-seq data and [HTSeq](https://htseq.readthedocs.io/en/latest/) for further analysis
2. Modeling:
- After pre-processing and normalization, machine learning models were trained to predict important miRNAs involved in ARDS.
- Models such as Support Vector Machines (SVM), Nearest Shrunken Centroids (NSC), and others were applied to classify and predict miRNAs.
- The trained models were compared using the MLSeq R package to assess their effectiveness in predicting important miRNAs.
- Evaluation metrics (accuracy, sensitivity, specificity) were used to select the best performing model.

**Requirements**
1. Snakemake: A workflow management system for data analysis.
2. STAR: RNA-seq aligner for genome-based alignment.
3. HTSeq: Python package for analysis of high-throughput sequencing data
4. R: For model training and evaluation, including the MLSeq package.
 
**Model Evaluation Results**

After evaluating several models using cross-validation, NSC (Nearest Shrunken Centroids) was identified as the best-performing model for predicting important miRNAs in ARDS.
Key Metrics:
- Accuracy: 94.12%
- Sensitivity: 81.82%
- Specificity: 100.00%

