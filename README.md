# Prediction Of Outcome & Location (POOL)

This is the github repo for code accompanying the paper: "Multi-centre Epilepsy Lesion Detection (MELD) Project: Predictors of lesion location and postsurgical seizure freedom in focal cortical dysplasia" (insert link to preprint / paper).

### To install the pool package:
Most packages will be automatically installed with
pip install -e .

### Preprocessing
We have scripts to run preprocessing on our demographic data and our lesion masks (scripts/preprocess_qc.sh). This was to deal with mistakes / n.a. in data submitted from individual sites and create preprocessed demographic data and lesion data. However, it is likely to only be useful for our MELD data. 

To run the notebooks on dummy data (dummy_demographics_qc.csv & dummy_lesions.npz) you do NOT need to do any preprocessing but you need to rename lesions_dummy.npz to lesions.npz & dummy_demographics_qc.csv to demographics_qc.csv & then move straight on to Notebook 1. 
Warning : dummy datasets are synthetic data generated to play with the notebooks. Results cannot be interpreted.

If you are using the notebooks on your own data, you need to QC your data and have a demographics.csv file and a lesions.npz file.
The demographics.csv file has a row for each patient and has the following columns: 
 - ID	
 - Age of onset	of epilepsy
 - Duration	of epilepsy
 - Age at preoperative scan
 - Sex	
 - Ever reported MRI negative	
 - Engel Outcome	Surgery	Lesion area	Hemisphere	Lobe	Histology	Seizure free	lesion_masked
The lesions.npz file has a row for each patient and a column for every vertex on fsaverage_sym. 0 = non lesional vertex. 1 = lesional vertex.

#### To run the preprocessing on the MELD data
cd scripts

bash preprocess_qc.sh

### Notebook 1: Demographics table
This notebook:
- Calculates data required for table 1 and 2 in the paper
- Calculates number of patients with LH or RH lesions & number of lesions in each cortical lobe
- Plots the proportion of lesions in each lobe for each histological subtype

### Notebook 2: Lesion locatioons
This notebook (notebooks/lesion_locations.ipynb):
- creates lesion map images for the whole cohort, operated and non-operated cohorts and histopathological subtypes
- estimates sample size needed for pattern consistency
These and script 1 make up figure 1 and supp fig 2

### script 1 make volumetric lesions gif
video still and gif of heat maps
scripts/make_lesions_gif.py

### script 2 demographic regression
scripts/run_all_regressions.sh
This call scripts/run_regressions.py with various arguments
methods for this described in supp fig 1
WARNING it is configured to run on an hpc and might be fiddly to get working

### Notebook 3: plotting regressions
This notebook generates figures 2 & 3A:
- Plots results from regression analyses
- Carries out spatial comparison with various explanatory maps
- Carries out spin testing

### Notebook 4: Seizure outcome plots
This notebook calculates seizure freedom scores for figure 3B.
- Plots per subject seizure freedom scores

### Notebook 5: Grid analysis
This notebook calculates the relationships between clinical variables
to create grid plots for Figure 4 and supp fig 3


### Notebook 6: Figure 4 plots
Plot remaining analyses of interest for Figure 4.
- Rainclouds plots for the distributions of age etc
- Plots for lesion size, age of onset and surgery relationships
- Seizure free duration raincloud plots


