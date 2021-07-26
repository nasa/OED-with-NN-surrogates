# Optimal Experimental Design with Neural Network Surrogates
Authors: Joshua Stuckner, Matt Piekenbrock

This code was used in the following papers:  <br>
 - Joshua Stuckner, Matthew Piekenbrock, Steven M Arnold, Trenton M Ricks. (2021) Optimal Experimental Design With Fast Neural Network Surrogate Models To be published in Computational Materials Science
 - Joshua Stuckner, Matthew Piekenbrock, Steven M Arnold, Trenton M Ricks. (2021) Optimal Experimental Design With Fast Neural Network Surrogate Models. NASA Technical Reports Server, TM-20205003868
 - Arnold, S. M., Piekenbrock, M., Ricks, T. M., & Stuckner, J. (2020). Multiscale Analysis of Composites Using Surrogate Modeling and Information Optimal Designs. In AIAA Scitech 2020 Forum (p. 1863).

## How to use
The entry point to run the four experiments in the Computational Materials Science paper are in the vignettes folder. Each experiment has its own notebook: PMC_class.Rmd, MMC_class.Rmd, CMC_class.Rmd, VF_experiment.Rmd. The data and trained models for these experiments are in the data folder.

The VF_experiment will be used as an example.

### Step 1 - Create training data (lines 30 - 121)
This step cannot be performed without MAC/GMC, a physics based composite modeling software. The parsed data from this step is included in the data folder and may be skipped.

### Step 2 - Load and clean the data (lines 124 - 213)

### Step 3 - Train the neural network surrogate (lines 214 - 372)
This step can be skipped by loading the trained models directly from the data folder. This step ended up being replaced for the VF experiment by the Hyperparameter optimization step.

### Step 3b - Hyperparameter optimization (lines 375 - 523)
This step can be skipped by loading the trained models directly from the data folder. 

### Step 4 - Perform OED (lines 525 - 801)
