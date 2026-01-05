# Shallow-Wetland-Modeling
![R](https://img.shields.io/badge/language-R-blue)

## Overview

This draft repository is a companion to two in-progress manuscripts. 

The first manuscript, "Mapping Shallow Wetlands for Wildlife Conservation", outlines an approach for identifying shallow freshwater wetlands at landscape scales.

The second manuscript, "Opportunistic Stopovers by Hudsonian Godwits Limosa haemaestica in a Midcontinental Agricultural Landscape," incorporates these and other predictions into habitat modeling for a migratory shorebird species.

<br/>

## Repository contents

- 'scripts/' contains python and R scripts for data processing
- 'data/' contains analysis-ready tabular data and GIS vector files

<br/>

## Scripts

- 'Shallow Wetlands/'
	- '01_rasters_to_superpixels.R' demonstrates conversion surface reflectance rasters to l*A*B-based superpixels 
	- '02_probability_forest_models.R' fits three competing models for shallow water detection
	- '03_climate_predictions.R' evaluates the influence of climate on shallow water area

- 'HUGO Habitat/'
	- '01_stopover_selection_model.R' fits a separate integrated point process model for stopover area selection
	- '02_survival_models.R' fits a Cormack-Jolly-Seber survival models for godwits

<br/>

## Acknowledgements

Code development was supported by a NASA FINESST award (#21-Earth210-0404) and a University of Massachusetts Amherst Agricultural Experiment Station Grant (#MAS00592). 

<br/>

## Contact

For questions, contact me at linscotj@email.sc.edu.

