# AffinitySurveyor

This database contains files used to analyze proteomics data from micro-mapping after the mass spectrometry data has already been processed by DIANN software. 

Polodsky.R is an R script used to generate volcano plots to understand the difference in fold difference in one condition compared to another.

The **R_affinity_mapping_preprocessing.R** contains code used to normalize and impute data generated from DIANN tsv files and output csv file used for full_plate_affinity_map.py.  (Previous version of preprocessing code in affinity_map_data_cleanup_from_DIANN.R and affinity data generating code in affinity_map_code.py)

**affinity_surveyor.py** contains code used to run a fully processed csv file for all identified and preprocessed proteins which contains apparent Kd values and sigmoidal curve fit values and errors. **full_plate_affinity_map.py** contains code used to generate the fits and apparent Kd values. 

**get_graph_for_protein.py** specifies the proteins interested in seeing the sigmoidal fit and generates sigmoidal fits and Kd apparent fits for these proteins, outputting graphs with EC50 fits and Cheng-Prusoff fits. 

**affinity_surveyor.py** and **get_graph_for_protein.py** use functions in affinity_map_functions.py for data processing and graph generation. 

Folder ./sample/ contains sample input and output files compatible as input for R_affinity_mapping_preprocessing_with_probe_conc_norm.R and subsequently affinity_map_code.py. 

