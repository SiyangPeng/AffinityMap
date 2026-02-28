"""
This script uses the R preprocessed csv file to measure binding affinities, outputting an Affinity Surveyor output table.
"""

# Processing options
# Outlier removal: remove or do not remove 3-sigma outlier data points
remove_outlier = True
# Save plot of affinities
savefig = True
# If this fraction of values are imputed or more, do not use a catalyst concentration for Cheng-Prusoff modeling
NA_impute_cutoff = 0.25 # More relaxed conditions are higher
# Define input and output files
input_file_path = 'input_path/data_preprocessed_separate.csv'
output_file_path = 'output_path/full.csv'

# Import functions used for data processing
import full_plate_affinity_map as full_plate

#### Example Setup ####

# Define experimental matrix based off of how experiment condition was set up. This is used to map columns in the R output file to conditions
catalyst_conc = [10,1,0.1,0.01] * 24 
off_compete_conc = [0, 10, 3.16, 1, .316, .1, .0316, .01, .00316, .001, .000316, .0001] * 8

# Perform Affinity Surveyor processing
sorted_data = full_plate.run_plate_affinity(input_file_path,output_file_path,remove_outlier,NA_impute_cutoff,catalyst_conc,off_compete_conc)


