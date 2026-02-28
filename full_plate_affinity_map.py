"""
This file contains functions used to model binding affinities, using the 
preprocessed csv file generated in R and outputting a fully processed csv 
file and a plot of Kd vs. combined p-value. 
"""

# Import functions 
import affinity_map_functions as affin_map
# Import libraries
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests

def run_plate_affinity(input_file_path,output_file_path,remove_outlier,NA_impute_cutoff,probe_conc,off_compete_conc,process_two = 'min_x', significant_row_p = 0.001):
    """Take relevant experimental information and R output .csv file containing mass spectrometry intensity data to determine affinities of proteins to the ligand of interest

    Args:
        input_file_path (str): absolute/relative path to input R script preprocessed csv file 
        output_file_path (str): absolute/relative path to output table
        remove_outlier (boolean): True = removes outlier input data points
        NA_impute_cutoff (float): if larger than this percentage of NA values for a protein in a probe concentration, 
                                  then do not fit the competitive model for this protein at this condition
        probe_conc (list): list of probe concentrations in same order as sample numbering
        off_compete_conc (list): list of off-compete concentrations in same order as sample numbering
        process_two (str, optional): method by which Kd is reported with the number of EC50 values available for protein. 
                                     Choose from {'min_x', 'mean_EC50', 'orig'}. Defaults to 'min_x'.
        significant_row_p (float, optional): Significance threshold cutoff for what is considered a valid fit used to determine EC50. Default value is 0.001.

    Returns:
        pandas.DataFrame: Output dataframe with affinity data and metadata related to data processing

    Output: 
        csv file with processed pandas.DataFrame 
    """    

    # Load data
    data = pd.read_csv(input_file_path, header=0)

    # Define column mapping 
    # Map experimental conditions to dataframe columns containing the signal intensities
    data_start_col_num = affin_map.find_start_column(data)
    column_mapping = affin_map.make_column_map(probe_conc,off_compete_conc,data,data_start_col_num)

    # Get the high and low bounds of the off-compete concentrations
    off_compete_conc_low = sorted(set(off_compete_conc))[1] # get the second smallest value since vehicle
    off_compete_conc_high = max(off_compete_conc)

    # Get column index of the percent NA of each probe concentration for each protein
    conc1_idx = data.columns.tolist().index("NA_percent_conc1")
    conc2_idx = data.columns.tolist().index("NA_percent_conc2")
    conc3_idx = data.columns.tolist().index("NA_percent_conc3")
    conc4_idx = data.columns.tolist().index("NA_percent_conc4")
    # NA column indexes, which maps to the probe concentration setup
    NA_idxs = [conc1_idx,conc2_idx,conc3_idx,conc4_idx] * 24

    # Determine global minimum intensity: used as lower bound limit in curve fitting optimization
    global_zmin = affin_map.find_global_zmin(data, column_mapping)

    # Organize data into slices:
    # List of lists
    # [[(0.001, 0.0, 467.926743499985),
    #  (0.001, 0.0, 5204.11421317932),...],
    # [(0.001, 0.0, 5602.10499679806),...]]
    # Each of the inner lists represent a protein at each experimental concentration
    # Data is in same order as experimental setup and protein sequence
    # Each entry in the list follows: (probe concentration, off-compete concentration, intensity of peak)
    datasets = []
    for index, row in data.iterrows():
        row_data = []
        for col_idx in column_mapping:
            if col_idx < len(row):
                x_val, y_val = column_mapping[col_idx]
                z_val = row[col_idx]
                # Filter x_val sets with too many missing values
                NA_perc_column = NA_idxs[probe_conc.index(x_val)]
                NA_perc_of_x = row[NA_perc_column]
                # Only include values if less than a certain number are imputed values
                if NA_perc_of_x < NA_impute_cutoff:
                    row_data.append((x_val, y_val, z_val))
        datasets.append(sorted(row_data))

    # First pass curve fitting. Bounds of measurable affinities set here. 
    # P values tested from F-statistics derived from a zero-slope linear fit as the null comparison.  

    # Store all fit results from testing Kd fit into fit_results
    fit_results = []
    # Go through data protein-by-protein
    for idx, dataset in enumerate(datasets):
        data_by_x = {}
        # x_val = concentration of probe
        # y_val = concentration of off-compete
        # z_val = MS intensity
        for x_val, y_val, z_val in dataset:
            # Sort into dictionaries based on probe concentration
            if x_val not in data_by_x:
                data_by_x[x_val] = {'y': [], 'z': []}

            data_by_x[x_val]['y'].append(y_val)
            data_by_x[x_val]['z'].append(z_val)

        # Fit sigmoidal curve
        # Get dictionary of fit by each probe cocentration
        fit_by_x = affin_map.fit_protein_sigmoid(data_by_x, global_zmin)

        # Remove outlier if it is more than 3 sigma away from the fit
        if remove_outlier:
            fit_by_x, data_by_x, protein_row_data = affin_map.outlier_detect_and_removal(fit_by_x,data_by_x,dataset,global_zmin)
            datasets[idx] = protein_row_data

        # If the sigmoid is significant, take out vehicle to test if it is a confounding factor 
        fit_by_x = affin_map.remove_vehicle_test(fit_by_x,data_by_x,dataset,global_zmin)
        
        # Store data related to the protein in fit_results_row
        fit_results_row = []
        for x_val in fit_by_x:
            (param_full, result_full, f_stat, p_value) = fit_by_x[x_val]
            fit_results_row.append((x_val, param_full, result_full, f_stat, p_value))

        # Store all output in the global list
        fit_results.append(fit_results_row)

    # Data organization used for Cheng-Prusoff linear fitting to deconvolute Kd of the unmodifeid ligand molecule   
    all_b_values = []
    all_b_errors = []
    all_x_values = []
    all_p_values = []

    # Process protein-by-protein
    for dataset_fits in fit_results:
        row_bs = []
        row_b_errors = []
        row_x_vals = []
        row_p_vals = []
        for x_val, params, fit_result, f_stat, p_val in dataset_fits:  # Unpacking sigmoidal fits
            row_x_vals.append(x_val)
            # Check if there are values in b (EC50)
            try:
                param_b = params['b']
            except IndexError as e:
                row_bs.append(np.nan)
                row_b_errors.append(np.nan)
                row_p_vals.append(np.nan)
                print(f'Index error for param, did not have index b. Eror message: {e}')
                continue

            if not np.isnan(param_b.value):  # Check if the b parameter (EC50) is valid
                row_bs.append(param_b.value)
                row_b_errors.append(param_b.stderr)  # Extract the standard error of b from the covariance matrix
                row_p_vals.append(p_val)

            else:
                row_bs.append(np.nan)
                row_b_errors.append(np.nan)
                row_p_vals.append(np.nan)

        all_b_values.append(row_bs)
        all_b_errors.append(row_b_errors)
        all_x_values.append(row_x_vals)
        all_p_values.append(row_p_vals)

    # Cheng-Prusoff modeling, Fisher p value combination
    all_row_fit_results = []

    for row_x_vals, row_bs, row_b_errors, row_p_vals in zip(all_x_values, all_b_values, all_b_errors, all_p_values):
        row_fit_results = affin_map.kd_by_lin(row_x_vals, row_bs, row_b_errors, row_p_vals, process_two, significant_row_p)
        all_row_fit_results.append(row_fit_results)

    # Put the new columns of kd fit at the front, calculate adjusted p value, reorder columns, and save .csv output file
    dataframe_processed = affin_map.fit_dataframe_for_output(data,all_row_fit_results)
    #BH FDR adjusted p value
    fdr = 0.0001
    com_p = dataframe_processed['combined_p'] 
    combined_p_list = com_p.apply(lambda x: 10**(-x) if isinstance(x, (int, float)) else 1) # Set all special instances to 1
    bh_adjust = multipletests(combined_p_list, alpha = fdr, method = 'fdr_bh')
    # Check if the Kd value is within off competition limits
    is_in_range = dataframe_processed['log10_kd'].between(off_compete_conc_low, off_compete_conc_high)
    dataframe_processed['Kd_within_range'] = is_in_range
    # Boolean: True or False
    dataframe_processed['passes_fdr'] = bh_adjust[0]
    # If there are no values reported for log10_kd, too many of the values were imputed; do not report row
    dataframe_processed = dataframe_processed[dataframe_processed['log10_kd'] != '']

    # Order rows that pass FDR by Kd
    pass_fdr = dataframe_processed.loc[dataframe_processed['passes_fdr']]
    not_pass_fdr = dataframe_processed[~dataframe_processed.index.isin(pass_fdr.index)]
    pass_fdr = pass_fdr.sort_values(by='combined_p',ascending=False) # sort by p value
    ordered_data = pd.concat([pass_fdr,not_pass_fdr])
    ordered_data.to_csv(output_file_path)    
    return ordered_data
