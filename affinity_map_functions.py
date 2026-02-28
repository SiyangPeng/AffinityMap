"""
This file contains ancillary functions used for affinity profiling
"""

# Import functions
import numpy as np
import pandas as pd
import lmfit
from scipy.optimize import curve_fit
from scipy.stats import f
from scipy.stats import t
import matplotlib.pyplot as plt
from scipy.stats import chi2
from matplotlib.gridspec import GridSpec
from scipy.stats import norm
import random
import regex as re
import matplotlib.lines as mlines

def model_func_sig(params, x):
    """Logistic model for competitive binding 
    Args:
        params (lmfit Parameter object): contains parameter 'a', 'b', and 'c' to be fit to data; 
                                         with boundary and initial conditions defined
        x (numpy.ndarray): off-compete concentrations
    Returns:
        numpy.ndarray or None: Modeled signal intensity value 
    """    
    a = params['a'].value
    b = params['b'].value 
    c = params['c'].value 

    d = a + c
    return (d / (1 + 10**(x-b))) + c

def model_reduced(params, x):
    """Reduced model: 1 degree of freedom for lmfit comparison

    Args:
        params (lmfit Parameter object): contains parameter 'a' to be fit to data; 
                                         with boundary and initial conditions defined
        x (numpy.ndarray): off-compete concentrations

    Returns:
        numpy.ndarray or None: Modeled signal intensity value 
    """ 
    a = params['a'].value
    return np.full_like(x, a)


def model_func_2(x, a, b):
    """Cheng-Prusoff model to measure Kd from EC50s 
    Args:
        x (numpy.ndarray): array of probe concentrations
        a (float): parameter variable with value determined by scipy.optimize.curve_fit
        b (float): parameter variable with value determined by scipy.optimize.curve_fit
    Returns:
        numpy.ndarray or None: array of fitted log10(EC50) values
    """    
    return np.log10(a + (b * x))

def model_reduced_curve(x, a):
    """Reduced model: 1 degree of freedom
    Args:
        x (numpy.ndarray): array of probe concentrations
        a (float): parameter variable with value determined by scipy.optimize.curve_fit
    Returns:
        numpy.ndarray or None: array of fitted log10(EC50) values
    """    
    return np.full_like(x, a)

def model_func_lin_fit(a, b, x):
    """Cheng-Prussof linear model for graphical depiction of fitted values.
       Produce schild plots for inspection. 
    Args:
        a (float): parameter variable with value determined by scipy.optimize.curve_fit
        b (float): parameter variable with value determined by scipy.optimize.curve_fit
        x (numpy.ndarray): array of probe concentrations
    Returns:
        numpy.ndarray or None: array of fitted EC50 values
    """    
    return (a + (b * x))

def objective(params, x, y, model):
    """Define the objective function for lmfit, which is sum of squares of deviation from the nonlinear function
       Goodness of fit for model defined by R squared of (fitted intensity) - (signal intensity) 
    Args:
        params (lmfit Parameter object): contains fitted model parameters (example 'a')
        x (array_like): array of input values used in model
        y (array_like): array of fitted values
        model (callable): model function used for fitting
    Returns:
        int or None: sum of squares of nonlinear function of fitted model
    """    
    fit = model(params, x)
    return (fit-y)**2

def find_start_column(df):
    """Find the column of the first sample from input data file using column annotation
    Args:
        df (pandas.DataFrame): input datafile with mass spectrometry signal and metadata
    Returns:
        int: column number of the first data column with mass spectrometry signal data
    """    
    # Define the regex pattern
    pattern = r'\d+_Slot\d\.\d+_'
    # Find column names that match the regex pattern
    matching_columns = df.filter(regex=pattern).columns

    if not matching_columns.empty:
        # Get the index of the first matching column
        first_matching_column = matching_columns[0]
        data_start_col_num = df.columns.get_loc(first_matching_column) - 1
    else:
        print('Did not find DIANN output data column marked by Slot1')
        return 
    
    return data_start_col_num

def make_column_map(catalyst_conc,off_compete_conc,df,data_start_col_num=9):
    """Matches the column in the preprocessed data to an experimental condition.
    Args:
        catalyst_conc (list): probe concentration in the same order as the experimental setup
        off_compete_conc (list): off-compete ligand concentration in the same order as experimental setup
        sample_number (int, optional): number of columns with samples. Defaults to 96.
        data_start_col_num (int, optional): number of metadata columns prior to columns with experimental values. Defaults to 9.
    Returns:
        dictionary: map dataframe column index to experimental condition
                    column_index: (ligand-catalyst concentration, off-compete concentration)
    """
    column_mapping = dict()
    count = 1
    pattern = r'(\d+)_Slot1\.\d+_'
    column_names = list(df.columns)
    for col in column_names:
        m = re.search(pattern, col)
        if m:
          i = int(m.group(1))
          if i % 2 == 1:
              column_mapping[count+data_start_col_num] = (catalyst_conc[i // 24],off_compete_conc[(i-1)//2])
          else:
              column_mapping[count+data_start_col_num] = column_mapping[count+data_start_col_num-1]
          count += 1
    return column_mapping

def find_global_zmin(data, column_mapping):
    """ Determine global minimum intensity
        Used as lower bound limit in curve fitting optimization
    Args:
        data (pandas.DataFrame): dataframe with all metadata and data of each sample with proteins as rows
        column_mapping (dictionary): dictionary with column index of each sample as key
    Returns:
        float: lowest signal intensity detected
    """
    z_min_determination = []
    for index, row in data.iterrows():
        # Use columns that have samples by the column indexes defined as having samples
        for col_idx in column_mapping.keys():
            if col_idx < len(row):
                # Put all values into a list
                z_min_determination.append(row.iloc[col_idx])
    z_min_array = np.array(z_min_determination)
    # Keep the smallest number in the array
    cutoff = np.percentile(z_min_array, 1)
    global_zmin = z_min_array[z_min_array <= cutoff].mean()
    return global_zmin

def sigmoid_curve_gen(x_val, y_values, z_values, global_zmin, b_seed):   
    """Generate competitive model (full) and null hypothesis model (reduced) for all probe concentrations of given protein
       Test random b value seed values and keep the most significant sigmoidal fit.
    Args:
        x_val (float): probe concentration for fitting
        y_values (numpy.ndarray): array of off-compete ligand concentrations 
        z_values (list): array of mass spectrometry normalized and imputed signal intensities in the same order as ligand array above
        global_zmin (float): lowest signal intensity detected
        b_seed (float): Seed value for parameter b (correponding to log10(EC50)) to be optimized by nonlinear fitting
    Returns:
       result_full (lmfit MinimizerResult): object containing optimized parameters and goodness-of-fit statistics for full logistic model after fitting
       result_reduced (lmfit MinimizerResult): object containing optimized parameters and goodness-of-fit statistics for reduced model after fitting
    """
    try:
        # Full model fitting
        # Create a Parameters object and set initial guesses
        params_full = lmfit.Parameters()
        params_full.add('a', value=max(z_values), min=0.001, max=3000000000000)
        params_full.add('b', value=b_seed, min=-10, max=4)
        params_full.add('c', value=global_zmin, min=global_zmin, max=3000000000000 - global_zmin - 1)
        # TODO: have a directional constraint for c to be less than a, so that the model is sigmoidal and not inverse sigmoidal
        # TODO: need to have constraint on b that it needs to be in the competition range
        # Perform the minimization of sum of squares of nonlinear function
        result_full = lmfit.minimize(objective, params_full, args=(y_values, z_values, model_func_sig))

        # Reduced model fitting
        # Create a Parameters object and set initial guesses
        params = lmfit.Parameters()
        params.add('a', value=np.mean(z_values))
        # Perform the minimization
        result_reduced = lmfit.minimize(objective, params, args=(y_values, z_values, model_reduced))
        return result_full, result_reduced 

    except RuntimeError as e:
        print(f"Fit did not converge for x={x_val}: {e}")
        return None, None
    
    except TypeError as e:
        print(f"Error: {e}")
        return None, None

# TODO: filter the ec50 values to only include values within competition range (if not done in the modeling itself)
def best_sigmoid(x_val, y_values, z_values, global_zmin):
    """Function that varies the starting value of b to try to find best sigmoidal fit to data, defined by lowest p value
    Args:
        x_val (float): probe concentration for fitting
        y_values (numpy.ndarray): array of off-compete ligand concentrations 
        z_values (list): array of mass spectrometry normalized and imputed signal intensities in the same order as ligand array above
        global_zmin (float): lowest signal intensity detected
    Returns:
        best_param_full (lmfit Parameters): best fitted parameters for model, lowest p value from F statistics test
        best_result_full (lmfit MinimizerResult): object containing best_param_full and goodness-of-fit statistics for sigmoidal model generated by lmfit
        best_f_stat (float): F statistic of best fitted sigmoidal model compared to the reduced model
        best_p_value (float): p value of the F statistic by model comparison
    """    
    # Save the best fit (lowest p_value)
    best_result_full = None
    best_param_full = None
    best_result_reduced = None
    best_f_stat = np.nan
    best_p_value = np.inf

    # Fit the sigmoidal curve
    for _ in range(30): 
        b_seed = random.uniform(-5,3)
        try:
            result_full, result_reduced = sigmoid_curve_gen(x_val, y_values, z_values, global_zmin, b_seed)
        except ValueError as e:
            print(f'No z values found for x_val {x_val}, with exception: {e}')
            continue

        # If fitting fails, then all parameters are nan
        if np.all(pd.isnull(result_full)):
            f_stat, p_value = np.nan, np.nan
        # F-statistic test
        else:
            f_stat, p_value = f_stat_test(y_values, z_values, result_full, result_reduced, objective, model_func_sig, model_reduced, 3)
            # Access the fitted parameters
            fitted_params = result_full.params
        
        # Ensure that the error for b is valid 
        if result_full is None:
            continue
        params = result_full.params
        error_b = params['b'].stderr 

        # Check if this p-value is the lowest value yet generated in the loop
        if (not p_value == np.nan) and (p_value < best_p_value) and error_b != None:
            best_p_value = p_value
            best_f_stat = f_stat
            best_param_full = fitted_params
            best_result_full = result_full
            best_result_reduced = result_reduced
    
    # Get the best fitted parameters
    return best_param_full, best_result_full, best_f_stat, best_p_value

def fit_protein_sigmoid(protein_x_data, global_zmin):
    """Generate sigmoidal fit for a given protein
    Args:
        protein_x_data (dict): dictionary with concentrations as keys containing nested dictionaries of two lists of values:
                               value of 'y' contains off-compete concentrations, 'z' contains signal intensities matching 'y' conditions
        global_zmin (float): lowest signal intensity detected
    Returns:
        dict: probe concentrations are the keys and tuple of best fitted competitive curve fitted lmfit Parameters, 
        lmfit MinimizerResult best fitted value, the F statistics of the nonlinear fit, and the p value of this fit compared to a reduced model
    """    
    prot_sig_fits = dict()
    for x_val, values in protein_x_data.items():
        # Since log of 0 is negative infinity, 0 concentration is defined as 0.00001
        y_values = np.log10([y if y > 0 else 0.00001 for y in values['y']])
        z_values = values['z']

        # Get best fit
        best_param_full, best_result_full, best_f_stat, best_p_value = best_sigmoid(x_val, y_values, z_values, global_zmin)

        # If fitting fails:
        if best_result_full is None:
            nan_array = np.full((3,), np.nan)
            best_p_value = np.nan
            best_f_stat = np.nan
            best_param_full =  nan_array 
            best_result_full = nan_array

        # Save fits for a given catalyst concentration
        prot_sig_fits[x_val] = (best_param_full, best_result_full, best_f_stat, best_p_value)

    return prot_sig_fits

def f_stat_test(orig_x_vals, orig_y_vals, result_full, result_reduced, rss, model_full, model_reduced, par_num):
    """F statistic testing for model parameters optimized with lmfit
       Indicates the goodness of fit of a model compared to a reduced model with fewer parameters
    Args:
        orig_x_vals (arry_like): array of off-compete concentrations
        orig_y_vals (array_like): array of signal intensity values the competitive sigmoidal model uses
        result_full (lmfit MinimizerResult): optimized parameters and goodness-of-fit statistics for full model
        result_reduced (lmfit MinimizerResult): optimized parameters and goodness-of-fit statistics for reduced model
        rss (callable): model to calculate residuals
        model_full (callable): model function with largest number of parameters used for fitting
        model_reduced (callable): model function with smaller numer of parameters used for fitting
        par_num (int): number of parameters in full model
    Returns:
        F (float): F statistic for full model compared to reduced model
        p_value (float): p value of the F statistic given the degrees of freedom of both the full and reduced model
    """    
    # Calculate residuals
    # Full model 
    rss_full = np.sum(rss(result_full.params, orig_x_vals, orig_y_vals, model_full))
    # Reduced model
    rss_reduced = np.sum(rss(result_reduced.params, orig_x_vals, orig_y_vals, model_reduced))

    # Degrees of freedom
    n = len(orig_y_vals) 
    df_full = n - par_num  # Degrees of freedom for full model
    df_reduced = n - 1  # Degrees of freedom for reduced model

    # F-statistic
    F = ((rss_reduced - rss_full) / (df_reduced - df_full)) / (rss_full / df_full)

    # p-value from the F-distribution
    p_value = 1 - f.cdf(F, df_reduced - df_full, df_full)

    return F, p_value

def f_stat_test_curvefit(orig_x_vals, orig_y_vals, popt_full, popt_reduced, func, func_reduced, par_num):
    """F Statistic testing for model parameters optimized with curve_fit
       Indicates the goodness of fit of a model compared to a reduced model with fewer parameters
    Args:
        orig_x_vals (array_like): array of probe concentrations
        orig_y_vals (array_like): array of log10(EC50) values used for modeling 
        popt_full (array): optimal values for parameters for full model by curve_fit
        popt_reduced (array): optimal values for parameters for reduced model by curve_fit
        func (callable): model function with largest number of parameters used for fitting
        func_reduced (callabe): model function with smaller number of parameters used for fitting
        par_num (int): number of parameters in the full model

    Returns:
        F (float): F statistic of full model compared to reduced model
        p_value (float): p value of the F statistic given the degrees of freedom of both the full and reduced model
    """    
    # Calculate residuals
    residuals_full = orig_y_vals - func(orig_x_vals, *popt_full)
    if par_num == 2:
        residuals_reduced = orig_y_vals - func_reduced(orig_x_vals, popt_reduced[0])
    if par_num == 3:
        residuals_reduced = orig_y_vals - func(orig_x_vals, popt_reduced[0])
    ss_res_full = np.sum(residuals_full**2)
    ss_res_reduced = np.sum(residuals_reduced**2)
    df_full = len(orig_y_vals) - par_num  # Degrees of freedom for full model
    df_reduced = len(orig_y_vals) - 1  # Degrees of freedom for reduced model
    
    # F-statistic
    f_stat = ((ss_res_reduced - ss_res_full) / (df_reduced - df_full)) / (ss_res_full / df_full)
    
    # p-value from the F-distribution
    p_value = 1 - f.cdf(f_stat, df_reduced - df_full, df_full)

    return f_stat, p_value

def outlier_detect_and_removal(fit_result,data_by_x,protein_data,global_zmin):
    """Handle ouliers, defined as being 3-sigma removed from the fit model
    Args:
        fit_result (dict): Contain EC50 values for each protein at each probe concentration,
                           probe concentrations are the keys and tuple of best fitted competitive curve fitted lmfit Parameters, 
                           lmfit MinimizerResult best fitted value, the F statistics of fit, and the p value of this fit
        data_by_x (dict): dictionary with catalyst concentration concentration as keys and nested dictionary containing two lists of values:
                          'y' contains a list of off-compete concentrations and 'z' contains a list of signal intensities matching 'y' conditions
        protein_data (list): contain tuples of 3 values, in the format of (probe concentration, off-compete concentration, signal intensity value)
        global_zmin (float): lowest signal intensity detected
    Returns:
        the following items with all values relating to the defined outliers removed
        fit_result (dict),data_by_x (dict),protein_data (list)     
    """    
    for x_value, values in data_by_x.items():
        # Model fitting parameters
        params, result, f_stat, p_value = fit_result[x_value]
        # Since log of 0 is negative infinity, 0 concentration is defined as 0.00001
        y_values = np.log10([y if y > 0 else 0.00001 for y in values['y']])
        z_values = values['z']
        # Check for valid fit
        if params is None:
            continue
        try:
            # Compute the fitted values
            fitted_values = model_func_sig(params, y_values)
        except IndexError as e:
            print(f'Error in params: {params}')
            continue

        # Calculate residuals
        residuals = z_values - fitted_values
        # Calculate the standard deviation of the residuals
        std_residuals = np.std(residuals)
        # Identify outliers as points where the residual is greater than 3 standard deviations away from the fit
        threshold = 3 * std_residuals
        outliers = np.abs(residuals) > threshold
        
        # Find indices of outliers and remove
        outlier_indices = np.where(outliers)[0]
        if outlier_indices.size > 0:
            remove_from_orig_data_idx = []
            ## Remove the outlier and rerun the fitting
            for index in sorted(outlier_indices, reverse = True):
                orig_y = values['y'][index]
                remove_value = (x_value, orig_y,z_values[index])
                remove_from_orig_data_idx.append(remove_value)
                y_values = np.delete(y_values,index)
                z_values = np.delete(z_values,index)
            new_y_values = np.power(10,y_values)
            data_by_x[x_value]['y'] = new_y_values
            data_by_x[x_value]['z'] = z_values
            try: 
                protein_data = [x for x in protein_data if x not in remove_from_orig_data_idx]
            except ValueError as e:
                print('Outlier removal failure:')
                print(f'protein data to remove from is {protein_data}')
                print(f'identified outlier to remove is {remove_from_orig_data_idx}')
        
    # Refit the data
    fit_result = fit_protein_sigmoid(data_by_x, global_zmin)
    return fit_result,data_by_x,protein_data

def remove_correspond_0(list1, list2):
    """Remove the values corresponding to 0 in list1 from list2 and return the altered list
    Args:
        list1 (list): The first list with 0 to be removed
        list2 (list): The second list with values at places corresponding to 0 to be removed
    Returns:
        list1_filtered (list): The first list with 0s removed
        list2_filtered (list): The second list with valuescorresponding to 0s in list1 removed
    """    
    # Find indices of 0 in list1
    indices_to_remove = [index for index, value in enumerate(list1) if value == 0 or value == -5] # because set the 0 to log10(0.00001) = -5
    # Remove corresponding values from list1 and list2
    list1_filtered = [value for index, value in enumerate(list1) if index not in indices_to_remove]
    list2_filtered = [value for index, value in enumerate(list2) if index not in indices_to_remove]

    return list1_filtered, list2_filtered

def remove_vehicle_test(fit_result,data_by_x,protein_data,global_zmin):
    """Remove the vehicle and check if the logistic fit is still significant.
       If it is still significant then keep the original. 
       If it is not significant anymore, then keep the removed version.
    Args:
        fit_result (dict): Contains EC50 values for proteins at each probe concentration,
                           probe concentrations are the keys and tuple of best fitted competitive curve fitted lmfit Parameters, 
                           lmfit MinimizerResult best fitted value, the F statistic of fit, and the p value of this fit
        data_by_x (dict): dictionary with catalyst concentration as keys and a nested dictionary containing two lists of values:
                          'y' contains a list of off-compete concentrations and 'z' contain a list of signal intensities matching 'y' conditions
        protein_data (list): contain tuples of 3 values, in the format of (probe concentration, off-compete concentration, signal intensity value)
        global_zmin (float): lowest signal intensity detected
    Returns:
         fit_result (dict): update the tuple of parameters if the p value is no longer significant (<= 0.001) compared to reduced model 
    """    
    for x_value, values in data_by_x.items():
        # Model fitting parameters
        params, result, f_stat, p_value = fit_result[x_value]

        # If the fit is significant (p < 0.001), remove vehicle and check again to test if that is still the case
        if p_value < 0.001:
            # Take out 0 values
            y_values, z_values = remove_correspond_0(values['y'], values['z'])

            # Fit the sigmoid again
            best_param_full, best_result_full, best_f_stat, best_p_value = best_sigmoid(x_value, y_values, z_values, global_zmin)

            # If the sigmoid is significant, then keep the original
            # If the sigmoid is no longer significant, change the keep the removed version
            if best_p_value > 0.001:
                data_by_x[x_value]['y'] = y_values
                data_by_x[x_value]['z'] = z_values
                protein_data = [entry for entry in protein_data if not (entry[0] == x_value and entry[1] == 0)]

                # If no sigmoidal fit:
                if best_result_full is None:
                    nan_array = np.full((3,), np.nan)
                    best_p_value = np.nan
                    best_f_stat = np.nan
                    best_param_full =  nan_array 
                    best_result_full = nan_array

                # Save the new fit for the given catalyst concentration
                fit_result[x_value] = (best_param_full, best_result_full, best_f_stat, best_p_value)
    
    return fit_result

def rearrange_sig_output_for_kd(fit_result):
    """Rearrange data from EC50 fitting for Kd fitting
    Args:
        fit_result (dict): Contain EC50 values for each protein at each probe concentration,
                           probe concentrations are the keys and tuple of best fitted competitive curve fitted lmfit Parameters, 
                           lmfit MinimizerResult best fitted value, the F statistics of fit, and the p value of this fit
    Returns:
        row_x_vals (list): list of probe concentrations
        row_bs (list): list of log10(EC50) values identified for each probe concentration
        row_b_errors (list): list of standard error found when fitting log10(EC50) values
        row_p_vals (list): list of p values defining significance of sigmoidal fits
    """    
    row_bs = []
    row_b_errors = []
    row_x_vals = []
    row_p_vals = []
    for x_val, values in fit_result.items(): 
        params = values[0]
        error_b = params['b'].stderr
        f_stat = values[2]
        p_val = values[3]

        row_x_vals.append(x_val)
        # Check if params has any values, if not append np.nan
        if len(params) == 0:
            row_bs.append(np.nan)
            row_b_errors.append(np.nan)
            row_p_vals.append(np.nan)

        # Check if the b parameter (EC50) is valid
        if (not np.isnan(params['b'].value)):  
            row_bs.append(params['b'].value)
            row_b_errors.append(error_b)
            row_p_vals.append(p_val)
        
        # append np.nan for failed fits 
        else:
            row_bs.append(np.nan)
            row_b_errors.append(np.nan)
            row_p_vals.append(np.nan)

    return row_x_vals,row_bs,row_b_errors,row_p_vals


def get_combined_p(filtered_p_vals, method):
    """Combine p values:
    Args:
        filtered_p_vals (numpy.ndarray): array of p values identified from tests
        method (str): method of p value combination, choose from {'stouffer','fisher'}
    Returns:
        float: combined p-value by defined method
    """    
    # Stouffer's method
    if method == 'stouffer':
        z_scores = norm.ppf(1 - np.array(filtered_p_vals))
        z_combined = np.sum(z_scores) / np.sqrt(len(filtered_p_vals))
        combined_p_value = 1 - norm.cdf(z_combined)

    # Fisher's method
    if method == 'fisher':
        valid_p = ((filtered_p_vals < 1)& (filtered_p_vals > 0)& ~np.isnan(filtered_p_vals))
        p_vals = np.array(filtered_p_vals[valid_p])
        chi_stat = -2 * np.sum(np.log(p_vals))
        combined_p_value = chi2.sf(chi_stat, 2 * len(p_vals))

    return combined_p_value

def kd_by_lin(row_x_vals, row_bs, row_b_errors, row_p_vals, process_two = 'min_x', significant_row_p = 0.001):
    """Determine unmodified ligand affinity (Kd) value for a given protein based on available data.
       Determine valid EC50 values, where EC50 value is a float, standard error for EC50 value fit is between 0 and 2, 
       and significance of fit less than significant_row_p.
       Fit predicted EC50s with Cheng-Prusoff derived equation whenever sufficient data is available. 
    Args:
        row_x_vals (list): list of probe concentrations
        row_bs (list): list of log10(EC50) values identified for each probe concentration
        row_b_errors (list): list of standard errors found when fitting log10(EC50) values
        row_p_vals (list): list of p values defining significance of sigmoidal fits
        process_two (str, optional): Define how to process cases where only two valid EC50 values are available. 
                                     Choose from {'min_x', 'mean_EC50', 'orig'}. Defaults to 'min_x'.
        significant_row_p (float, optional): Significance threshold cutoff for valid fits used to determine EC50. Defaults to 0.001.
    Returns:
        dict: values determined for protien of interest
                'a': determined log10(Kd) value
                'b': determined Kd/Kc value from model_func_2 fitting
                'error_a': standard error of determined log10(Kd) value
                'error_b': standard error of determined Kd/Kc value from model_func_2 fitting
                'combined_p': combined p value from Fisher combination test of sigmoidal model p value
                'index_count': cnumber of EC50 values considered valid for a given protein
                'fit_p_value': p value from F statistic test of model_func_2 fitting
    """    
    # Transform lists to numpy float arrays
    row_fit_results = {}
    row_bs = np.array(row_bs)
    row_b_errors = np.array(row_b_errors)
    row_p_vals = np.array(row_p_vals)
    row_x_vals = np.array(row_x_vals)
    # Convert array to float64, replacing None with np.nan
    row_bs = np.array(row_bs, dtype=np.float64)
    row_b_errors = np.array(row_b_errors, dtype=np.float64)
    np.nan_to_num(row_b_errors, copy=False) 
    row_p_vals = np.array(row_p_vals, dtype=np.float64)

    # Identify valid EC50 fits
    # Only keep indices meeting the following conditions
    valid_indices = (~np.isnan(row_bs) & (row_p_vals < significant_row_p) & (row_b_errors < 2)& (row_b_errors > 0))

    #Slice filtered data 
    filtered_x_vals = row_x_vals[valid_indices]
    filtered_bs_log = row_bs[valid_indices]
    filtered_b_errors_log = row_b_errors[valid_indices]

    # Assuming a Gaussian distribution, report EC50 with maximum likelihood function
    # Define error-weighted combined p values and errors
    error_weights = 1 / filtered_b_errors_log**2
    error_weights[filtered_b_errors_log == 0] = 0
    # Avoid division by zero errors
    error_weights_sum = np.sum(error_weights)
    if error_weights_sum == 0: 
        weighted_mean_b = np.mean(filtered_bs_log) # If weights cannot be determined, take the mean
    else:    
        weighted_mean_b = np.sum(error_weights * filtered_bs_log) / error_weights_sum
    weighted_error_b = np.sqrt(1 / error_weights_sum)

    # Get combined p value
    combined_p_value = get_combined_p(row_p_vals, 'fisher')
    try:
        com_p = (-np.log10(combined_p_value))
    except ValueError as e:
        com_p = 'no_p_value_in_fit'
        print(f"No p values for fit: {e}")
    
    #Define reporter for missing values
    null_value = ''

    if np.sum(valid_indices) == 2:
        # Method: If there are exactly two EC50 values corresponding to lowest catalyst concentration, report value at minimum catalyst concentration
        if process_two == 'min_x':
            lowest_x_idx = np.argmin(filtered_x_vals)
            # Since missing b error values may have been converted to 0 for calculations, correctly report as None in output
            b_error = filtered_b_errors_log[lowest_x_idx]
            if b_error == 0:
                b_error = null_value
            row_fit_results.update({
                'a': filtered_bs_log[lowest_x_idx], 'b': null_value, 'error_a': b_error, 
                'error_b': null_value, 'combined_p': com_p, 'index_count': np.sum(valid_indices), 'fit_p_value': null_value
            })
            return row_fit_results
        
        # Method: report weighted mean EC50 value
        if process_two == 'mean_EC50':
            row_fit_results.update({
                    'a': weighted_mean_b, 'b': null_value, 'error_a': weighted_error_b, 'error_b': null_value, 
                    'combined_p': com_p, 'index_count': np.sum(valid_indices), 'fit_p_value': null_value
                })
            return row_fit_results

    if np.sum(valid_indices) > 1:
        # If at least three points fit data using model_func_2
        try:
            #Set bounds for fitting and define data quality requirements for using points for fitting
            bounds = ([0, 0], [10, 10])
            p0s = [0.1, 0.1]
            popt, pcov = curve_fit(model_func_2, filtered_x_vals, filtered_bs_log, bounds=bounds, p0=p0s, maxfev=100)
            a = np.log10(popt[0])
            perr = np.sqrt(np.diag(pcov))
            std_a = perr[0]
            b = popt[1]
            std_b = perr[1]
            # Assume a normal distribution for the fitted parameter
            CI_95_slope = [round(b-1.96*std_b, 3), round(b+1.96*std_b, 3)]

            # Fit reduced function
            bounds = ([-10], [10])
            p0s = [weighted_mean_b]
            popt_reduced, pcov_reduced = curve_fit(model_reduced_curve, filtered_x_vals, filtered_bs_log, bounds=bounds, p0=p0s, maxfev=100)

        except RuntimeError as e:
            print(f"Fit did not converge: {e}")
            row_fit_results.update({
                'a': weighted_mean_b, 'b': null_value, 'error_a': weighted_error_b, 'error_b': "fit failed RuntimeError",
                'combined_p': com_p, 'index_count': np.sum(valid_indices), 'fit_p_value': null_value
            })
        except TypeError as e:
            print(f'Input vector issue. Input x: {filtered_x_vals}, input b: {filtered_bs_log}. Error: {e}')
            row_fit_results.update({
                'a': weighted_mean_b, 'b': null_value, 'error_a': weighted_error_b, 'error_b': "fit failed TypeError",
                'combined_p': com_p, 'index_count': np.sum(valid_indices), 'fit_p_value': null_value
            })
        except ValueError as e:
            print(f'ValueError: Input x: {filtered_x_vals}, input b: {filtered_bs_log}. Error: {e}')
            row_fit_results.update({
                'a': weighted_mean_b, 'b': null_value, 'error_a': weighted_error_b, 'error_b': "fit failed ValueError",
                'combined_p': com_p, 'index_count': np.sum(valid_indices), 'fit_p_value': null_value
            })
        
        # Determine if the slope is significant
        try:
            # If curve fitting fails:
            if np.all(np.isnan(popt)):
                f_stat, p_value = np.nan, np.nan
            # F-statistic test
            else:
                # Test against the reduced model with 0 = kd/kc
                f_stat, p_value = f_stat_test_curvefit(filtered_x_vals, filtered_bs_log, popt, popt_reduced, model_func_2, model_reduced_curve, 2)
            
            # Slope is significant if p < 0.1
            bad_fit = p_value > 0.1 

            #If the slope is significant, report all relevant parameters
            if bad_fit == False:
                row_fit_results.update({
                    'a': a, 'b': b, 'error_a': std_a, 'error_b': std_b,
                    'combined_p': com_p, 'index_count': np.sum(valid_indices), 'fit_p_value': p_value, 'slope_CI_95%':CI_95_slope
                })
                
            #If the slope is not significant, the affinity of the catalyst conjugate is too small to measure or error is too large. We report the weighted mean average of valid EC50s in such cases.
            if bad_fit == True:
                row_fit_results.update({
                    'a': weighted_mean_b, 'b': null_value, 'error_a': weighted_error_b, 'error_b': std_b, 
                    'combined_p':  com_p, 'index_count': np.sum(valid_indices), 'fit_p_value': p_value, 'slope_CI_95%':CI_95_slope
                })

        except TypeError as e:
            print(f'Significance determination error, b: {b} , standard_error: {std_b}.')
            print(f'error message: {e}')
            row_fit_results.update({
                'a': weighted_mean_b, 'b': null_value, 'error_a': weighted_error_b, 'error_b': null_value,
                'combined_p': com_p, 'index_count': np.sum(valid_indices), 'fit_p_value': null_value
            })

    elif np.sum(valid_indices) == 1:
        # Report the EC50 value as the affinity if it is the only data available. This is correct in the zero slope case, which holds for most protein-ligand pairs.
        if filtered_b_errors_log == 0:
            filtered_b_errors_log = null_value 
        row_fit_results.update({
            'a': (filtered_bs_log), 'b': null_value, 'error_a': filtered_b_errors_log, 'error_b': null_value,
            'combined_p': com_p, 'index_count': np.sum(valid_indices), 'fit_p_value': null_value
        })

    elif np.sum(valid_indices) == 0:
        #If there are no valid points for measurement, report the EC50 with the best associated p value. 
        try:
            index_lowest_p = np.argmin(row_p_vals)
            best_b = row_bs[index_lowest_p] # Report EC50 with best p value

        # When logistic fitting fails:
        except ValueError as e: 
            best_b = null_value 
            print(f"No p values for fit: {e}")
        
        if np.isnan(com_p):
            best_b = null_value 
            print(f" Low quality p values for fits: {row_p_vals}")

        row_fit_results.update({
            'a': best_b, 'b': null_value, 'error_a': null_value, 'error_b': null_value,
            'combined_p': com_p, 'index_count': np.sum(valid_indices), 'fit_p_value': null_value
        })

    return row_fit_results

def fit_dataframe_for_output(data,all_row_fit_results):
    """Append summary data of determined Kd values and metadata during data processing into an output dataframe:
       'log10_kd', 'kd/kc', 'log10_kd_error', 'kd/kc_error', 'combined_p', 'number_of_sigmoids', 'CP_p_value', and 'slope_CI_95%'
    Args:
        data (pandas.DataFrame): input datafile with mass spectrometry signals and metadata
        all_row_fit_results (list): list of all output results for Kd determination and associated metadata 
                                    the list index of genes correpond to the row index in the input dataframe
    Returns:
        pandas.DataFrame: Output dataframe with processed data from all_row_fit_results for each protein appended
                          Add columns: 'log10_kd', 'kd/kc', 'log10_kd_error', 'kd/kc_error', 'combined_p', 'number_of_sigmoids', 'CP_p_value', 'slope_CI_95%'
    """    
    #Copy dataframe for annotation.
    dataframe_processed = data.copy()
    #Create new dataframe containing reported values
    new_columns = ['log10_kd', 'kd/kc', 'log10_kd_error', 'kd/kc_error', 'combined_p', 'number_of_sigmoids', 'CP_p_value', 'slope_CI_95%']
    for col in new_columns:
        dataframe_processed[col] = np.nan
    for idx, fit_result in enumerate(all_row_fit_results):
        if 'a' in fit_result:
            dataframe_processed.at[idx, 'log10_kd'] = fit_result['a']
            dataframe_processed.at[idx, 'log10_kd_error'] = fit_result['error_a']

        if 'b' in fit_result:
            dataframe_processed.at[idx, 'kd/kc'] = fit_result['b']

        if 'error_b' in fit_result:
            dataframe_processed.at[idx, 'kd/kc_error'] = fit_result['error_b']

        if 'combined_p' in fit_result:
            dataframe_processed.at[idx, 'combined_p'] = fit_result['combined_p']

        if 'index_count' in fit_result: 
            dataframe_processed.at[idx, 'number_of_sigmoids'] = fit_result['index_count']
        
        if 'fit_p_value' in fit_result: 
            dataframe_processed.at[idx, 'CP_p_value'] = fit_result['fit_p_value']

        if 'slope_CI_95%' in fit_result:
            dataframe_processed.at[idx, 'slope_CI_95%'] = fit_result['slope_CI_95%']
        else: 
            dataframe_processed.at[idx, 'slope_CI_95%'] = ''

    # Reorder columns
    column_order = new_columns + [col for col in dataframe_processed.columns if col not in new_columns]
    dataframe_processed = dataframe_processed[column_order]

    return dataframe_processed


# Plot -log(p) vs. Kd; , highlighting by adjusted p value (conservative estimate)
highlight_p_cutoff = 1.3 # Define cutoff for highlighting combined_p value of sigmoidal fit
# Define the pastel color mapping
pastel_color_map = {
     'No EC50 pass filter': '#87ae73',  
     '1 EC50 found': '#00fbb0',  
     '2 EC50s found': '#a2bffe',  
     '3 or more EC50s': '#380282', 
 }
# Need to have dataframe with all gene fits values contained within
def plot_logp_vs_kd(sorted_data, x_low = -7.9, x_high =2.9,savefig = False, figname = ''):
    """Plot log(Kd) vs. -log(p) of sigmoidal fit, highlighting points with adjusted p value of sigmoidal fit above defined cutoff.

    Args:
        sorted_data (dict): DataFrame containing all data for plotting, including columns 'log10_kd', 'combined_p', 'log10_kd_error', 'number_of_sigmoids', and 'passes_fdr'
        x_low (float, optional): _description_. Defaults to -7.9.
        x_high (float, optional): _description_. Defaults to 2.9.
        savefig (bool, optional): _description_. Defaults to False.
        figname (str, optional): _description_. Defaults to ''.
    
    Returns:     
        None: Displays the plot of log(Kd) vs. -log(p) of sigmoidal fit, highlighting points with adjusted p value of sigmoidal fit above defined cutoff. 
        Saves the figure if savefig is True.
    """
    # only plot the data if the kd is within the sensible range of detection
    sorted_data = sorted_data[sorted_data['valid_b'] == True]

    # BH adjust of FDR
    highlighted_data = sorted_data[sorted_data['passes_fdr'] == True]
    
    # Replace empty strings with 0 in 'log10_kd_error' so would plot with only numericals
    highlighted_data['log10_kd_error'] = highlighted_data['log10_kd_error'].replace('', 0)
    # Replace all string in column 'combined_p' with 0   
    sorted_data['combined_p'] = sorted_data['combined_p'].replace(to_replace=r'.*', value=0, regex=True)
    # Map the pastel colors to the DataFrame
    highlighted_data['number_of_sigmoids'] = highlighted_data['number_of_sigmoids'].apply(lambda x: '3 or more EC50s' if x >= 3 else x)
    highlighted_data['number_of_sigmoids'] = highlighted_data['number_of_sigmoids'].apply(lambda x: 'No EC50 pass filter' if x == 0 else x)
    highlighted_data['number_of_sigmoids'] = highlighted_data['number_of_sigmoids'].apply(lambda x: '1 EC50 found' if x == 1 else x)
    highlighted_data['number_of_sigmoids'] = highlighted_data['number_of_sigmoids'].apply(lambda x: '2 EC50s found' if x == 2 else x)
    highlighted_data['color'] = highlighted_data['number_of_sigmoids'].map(pastel_color_map)

    # Plot data
    fig, ax = plt.subplots(figsize=(10, 8))
    al = ax.scatter(sorted_data['log10_kd'], sorted_data['combined_p'], color='grey', alpha=0.3, label='All Data')
    # Scatter plot
    sc = ax.scatter(
        highlighted_data['log10_kd'],
        highlighted_data['combined_p'],
        c=highlighted_data['color'],label='Adj P > 0.05'
    )

    # Error bars
    ax.errorbar(
        highlighted_data['log10_kd'],
        highlighted_data['combined_p'],
        xerr=highlighted_data['log10_kd_error'],
        fmt='o',
        ecolor='black',
        alpha=0.5
    )

    if not highlighted_data.empty:
        # Standardize the axis range
        ax.set_xlim(x_low, x_high)  # Add some padding

    ax.set_xlabel('log(Kd)')
    ax.set_ylabel('-Log(p) of Sigmoidal Fit')
    ax.set_title('Plot of log(Kd) and its error vs. -log(P) Value Highlighting Proteins Passing 5% FDR')

    if savefig:
        fig.savefig(figname, bbox_inches='tight',  dpi = 400)