# Importing necessary libraries
import pandas as pd
import numpy as np
import os
import glob
import scipy as sp
import scipy.optimize
import matplotlib.pyplot as plt
from scipy import stats

# Suppressing division and invalid warnings
np.seterr(divide='ignore', invalid='ignore')

# Global variable to decide whether to read protein fit curves from file
Read_Protein_Fit_Curves = False   


def model_func(t, A, K, C):
    """
    Model function for fitting the data.
    
    Parameters:
    t : independent variable
    A, K, C : parameters of the model
    
    Returns:
    Value of the model function for given parameters
    """
    return A * np.exp(K * t) + C

def prot_standars(data):
    """
    Process protein standards and fit them to the model.
    
    Parameters:
    data : DataFrame containing protein data
    
    Returns:
    DataFrame with fitted curves
    """
    # Subtracting background to get the Fprot value
    data["Fprot"] = data["F670 Mean"] - data["B670 Mean"]
    
    # Setting up directory for saving protein fits
    current_directory = os.getcwd()
    final_directory = os.path.join(current_directory, 'Protein fits')
    
    # If the directory doesn't exist, create it, otherwise clean it
    if not os.path.exists(final_directory):
        os.makedirs(final_directory)
    else:
        fig_files = glob.glob(final_directory + "\\*")
        for f in fig_files:
            os.remove(f)
    
    lysates = data["Lysate code"].dropna().unique()
    curves = pd.DataFrame()
    
    for lysate in lysates:
        cut = data[data["Lysate code"] == lysate][["Dilution ratio", "Fprot"]]
        cut_temp = pd.DataFrame()
        
        for dil in cut["Dilution ratio"].unique():
            valid_data = cut[cut["Dilution ratio"] == dil]
            filtered_data = valid_data[np.abs(stats.zscore(valid_data)["Fprot"]) < 2]
            cut_temp = pd.concat([cut_temp, filtered_data])
        
        cut = cut_temp[cut_temp["Fprot"] > 3000]
        
        # Fitting the data to the model function
        try:
            opt_parms, parm_cov = sp.optimize.curve_fit(model_func, cut["Dilution ratio"], cut["Fprot"], p0 = [30000, -0.5, 15000], maxfev=10000)
            A, K, C = opt_parms
            residuals = cut["Fprot"] - model_func(cut["Dilution ratio"], *opt_parms)
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((cut["Fprot"] - np.mean(cut["Fprot"]))**2)
            r_squared = 1 - (ss_res / ss_tot)
            
            # Plotting the data and the fit
            ax = cut.plot.scatter("Dilution ratio", "Fprot")
            line_x = np.linspace(cut["Dilution ratio"].min(), cut["Dilution ratio"].max(), 100)
            line_y = model_func(line_x, *opt_parms)
            ax.plot(line_x, line_y, color='red')
            ax.set_title('Lysate: %s, R^2: %.2f' % (lysate, r_squared))
            plt.savefig(final_directory + "//" + lysate + ".png")
            plt.close()
            
            # Storing the fit parameters in the curves DataFrame
            curves = pd.concat([curves, pd.DataFrame([{"Lysate": lysate, "A": A, "K": K, "C": C, "R2": r_squared}])], ignore_index=True)
        
        except Exception as e:
            print(f"Error fitting lysate {lysate}: {e}")
    
    curves.index = curves["Lysate"]
    curves.to_csv(final_directory + "//Protein_curves.csv")
    return curves

def normalize_reads(data, standards):
    # Initialize an empty DataFrame for calculations
    calc = pd.DataFrame()
    calcs = []
    
    # Reset the index of the input data
    data = data.reset_index()
    
    # Calculate the Fantibody value by subtracting background
    data["Fantibody"] = data["F785 Mean"] - data["B785 Mean"]
    
    # Get the unique lysate codes
    lysates = data["Lysate code"].dropna().unique()
    
    # Process each lysate
    for lysate in lysates:
        try:
            # Filter the data by lysate code and relevant columns
            cut = data[data["Lysate code"] == lysate][["Dilution ratio", "Fantibody"]]
            
            # Further filter data by Fantibody values
            cut = cut[cut["Fantibody"] > 3000]
            
            # Calculate Fprot values using the model function
            cut["Fprot_calc"] = cut["Dilution ratio"].apply(model_func, args=(standards.loc[lysate]["A"], standards.loc[lysate]["K"], standards.loc[lysate]["C"]))
            
            # Normalize Fantibody values using Fprot_calc
            cut["Fantibody_norm"] = cut["Fantibody"] / cut["Fprot_calc"]
            
            # Filter normalized Fantibody values using z-score
            cut_temp = cut[(np.abs(stats.zscore(cut["Fantibody_norm"])) < 3)]
            cut["Fantibody_norm_filt"] = cut_temp["Fantibody_norm"]
            
            # Drop the dilution ratio column
            cut = cut.drop("Dilution ratio", axis=1)
            
            # Append the filtered and processed data to the calc DataFrame
            calcs.append(cut)
        except Exception:
            print(f"{lysate} is not in protein file.")
    calc = pd.concat(calcs, axis=1)
    
    # Concatenate the original data with the calculated data
    return pd.concat([data, calc], axis=1)

def graph_data(file, graph, data, groups):
    
    # Setting up directory for saving Figures
    current_directory = os.getcwd()
    final_directory = os.path.join(current_directory, 'Figures')
    
    # If the directory doesn't exist, create it
    if not os.path.exists(final_directory):
        os.makedirs(final_directory)
        
    # Define plot colors
    pcolors =  ['r', 'r', 'g', 'g', 'b', 'b', 'y', 'y']
    
    # Create an empty DataFrame for plotting
    plotdata = pd.DataFrame()
    
    # Initialize a new plot
    fig, ax = plt.subplots()
    c = 0
    
    # Process data for each group
    for group in groups.columns:
        cut = data[data["Lysate code"].isin(groups[group])]
        lysates = cut["Lysate code"].dropna().unique()
        cut_mean = pd.DataFrame(index=range(2))
        all_data = []
        # Compute mean and standard deviation for each lysate in the group
        for lysate in lysates:
            mean = cut[cut["Lysate code"] == lysate][["Fantibody_norm_filt"]].mean().values[0]
            std = cut[cut["Lysate code"] == lysate][["Fantibody_norm_filt"]].std().values[0]
            all_data.append(pd.DataFrame({lysate: [mean, std]}))
        cut_mean = pd.DataFrame(index=range(2))
        try:
            cut_mean = pd.concat(all_data, axis=1)
        except Exception:
            pass
        
        # Reformat the DataFrame for plotting
        cut_mean = cut_mean.T.reset_index()
        cut_mean.columns = [group + "_Lysates", group + "_Mean", group + "_std"]
        cut_mean = cut_mean.set_index(group + "_Lysates").reindex(index=groups[group]).reset_index()
        cut_mean.columns = [group + "_Lysates", group + "_Mean", group + "_std"]
        
        # Append the group data to the plot data
        plotdata = pd.concat([plotdata, cut_mean], axis=1)
        plotdata[group + "_x"] = c
        plotdata[group + "_x"] += (np.random.sample(plotdata[group + "_x"].shape[0]) - 0.5) * 0.25
                
        # Plot the group data with error bars
        ax.scatter(plotdata[group + "_x"], plotdata[group + "_Mean"], label=group, alpha=0.5, color=pcolors[c % len(pcolors)])
        ax.errorbar(c, plotdata[group + "_Mean"].mean(), plotdata[group + "_Mean"].std(), color="black", linewidth=3)
        ax.hlines(plotdata[group + "_Mean"].mean(), c - 0.15, c + 0.15, color="black", linewidth=3)
        
        # Connect data points between adjacent groups
        if c % 2 == 1:
            g1 = groups.columns[c - 1]
            g2 = groups.columns[c]
            l1 = plotdata[g1 + "_Lysates"].str.split("_", expand=True)[0].dropna()
            l2 = plotdata[g2 + "_Lysates"].str.split("_", expand=True)[0].dropna()
            for ix1 in l1.index:
                if l1[ix1] in set(l2):
                    ix2 = l2[l2 == l1[ix1]].index[0]
                    ax.plot([plotdata[g1 + "_x"][ix1], plotdata[g2 + "_x"][ix2]], 
                            [plotdata[g1 + "_Mean"][ix1], plotdata[g2 + "_Mean"][ix2]], 
                            color=pcolors[c % len(pcolors)], alpha=0.5, linewidth=1)
        c += 1
    
    # Set plot parameters and save
    plt.xticks(np.arange(c), groups.columns, rotation=45)
    plt.ylabel("fprot")
    plt.savefig("Figures//{}_{}.jpeg".format(file.split(".xlsx")[0], graph.split(".xlsx")[0]), bbox_inches='tight', dpi=600)
    plotdata.to_csv("Figures//{}_{}.csv".format(file.split(".xlsx")[0], graph.split(".xlsx")[0]))


def flip_wells(data):
    data["New well sign"] = (65+80-data["Well sign"].str[:1].apply(ord)).apply(chr) + (25-data["Well sign"].str[1:].astype(int)).astype(str)
    copy = data[data.columns[:8].append(data.columns[22:-1])]

    data = data.set_index(["New well sign", "Well plate barcode"])
    data = data.sort_index()
    data = data[data.columns[8:22]]
    copy = copy.set_index(["Well sign", "Well plate barcode"])
    copy = copy.sort_index()
    

    data = pd.concat([data,copy], axis=1)
    #data.to_excel(file, merge_cells=False)
    return data

def Main():
    # Check if protein fit curves should be read from a file
    if Read_Protein_Fit_Curves:
        standards = pd.read_csv("Protein fits//Protein_curves.csv", index_col=0)
    else:
        # Fit calibration curves if not reading from file
        print("Fitting calibration curves...")
        
        # List all files in the current directory
        files = os.listdir()
        
        # Filter for files related to protein staining
        proteins = [x for x in files if "PROTEIN STAIN.xlsx" in x]
        data = pd.DataFrame()
        
        # Load protein data from each protein file and concatenate
        for protein in proteins:
            pdata = pd.io.excel.read_excel(protein, "Sheet1")
            data = pd.concat([data, pdata])
        
        # Reset index for the combined data
        data.reset_index(drop=True, inplace=True)
        
        # Get standard protein data
        standards = prot_standars(data)
        print("Done!")
    
    # List all files again and filter for specific file types
    files = os.listdir()
    files = [x for x in files if ".xlsx" in x and "PROTEIN STAIN.xlsx" not in x and "_norm.xlsx" not in x and "GROUP.xlsx" not in x]
    
    # Extract unique antibodies from the filenames
    antibodies = list(set([x.split(" ")[-1].split(".xlsx")[0] for x in files]))
    
    # Process data for each antibody
    for ab in antibodies:
        data = pd.DataFrame()
        
        # Get files specific to the current antibody
        abfiles = [x for x in files if ab in x]
        
        # Load data from each file and concatenate
        for file in abfiles:
            print("Processing %s" % file)
            fdata = pd.io.excel.read_excel(file, "Sheet1")
            data = pd.concat([data, fdata])
        
        # Normalize reads using standards
        data = normalize_reads(data, standards)
        
        # Get list of graph files
        graphs = [x for x in os.listdir() if "GROUP.xlsx" in x]
        
        # Process data for each graph
        for graph in graphs:
            groups = pd.io.excel.read_excel(graph, "Sheet1")
            graph_data(file, graph, data, groups)
        
        # Save normalized data to a file
        data.to_excel(file.split(".xlsx")[0] + "_norm.xlsx", merge_cells=False)
    
    print("Done")
    
Main()