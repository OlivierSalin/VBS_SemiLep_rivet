import glob
import os
from array import array
import ROOT
import shutil
import matplotlib.pyplot as plt
import numpy as np
import uproot
from matplotlib.colors import LinearSegmentedColormap

plt.rcParams['text.usetex'] = True


from optparse import OptionParser
parser = OptionParser()
parser.add_option("--Ana", default = "WpZ_llqq")
parser.add_option("--DOCUT", default = "YES")
parser.add_option("--SM", default = "YES")
parser.add_option("--QUAD1", default = "FS1")
parser.add_option("--QUAD2", default = "FT1")
parser.add_option("--bins", default = 50)
opts, _ = parser.parse_args()

all_ops_SM = ["SM","FM0","FM1","FM2","FM3","FM4","FM5","FM7",
        "FS02","FS1",
        "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]

base_dir_plot = "/exp/atlas/salin/ATLAS/VBS_mc/EFT_files_AMI/Plots_op/python/"
# Define the color gradient for each category
categories = {
    "SM": ["violet", "magenta"], 
    "FM": ["red", "orange"], 
    "FS": ["blue", "cyan"],  
    "FT": ["green", "teal"]  
}

# Function to generate a gradient of colors
def get_gradient_color(color1, color2, steps, step):
    cmap = LinearSegmentedColormap.from_list("", [color1, color2])
    return cmap(step/steps)

def plot_histograms(desired_num_bins, root_files_info, output_plot):
    # Create a dictionary to map each operator to a color
    operator_colors = {}
    for operator in all_ops_SM:
        category = operator[:2]  # Get the category from the operator name
        operators_in_category = [op for op in all_ops_SM if op.startswith(category)]  # Get all operators in this category
        color_index = operators_in_category.index(operator)  # Get the index of the operator within its category
        color1, color2 = categories[category]  # Get the colors for this category
        color = get_gradient_color(color1, color2, len(operators_in_category), color_index)  # Create a gradient color
        operator_colors[operator] = color  # Store the color

    # Open the first ROOT file to retrieve the list of parameters (branches)
    first_root_file_path = next(iter(root_files_info.values()))
    first_root_file = uproot.open(first_root_file_path)
    keys = first_root_file.keys()

    for parameter_to_plot in keys:    
        # Create a new figure
        plt.figure(figsize=(10, 10))

        # Loop over the files to retrieve and stack the histograms
        for legend_name, file_path in root_files_info.items():
            # Open the ROOT file
            root_file = uproot.open(file_path)
        
            # Retrieve the histogram
            histogram = root_file[parameter_to_plot]
            counts, bins = histogram.to_numpy()

            # Rebin the histogram if necessary
            if len(bins) > desired_num_bins:
                counts, bins = np.histogram(bins[:-1], bins=desired_num_bins, weights=counts)
                
            op = legend_name.split("_")[-1]
            color = operator_colors[op]  # Get the color for this operator

            # Plot the histogram with larger markers
            bin_centers = 0.5 * (bins[:-1] + bins[1:])
            plt.plot(bin_centers, counts, marker='o', markersize=8, linestyle='', color=color, label=legend_name)

        # Set log scale for y-axis
        plt.yscale('log')

        # Set labels
        plt.xlabel(parameter_to_plot, fontsize=14)
        plt.ylabel("Events / bin", fontsize=14)

        # Set title
        plt.title("Histogram of " + parameter_to_plot, fontsize=16)

        # Add gridlines
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)

        # Increase legend font size
        plt.legend(fontsize=12)

        # Save the figure to a file
        plt.savefig(output_plot + parameter_to_plot + "_hist.png")
        plt.close()

# Call the function with Root_paths as argument

base_dir = "/exp/atlas/salin/ATLAS/VBS_mc/EFT_files_AMI/"

all_ops_SM = ["SM","FM0","FM1","FM2","FM3","FM4","FM5","FM7",
        "FS02","FS1",
        "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]



for process in ["WmZ"]:
    Root_paths = {}

    for decay in ["llqq"]:
        for op in all_ops_SM:
            if op=="SM":
                path = os.path.join(base_dir, f"{process}_{decay}/mc16_13TeV.*.MGPy8EG_aQGCFM0_{op}_1_{process}_{decay}.merge.EVNT.*/DOCUT_YES/Tables/Tables_file_first/hists.root")
            else:
                path = os.path.join(base_dir, f"{process}_{decay}/mc16_13TeV.*.MGPy8EG_aQGC{op}_QUAD_1_{process}_{decay}.merge.EVNT.*/DOCUT_YES/Tables/Tables_file_first/hists.root")
            matches = glob.glob(path)
            #print(matches)
            if not matches:
                print(f"No match found for operator {op} and process {process}")
                continue
            Root_paths[f"{process}_{decay}_{op}"] = matches[0]
    #print(Cutflow_paths)



    # Call the function with the desired parameters
    desired_num_bins = int(opts.bins)  # Replace with your desired number of bins

    output_plot = f"{base_dir_plot}/{process}/"

    # Check if the directory exists, if not, create it
    if not os.path.exists(output_plot):
        os.makedirs(output_plot)

    plot_canvas = plot_histograms(desired_num_bins, Root_paths,output_plot)



    