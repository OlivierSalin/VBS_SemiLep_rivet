import glob
import os
from array import array
import ROOT
ROOT.gROOT.SetBatch(True)
import shutil
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import pandas as pd
import numpy as np
import math
import json
import re
import mplhep as hep
hep.style.use([hep.style.ATLAS])
import AtlasStyle
import numpy as np




from optparse import OptionParser
parser = OptionParser()
parser.add_option("--DOCUT", default = "YES")
parser.add_option("--Full_op", default =False)
parser.add_option("--linear", default =True)
parser.add_option("--Name", default = "test")
parser.add_option("--type_MC", default = "Run3")

parser.add_option("--bins", default = 25)
opts, _ = parser.parse_args()





all_ops_SM = ["SM","FM0","FM1","FM2","FM3","FM4","FM5","FM7",
        "FS02","FS1",
        "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]

all_ops_cat = ["FM0","FM2","FS1","FT1","FT5"]
all_ops_cat = ["FM0","FS1"]
processes=["WpZ"]
decays=["llqq"]
order=["QUAD"]
#all_ops_cat = ["SM"]

Complement_path="/Test_script_/"
Complement_path="/Test_bis/test_Run3/"

if not opts.linear:
    base_dir_plot = f"/exp/atlas/salin/ATLAS/VBS_mc/Plots/Plot_aqgc_test/{Complement_path}/bins_{opts.bins}/"
else:
    base_dir_plot = f"/exp/atlas/salin/ATLAS/VBS_mc/Plots/Plot_aqgc_test/{Complement_path}/bins_{opts.bins}/"

if "aqgc" in opts.type_MC or "model" in opts.type_MC:
    all_ops_cat = ["FM0","FM1","FM2","FM8","FM9",
            "FS2","FS1",
            "FT0","FT1","FT5","FT6","FT8",
            "FM1odd","FM2odd","FM3odd","FM6odd","FT2odd","FT4odd","FT6odd"]

all_ops_both= ["FM0","FM2","FS1","FT1","FT5"]
all_ops_both= ["FM0","FM2"]

base_dir = "/exp/atlas/salin/ATLAS/VBS_mc/eft_files/Tables/"

Special_name = "UFO_aqgc"
aQGC_models_name= ["Aqgc_Eboli","aqgc_Eboli_Run3"]




variables_plot= ['merged_Vhad_DeltaEta_Vlep', 'merged_Vhad_DeltaPhi_Vlep', 'merged_Vhad_DeltaR_Vlep', 'merged_VlepVhad_eta', 'merged_VlepVhad_mass', 'merged_VlepVhad_pt']
variables_plot= [ 'merged_VlepVhad_eta', 'merged_VlepVhad_mass', 'merged_VlepVhad_pt']

def plot_with_ratio_from_root_paths(
    root_paths,
    output_dir,
    variable,
    aQGC_models_name,
    desired_num_bins=50,
    tree_name="Merged",
    xlabel="$p_T$(VV) [GeV]",
    ylabel="Events",
    ratio_ylabel="newModel / baseline",
):
    # Prepare histograms
    histograms = {}
    for label, file_path in root_paths.items():
        root_file = ROOT.TFile(file_path, "READ")
        tree = root_file.Get(tree_name)
        if not tree:
            print(f"Error: Tree '{tree_name}' not found in {file_path}. Skipping...")
            continue

        # Determine min/max dynamically
        tree.Draw(f"{variable}>>dummyHist", "", "goff")
        dummy_hist = ROOT.gDirectory.Get("dummyHist")
        min_x = dummy_hist.GetXaxis().GetXmin()
        max_x = dummy_hist.GetXaxis().GetXmax()
        dummy_hist.Delete()

        # Create the histogram
        print(f"Creating histogram for {label} with {desired_num_bins} bins in range ({min_x}, {max_x})")
        print(f"Drawing {variable} from {file_path}")
        hist = ROOT.TH1F(label, label, desired_num_bins, min_x, max_x)
        tree.Draw(f"{variable}>>{label}", "", "goff")
        hist.Scale(1.0 / hist.Integral())  # Normalize histogram to 1
        histograms[label] = hist.Clone()
        histograms[label].SetDirectory(0)
        root_file.Close()

    # Extract data for matplotlib
    bin_edges = np.array([histograms[next(iter(histograms))].GetBinLowEdge(i+1) for i in range(desired_num_bins+1)])
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    data = {label: [histograms[label].GetBinContent(i+1) for i in range(desired_num_bins)] for label in histograms}

    # Main plot and ratio subplot
    fig, (ax_main, ax_ratio) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, figsize=(8, 6), sharex=True)
    #print("data", data)
    # Main plot
    baseline_key = aQGC_models_name[0]
    for label, y_vals in data.items():
        ax_main.errorbar(bin_centers, y_vals, yerr=0, fmt='o', label=label)  # Data points with error bars
        ax_main.step(bin_edges, np.append(y_vals, y_vals[-1]), where='post', label=f"{label} hist")  # Histogram lines
    ax_main.set_yscale("linear")  # Change to linear scale
    ax_main.set_ylabel(ylabel)
    ax_main.set_ylim(bottom=0)  # Ensure y-axis starts at 0
    hep.atlas.text("Internal Simulation", ax=ax_main, loc=1)

    print("baseline_key", baseline_key)
    # Ratio subplot
    baseline = np.array(data[baseline_key])
    for label, y_vals in data.items():
        if label != baseline_key:
            ratio = np.divide(y_vals, baseline, out=np.zeros_like(y_vals), where=baseline!=0)  # Avoid division by zero
            ax_ratio.errorbar(bin_centers, ratio, yerr=0, fmt='o', label=f"{label} / {baseline_key}")
            #ax_ratio.step(bin_edges, np.append(ratio, ratio[-1]), where='post', label=f"{label} / {baseline_key} hist")
    ax_ratio.set_ylim(0.5, 1.5)
    ax_ratio.axhline(1, color="k", linestyle="dashed", linewidth=0.8)
    ax_ratio.set_xlabel(xlabel)
    ax_ratio.set_ylabel(ratio_ylabel)

    # Legends and save
    ax_main.legend()
    ax_ratio.legend(loc="upper right")
    fig.tight_layout()

    # Save plot
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f"{variable}_plot.png")
    fig.savefig(output_path)
    print(f"Plot saved to {output_path}")

def plot_with_ratio_per_operator(
    order_EFT,
    processes,
    decays,
    variables_plot,
    base_path,
    base_dir_plot,
    tree_name="Merged",
    aQGC_models_name=["Aqgc_Eboli", "aqgc_new"],
    special_name="UFO_aqgc",
    bins=50,
):

    for order in order_EFT:  # Loop over operators
        for process in processes:  # Loop over processes
            for decay in decays:
                for op in all_ops_both:   # Loop over decay channels
                    root_paths = {}

                    for model in aQGC_models_name:  # Loop over models
                        folder_name = f"{model}/aqgc/" if model == "aqgc_new" else f"{model}/"
                        if op == "SM":
                            path = f"{base_path}/{folder_name}/2Lepton/{process}_{decay}/{op}_SM/ntuple_rivet.root"
                        else:
                            path = f"{base_path}/{folder_name}/2Lepton/{process}_{decay}/{op}_{order}/ntuple_rivet.root"
                        print(f"Searching for ROOT files in {path}")
                        matches = glob.glob(path)
                        print("Matches:", matches)
                        if matches:
                            root_paths[f"{model}"] = matches[0]
                        else:
                            print(f"No match for operator: {op}, process: {process}, decay: {decay}, model: {model}")

                    if not root_paths:
                        print(f"No ROOT files found for operator: {op}, process: {process}, decay: {decay}. Skipping...")
                        continue

                    # Output directory for the operator
                    output_dir = os.path.join(base_dir_plot, f"{process}_{decay}/{op}/")
                    os.makedirs(output_dir, exist_ok=True)

                    for variable in variables_plot:  # Loop over variables
                        plot_with_ratio_from_root_paths(
                            root_paths=root_paths,
                            output_dir=output_dir,
                            variable=variable,
                            aQGC_models_name=aQGC_models_name,
                            desired_num_bins=bins,
                            tree_name=tree_name,
                            xlabel=f"{variable} [Units]",  # Customize units if needed
                            ylabel="Events",
                            ratio_ylabel="newModel / baseline",
                        )
                        print(f"Completed plot for variable {variable} under operator {op}, process {process}, decay {decay}.")


processes=["WpZ"]
decays=["llqq"]
order_EFT=["QUAD"]
variables_plot = [
    'merged_Vhad_DeltaEta_Vlep',
    'merged_Vhad_DeltaPhi_Vlep',
    'merged_Vhad_DeltaR_Vlep',
    'merged_VlepVhad_eta',
    'merged_VlepVhad_mass',
    'merged_VlepVhad_pt'
]

base_dir = "/exp/atlas/salin/ATLAS/VBS_mc/eft_files/Tables/"
base_path= base_dir + f"{Special_name}/"


plot_with_ratio_per_operator(
    order_EFT=order_EFT,
    processes=processes,
    decays=decays,
    variables_plot=variables_plot,
    base_path=base_path,
    base_dir_plot=base_dir_plot,
    aQGC_models_name=aQGC_models_name,
    bins=20,
)