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
Complement_path="/Test_bis/df/test_Run3_bis/"

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
all_ops_both= ["FM0"]

base_dir = "/exp/atlas/salin/ATLAS/VBS_mc/eft_files/Tables/"

Special_name = "UFO_aqgc"
aQGC_models_name= ["aqgc_Eboli_Run3","aqgc_new"]




variables_plot= ['merged_Vhad_DeltaEta_Vlep', 'merged_Vhad_DeltaPhi_Vlep', 'merged_Vhad_DeltaR_Vlep', 'merged_VlepVhad_eta', 'merged_VlepVhad_mass', 'merged_VlepVhad_pt']
variables_plot= [ 'merged_VlepVhad_eta', 'merged_VlepVhad_mass', 'merged_VlepVhad_pt']

def root_to_dataframe(file_path, tree_name):
    root_file = ROOT.TFile(file_path, "READ")
    tree = root_file.Get(tree_name)
    if not tree:
        print(f"Error: Tree '{tree_name}' not found in {file_path}.")
        return None

    data = {branch.GetName(): [] for branch in tree.GetListOfBranches()}
    for event in tree:
        for branch in data:
            data[branch].append(getattr(event, branch))
    
    root_file.Close()
    return pd.DataFrame(data)

def plot_with_ratio_from_dataframe(
    dataframes,
    output_dir,
    variable,
    aQGC_models_name,
    desired_num_bins=50,
    xlabel="$p_T$(VV) [GeV]",
    ylabel="Events",
    ratio_ylabel="newModel / baseline",
):
    # Prepare histograms
    histograms = {}
    bin_edges = None
    for label, df in dataframes.items():
        hist, bin_edges = np.histogram(df[variable], bins=desired_num_bins)
        hist = hist / np.sum(hist)  # Normalize histogram to 1
        histograms[label] = hist

    bin_centers = 0.5 *(bin_edges[:-1] + bin_edges[1:])

    # Main plot and ratio subplot
    fig, (ax_main, ax_ratio) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, figsize=(8, 6), sharex=True)

    # Main plot
    baseline_key = aQGC_models_name[0]
    for label, hist in histograms.items():
        ax_main.errorbar(bin_centers, hist, yerr=0, fmt='o', label=label)
        ax_main.step(bin_edges, np.append(hist, hist[-1]), where='post', label=f"{label} hist")
    ax_main.set_yscale("linear")
    ax_main.set_ylabel(ylabel)
    ax_main.set_ylim(bottom=0)
    hep.atlas.text("Internal Simulation", ax=ax_main, loc=1)

    # Ratio subplot
    baseline = histograms[baseline_key]
    for label, hist in histograms.items():
        if label != baseline_key:
            ratio = np.divide(hist, baseline, out=np.zeros_like(hist), where=baseline!=0)
            ax_ratio.errorbar(bin_centers, ratio, yerr=0, fmt='o', label=f"{label} / {baseline_key}")
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
    aQGC_models_name=["Aqgc_Eboli", "aqgc_new"],
    special_name="UFO_aqgc",
    bins=50,
):
    for order in order_EFT:
        for process in processes:
            for decay in decays:
                for op in all_ops_both:
                    dataframes = {}
                    i=0
                    for model in aQGC_models_name:
                        folder_name = f"{model}/aqgc/" if model == "aqgc_new" else f"{model}/"
                        if op == "SM":
                            path = f"{base_path}/{folder_name}/2Lepton/{process}_{decay}/{op}_SM/ntuple_rivet.root"
                        else:
                            path = f"{base_path}/{folder_name}/2Lepton/{process}_{decay}/{op}_{order}/ntuple_rivet.root"
                        print(f"Searching for ROOT files in {path}")
                        matches = glob.glob(path)
                        print("Matches:", matches)
                        if matches:
                            df = root_to_dataframe(matches[0], "Merged")
                            print(f"Dataframe for {model} created with {len(df)} entries.")
                            #print(df)
                            #print(df.head())
                            if df is not None:
                                dataframes[f"{model}"] = df
                        else:
                            print(f"No match for operator: {op}, process: {process}, decay: {decay}, model: {model}")
                        i+=1
                    if not dataframes:
                        print(f"No dataframes found for operator: {op}, process: {process}, decay: {decay}. Skipping...")
                        continue

                    output_dir = os.path.join(base_dir_plot, f"{process}_{decay}/{op}/")
                    os.makedirs(output_dir, exist_ok=True)

                    for variable in variables_plot:
                        plot_with_ratio_from_dataframe(
                            dataframes=dataframes,
                            output_dir=output_dir,
                            variable=variable,
                            aQGC_models_name=aQGC_models_name,
                            desired_num_bins=bins,
                            xlabel=f"{variable} [Units]",
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