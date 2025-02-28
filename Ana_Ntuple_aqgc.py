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
import AtlasStyle
import mplhep as hep

#AtlasStyle.SetAtlasStyle()


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
Complement_path="/Compare_model/Run3_agqc_bis/"

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
#all_ops_both= ["FM0","FM2"]

base_dir = "/exp/atlas/salin/ATLAS/VBS_mc/eft_files/Tables/"

Special_name = "UFO_aqgc"
aQGC_models_name= ["Aqgc_Eboli","aqgc_Eboli_Run3"]
aQGC_models_name= ["aqgc_Eboli_Run3","aqgc_new"]

if "aqgc" in opts.type_MC or "model" in opts.type_MC:
    aQGC_models_name= ["aqgc_new"]




variables_plot= ['merged_CS_V_cos_theta','merged_Vhad_DeltaEta_Vlep', 'merged_Vhad_DeltaPhi_Vlep', 
                 'merged_Vhad_DeltaR_Vlep', 'merged_VlepVhad_eta', 
                 'merged_VlepVhad_mass', 'merged_VlepVhad_pt',
                 'merged_tagjets_delta_eta', 'merged_tagjets_mass','merged_tagjets_pt',
                 'merged_Full_pt','merged_fjet_pt',]

categories = {
    "cat1": ["FM0", "FM1", "FM7", ROOT.kRed-2, ROOT.kPink+6],
    "cat2": ["FM2", "FM3", "FM4", "FM5", ROOT.kOrange+10, ROOT.kYellow-4],
    "cat3": ["FT0", "FT1", "FT2", ROOT.kBlue-9, ROOT.kViolet-6],
    "cat4": ["FT5", "FT6", "FT7", ROOT.kAzure, ROOT.kCyan],
    "cat5": ["FT8", "FT9", ROOT.kCyan-10, ROOT.kTeal-6],
    "cat6": ["FS1","FS02", ROOT.kGreen+3, ROOT.kSpring-9],
    "cat7": ["SM", ROOT.kBlack, ROOT.kBlack]
}

color_root=[ROOT.kRed,ROOT.kBlue,ROOT.kGreen,ROOT.kYellow,ROOT.kOrange,ROOT.kPink,ROOT.kViolet,ROOT.kCyan,ROOT.kTeal,ROOT.kAzure,ROOT.kSpring,ROOT.kMagenta,ROOT.kGray,ROOT.kBlack,ROOT.kWhite]

categories_bis = {
    "FM0": [ROOT.kMagenta-4,1],
    "FM0": [ROOT.kRed+1,1],
    "FM1": [ROOT.kPink+1,9],
    "FM7": [ROOT.kPink+6,7],
    "FM2": [ROOT.kRed+1,1],
    "FM3": [ROOT.kOrange+10,9],
    "FM4": [ROOT.kOrange-3,7],
    "FM5": [ROOT.kYellow-9,4],
    "FT0": [ROOT.kBlue-9,1],
    "FT1": [ROOT.kViolet-9,9],
    #"FT1": [ROOT.kAzure,1],
    "FT2": [ROOT.kViolet-6,7],
    "FT5": [ROOT.kAzure,1],
    "FT6": [ROOT.kAzure+7,9],
    "FT7": [ROOT.kCyan-3,7],
    "FT8": [ROOT.kCyan-10,1],
    "FT9": [ROOT.kTeal-6,9],
    "FS1": [ROOT.kGreen+3,1],
    "FS02": [ROOT.kSpring-9,9],
    "SM": [ROOT.kBlack,1]
}


# Function to generate a gradient of colors
def get_gradient_color(color1, color2, steps, step):
    r1 = ROOT.gROOT.GetColor(color1).GetRed()
    g1 = ROOT.gROOT.GetColor(color1).GetGreen()
    b1 = ROOT.gROOT.GetColor(color1).GetBlue()
    r2 = ROOT.gROOT.GetColor(color2).GetRed()
    g2 = ROOT.gROOT.GetColor(color2).GetGreen()
    b2 = ROOT.gROOT.GetColor(color2).GetBlue()
    r = r1 + (r2 - r1) * step / steps
    g = g1 + (g2 - g1) * step / steps
    b = b1 + (b2 - b1) * step / steps
    return ROOT.TColor.GetColor(r, g, b)



def plot_histograms2(desired_num_bins, root_files_info, output_plot,aQGC_models_name,key_name,all_op_plot, tree_name="Merged"):
    # Open the first ROOT file to retrieve the list of parameters (branches)
    first_root_file_path = next(iter(root_files_info.values()))
    first_root_file = ROOT.TFile(first_root_file_path, "READ")
    tree = first_root_file.Get(tree_name)  # Get the TTree named "Merged"
    if not tree:
        print(f"Error: No TTree named 'Merged' found in file {first_root_file_path}.")
        return
    keys = [branch.GetName() for branch in tree.GetListOfBranches() if branch.GetName() != "Event Number"]  # Get the list of branches
    print(keys)
    first_root_file.Close()

    output_plot = output_plot + f'{tree_name}/'

    if not os.path.exists(output_plot):
        os.makedirs(output_plot)

    operator_colors_bis = {}
    operator_style_bis = {}
    for operator, values in categories_bis.items():
        operator_colors_bis[operator] = values[0]
        operator_style_bis[operator] = values[1] 

    for parameter_to_plot in keys:
        if parameter_to_plot in variables_plot:
            for op in all_op_plot:
                # Create a new TCanvas
                canvas_name = f"canvas_{parameter_to_plot}_{op}"
                canvas = ROOT.TCanvas(canvas_name, "Stacked Histograms", 2000, 2000)
                ROOT.gPad.SetRightMargin(0.2)
                param_name= parameter_to_plot.lstrip("merged_")
                process, decay, op__, order_EFT = key_name.split('_')[0:4]
                key_name_=f"{process} {decay} {op} {order_EFT}"
                key_name_title= key_name_
                hs = ROOT.THStack("hs", f"{param_name} for {key_name_title}")
                legend = ROOT.TLegend(0.55, 0.9, 0.85, 0.8)  # Move legend to be more central but within bounds
                legend.SetTextSize(0.03)  # Decrease the font size for the legend
                legend.SetBorderSize(0)
                legend.SetFillStyle(0)
                min_hist = 0
                max_hist = float('-inf')
                max_bin_content = 0
                for legend_name, file_path in root_files_info.items():
                    if op in legend_name:
                        root_file = ROOT.TFile(file_path, "READ")
                        tree = root_file.Get(tree_name)
                        if not tree:
                            print(f"Error: No TTree named 'Merged' found in file {file_path}. Skipping this file.")
                            root_file.Close()
                            continue
                        min_hist = min(min_hist, tree.GetMinimum(parameter_to_plot))
                        max_hist = max(max_hist, tree.GetMaximum(parameter_to_plot))
                        histogram = ROOT.TH1D("histogram_temp", "title", desired_num_bins, tree.GetMinimum(parameter_to_plot), tree.GetMaximum(parameter_to_plot))
                        tree.Draw(parameter_to_plot + ">>histogram_temp", "", "norm")
                        max_bin_content = max(max_bin_content, histogram.GetMaximum())*1.1
                i=0
                for idx, (legend_name, file_path) in enumerate(root_files_info.items()):
                    
                    if op in legend_name:
                        root_file = ROOT.TFile(file_path, "READ")
                        tree = root_file.Get(tree_name)
                        if not tree:
                            print(f"Error: No TTree named 'Merged' found in file {file_path}. Skipping this file.")
                            root_file.Close()
                            continue
                        tolerance = 0.1
                        min_hist_ = min_hist - tolerance * abs(min_hist)
                        max_hist_ = max_hist + tolerance * abs(max_hist)
                        if parameter_to_plot == "merged_ll_mass":
                            min_hist_ = 88
                            max_hist_ = 94
                        histogram = ROOT.TH1D("histogram", "title", desired_num_bins, min_hist_, max_hist_)
                        tree.Draw(parameter_to_plot + ">>histogram", "", "norm")
                        histogram.SetDirectory(0)
                        color = color_root[i]
                        histogram.SetLineColor(color)
                        histogram.SetLineWidth(3)
                        style = operator_style_bis[op]
                        hs.Add(histogram)
                        print(legend_name)
                        process, decay, op, order_EFT = legend_name.split('_')[0:4]
                        if aQGC_models_name[i] in legend_name:
                            print("we are in")
                            print(aQGC_models_name[i])
                            legend.AddEntry(histogram, f"{aQGC_models_name[i]}", "l")
                        else:
                            print(aQGC_models_name[i])
                            print(legend_name)
                            legend.AddEntry(histogram, f"{legend_name}", "l")
                        root_file.Close()
                        i+=1
                if hs.GetNhists() > 0:
                    hs.Draw("nostack HIST")
                    X_axis_param = parameter_to_plot.lstrip("merged_")
                    if "VhadVlep" in X_axis_param:
                        X_axis_param = X_axis_param.replace("VhadVlep", "VV")
                    hs.GetXaxis().SetTitle(X_axis_param)
                    hs.SetTitle(f"{X_axis_param} for {key_name_title}")
                    #hs.SetMaximum(max_bin_content * 1.1)  # Increase y-axis upper limit by 15%
                    hs.GetXaxis().SetNdivisions(505)  # Reduce the number of ticks on the x-axis
                if not opts.linear:
                    canvas.SetLogy()
                    ROOT.gPad.SetLogy()
                legend.Draw()
                canvas.SetCanvasSize(1000, 1000)
                canvas.Update()
                outdir = output_plot + f"/{op}/"
                os.makedirs(outdir, exist_ok=True)
                canvas.SaveAs(outdir + f"{parameter_to_plot}_hist.png")

base_dir = "/exp/atlas/salin/ATLAS/VBS_mc/eft_files/Tables/"
base_path= base_dir + f"{Special_name}/"
for order_EFT in order:
    for process in processes:
        Root_paths = {}

        for decay in decays:
            if opts.Full_op:
                all_op_plot = all_ops_SM
            else:
                all_op_plot = all_ops_both
            for op in all_op_plot:
                for name_model in aQGC_models_name:
                    if name_model == "aqgc_new":
                        Folder_name= f"{name_model}/aqgc/"
                    else:
                        Folder_name= f"{name_model}/"
                    if op=="SM":
                        #path = os.path.join(base_dir, f"{process}_{decay}/mc16_13TeV.*.MGPy8EG_aQGCFM0_{op}_1_{process}_{decay}.merge.EVNT.*/DOCUT_YES/Tables/All_channel/Stats/combined/ntuple_rivet.root")
                        path = base_path+ f"/{Folder_name}/2Lepton/{process}_{decay}/{op}_SM/ntuple_rivet.root"
                    else:
                        #path = os.path.join(base_dir, f"{process}_{decay}/mc16_13TeV.*.MGPy8EG_aQGC{op}_QUAD_1_{process}_{decay}.merge.EVNT.*/DOCUT_YES/Tables/All_channel/Stats/ntuple_rivet.root")
                        path = base_path+ f"/{Folder_name}/2Lepton/{process}_{decay}/{op}_{order_EFT}/ntuple_rivet.root"
                    print("Path: ",path)
                    matches = glob.glob(path)
                    print("Match: ",matches)
                    if not matches:
                        print(f"No match found for operator {op} and process {process}")
                        continue
                    Root_paths[f"{process}_{decay}_{op}_{order_EFT}_{name_model}"] = matches[0]
        #print(Cutflow_pathsn A


        
            # Call the function with the desired parameters
            desired_num_bins = int(opts.bins)  # Replace with your desired number of bins
            key_name = f"{process}_{decay}_{op}_{order_EFT}"
            print("Key name: ",key_name)
            if opts.Full_op:
                output_plot = f"{base_dir_plot}/full_op/{process}_{decay}/"
            else:   
                output_plot = f"{base_dir_plot}/{process}_{decay}/"


            # Check if the directory exists, if not, create it
            if not os.path.exists(output_plot):
                os.makedirs(output_plot)

            #plot_canvas = plot_histograms2(desired_num_bins, Root_paths,output_plot,tree_name="Angle")
            plot_canvas = plot_histograms2(desired_num_bins, Root_paths,output_plot,aQGC_models_name,key_name,all_op_plot,tree_name="Merged")