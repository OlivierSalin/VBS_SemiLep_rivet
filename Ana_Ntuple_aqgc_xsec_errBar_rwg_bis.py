import glob
import os
from array import array
import ROOT
import utils_func as uf
ROOT.gROOT.SetBatch(True)
import shutil
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import pandas as pd
import numpy as np
import math
import json
import re
#import AtlasStyle
import mplhep as hep

#AtlasStyle.SetAtlasStyle(0)
ROOT.gROOT.LoadMacro("~/ATLAS/atlasrootstyle/AtlasStyle.C")
ROOT.gROOT.LoadMacro("~/ATLAS/atlasrootstyle/AtlasLabels.C")
ROOT.gROOT.LoadMacro("~/ATLAS/atlasrootstyle/AtlasUtils.C")
ROOT.SetAtlasStyle()

ROOT.gStyle.SetLegendBorderSize(0)
ROOT.gStyle.SetLegendFillColor(0)
# ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetLegendTextSize(0.03)

from optparse import OptionParser
parser = OptionParser()
parser.add_option("--DOCUT", default = "YES")
parser.add_option("--All_channel", default = False)
parser.add_option("--Full_op", default =False)
parser.add_option("--linear", default =True)
parser.add_option("--Name", default = "test")
parser.add_option("--type_MC", default = "")
parser.add_option("--EFT_order", default = "QUAD")
parser.add_option("--Channel", default = "")
parser.add_option("--nb_lep", default = 2)
parser.add_option("--bins", default = 25)
parser.add_option("--Xsec", default =False)    
opts, _ = parser.parse_args()

Processes = ["WmZ","WpZ","ZZ","WmWm","WpWm","WpWp"]
Decays = ["llqq",'lvqq','vvqq']

all_ops_SM = ["SM","FM0","FM1","FM2","FM3","FM4","FM5","FM7",
        "FS02","FS1",
        "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]

order_=opts.EFT_order
nb_lepton= int(opts.nb_lep)
name_plt=opts.Name
norm_xsec = bool(opts.Xsec)
print(f'norm_xsec: {norm_xsec}')
#norm_xsec = True if norm_xsec == "True" else False
if norm_xsec==True:
    print("Normalizing to cross section")
else:   
    print("Not normalizing to cross section")

all_ops_cat = ["FM0","FM2","FS1","FT1","FT5"]
all_ops_cat = ["FM0","FS1"]
processes=["WpZ"]
decays=["llqq"]
order=[order_]

if opts.All_channel:
    processes = ["WmZ","WpZ","ZZ","WmWm","WpWm","WpWp"]
    decays = ['lvqq','llqq','vvqq']
else:
    if nb_lepton == 1:
        print("1 lepton")
        processes = ["WmZ","WpZ","WmWm","WpWm","WpWp"]
        decays = ['lvqq']
    elif nb_lepton == 2:
        processes=["WpZ","WmZ","ZZ"]
        decays = ['llqq']
    elif nb_lepton == 0:
        processes=["WpZ","WmZ","ZZ"]
        decays = ['vvqq']
  
channel_name = opts.Channel
channel_split = channel_name.split("_")    
# Generate all possible process_decay combinations
process_decay_combinations = [f"{process}_{decay}" for process in Processes for decay in Decays]

# Check if channel_name is in the combinations
if channel_name in process_decay_combinations:
    print(f"Channel name '{channel_name}' is a valid process_decay combination.")
    processes, decays = [channel_name.split("_")[0]], [channel_name.split("_")[1]]
else:
    print(f"Channel name '{channel_name}' is NOT a valid process_decay combination.")
#all_ops_cat = ["SM"]

Complement_path="/Test_script_/"
Complement_path="Compare_shape_aqgc_noweight/"
#Complement_path="/Models/Stats_100k/Eboli_aqgc_new/"

if not opts.linear:
    base_dir_plot = f"/exp/atlas/salin/ATLAS/VBS_mc/Plots/Plot_aqgc_Model_test/Rwg/Test06/{Complement_path}/bins_{opts.bins}/"
else:
    base_dir_plot = f"/exp/atlas/salin/ATLAS/VBS_mc/Plots/Plot_aqgc_Model_test/Rwg/Rwg_test/{Complement_path}/bins_{opts.bins}/"



all_ops_both= ["FM0","FM1","FM2","FM3","FM4","FM5","FM7","FM8","FM9",]
#all_ops_both= ["FM0"]
#all_ops_both= ["FM0vsFM1"]

if "aqgc" in opts.type_MC or "model" in opts.type_MC:
    all_ops_both = ["FM0","FM1","FM2","FM8","FM9",
            "FS2","FS1",
            "FT0","FT1","FT5","FT6",
            "FM1odd","FM2odd","FM3odd","FM6odd","FT2odd","FT4odd","FT6odd"]

base_dir = "/exp/atlas/salin/ATLAS/VBS_mc/eft_files/Tables/"

Special_name = f"UFO_new_aQGC_test/Cross_term/"
Special_name = f"UFO_aQGC_Model/Stats_100k/"
Special_name = f"Rwg_fm/"

aQGC_models_name= ["Aqgc_Eboli","aqgc_Eboli_Run3"]
aQGC_models_name= ["Eboli_run3","aqgc_model"]
aQGC_models_name= ["rwg_fm_Nohel","aqgc_new"]
aQGC_models_name= ["rwg_fm_Nohel","aqgc_new"]
#aQGC_models_name= ["aqgc_Eboli_Run3","aqgc_Eboli_Run3","aqgc_new"]
Folder_MC= ["/Cross_term/","/Cross_term/",]
Folder_MC= ["Rwg_fm/QUAD/","Rwg_fm/QUAD/"]

if "aqgc" in opts.type_MC or "model" in opts.type_MC:
    aQGC_models_name= ["aqgc_new"]




variables_plot= ['EventWeight','merged_cos_theta_star','merged_Vhad_DeltaEta_Vlep', 'merged_Vhad_DeltaPhi_Vlep', 
                 'merged_Vhad_DeltaR_Vlep', 'merged_VlepVhad_eta', 
                 'merged_VlepVhad_mass', 'merged_VlepVhad_pt',
                 'merged_tagjets_delta_eta', 'merged_tagjets_mass','merged_tagjets_pt',
                 'merged_Full_pt','merged_fjet_pt']

#variables_plot= ['merged_VlepVhad_mass']

categories = {
    "cat1": ["FM0", "FM1", "FM7", ROOT.kRed-2, ROOT.kPink+6],
    "cat2": ["FM2", "FM3", "FM4", "FM5", ROOT.kOrange+10, ROOT.kYellow-4],
    "cat3": ["FT0", "FT1", "FT2", ROOT.kBlue-9, ROOT.kViolet-6],
    "cat4": ["FT5", "FT6", "FT7", ROOT.kAzure, ROOT.kCyan],
    "cat5": ["FT8", "FT9", ROOT.kCyan-10, ROOT.kTeal-6],
    "cat6": ["FS1","FS02", ROOT.kGreen+3, ROOT.kSpring-9],
    "cat7": ["SM", ROOT.kBlack, ROOT.kBlack]
}

def get_cross_section(EFT_op, EFT_type, proc, decay, name_model_):
    keyy = f"{EFT_op}_{EFT_type}_{proc}_{decay}"
    if "aqgc_new" in name_model_ or "model" in name_model_ or "new" in name_model_:
        print("New aQGC model")
        VBS_txt = 'VBS_cross_section_aqgc.txt'
        xsection_fb = uf.take_xsec_fb_aqgc(VBS_txt, EFT_op, EFT_type, proc, decay)
    elif "Eboli_Run3" in name_model_ or  "Run3" in name_model_ or "run3" in name_model_:
        print("Eboli Run 3")
        VBS_txt = 'VBS_cross_section_Eboli_run3.txt'
        xsection_fb = uf.take_xsec_fb_aqgc(VBS_txt, EFT_op, EFT_type, proc, decay)
    else:
        print("Eboli run 2")
        VBS_txt = 'VBS_xsection_test.txt'
        xsection_fb = uf.take_xsec_fb(VBS_txt, EFT_op, EFT_type, proc, decay)
    return xsection_fb

color_root=[ROOT.kRed,ROOT.kBlue,ROOT.kGreen+3,ROOT.kYellow,ROOT.kOrange,ROOT.kPink,ROOT.kViolet,ROOT.kCyan,ROOT.kTeal,ROOT.kAzure,ROOT.kSpring,ROOT.kMagenta,ROOT.kGray,ROOT.kBlack,ROOT.kWhite]


def plot_histograms2(desired_num_bins, root_files_info, output_plot, aQGC_models_name, Process_name, Xsec_model,Xsec_op, normalize_to_xsec=False, tree_name="Merged"):
    def get_histogram(tree, parameter, bins, weight_branch_name, min_val, max_val, scale_factor=1.0):
        print(f"Weight branch name: {weight_branch_name}")
        print(f"Scale factor: {scale_factor}")
        histogram = ROOT.TH1D("histogram_temp", "title", bins, min_val, max_val)
        if weight_branch_name in [branch.GetName() for branch in tree.GetListOfBranches()]:
            tree.Draw(f"{parameter}>>histogram_temp", weight_branch_name, "norm")
        else:
            print(f"Warning: Branch '{weight_branch_name}' not found in the tree. Using default weight.")
            tree.Draw(f"{parameter}>>histogram_temp")
        histogram.Scale(scale_factor)
        histogram.Sumw2()  # Add this line to calculate the sum of squares of weights
        return histogram

    def setup_canvas(parameter, op):
        canvas_name = f"canvas_{parameter}_{op}"
        canvas = ROOT.TCanvas(canvas_name, "Stacked Histograms", 3000, 3000)
        return canvas

    def setup_legend():
        legend = ROOT.TLegend(0.60, 0.78, 0.65, 0.85)
        legend.SetTextSize(0.03)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        return legend

    def get_histogram_range(tree, parameter):
        min_val = tree.GetMinimum(parameter)
        max_val = tree.GetMaximum(parameter)
        return min_val, max_val

    def adjust_histogram_range(min_hist, max_hist, parameter, tolerance=0.1):
        min_hist -= tolerance * abs(min_hist)
        max_hist += tolerance * abs(max_hist)
        return min_hist, max_hist

    def save_canvas(canvas, output_dir, parameter):
        os.makedirs(output_dir, exist_ok=True)
        canvas.SaveAs(output_dir + f"{parameter}_hist.png")
        canvas.SaveAs(output_dir + f"{parameter}_hist.pdf")

    def format_model_name(model_name):
        if "aqgc" in model_name:
            model_name = model_name.replace("aqgc", "aQGC")
        if "Aqgc" in model_name:
            model_name = model_name.replace("Aqgc", "aQGC")
        if "new" in model_name:
            model_name = model_name.replace("new", "new basis")
        return model_name.replace("_", " ")
    
    def format_title_name(title_name):
        if "WpZ" in title_name:
            title_name = title_name.replace("WpZ", "W^{+}Z")
        if "WmZ" in title_name:
            title_name = title_name.replace("WmZ", "W^{-}Z")
        if "llqq" in title_name:
            title_name = title_name.replace("llqq", "#rightarrow llqq")
        return title_name
    
    def format_param_name(param_name):
        param_name = param_name.lstrip("merged_").replace("VlepVhad", "VV").replace("_", " ")
        if param_name == "VV mass":
            param_name = "m_{VV} (GeV)"
        elif param_name == "CS V cos theta":
            param_name = "cos#theta_{CS}"
        elif "pt" in param_name:
            param_name = "p_{T}(" + param_name.replace(" pt", "") + ")"
        elif "eta" in param_name and "Delta" not in param_name:
            param_name = "#eta(" + param_name.replace("eta", "") + ")"
        elif "DeltaEta" in param_name:
            param_name = "#Delta#eta(" + param_name.replace("DeltaEta", "") + ")"
        elif "DeltaPhi" in param_name:
            param_name = "#Delta#phi(" + param_name.replace("DeltaPhi", "") + ")"
        return param_name
    
    def format_Y_axis(X_param, nb_bins, max_hist):
        per_bins = max_hist / nb_bins
        rounded_per_bins = round(per_bins, 2)
        if "pt" in X_param or "mass" in X_param:
            return f"Normalized Entry / {rounded_per_bins} (GeV)"
        else:
            return f"Normalized Entry / {rounded_per_bins}"

    first_root_file_path = next(iter(root_files_info.values()))
    first_root_file = ROOT.TFile(first_root_file_path, "READ")
    tree = first_root_file.Get(tree_name)
    if not tree:
        print(f"Error: No TTree named '{tree_name}' found in file {first_root_file_path}.")
        return
    keys = [branch.GetName() for branch in tree.GetListOfBranches() if branch.GetName() != "Event Number"]
    print(f"Branches in the tree: {keys}")
    first_root_file.Close()

    #output_plot = output_plot + f'{tree_name}/'
    os.makedirs(output_plot, exist_ok=True)

    for parameter_to_plot in keys:
        if parameter_to_plot not in variables_plot:
            continue
        process, decay, op__, order_EFT = Process_name.split('_')[0:4]    
        canvas = setup_canvas(parameter_to_plot, op__)
        param_name = format_param_name(parameter_to_plot)
        Process_name_title = f"{process} {decay} {op__} {order_EFT}"
        ecm = "13.6 TeV"
        sample_legend = f"{process} {decay} {op__} {order_EFT} {ecm} (Merged)"
        legend_title= format_title_name(sample_legend)
        print(f"Plotting {param_name} for {Process_name_title}")
        hs = ROOT.THStack("hs", f"{param_name} for {Process_name_title}")
        legend = setup_legend()
        legend_bis = ROOT.TLegend(0.20, 0.85, 0.55, 0.90)
        legend_bis.SetHeader(legend_title)
        min_hist = 0
        max_hist = float('-inf')
        max_bin_content = 0

        for legend_name, file_path in root_files_info.items():
            root_file = ROOT.TFile(file_path, "READ")
            model_naming_=legend_name.split('_')[4:]
            tree = root_file.Get(tree_name)
            if not tree:
                print(f"Error: No TTree named '{tree_name}' found in file {file_path}. Skipping this file.")
                root_file.Close()
                continue
            min_val, max_val = get_histogram_range(tree, parameter_to_plot)
            min_hist = min(min_hist, min_val)
            max_hist = max(max_hist, max_val)            

            scale_factor = Xsec_model[legend_name] if normalize_to_xsec else 1.0
            weight_branch_name = "EventWeight"

            print(f"Model naming: {model_naming_}")
            if 'rwg' in model_naming_ or 'reweight' in model_naming_ or "Rwg" in model_naming_:
                scale_factor = Xsec_op[Process_name] if normalize_to_xsec else 1.0
                op_order_strg = op__ + "_" + order_EFT
                op_order_strg_rwg = op_order_strg.lower().replace("vs", "_")
                #weight_branch_name = f"EventWeight_{op_order_strg_rwg}"

            histogram = get_histogram(tree, parameter_to_plot, desired_num_bins,weight_branch_name, min_val, max_val, scale_factor)
            print(f"Max bin content: {histogram.GetMaximum()}")
            if histogram.GetMaximum() < 1e-6:
                max_bin_content = max(max_bin_content, histogram.GetMaximum()) + abs(histogram.GetMinimum()) * 0.2

            if histogram.GetMaximum() > 0.4:
                
                max_bin_content = max(max_bin_content, histogram.GetMaximum()) * 1.25
                
                
            else:
                max_bin_content = max(max_bin_content, histogram.GetMaximum()) * 1.3
            #print(f"Max bin content: {max_bin_content}")
            root_file.Close()

        i = 0
        for idx, (legend_name, file_path) in enumerate(root_files_info.items()):
            process, decay, op, order_EFT = legend_name.split('_')[0:4]
            model_naming_=legend_name.split('_')[4:]    
            #print(f"Legend: Plotting {param_name} for {process} {decay} {op} {order_EFT}")
            root_file = ROOT.TFile(file_path, "READ")
            tree = root_file.Get(tree_name)
            if not tree:
                #print(f"Error: No TTree named '{tree_name}' found in file {file_path}. Skipping this file.")
                root_file.Close()
                continue
            min_hist_, max_hist_ = adjust_histogram_range(min_hist, max_hist, parameter_to_plot)
            if parameter_to_plot == "merged_ll_mass":
                min_hist_, max_hist_ = 88, 94
            

            scale_factor = Xsec_model[legend_name] if normalize_to_xsec else 1.0
            print(f"Scale factor: {scale_factor}")
            histogram = ROOT.TH1D("histogram", "title", desired_num_bins, min_hist_, max_hist_)
            weight_branch_name_bis = "EventWeight"

            if 'rwg' in model_naming_ or 'reweight' in model_naming_ or "Rwg" in model_naming_:
                print("Rwg model")
                scale_factor = Xsec_op[Process_name] if normalize_to_xsec else 1.0
                op_order_strg = op__ + "_" + order_EFT
                op_order_strg_rwg = op_order_strg.lower().replace("vs", "_")
                #weight_branch_name_bis = f"EventWeight_{op_order_strg_rwg}"
        
            if weight_branch_name_bis in [branch.GetName() for branch in tree.GetListOfBranches()]:
                print(f"Model name: {model_naming_}")
                print(f"Weight branch name in if: {weight_branch_name_bis}")
                tree.Draw(f"{parameter_to_plot}>>histogram", weight_branch_name_bis, "norm")
            else:
                print(f"Model name in else: {model_naming_}")
                print(f"Weight branch name in else: {weight_branch_name_bis}")
                print(f"Warning: Branch '{weight_branch_name_bis}' not found in the tree. Using default weight.")
                tree.Draw(f"{parameter_to_plot}>>histogram")
            histogram.Scale(scale_factor)
            histogram.Sumw2()  # Add this line to calculate the sum of squares of weights
            histogram.SetDirectory(0)
            color = color_root[i]
            histogram.SetLineColor(color)
            histogram.SetLineWidth(3)
            hs.Add(histogram)
            formatted_model_name = format_model_name(aQGC_models_name[i])
            if i%2==0:
                legend.AddEntry(histogram, f"{formatted_model_name}", "l")
            else:
                legend.AddEntry(histogram, f"{formatted_model_name}", "l")
            root_file.Close()
            i += 1

        if hs.GetNhists() > 0:
            hs.Draw("nostack HIST E")  # Add "E" to draw error bars
            X_axis_param = format_param_name(parameter_to_plot)
            #print(f"X axis: {X_axis_param}")
            hs.GetXaxis().SetTitle(X_axis_param)
            hs.SetTitle(f"{X_axis_param} for {Process_name_title}")
            hs.SetMaximum(max_bin_content)
            hs.GetXaxis().SetNdivisions(505)
            y_axis_label = format_Y_axis(X_axis_param, desired_num_bins, max_hist)
            y_axis_label = ""
            hs.GetYaxis().SetTitle(y_axis_label)
            if not opts.linear:
                canvas.SetLogy()
                ROOT.gPad.SetLogy()
            legend_bis.Draw()
            legend.Draw()
            ROOT.ATLASLabel(0.20, 0.9, "Work in Progress")
            canvas.SetCanvasSize(1000, 1000)
            canvas.Update()
            os.makedirs(output_plot + f"/{op}_{order_EFT}/", exist_ok=True)
            save_canvas(canvas, output_plot + f"/{op}_{order_EFT}/", parameter_to_plot)

base_dir = "/exp/atlas/salin/ATLAS/VBS_mc/eft_files/Tables/"
base_path = base_dir + f"{Special_name}/"


Root_paths = {}
X_section= {}
for process in processes:
    for decay in decays:
        all_op_plot = all_ops_SM if opts.Full_op else all_ops_both
        for op in all_op_plot:
            Root_paths_model = {}
            X_section_model = {}
            for order_EFT in order:
                for index, name_model in enumerate(aQGC_models_name):
                    name_model_ = name_model
                    if aQGC_models_name[0] == aQGC_models_name[1]:
                        base_path = base_dir + f"{Folder_MC[index]}/"
                        name_model_= f"{name_model}_{index}"


                    
                    Folder_name = f"{name_model}/" if name_model == "aqgc_new" else f"{name_model}/"
                    base_path_ = base_path + f"{order_EFT}/"
                    #base_path_ = base_path 
                    path = (base_path_ + f"/{Folder_name}/2Lepton/{process}_{decay}/{op}_SM/ntuple_rivet.root"
                            if op == "SM" else
                            #base_path_ + f"/{Folder_name}/{op}_{order_EFT}/ntuple_rivet.root")
                            base_path_ + f"/{Folder_name}/2Lepton/{process}_{decay}/{op}_{order_EFT}/ntuple_rivet.root")
                    if 'rwg' in name_model_ or 'reweight' in name_model_ or 'wg' in name_model_:
                        path = base_path_ + f"/{Folder_name}/2Lepton/{process}_{decay}/ntuple_rivet.root"
                    print("Path: ", path)
                    matches = glob.glob(path)
                    print("Match: ", matches)
                    if not matches:
                        print(f"No match found for operator {op} and process {process}")
                        continue
                    print(f"Match found for operator {op} and process {process}: {matches[0]}")
                    
                    xsec = get_cross_section(op, order_EFT, process, decay, "aqgc_new")
                    print(f"Cross section for {op} and {process}: {xsec} fb")
                    if xsec is not None:
                        Root_paths[f"{process}_{decay}_{op}_{order_EFT}_{name_model_}"] = matches[0]
                        X_section[f"{process}_{decay}_{op}_{order_EFT}"] = xsec
                        Root_paths_model[f"{process}_{decay}_{op}_{order_EFT}_{name_model_}"] = matches[0]
                        X_section_model[f"{process}_{decay}_{op}_{order_EFT}_{name_model_}"] = xsec

                num_bins = int(opts.bins)
                Process_name = f"{process}_{decay}_{op}_{order_EFT}"
                print("Key name: ", Process_name)
                outPlot = (f"{base_dir_plot}/full_op/{process}_{decay}/"
                                if opts.Full_op else
                                f"{base_dir_plot}/{process}_{decay}_rwg/")
                os.makedirs(outPlot, exist_ok=True)
                plot_histograms2(num_bins, Root_paths_model, outPlot, aQGC_models_name, Process_name, X_section_model,X_section, norm_xsec, tree_name="Merged")

