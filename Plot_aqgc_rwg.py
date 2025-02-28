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
parser.add_option("--Xsec", default = False)    
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
if norm_xsec:
    print("Normalizing to cross section")
else:   
    print("Not normalizing to cross section")
norm_xsec = False


def set_processes_and_decays(opts, Processes, Decays):
    if opts.All_channel:
        return Processes, Decays

    nb_lepton = int(opts.nb_lep)
    if nb_lepton == 1:
        print("1 lepton")
        return ["WmZ", "WpZ", "WmWm", "WpWm", "WpWp"], ['lvqq']
    elif nb_lepton == 2:
        return ["WpZ", "WmZ", "ZZ"], ['llqq']
    elif nb_lepton == 0:
        return ["WpZ", "WmZ", "ZZ"], ['vvqq']
    return Processes, Decays

def validate_channel_name(channel_name, Processes, Decays):
    process_decay_combinations = [f"{process}_{decay}" for process in Processes for decay in Decays]
    if channel_name in process_decay_combinations:
        print(f"Channel name '{channel_name}' is a valid process_decay combination.")
        return [channel_name.split("_")[0]], [channel_name.split("_")[1]]
    else:
        print(f"Channel name '{channel_name}' is NOT a valid process_decay combination.")
        return None, None

processes, decays = set_processes_and_decays(opts, Processes, Decays)
channel_name = opts.Channel
channel_processes, channel_decays = validate_channel_name(channel_name, Processes, Decays)
if channel_processes and channel_decays:
    processes, decays = channel_processes, channel_decays


###########################################


processes=["WpZ"]
decays=["llqq"]
order=["QUAD","CROSS"]


Complement_path="Plot_rweighting_study/S_M_op_hel_aware/Negative_weight_cut/"

base_dir_plot = (f"/exp/atlas/salin/ATLAS/VBS_mc/Plots/Plot_Reweighting/{Complement_path}/bins_{opts.bins}/"
                 if not opts.linear else
                 f"/exp/atlas/salin/ATLAS/VBS_mc/Plots/Plot_Reweighting/{Complement_path}/bins_{opts.bins}/")

all_ops_both= ["FM0","FM1","FM2","FM3","FM4","FM5","FM7","FM8","FM9","FS0","FS1","FS2"]
#all_ops_both= ["FT0","FT1","FT2","FT5","FT6","FT7","FS0","FS1","FS2"]
#all_ops_both= ["FT0"]
#all_ops_both= ["FM0","FS1","FM2"]

base_dir = "/exp/atlas/salin/ATLAS/VBS_mc/eft_files/Histograms/"
Special_name = "/Reweighting/STM_op_50k/"

aQGC_models_name = ["rwg_fmfs_Nohel","rwg_fmfs_hel_aware","aqgc"]

if "aqgc" in opts.type_MC or "model" in opts.type_MC:
    aQGC_models_name = ["aqgc_new"]

variables_plot = ['EventWeight','merged_cos_theta_star','merged_Vhad_DeltaEta_Vlep', 'merged_Vhad_DeltaPhi_Vlep', 
                 'merged_Vhad_DeltaR_Vlep', 'merged_VlepVhad_eta', 
                 'merged_VlepVhad_mass', 'merged_VlepVhad_pt',
                 'merged_tagjets_delta_eta', 'merged_tagjets_mass','merged_tagjets_pt',
                 'merged_Full_pt','merged_fjet_pt']

#variables_plot= ['merged_VlepVhad_mass','merged_fjet_pt']




###############################################
weight_cut = {
    "fm0_quad":0.00001,
    "fm0_fm1_cross":0.0006,
    "fm1_quad":0.0005,
    "fm2_quad":0.000005,
    "fm8_quad":0.00001,
    "fs0_quad":2e-7,
    "fs1_quad":5e-7
    }

categories_bis = {
    "FM0": [ROOT.kMagenta-4,1],
    #"FM0": [ROOT.kRed+1,1],
    "FM1": [ROOT.kPink+1,9],
    "FM7": [ROOT.kPink+6,7],
    "FM2": [ROOT.kRed+1,1],
    "FM3": [ROOT.kOrange+10,9],
    "FM4": [ROOT.kOrange-3,7],
    "FM5": [ROOT.kYellow-9,4],
    "FM8": [ROOT.kYellow-6,4],
    "FM9": [ROOT.kYellow-3,4],
    "FT0": [ROOT.kBlue-9,1],
    "FT1": [ROOT.kViolet-9,9],
    "FM0vsFM1": [ROOT.kAzure,1],
    #"FT1": [ROOT.kAzure,1],
    "FT2": [ROOT.kViolet-6,7],
    "FT5": [ROOT.kAzure,1],
    "FT6": [ROOT.kAzure+7,9],
    "FT7": [ROOT.kCyan-3,7],
    "FT8": [ROOT.kCyan-10,1],
    "FT9": [ROOT.kTeal-6,9],
    "FS1": [ROOT.kGreen+3,1],
    "FS0": [ROOT.kSpring-9,9],
    "FS2": [ROOT.kGreen+7,9],
    "SM": [ROOT.kBlack,1]
}

def get_cross_section(EFT_op, EFT_type, proc, decay, name_model_):
    if "aqgc_new" in name_model_ or "model" in name_model_ or "new" in name_model_:
        VBS_txt = 'VBS_cross_section_aqgc.txt'

    elif "Eboli_Run3" in name_model_ or  "Run3" in name_model_ or "run3" in name_model_:
        VBS_txt = 'VBS_cross_section_Eboli_run3.txt'

    elif "run2" in name_model_ or "Run2" in name_model_:
        VBS_txt = 'VBS_xsection_test.txt'

    else:
        VBS_txt = 'VBS_cross_section_aqgc.txt'
    xsection_fb = uf.take_xsec_fb_aqgc(VBS_txt, EFT_op, EFT_type, proc, decay)
    return xsection_fb

color_root = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+3, ROOT.kYellow, ROOT.kOrange, ROOT.kPink, ROOT.kViolet, ROOT.kCyan, ROOT.kTeal, ROOT.kAzure, ROOT.kSpring, ROOT.kMagenta, ROOT.kGray, ROOT.kBlack, ROOT.kWhite]

def get_histogram_limits(root_files_info, parameter_to_plot, tree_name, Xsec_model, Process_name, norm_to_xsec):
    min_hist = 0
    max_hist = float('-inf')
    max_bin_content = 0

    for legend_name, file_path in root_files_info.items():
        process, decay, op, order_EFT = legend_name.split('_')[0:4]
        Process_op_order = f"{process}_{decay}_{op}_{order_EFT}"
        root_file = ROOT.TFile(file_path, "READ")
        tree = root_file.Get(tree_name)
        if not tree:
            root_file.Close()
            continue
        min_val, max_val = tree.GetMinimum(parameter_to_plot), tree.GetMaximum(parameter_to_plot)
        min_hist, max_hist = min(min_hist, min_val), max(max_hist, max_val)

        scale_factor = Xsec_model[legend_name] if norm_to_xsec else 1.0
        weight_branch_name = "EventWeight"

        if any(keyword in legend_name.split("_")[4:] for keyword in ['rwg', 'reweight', 'wg']):
            scale_factor = Xsec_model[Process_op_order] if norm_to_xsec else 1.0
            op_order_strg = op + "_" + order_EFT
            op_order_strg_rwg = op_order_strg.lower().replace("vs", "_")
            weight_branch_name = f"EventWeight_{op_order_strg_rwg}"

        histogram = ROOT.TH1D("histogram_temp", "title", 100, min_val, max_val)
        draw_option = weight_branch_name if weight_branch_name in [branch.GetName() for branch in tree.GetListOfBranches()] else ""
        tree.Draw(f"{parameter_to_plot}>>histogram_temp", draw_option)
        histogram.Scale(scale_factor)
        if histogram.GetEntries() > 0:
            histogram.Sumw2()
        max_bin_content = max(max_bin_content, histogram.GetMaximum() * 1.3)
        root_file.Close()

    return min_hist, max_hist

def calculate_event_loss_fraction(tree, weight_branch_name, weight_cut=None,log_fraction_lost=None):
    total_events = tree.GetEntries()
    if weight_cut:
        cut_events = tree.GetEntries(f"{weight_branch_name} < {weight_cut}")
    else:
        cut_events = total_events
    fraction_lost = (total_events - cut_events) / total_events if total_events > 0 else 0
    log_fraction_lost += f"\nWeight branch name: {weight_branch_name}, Weight cut: {weight_cut}, Total events: {total_events}, Events after cut: {cut_events}, Fraction lost: {fraction_lost:.2%}\n"
    print(log_fraction_lost)
    print(f"Total events: {total_events}, Events after cut: {cut_events}, Fraction lost: {fraction_lost:.2%}")

    return fraction_lost

def get_histogram(tree,model_name, parameter, bins, weight_branch_name, min_val, max_val, scale_factor=1.0, weight_cut=None,log_fraction_lost=None):
    print(parameter)
    if parameter == "merged_VlepVhad_mass":
        calculate_event_loss_fraction(tree, weight_branch_name, weight_cut,log_fraction_lost)
    
    histogram = ROOT.TH1D("histogram", "title", bins, min_val, max_val)
    draw_option = weight_branch_name if weight_branch_name in [branch.GetName() for branch in tree.GetListOfBranches()] else ""
    cut_condition = f"{weight_branch_name} < {weight_cut}" if weight_cut else f"{weight_branch_name}"
    if "hel_aware" in model_name:
        tree.Draw(f"{parameter}>>histogram", f"(({weight_branch_name} < {weight_cut})&&({weight_branch_name}>0))*{weight_branch_name}", "norm")
        #tree.Draw(f"{parameter}>>histogram", draw_option, "norm")
    else:
        tree.Draw(f"{parameter}>>histogram", draw_option, "norm")
    histogram.Scale(scale_factor)
    if histogram.GetEntries() > 0:
        histogram.Sumw2()
    return histogram

def plot_histograms2(desired_num_bins, root_files_info, output_plot, Process_name, Xsec_model, norm_to_xsec=False, tree_name="Merged", weight_cut={"fm0_quad":0.002}):
    
    def setup_canvas(parameter):
        canvas = ROOT.TCanvas(f"canvas_{parameter}", "Stacked Histograms", 3000, 3000)
        #canvas.SetRightMargin(0.25)  # Adjust right margin to make space for the legend
        return canvas

    def setup_legend(num_entry):
        length = 0.01 * num_entry

        legend = ROOT.TLegend(0.60, 0.78 - length, 0.90, 0.85)
        legend.SetTextSize(0.017)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        return legend

    def adjust_histogram_range(min_hist, max_hist, tolerance=0.1):
        return min_hist - tolerance * abs(min_hist), max_hist + tolerance * abs(max_hist)

    def save_canvas(canvas, output_dir, parameter):
        os.makedirs(output_dir, exist_ok=True)
        canvas.SaveAs(output_dir + f"{parameter}_hist.png")
        canvas.SaveAs(output_dir + f"{parameter}_hist.pdf")

    def format_model_name(model_name):
        if "aqgc" in model_name:
            model_name= model_name.replace("aqgc", "aQGC").replace("Aqgc", "aQGC").replace("new", "new basis").replace("_", " ")
        elif "rwg" in model_name or "reweight" in model_name or "rwg" in model_name:
            model_name= model_name.replace("rwg", "Rwg").replace("_", " ").replace("fs", "").replace("fm", "").replace("ft", "")
            model_name= model_name.replace("Nohel", "Hel ignorant")
            model_name= model_name.replace("hel_aware", "Hel aware")
        return model_name

    def format_title_name(title_name):
        return title_name.replace("WpZ", "W^{+}Z").replace("WmZ", "W^{-}Z").replace("llqq", "#rightarrow llqq")

    def format_param_name(param_name):
        param_name = param_name.lstrip("merged_").replace("VlepVhad", "VV").replace("_", " ")
        if param_name == "VV mass":
            return "m_{VV} (GeV)"
        elif param_name == "CS V cos theta":
            return "cos#theta_{CS}"
        elif "phi" in param_name:
            return param_name.replace("phi", "#phi")
        elif param_name == "cos theta star":
            return "cos#theta*"
        elif "pt" in param_name:
            return "p_{T}(" + param_name.replace(" pt", "") + ")"
        elif "DeltaEta" in param_name or "delta eta" in param_name:
            return "#Delta#eta(" + param_name.replace("DeltaEta", "").replace("delta eta","") + ")"
        elif "DeltaPhi" in param_name:
            return "#Delta#phi(" + param_name.replace("DeltaPhi", "") + ")"
        return param_name

    def format_Y_axis(X_param, nb_bins, max_hist):
        per_bins = round(max_hist / nb_bins, 2)
        return f"Normalized Entry / {per_bins} (GeV)" if "pt" in X_param or "mass" in X_param else f"Normalized Entry / {per_bins}"

    first_root_file_path = next(iter(root_files_info.values()))
    first_root_file = ROOT.TFile(first_root_file_path, "READ")
    tree = first_root_file.Get(tree_name)
    if not tree:
        print(f"Error: No TTree named '{tree_name}' found in file {first_root_file_path}.")
        return
    keys = [branch.GetName() for branch in tree.GetListOfBranches() if branch.GetName() != "Event Number"]
    first_root_file.Close()

    os.makedirs(output_plot, exist_ok=True)
    Fraction_txt=""

    for parameter_to_plot in keys:
        if parameter_to_plot not in variables_plot:
            continue
        process, decay= Process_name.split('_')[0:2]
        canvas = setup_canvas(parameter_to_plot)
        param_name = format_param_name(parameter_to_plot)
        Process_name_title = f"{process} {decay}"
        ecm = "13.6 TeV"
        sample_legend = f"{process} {decay} {ecm} (Merged)"
        legend_title = format_title_name(sample_legend)
        hs = ROOT.THStack("hs", f"{param_name} for {Process_name_title}")
        legend = setup_legend(len(root_files_info))
        legend_bis = ROOT.TLegend(0.20, 0.85, 0.55, 0.90)
        legend_bis.SetHeader(legend_title)

        min_hist, max_hist = get_histogram_limits(root_files_info, parameter_to_plot, tree_name, Xsec_model, Process_name, norm_to_xsec)
        #print(f"min_hist: {min_hist}, max_hist: {max_hist}")
        i = 0
        max_bin_content = float('-inf')
        min_bin_content = float('inf')
        for legend_name, file_path in root_files_info.items():
            process, decay, op, order_EFT = legend_name.split('_')[0:4]
            Process_op_order = f"{process}_{decay}_{op}_{order_EFT}"
            model_element = legend_name.split('_')[4:]
            model_name = "_".join(model_element)
            root_file = ROOT.TFile(file_path, "READ")
            tree = root_file.Get(tree_name)
            if not tree:
                root_file.Close()
                continue
            min_hist_, max_hist_ = adjust_histogram_range(min_hist, max_hist)
            if parameter_to_plot == "merged_ll_mass":
                min_hist_, max_hist_ = 88, 94
            scale_factor = Xsec_model[legend_name] if norm_to_xsec else 1.0
            weight_branch_name = "EventWeight"

            if any(keyword in model_name for keyword in ['rwg', 'reweight', 'wg']):
                scale_factor = Xsec_model[Process_op_order] if norm_to_xsec else 1.0
                op_order_strg = op + "_" + order_EFT
                op_order_strg_rwg = op_order_strg.lower().replace("vs", "_")
                weight_branch_name = f"EventWeight_{op_order_strg_rwg}"
                print(f"op_order_strg_rwg: {op_order_strg_rwg}")
                print(f"weight_cut keys: {weight_cut.keys()}")
                weight_cut_ = weight_cut[op_order_strg_rwg] if op_order_strg_rwg in weight_cut else 100
            else:
                weight_cut_ = 100

            print(f"Model name: {model_name}, weight_branch_name: {weight_branch_name}, weight_cut: {weight_cut_}")
            histogram = get_histogram(tree, model_name,parameter_to_plot, desired_num_bins, weight_branch_name, min_hist_, max_hist_, scale_factor, weight_cut_,Fraction_txt)
            print("fraction",Fraction_txt)
            histogram.SetDirectory(0)
            if "hel_aware" in model_name:
                histogram.SetLineColor(ROOT.kBlack)
            else:
                histogram.SetLineColor(categories_bis[op][0])
            histogram.SetLineWidth(3)
            if any(keyword in model_name for keyword in ['rwg', 'reweight', 'wg']):
                histogram.SetLineStyle(1)
            else:
                histogram.SetLineStyle(2)
            max_bin_content = max(max_bin_content, histogram.GetMaximum())
            min_bin_content = min(min_bin_content, histogram.GetMinimum())
            hs.Add(histogram)

            print('Model name: ', model_name)
            formatted_model_name = format_model_name(model_name)

            legend.AddEntry(histogram, f"{formatted_model_name} {op} {order_EFT}", "l")
            root_file.Close()
            i += 1

        if hs.GetNhists() > 0:
            hs.Draw("nostack HIST E")
            hs.GetXaxis().SetTitle(format_param_name(parameter_to_plot))
            hs.SetTitle(f"{param_name} for {Process_name_title}")
            print(f"min_bin_content: {min_bin_content}, max_bin_content: {max_bin_content}")
            hs.SetMaximum(max_bin_content * 1.5)
            if max_bin_content ==0:
                print(f"max_bin_content is 0")
                hs.SetMaximum(abs(min_bin_content)*0.5)
            hs.GetXaxis().SetNdivisions(505)
            hs.GetYaxis().SetTitle("")
            if not opts.linear:
                canvas.SetLogy()
                ROOT.gPad.SetLogy()
            legend_bis.Draw()
            legend.Draw()
            ROOT.ATLASLabel(0.20, 0.9, "Internal")
            canvas.SetCanvasSize(1000, 1000)
            canvas.Update()
            save_canvas(canvas, output_plot + f"/{op}_{order_EFT}/", parameter_to_plot)
            open(f"{output_plot}/Fraction_lost.txt", "a").write(Fraction_txt)
            outdir_for_plot = output_plot + f"/S_T_op/"    
            outdir_for_plot = output_plot + f"/{op}_{order_EFT}/"    
            os.makedirs(outdir_for_plot, exist_ok=True)
            save_canvas(canvas, outdir_for_plot , parameter_to_plot)

def construct_path(base_path, folder_name, process, decay, op, order_EFT, name_model_):
    base_dir = f"{base_path}/{process}_{decay}/{folder_name}"

    if op == "SM":
        return f"{base_dir}/{op}_SM/ntuple_rivet.root"

    if any(keyword in name_model_ for keyword in ['rwg', 'reweight', 'wg']):
        return f"{base_dir}/ntuple_rivet.root"

    return f"{base_dir}/{op}_{order_EFT}/ntuple_rivet.root"

def find_files_and_cross_sections(processes, decays, all_op_plot, order, aQGC_models_name, base_dir, Special_name):
    Root_paths = {}
    X_section = {}
    base_path = base_dir + f"{Special_name}/"

    for process, decay, op in [(p, d, o) for p in processes for d in decays for o in all_op_plot]:
        if "vs" in op:
            order_EFT = "CROSS"
        else:
            order_EFT = order[0]

        xsec_noModel = get_cross_section(op, order_EFT, process, decay, "aqgc_new")
        if xsec_noModel is not None:
            key_noModel = f"{process}_{decay}_{op}_{order_EFT}"
            X_section[key_noModel] = xsec_noModel

        for index, name_model_ in enumerate(aQGC_models_name):
            base_path = base_path
            Folder_name = f"{name_model_}/"
            path = construct_path(base_path, Folder_name, process, decay, op, order_EFT, name_model_)
            print(f"Model: {name_model_}")
            print(f"Looking for file at path: {path}")
            matches = glob.glob(path)
            print(f"Found {len(matches)} matches for operator {op} and process {process}")
            print(f"Matches: {matches}")
            if not matches:
                print(f"No match found for operator {op} and process {process}")
                continue
            xsec = get_cross_section(op, order_EFT, process, decay, name_model_)
            if xsec is not None:
                key = f"{process}_{decay}_{op}_{order_EFT}_{name_model_}"
                Root_paths[key] = matches[0]
                X_section[key] = xsec

    return Root_paths, X_section

Root_paths, X_section = find_files_and_cross_sections(processes, decays, all_ops_both, order, aQGC_models_name, base_dir, Special_name)

for process in processes:
    for decay in decays:
        all_op_plot = all_ops_SM if opts.Full_op else all_ops_both
        Root_paths_model = {k: v for k, v in Root_paths.items() if f"{process}_{decay}" in k}
        Process_name = f"{process}_{decay}"

        X_section_model = {k: v for k, v in X_section.items() if f"{process}_{decay}" in k}
        for op in all_op_plot:
            


            num_bins = int(opts.bins)
            order_eft_=order[0]
            if "vs" in op:
                order_eft_ = "CROSS" 
            Process_name = f"{process}_{decay}_{op}_{order_eft_}"
            Root_paths_model = {k: v for k, v in Root_paths.items() if f"{process}_{decay}_{op}" in k}

            outPlot = (f"{base_dir_plot}/{process}_{decay}/")
            os.makedirs(outPlot, exist_ok=True)
            plot_histograms2(num_bins, Root_paths_model, outPlot, Process_name, X_section, norm_xsec, tree_name="Merged",weight_cut=weight_cut)