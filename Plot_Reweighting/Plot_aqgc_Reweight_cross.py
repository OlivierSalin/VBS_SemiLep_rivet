import glob
import os
from array import array
import ROOT
import utils_func as uf
import shutil
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math
import json
import re
import mplhep as hep


from format_utils import (
    format_model_name,
    format_title_name,
    format_param_name,
    format_Y_axis,
    get_cross_section,
)

# --- ROOT and plotting configuration ---
ROOT.gROOT.SetBatch(True)
plt.rcParams['text.usetex'] = True
ROOT.gROOT.LoadMacro("~/ATLAS/atlasrootstyle/AtlasStyle.C")
ROOT.gROOT.LoadMacro("~/ATLAS/atlasrootstyle/AtlasLabels.C")
ROOT.gROOT.LoadMacro("~/ATLAS/atlasrootstyle/AtlasUtils.C")
ROOT.SetAtlasStyle()
ROOT.gStyle.SetLegendBorderSize(0)
ROOT.gStyle.SetLegendFillColor(0)
ROOT.gStyle.SetLegendTextSize(0.03)

# --- Option parsing and global setup ---
def parse_options():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--DOCUT", default="YES")
    parser.add_option("--All_channel", default=False)
    parser.add_option("--linear", default=True)
    parser.add_option("--Name", default="test")
    parser.add_option("--type_MC", default="")
    parser.add_option("--EFT_order", default="QUAD")
    parser.add_option("--Channel", default="")
    parser.add_option("--nb_lep", default=2)
    parser.add_option("--bins", default=25)
    parser.add_option("--Xsec", default=False)
    opts, _ = parser.parse_args()
    return opts

opts = parse_options()
# --- Constants and configuration ---
processes = ["WpZ"]
decays = ["llqq"]
order = ["QUAD"]


weight_cut = {
    "fm0_quad": 0.00001, "fm0_fm1_cross": 0.0006, "fm1_quad": 0.0005, "fm2_quad": 0.000005,
    "fm8_quad": 0.00001, "fs0_quad": 2e-7, "fs1_quad": 5e-7
}

categories_bis = {
    "FM0": [ROOT.kMagenta-4, 1], "FM1": [ROOT.kPink+1, 9], "FM7": [ROOT.kPink+6, 7], "FM2": [ROOT.kRed+1, 1],
    "FM3": [ROOT.kOrange+10, 9], "FM4": [ROOT.kOrange-3, 7], "FM5": [ROOT.kYellow-9, 4], "FM8": [ROOT.kYellow-6, 4],
    "FM9": [ROOT.kYellow-3, 4], "FT0": [ROOT.kBlue-9, 1], "FT1": [ROOT.kViolet-9, 9], "FM0vsFM1": [ROOT.kAzure, 1],
    "FT2": [ROOT.kViolet-6, 7], "FT3": [ROOT.kViolet-3, 4], "FT4": [ROOT.kAzure-3, 1], "FT5": [ROOT.kAzure, 1],
    "FT6": [ROOT.kAzure+7, 9], "FT7": [ROOT.kCyan-3, 7], "FT8": [ROOT.kCyan-10, 1], "FT9": [ROOT.kTeal-6, 9],
    "FS1": [ROOT.kGreen+3, 1], "FS0": [ROOT.kSpring-9, 9], "FS2": [ROOT.kGreen-4, 5], "SM": [ROOT.kBlack, 1]
}

# --- Utility functions ---

def get_histogram_limits(root_files_info, parameter_to_plot, tree_name, Xsec_model, Process_name, norm_to_xsec):
    min_hist, max_hist = 0, float('-inf')
    
    for legend_name, file_path in root_files_info.items():
        process, decay, op, order_EFT = legend_name.split('_')[0:4]
        Process_op_order = f"{process}_{decay}_{op}_{order_EFT}"
        root_file = ROOT.TFile(file_path, "READ")
        tree = root_file.Get(tree_name)
        
        if not tree:
            root_file.Close()
            continue
        
        min_val,max_val  = tree.GetMinimum(parameter_to_plot), tree.GetMaximum(parameter_to_plot)
        min_hist, max_hist = min(min_hist, min_val), max(max_hist, max_val)

    return min_hist, max_hist


def get_histogram(tree, model_name, parameter, bins, weight_branch_name, min_val, max_val, scale_factor=1.0, weight_cut=None):
    histogram = ROOT.TH1D("histogram", "title", bins, min_val, max_val)
    draw_option = weight_branch_name if weight_branch_name in [branch.GetName() for branch in tree.GetListOfBranches()] else ""
    if "hel_aware" in model_name:
        tree.Draw(f"{parameter}>>histogram", f"(({weight_branch_name} < {weight_cut})&&({weight_branch_name}>0))*{weight_branch_name}", "norm")
    else:
        tree.Draw(f"{parameter}>>histogram", draw_option, "norm")
    histogram.Scale(scale_factor)
    if histogram.GetEntries() > 0:
        histogram.Sumw2()
    return histogram


def construct_path(base_path, folder_name, process, decay, op, op_rwg, order_EFT, name_model_):
    base_dir = f"{base_path}/{folder_name}/{process}_{decay}/"
    if op == "SM":
        return f"{base_dir}/{op}_SM/"
    if any(keyword in name_model_ for keyword in ["Reweight",'rwg', 'reweight', 'wg']):
        return f"{base_dir}/{op_rwg}_{order_EFT}/"
    return f"{base_dir}/{op}_{order_EFT}/"


def find_files_and_Xsections_sumW(processes, decays, all_op_plot, op_rwg, order, Run_model_name, base_dir, Special_name):
    Root_paths = {}
    X_section = {}
    SumWeights = {}
    base_path = base_dir + f"{Special_name}/"

    for process, decay, op in [(p, d, o) for p in processes for d in decays for o in all_op_plot]:
        order_EFT = "CROSS" if "vs" in op else order[0]
        for name_model_ in Run_model_name:
            Folder_name = f"{name_model_}/"
            path = construct_path(base_path, Folder_name, process, decay, op, op_rwg, order_EFT, name_model_)
            print(f"Model: {name_model_}")
            print(f"Looking for file at path: {path}")
            matches = glob.glob(path)
            print(f"Found {len(matches)} matches for operator {op} and process {process}")
            print(f"Matches: {matches}")
            if not matches:
                print(f"No match found for operator {op} and process {process}")
                continue
            good_path = matches[0]
            print(f"Good path: {good_path}")
            key = f"{process}_{decay}_{op}_{order_EFT}_{name_model_}"
            Root_paths[key] = good_path + '/ntuple_rivet.root'
            xsec_path = good_path + '/X_section_fb.txt'
            sumW_path = good_path + '/SumW.txt'
            if glob.glob(xsec_path):
                xsec = float(open(xsec_path, "r").readline())
            else:
                xsec = get_cross_section(op, order_EFT, process, decay, name_model_)
            if glob.glob(sumW_path):
                sumW = float(open(sumW_path, "r").readline())
            else:
                sumW = None
            print(f"Cross section: {xsec} fb")
            print(f"SumW: {sumW}")
            if xsec is not None:
                X_section[key] = xsec
                SumWeights[key] = sumW
    return Root_paths, X_section, SumWeights


# --- Plotting function ---
def plot_histograms2(desired_num_bins,variables_plot, root_files_info, output_plot, Process_name, Xsec_model,
                     norm_to_xsec=False, tree_name="Merged", weight_cut=None,perOp=True, op_rwg=None):
    def setup_canvas(parameter):
        return ROOT.TCanvas(f"canvas_{parameter}", "Stacked Histograms", 3000, 3000)

    def setup_legend(num_entry):
        length = 0.01 * num_entry
        legend = ROOT.TLegend(0.50, 0.78 - length, 0.90, 0.85)
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

    first_root_file_path = next(iter(root_files_info.values()))
    first_root_file = ROOT.TFile(first_root_file_path, "READ")
    tree = first_root_file.Get(tree_name)
    if not tree:
        print(f"Error: No TTree named '{tree_name}' found in file {first_root_file_path}.")
        return
    keys = [branch.GetName() for branch in tree.GetListOfBranches() if branch.GetName() != "Event Number"]
    first_root_file.Close()

    os.makedirs(output_plot, exist_ok=True)

    for parameter_to_plot in keys:
        if parameter_to_plot not in variables_plot:
            continue
        process, decay = Process_name.split('_')[0:2]
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
        i = 0
        max_bin_content, min_bin_content = float('-inf'), float('inf')
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
            
            if parameter_to_plot == "merged_ll_mass" or parameter_to_plot == "merged_Vlep_mass":
                min_hist_, max_hist_ = 80, 100
                
            scale_factor = Xsec_model[legend_name] if norm_to_xsec else 1.0
            weight_branch_name = "EventWeight"
            if any(keyword in model_name for keyword in ["Reweight",'rwg', 'reweight', 'wg']):
                scale_factor = Xsec_model[Process_op_order] if norm_to_xsec else 1.0
                op_order_strg = op + "_" + order_EFT
                op_order_strg_rwg = op_order_strg.lower().replace("vs", "_")
                weight_branch_name = f"EventWeight_{op_order_strg_rwg}"
                weight_cut_ = weight_cut[op_order_strg_rwg] if weight_cut and op_order_strg_rwg in weight_cut else 100
            else:
                weight_cut_ = 100
            histogram = get_histogram(tree, model_name, parameter_to_plot, desired_num_bins, weight_branch_name, min_hist_, max_hist_, scale_factor, weight_cut_)
            histogram.SetDirectory(0)
            if "hel_aware" in model_name:
                histogram.SetLineColor(ROOT.kBlack)
            else:
                histogram.SetLineColor(categories_bis[op][0])
            histogram.SetLineWidth(3)
            histogram.SetLineStyle(1 if any(keyword in model_name for keyword in ["Reweight",'rwg', 'reweight', 'wg']) else 2)
            max_bin_content = max(max_bin_content, histogram.GetMaximum())
            min_bin_content = min(min_bin_content, histogram.GetMinimum())
            hs.Add(histogram)
            print('Model name: ', model_name)
            formatted_model_name = format_model_name(model_name, op_rwg)
            legend.AddEntry(histogram, f"{formatted_model_name} {op} {order_EFT}", "l")
            root_file.Close()
            i += 1
        if hs.GetNhists() > 0:
            hs.Draw("nostack HIST E")
            hs.GetXaxis().SetTitle(format_param_name(parameter_to_plot))
            hs.SetTitle(f"{param_name} for {Process_name_title}")
            print(f"min_bin_content: {min_bin_content}, max_bin_content: {max_bin_content}")
            hs.SetMaximum(max_bin_content * 1.7)
            if max_bin_content == 0:
                print(f"max_bin_content is 0")
                hs.SetMaximum(abs(min_bin_content) * 0.5)
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
            if perOp:
                outdir_for_plot = output_plot + f"/{op}_{order_EFT}/"
            else:
                outdir_for_plot = output_plot + f"/All_op/"
            os.makedirs(outdir_for_plot, exist_ok=True)
            save_canvas(canvas, outdir_for_plot, parameter_to_plot)



# --- Main execution ---
def main():
    opts = parse_options()
    order_ = opts.EFT_order
    nb_lepton = int(opts.nb_lep)
    num_bins = int(opts.bins)
    name_plt = opts.Name
    norm_xsec = False  # Always false as in original code


    types = ["Madspin", "NoSpin"]
    #types = ["Madspin"]

    all_operators=["FM0","FM1","FM2","FM3","FM4","FM5","FM7","FM8","FM9",
                "FS0","FS1","FS2",
                "FT0","FT1","FT2","FT3","FT4","FT5","FT6","FT7"]
    
    all_operators=["FM0","FM2","FM7","FM8",
                "FS0","FS1",
                "FT0","FT2","FT5"]
    operators_rwg= ["FM","FS","FT"]

    #all_ops_both = ["FS0", "FS1", "FS2"]

    variables_plot = [
        'EventWeight', 'merged_cos_theta_star', 'merged_Vhad_DeltaEta_Vlep', 'merged_Vhad_DeltaPhi_Vlep',
        'merged_Vhad_DeltaR_Vlep', 'merged_VlepVhad_eta', 'merged_VlepVhad_mass', 'merged_VlepVhad_pt',
        'merged_tagjets_delta_eta', 'merged_tagjets_mass', 'merged_tagjets_pt', 'merged_Full_pt',"merged_fjet_mass", 'merged_fjet_pt',
        "merged_Vlep_pt","merged_Vlep_mass", "merged_lepton*","merged_Centrality*", "merged_CS_V_cos_theta",
    ]
    #variables_plot = ["merged_VlepVhad_mass"]



    base_dir = "/exp/atlas/salin/ATLAS/VBS_mc/eft_files/Histograms/"
    Special_name = "/Reweighting/Rwg_test/"

    for type_ in types:

        Complement_path = f"Study_SpinCorrelation/Type_lessOP/{type_}/"
        base_dir_plot = f"/exp/atlas/salin/ATLAS/VBS_mc/Plots/Plot_Reweighting/{Complement_path}/bins_{opts.bins}/"

        reweighting_cat = f"Reweighting_{type_}"
        eftdec_cat = f"EFTDec_{type_}"
        aQGC_models = [reweighting_cat, eftdec_cat]
        for op_rwg in operators_rwg:
            if op_rwg in ["FM", "FS", "FT"]:
                eftdec_ops = [op for op in all_operators if op.startswith(op_rwg)]
            else:
                eftdec_ops = all_operators
                
            all_ops_both = eftdec_ops

            Root_paths, X_section, SumW = find_files_and_Xsections_sumW(
                processes, decays, all_ops_both, op_rwg, order, aQGC_models, base_dir, Special_name
            )

            for process in processes:
                for decay in decays:
                    all_op_plot = all_ops_both
                    Root_paths_model = {k: v for k, v in Root_paths.items() if f"{process}_{decay}" in k}
                    
                    Process_name = f"{process}_{decay}"
                    
                    outPlot = f"{base_dir_plot}/{process}_{decay}/{op_rwg}/"
                    os.makedirs(outPlot, exist_ok=True)
                    
                    plot_histograms2(num_bins,variables_plot , Root_paths_model, outPlot, Process_name, X_section,
                                     norm_xsec, tree_name="Merged", weight_cut=weight_cut, perOp=False, op_rwg=op_rwg)
                    for op in all_op_plot:
                        order_eft_ = "CROSS" if "vs" in op else order[0]
                        Process_name = f"{process}_{decay}_{op}_{order_eft_}"
                        Root_paths_model = {k: v for k, v in Root_paths.items() if f"{process}_{decay}_{op}" in k}
                        outPlot = f"{base_dir_plot}/{process}_{decay}/{op_rwg}/"
                        os.makedirs(outPlot, exist_ok=True)
                        plot_histograms2(num_bins,variables_plot, Root_paths_model, outPlot, Process_name, X_section
                                         , norm_xsec, tree_name="Merged", weight_cut=weight_cut, perOp=True,op_rwg=op_rwg)
            # Add the new comparison logic
  
if __name__ == "__main__":
    main() 