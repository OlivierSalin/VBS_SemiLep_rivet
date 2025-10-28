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
import itertools


from format_utils import (
    format_model_name,
    format_title_name,
    format_param_name,
    format_Y_axis,
    get_cross_section,
)

import array
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
    parser.add_option("--CROSS", default="True")
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
    if "hel_awaredihfihfdi" in model_name:
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
        order_EFT_=order_EFT
        if order_EFT == "CROSS":
            order_EFT_ = "QUAD"
        return f"{base_dir}/{op_rwg}_{order_EFT_}/"
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
            if any(keyword in name_model_ for keyword in ["Reweight", "reweight", "rwg","Rwg"]):
                xsec_path= None
                sumW = None
            else:
                xsec_path = good_path + '/cross_section_fb.txt'
                sumW_path = good_path + '/SumW.txt'
            if xsec_path and glob.glob(xsec_path):
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




def calculate_xsec_unc_per_typeMC(Root_paths, X_section, op_rwg, all_ops_both, type_MC, output_dir, eff_sr_dict=None):

    xsec_lines = []
    unc_lines = []
    if eff_sr_dict is None:
        eff_sr_dict = {}
    #print(f"Root paths: {Root_paths}")
    print(f"\n--- Calculating cross-section and uncertainty for {type_MC} (op_rwg={op_rwg}) ---")
    for op in all_ops_both:
        order_EFT = "CROSS" if "vs" in op else "QUAD"
        is_reweight = any(keyword in type_MC for keyword in ["Reweight", "reweight", "rwg","Rwg"])
        if is_reweight:
            key = f"WpZ_llqq_{op}_{order_EFT}_{type_MC}"
            op_order_strg = op + "_" + order_EFT
            op_order_strg_rwg = op_order_strg.lower().replace("vs", "_")
            weight_branch_name = f"EventWeight_{op_order_strg_rwg}"
            #print(f"Using weight branch: {weight_branch_name} for reweighting model {type_MC}")
        else:
            key = f"WpZ_llqq_{op}_{order_EFT}_{type_MC}"
            weight_branch_name = "EventWeight"

        file_path = Root_paths.get(key)
        if not file_path or not os.path.exists(file_path):
            print(f"File not found for key: {key}")
            continue

        root_file = ROOT.TFile(file_path, "READ")
        tree = root_file.Get("Merged")
        tree_before_cut = root_file.Get("Bef_cut")
        if not tree or not tree_before_cut:
            print(f"Tree missing in {file_path} for {key}")
            root_file.Close()
            continue

        # For EFTDec: calculate and store efficiency only if not already stored
        if not is_reweight:
            if op not in eff_sr_dict:
                sum_weights_before = 0.0
                for event in tree_before_cut:
                    w = getattr(event, weight_branch_name, getattr(event, "EventWeight", 0.0))
                    sum_weights_before += w
                nbins = 50
                min_mass = 0
                max_mass = 5000
                hist_sr = ROOT.TH1D("hist_sr", "hist_sr", nbins, min_mass, max_mass)
                for event in tree:
                    w = getattr(event, weight_branch_name, getattr(event, f"{weight_branch_name}", 0.0))
                    mass = getattr(event, "merged_VlepVhad_mass", 0.0)
                    hist_sr.Fill(mass, w)
                sum_weights_sr = hist_sr.Integral(0, nbins+1)
                eff_sr = sum_weights_sr / sum_weights_before if sum_weights_before != 0 else 0.0
                eff_sr_dict[op] = eff_sr
            else:
                eff_sr = eff_sr_dict[op]
        else:
            # For reweighting, reuse efficiency from EFTDec
            eff_sr = eff_sr_dict.get(op, 1.0)
            nbins = 50
            min_mass = 0
            max_mass = 5000
            hist_sr = ROOT.TH1D("hist_sr", "hist_sr", nbins, min_mass, max_mass)
            for event in tree:
                w = getattr(event, weight_branch_name, getattr(event, f"{weight_branch_name}", 0.0))
                mass = getattr(event, "merged_VlepVhad_mass", 0.0)
                hist_sr.Fill(mass, w)
            sum_weights_sr = hist_sr.Integral(0, nbins+1)
        # Error on sum of weights in SR using ROOT histogram IntegralAndError
        err_sr = array.array('d', [0.])
        integral = hist_sr.IntegralAndError(0, nbins+1, err_sr)
        error = err_sr[0]
        # Scale by cross-section and efficiency
        xsec = X_section.get(key, 1.0)
        xsec_final = xsec * eff_sr
        unc_final = xsec_final * (error / sum_weights_sr) if sum_weights_sr != 0 else 0.0
        print(f"{key}: xsec={xsec_final:.6f} fb, unc={unc_final:.6f} fb, eff={eff_sr:.6f}, xsec_raw={xsec:.6f}")
        xsec_lines.append(f"WpZ_llqq_{op}_{order_EFT}: {xsec_final:.6f}\n")
        unc_lines.append(f"WpZ_llqq_{op}_{order_EFT}: {unc_final:.6f}\n")
        root_file.Close()
    # Write to txt files (append mode)
    os.makedirs(output_dir, exist_ok=True)
    with open(os.path.join(output_dir, f"VBS_cross_section_run2_{type_MC}.txt"), "a") as fx:
        fx.writelines(xsec_lines)
    with open(os.path.join(output_dir, f"VBS_Uncertainty_cross_section_run2_{type_MC}.txt"), "a") as fu:
        fu.writelines(unc_lines)
    #print(f"Written: {output_dir}/VBS_fromPlot_xsec_run2_{type_MC}.txt and ...unc...")
    

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
    
    all_operators=["FM0","FM2","FM3","FM7","FM9",
                "FS0","FS1",
                "FT0","FT2","FT5"]
    all_operators=["FM0","FM1","FM2","FM3","FM4","FM5","FM7","FM8","FM9",
                "FS0","FS1","FS2",
                "FT0","FT1","FT2","FT3","FT4","FT5","FT6"]
    
    
    all_operators=["FM0"]
    cross_terms_ = [f"{op1}vs{op2}" for op1, op2 in itertools.combinations(all_operators, 2) if op1[:2] == op2[:2]]

    if "True" in opts.CROSS:
        all_operators = cross_terms_
    if "both" in opts.CROSS:
        all_operators = all_operators + cross_terms_
    
    #all_operators = cross_terms_
    print(f"All operators: {all_operators}")
    operators_rwg= ["FM","FS","FT"]

    #all_ops_both = ["FS0", "FS1", "FS2"]

    variables_plot = ["merged_VlepVhad_mass",'merged_cos_theta_star','merged_tagjets_delta_eta', "merged_CS_V_cos_theta","merged_fjet_mass","merged_Vlep"]
    variables_plot = ["merged_VlepVhad_mass",'merged_cos_theta_star','merged_tagjets_delta_eta']



    base_dir = "/exp/atlas/salin/ATLAS/VBS_mc/eft_files/Histograms/"
    Special_name = "/Reweighting/Rwg_Unc/"


    Complement_path = f"Polarisation/TEST/Cross_section/Uncertainty/test03/"
    base_dir_plot = f"/exp/atlas/salin/ATLAS/VBS_mc/Plots/Plot_Reweighting/{Complement_path}/"


    aQGC_models = ["EFTDec_Madspin","Rwg_pol_50k","Rwg_pol_100k", "Rwg_InvSqrtXsec_50k", "Rwg_InvXsec_50k"]
    #aQGC_models = ["EFTDec_Madspin", "Rwg_InvSqrtXsec_50k", "Rwg_InvXsec_50k"]
    #aQGC_models = ["EFTDec_Madspin","EFTDec_Polarisation"]
 
    all_operators=["FM0","FM1","FM2","FM3","FM4","FM5","FM7","FM8","FM9",
                "FS0","FS1","FS2",
                "FT0","FT1","FT2","FT3","FT4","FT5","FT6"]
    
    operators_rwg= ["FM","FS","FT"]
    all_operators=["FM0","FM1","FM2","FM3","FM4","FM5","FM7","FM8","FM9",
                "FS0","FS1","FS2",
                "FT0","FT1","FT2","FT3","FT4","FT5","FT6"]

    #all_operators=["FS0","FS1",]
    cross_terms_ = [f"{op1}vs{op2}" for op1, op2 in itertools.combinations(all_operators, 2) if op1[:2] == op2[:2]]

    if "True" in opts.CROSS:
        all_operators = cross_terms_
    if "both" in opts.CROSS:
        all_operators = all_operators + cross_terms_
            
    eff_sr_dict = {}
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
                #plot_histograms2(num_bins,variables_plot , Root_paths_model, outPlot, Process_name, Xsec_model,
                #                    norm_xsec, tree_name="Merged", weight_cut=weight_cut, perOp=False, op_rwg=op_rwg)
                for op in all_op_plot:
                    order_eft_ = "CROSS" if "vs" in op else order[0]
                    Process_name = f"{process}_{decay}_{op}_{order_eft_}"
                    Root_paths_model = {k: v for k, v in Root_paths.items() if f"{process}_{decay}_{op}" in k}
                    outPlot = f"{base_dir_plot}/{process}_{decay}/{op_rwg}/"
                    os.makedirs(outPlot, exist_ok=True)
                    #plot_histograms2(num_bins,variables_plot, Root_paths_model, outPlot, Process_name, Xsec_model
                    #                    , norm_xsec, tree_name="Merged", weight_cut=weight_cut, perOp=True,op_rwg=op_rwg)
            # Call the cross-section/uncertainty calculation for each type_MC
            outDir_Unc = "/exp/atlas/salin/ATLAS/VBS_mc/plotting/Plot_Reweighting/Tables/Tables/Uncertainty/Xsec/FromPlot/test04/"
            for type_MC in aQGC_models:
                calculate_xsec_unc_per_typeMC(Root_paths, X_section, op_rwg, all_ops_both, type_MC, outDir_Unc, eff_sr_dict)
                
if __name__ == "__main__":
    main()