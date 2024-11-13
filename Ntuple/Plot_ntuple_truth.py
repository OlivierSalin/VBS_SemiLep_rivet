import ROOT
import os
from itertools import product

# Open the ROOT file
file_path = "Truth_eta/Truth_aQGC/WmZ_llqq_ntuple_truth.root"

def plot_histograms1d(root_files_path, output_plot,tree_name,desired_num_bins,legend_name="WmZ_llqq",eta_cut="",linear=False):
    root_file = ROOT.TFile(root_files_path, "READ")
    tree = root_file.Get(tree_name)  # Get the TTree named "Merged"
    if not tree:
        print(f"Error: No TTree named 'Truth' found in file {root_files_path}.")
        return
    keys = [branch.GetName() for branch in tree.GetListOfBranches() if branch.GetName() != "Event Number"]  # Get the list of branches
    print(keys)
    if not os.path.exists(output_plot):
        os.makedirs(output_plot)

    for parameter_to_plot in keys:    
        if tree.GetBranch(parameter_to_plot) is None:
            print(f"Error: No branch named '{parameter_to_plot}' found in the tree.")
            continue

        canvas = ROOT.TCanvas(parameter_to_plot, parameter_to_plot, 1000, 600)
        ROOT.gStyle.SetOptStat(0)
        histogram = ROOT.TH1F(parameter_to_plot, parameter_to_plot, desired_num_bins, tree.GetMinimum(parameter_to_plot), tree.GetMaximum(parameter_to_plot))

        tree.Draw(parameter_to_plot + ">>" + parameter_to_plot)
        
        # Normalize the histogram
        integral = histogram.Integral(0,-1)
        if integral > 0:  # Avoid division by zero
            histogram.Scale(1.0 / integral)
            
        histogram.SetLineWidth(3)    
        histogram.Draw("HIST")
        canvas.Update()

        # Adjust legend position
        legend = ROOT.TLegend(0.85, 0.85, 1.05, 0.9)  # Move legend to the right of the plot        
       # Adjust these coordinates as needed
        legend.AddEntry(histogram, legend_name, "l")
        legend.Draw()

        if not linear:
            ROOT.gPad.SetLogy()

        canvas.SaveAs(os.path.join(output_plot, eta_cut + "_"+ parameter_to_plot + ".png"))

    root_file.Close()
    
def plot_histograms2d(root_files_path, output_plot,tree_name, desired_num_bins, pairs, legend_name="WmZ_llqq",eta_cut="", linear=False):
    root_file = ROOT.TFile(root_files_path, "READ")
    tree = root_file.Get(tree_name)  
    if not tree:
        print(f"Error: No TTree named 'Truth' found in file {root_files_path}.")
        return

    if not os.path.exists(output_plot):
        os.makedirs(output_plot)

    for pair in pairs:
        x_parameter, y_parameter = pair

        if tree.GetBranch(x_parameter) is None or tree.GetBranch(y_parameter) is None:
            print(f"Error: No branch named '{x_parameter}' or '{y_parameter}' found in the tree.")
            continue

        canvas = ROOT.TCanvas(f"{x_parameter}_{y_parameter}", f"{legend_name}: {x_parameter} vs {y_parameter}", 1000, 600)
        histogram = ROOT.TH2F(f"{x_parameter}_{y_parameter}", f"{legend_name}: {x_parameter} vs {y_parameter}", 
                               desired_num_bins, tree.GetMinimum(x_parameter), tree.GetMaximum(x_parameter),
                               desired_num_bins, tree.GetMinimum(y_parameter), tree.GetMaximum(y_parameter))

        histogram.SetXTitle(x_parameter)
        histogram.SetYTitle(y_parameter)
        tree.Draw(f"{y_parameter}:{x_parameter}>>{x_parameter}_{y_parameter}", "", "COLZ")

        # Normalize the histogram
        integral = histogram.Integral()
        if integral > 0:  # Avoid division by zero
            histogram.Scale(1.0 / integral)
            
        canvas.Update()  # Force canvas update to ensure the stats box is created

        # Attempt to move the statistics box



        if not linear:
            ROOT.gPad.SetLogz()

        canvas.SaveAs(os.path.join(output_plot, f"{eta_cut}_{x_parameter}_{y_parameter}.png"))

    root_file.Close()


Keys=[
    "VBS_Quarks_truth_Delta_Eta", "VBS_Quark1_truth_pT", "VBS_Quark1_truth_eta", "VBS_Quark1_truth_mass",
    "VBS_Quark2_truth_pT", "VBS_Quark2_truth_eta", "VBS_Quark2_truth_mass", "VBS_Quarks_truth_pT",
    "VBS_Quarks_truth_eta", "VBS_Quarks_truth_mass", "Tagging_Jet1_pT", "Tagging_Jet1_eta", "Tagging_Jet1_mass",
    "Tagging_Jet2_pT", "Tagging_Jet2_eta", "Tagging_Jet2_mass", "Tagging_Jets_mass", "Tagging_Jets_delta_Eta",
    "Delta_Eta_Tagging_jets_VBS_q_truth", "Delta_Mass_Tagging_jets_VBS_q_truth", "Delta_pT_Tagging_jets_VBS_q_truth",
    "Delta_eta_Tagging_Jet1_VBS_q", "Delta_eta_Tagging_Jet2_VBS_q", "tag_jets_matched", "tag_jet1_matched",
    "tag_jet2_matched", "Pass_VBS_jet", "W_Quark1_pT", "W_pT", "W_Quark1_eta", "W_Quark1_mass", "W_Quark2_pT",
    "W_Quark2_eta", "W_Quark2_mass", "W_Quarks_delta_Eta", "Delta_Eta_TaggingQuark1_Wquarks", "Delta_Eta_TaggingQuark2_Wquarks"
]

Keys_2D=[
    "VBS_Quarks_truth_Delta_Eta", "VBS_Quark1_truth_pT", "VBS_Quarks_truth_pT",
    "VBS_Quarks_truth_eta", "VBS_Quarks_truth_mass", "Tagging_Jet1_pT", "Tagging_Jet1_eta", "Tagging_Jet1_mass",
    "Tagging_Jet2_pT", "Tagging_Jet2_eta", "Tagging_Jet2_mass", "Tagging_Jets_mass", "Tagging_Jets_delta_Eta",
    "Delta_Eta_Tagging_jets_VBS_q_truth", "Delta_Mass_Tagging_jets_VBS_q_truth", "Delta_pT_Tagging_jets_VBS_q_truth",
    "Delta_eta_Tagging_Jet1_VBS_q", "Delta_eta_Tagging_Jet2_VBS_q", "tag_jets_matched", "tag_jet1_matched",
    "tag_jet2_matched",  "W_pT", "W_Quarks_delta_Eta", "Delta_Eta_TaggingQuark1_Wquarks", "Delta_Eta_TaggingQuark2_Wquarks"
]


detail_path=f"/Tables/VBS_jet_algo_plot/Highest_eta/2Lepton/"
f"/exp/atlas/salin/ATLAS/VBS_mc/eft_files/{detail_path}/Combined/"

Eta_cut=["CutEta_No","CutEta_25","CutEta_35","CutEta_40"]
option_SM_aQGC=["SM","aQGC"]
for eta_cut in Eta_cut:
    for op in option_SM_aQGC:
        Cut_folder= "Truth_bis/" + eta_cut +"/Truth_" + op +"/"
        file_path = Cut_folder+ "/WmZ_llqq_ntuple_truth.root"
        output_plot_1D="Truth_bis/Plot/Plot_1D/"+Cut_folder

        desired_num_bins=100
        process_op="WmZ_llqq" + "_"+ op 
        plot_histograms1d(file_path, output_plot_1D,desired_num_bins,legend_name=process_op,eta_cut=eta_cut,linear=False)

        output_plot_2D="Truth_bis/Plot/Plot_2D/"+Cut_folder
        desired_num_bins_2D=20



        parameters = ['Delta_Eta_Tagging_jets_VBS_q_truth', 'Delta_Mass_Tagging_jets_VBS_q_truth', "tag_jets_matched","Tagging_Jets_mass", "Tagging_Jets_delta_Eta"]
        targets = ['VBS_Quarks_truth_mass', 'VBS_Quarks_truth_Delta_Eta', 'Tagging_Jets_mass', 'W_pT']

        pairs = [list(pair) for pair in product(parameters, targets) if pair[0] != pair[1]]
        
        plot_histograms2d(file_path, output_plot_2D, desired_num_bins_2D, pairs, legend_name=process_op,eta_cut=eta_cut, linear=False)
