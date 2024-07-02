import glob
import os
from array import array
import ROOT
import shutil
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

import lib_utils as lu

from optparse import OptionParser
parser = OptionParser()
parser.add_option("--Ana", default = "WpZ_llqq")
parser.add_option("--DOCUT", default = "YES")
parser.add_option("--SM", default = "YES")
parser.add_option("--QUAD1", default = "FS1")
parser.add_option("--QUAD2", default = "FT1")
parser.add_option("--bins", default = 50)
opts, _ = parser.parse_args()


base_dir = "/exp/atlas/salin/ATLAS/VBS_mc/eft_files"

root_files_info = {
    f"{base_dir}/{opts.Ana}/user.osalin.MadGraph_{opts.Ana}_FM0_SM_EXT0/DOCUT_{opts.DOCUT}/hists.root": {'color': ROOT.kRed, 'label': f'{opts.Ana}_FM0_SM'},
    f"{base_dir}/{opts.Ana}/user.osalin.MadGraph_{opts.Ana}_{opts.QUAD1}_QUAD_EXT0/DOCUT_{opts.DOCUT}/hists.root": {'color': ROOT.kBlue, 'label': f'{opts.Ana}_{opts.QUAD1}_QUAD'},
    f"{base_dir}/{opts.Ana}/user.osalin.MadGraph_{opts.Ana}_{opts.QUAD2}_QUAD_EXT0/DOCUT_{opts.DOCUT}/hists.root": {'color': ROOT.kGreen, 'label': f'{opts.Ana}_{opts.QUAD2}_QUAD'},
}
# Call the function with the desired parameters
desired_num_bins = int(opts.bins)  # Replace with your desired number of bins

output_plot = f"{base_dir}/{opts.Ana}/Plots/{opts.QUAD1}_{opts.QUAD2}/DOCUT_{opts.DOCUT}/"

# Check if the directory exists, if not, create it
if not os.path.exists(output_plot):
    os.makedirs(output_plot)
    
output_plot_cutflow = f"{base_dir}/{opts.Ana}/Plots/{opts.QUAD1}_{opts.QUAD2}/DOCUT_{opts.DOCUT}/Cutflow_comp/"
# Check if the directory exists, if not, create it
if not os.path.exists(output_plot_cutflow):
    os.makedirs(output_plot_cutflow)
    
Cutflow_paths = {
    "SM": f"{base_dir}/{opts.Ana}/user.osalin.MadGraph_{opts.Ana}_FM0_SM_EXT0/DOCUT_{opts.DOCUT}/",
    opts.QUAD1: f"{base_dir}/{opts.Ana}/user.osalin.MadGraph_{opts.Ana}_{opts.QUAD1}_QUAD_EXT0/DOCUT_{opts.DOCUT}/",
    opts.QUAD2: f"{base_dir}/{opts.Ana}/user.osalin.MadGraph_{opts.Ana}_{opts.QUAD2}_QUAD_EXT0/DOCUT_{opts.DOCUT}/",
}
    
    
def plot_histograms(desired_num_bins, root_files_info,output_plot):
    # Open the first ROOT file to retrieve the list of parameters (branches)
    first_root_file_path = next(iter(root_files_info))
    first_root_file = ROOT.TFile(first_root_file_path, "READ")
    keys = [key.GetName() for key in first_root_file.GetListOfKeys()]
    print(keys)
    first_root_file.Close()

    for parameter_to_plot in keys:    
        # Create a new TCanvas
        canvas_name = "canvas_" + parameter_to_plot
        canvas = ROOT.TCanvas(canvas_name, "Stacked Histograms", 2000, 2000)
        canvas.Divide(1, 2)  # Divide canvas into 2 pads: one for the main plot and one for the ratio plot

        # Create a THStack
        hs = ROOT.THStack("hs", "Distribution of " + parameter_to_plot+ ' \nnb of bins:'+str(desired_num_bins))
        
        # Create a legend
        legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
        
        # List to store histograms for ratio plot
        histograms = []

        # Loop over the files to retrieve and stack the histograms
        for file_path, info in root_files_info.items():
            # Open the ROOT file
            root_file = ROOT.TFile(file_path, "READ")
        
            # Retrieve the histogram
            histogram = root_file.Get(parameter_to_plot)
            histogram.SetDirectory(0)  
            current_num_bins = histogram.GetNbinsX()
            
            rebin_factor = int(current_num_bins / desired_num_bins)
            if rebin_factor > 0:
                histogram.Rebin(rebin_factor)
                
            # Style the histogram
            histogram.SetFillColor(info['color'])
            histogram.SetMarkerStyle(21)
            histogram.SetMarkerColor(info['color'])
            histogram.SetMarkerSize(1.2)
            histogram.SetLineWidth(1)
            legend2 = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
            
            # Add the histogram to the stack
            hs.Add(histogram)
        
            # Add entry to the legend
            legend.AddEntry(histogram, info['label'], "f")
        
            # Close the ROOT file
            root_file.Close()

            # Add histogram to list for ratio plot
            histograms.append(histogram)

        # Draw the stacked histograms
        canvas.cd(1)  # Switch to first pad
        hs.Draw("nostack")
        hs.GetXaxis().SetTitle(parameter_to_plot )
        hs.GetYaxis().SetTitle("Events / bin")
        
        # Set log scale for y-axis
        canvas.SetLogy()
        ROOT.gPad.SetLogy()
        
        # Draw the legend
        legend.Draw()

        canvas.cd(2)  # Switch to second pad

        # Check if histograms list has at least 3 elements
        if len(histograms) < 3:
            print("Error: Less than 3 histograms")
        else:
            ratio_histogram = histograms[1].Clone()  # Clone first histogram to create ratio histogram

            # Check if histograms have the same number of bins
            if histograms[1].GetNbinsX() != histograms[2].GetNbinsX():
                print("Error: Histograms have different number of bins")
            else:
                # Check for empty bins
                for bin in range(1, histograms[1].GetNbinsX() + 1):
                    if histograms[1].GetBinContent(bin) == 0 or histograms[2].GetBinContent(bin) == 0:
                        print(f"Warning: Empty bin at index {bin}")
                    else:
                        # Divide bin by bin
                        ratio_histogram.SetBinContent(bin, histograms[1].GetBinContent(bin) / histograms[2].GetBinContent(bin))

            ratio_histogram.SetTitle(f"{opts.QUAD1}/{opts.QUAD2}" )  # Set title
            ratio_histogram.Draw("ep")  # Draw ratio histogram with

            # Create a legend
            legend2 = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
            
            legend2.SetBorderSize(0)  
            legend2.SetTextSize(0.03) 
            legend2.AddEntry(ratio_histogram, f"{opts.QUAD1}/{opts.QUAD2}", "l")
            legend2.Draw()

        canvas.SetCanvasSize(1000, 1000)
        ROOT.gPad.SetLogy()
        # Update the canvas
        
        canvas.Update()
        
        
        # Save the canvas to a file
        canvas.SaveAs(output_plot+parameter_to_plot + "_hist.png")
        canvas.SaveAs(output_plot+parameter_to_plot + "_hist.pdf")


plot_canvas = plot_histograms(desired_num_bins, root_files_info,output_plot)


def Cutflow_comp(cutflow_paths):
    # Initialize lists to store cutflow data for each case
    cut_names_list_merged = []
    cut_cumu_list_merged = []
    cut_incr_list_merged = []
    cut_names_list_resolved = []
    cut_cumu_list_resolved = []
    cut_incr_list_resolved = []

    for case in ["SM", opts.QUAD1, opts.QUAD2]:
        cutflow_merged_file = cutflow_paths[case] + "cutflow_merged.txt"
        if os.path.exists(cutflow_merged_file):
            shutil.copyfile(cutflow_merged_file, output_plot_cutflow + f"cutflow_{case}_merged.txt")
            cut_merged_names, cut_merged_cumu, cut_merged_incr = lu.get_cutflow_arrays(cutflow_merged_file)
            cut_names_list_merged.append(cut_merged_names)
            cut_cumu_list_merged.append(cut_merged_cumu)
            cut_incr_list_merged.append(cut_merged_incr)

        cutflow_resolved_file = cutflow_paths[case] + "cutflow_resolved.txt"
        if os.path.exists(cutflow_resolved_file):
            shutil.copyfile(cutflow_resolved_file, output_plot_cutflow + f"cutflow_{case}_resolved.txt")
            cut_resolved_names, cut_resolved_cumu, cut_resolved_incr = lu.get_cutflow_arrays(cutflow_resolved_file)
            cut_names_list_resolved.append(cut_resolved_names)
            cut_cumu_list_resolved.append(cut_resolved_cumu)
            cut_incr_list_resolved.append(cut_resolved_incr)

    # Draw cumulative cutflows in one file for merged
    lu.draw_cutflows(cut_names_list_merged[0], cut_cumu_list_merged, ["SM", {opts.QUAD1}, {opts.QUAD2}],
                     output_plot_cutflow + "/cutflow_cumu_merged.png", opts.Ana)

    # Draw incremental cutflows in another file for merged
    lu.draw_cutflows(cut_names_list_merged[0], cut_incr_list_merged, ["SM", {opts.QUAD1}, {opts.QUAD2}],
                     output_plot_cutflow + "/cutflow_incr_merged.png", opts.Ana)

    # Draw cumulative cutflows in one file for resolved
    lu.draw_cutflows(cut_names_list_resolved[0], cut_cumu_list_resolved, ["SM", {opts.QUAD1}, {opts.QUAD2}],
                     output_plot_cutflow + "/cutflow_cumu_resolved.png", opts.Ana)

    # Draw incremental cutflows in another file for resolved
    lu.draw_cutflows(cut_names_list_resolved[0], cut_incr_list_resolved, ["SM", {opts.QUAD1}, {opts.QUAD2}],
                     output_plot_cutflow + "/cutflow_incr_resolved.png", opts.Ana)
    
Cutflow_comp(Cutflow_paths)

    