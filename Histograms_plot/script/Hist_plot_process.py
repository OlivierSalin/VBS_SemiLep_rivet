import glob
import os
from array import array
import ROOT
import shutil
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True


from optparse import OptionParser
parser = OptionParser()
parser.add_option("--Ana", default = "WpZ_llqq")
parser.add_option("--DOCUT", default = "YES")
parser.add_option("--bins", default = 25)
opts, _ = parser.parse_args()

all_ops_SM = ["SM","FM0","FM1","FM2","FM3","FM4","FM5","FM7",
        "FS02","FS1",
        "FT0","FT1","FT2","FT5","FT6","FT7"]

base_dir_plot = "/exp/atlas/salin/ATLAS/VBS_mc/EFT_files_AMI/Plots_op2/"

categories = {
    "cat1": ["FM0", "FM1", "FM7", ROOT.kOrange+10, ROOT.kPink+10],
    "cat2": ["FM2", "FM3", "FM4", "FM5", ROOT.kRed+2, ROOT.kOrange-4],
    "cat3": ["FT0", "FT1", "FT2", ROOT.kBlue-9, ROOT.kViolet-3],
    "cat4": ["FT5", "FT6", "FT7", ROOT.kAzure, ROOT.kCyan],
    "cat5": ["FT8", "FT9", ROOT.kCyan-10, ROOT.kTeal-6],
    "cat6": ["FS1","FS02", ROOT.kGreen+3, ROOT.kSpring-9],
    "cat7": ["SM", ROOT.kBlack, ROOT.kBlack]
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

process_colors = {
    "WpZ_llqq": ROOT.kRed,
    "WmZ_llqq": ROOT.kBlue,
    "ZZ_llqq": ROOT.kGreen
}


def plot_histograms(desired_num_bins, root_files_info, output_plot):
    # Open the first ROOT file to retrieve the list of parameters (branches)
    first_root_file_path = next(iter(root_files_info.values()))
    first_root_file = ROOT.TFile(first_root_file_path, "READ")
    keys = [key.GetName() for key in first_root_file.GetListOfKeys()]
    print(keys)
    first_root_file.Close()


    for parameter_to_plot in keys:    
        # Create a new TCanvas
        canvas_name = "canvas_" + parameter_to_plot
        canvas = ROOT.TCanvas(canvas_name, "Stacked Histograms", 2000, 2000)
        ROOT.gPad.SetRightMargin(0.2)
        # Create a THStack
        hs = ROOT.THStack("hs", "Distribution of " + parameter_to_plot+ ' \n\nnb of bins:'+str(desired_num_bins))
        
        # Create a legend
        legend = ROOT.TLegend(0.82, 0.75, 1.0, 0.9)
        
        # Loop over the files to retrieve and stack the histograms
        for legend_name, file_path in root_files_info.items():
            # Open the ROOT file
            root_file = ROOT.TFile(file_path, "READ")
        
            # Retrieve the histogram
            histogram = root_file.Get(parameter_to_plot)
            histogram.SetDirectory(0)  
            current_num_bins = histogram.GetNbinsX()
            
            rebin_factor = int(current_num_bins / desired_num_bins)
            if rebin_factor > 0:
                histogram.Rebin(rebin_factor)
                
            proce = legend_name.split("_")[0]
            decay= legend_name.split("_")[1]
            process = proce+"_"+decay
            color = process_colors[process]  # Get the color for this process 
            histogram.SetFillColor(color)
            histogram.SetMarkerStyle(20)
            histogram.SetMarkerColor(color)
            histogram.SetMarkerSize(1.5)
            histogram.SetLineColorAlpha(color, 1.0)
            histogram.SetLineWidth(1)
            histogram.SetLineStyle(3)
            
            # Add the histogram to the stack
            hs.Add(histogram)
        
            # Add entry to the legend
            legend.AddEntry(histogram, legend_name, "lep")
        
            # Close the ROOT file
            root_file.Close()

        # Draw the stacked histograms
        hs.Draw("nostack")
        hs.GetXaxis().SetTitle(parameter_to_plot )
        hs.GetYaxis().SetTitle("Events / bin")
        
        # Set log scale for y-axis
        canvas.SetLogy()
        ROOT.gPad.SetLogy()
        
        # Draw the legend
        legend.Draw()
        canvas.SetCanvasSize(1000, 1000)
        ROOT.gPad.SetLogy()
        # Update the canvas
        canvas.Update()
        
        # Save the canvas to a file
        canvas.SaveAs(output_plot+parameter_to_plot + "_hist.png")
        #canvas.SaveAs(output_plot+parameter_to_plot + "_hist.pdf")

# Call the function with Root_paths as argument

base_dir = "/exp/atlas/salin/ATLAS/VBS_mc/EFT_files_AMI/"

all_ops_SM = ["SM","FM0","FM1","FM2","FM3","FM4","FM5","FM7",
        "FS02","FS1",
        "FT0","FT1","FT2","FT5","FT6","FT7"]




for op in all_ops_SM:
    Root_paths = {}
    for decay in ["llqq"]:
        for process in ["WmZ", "WpZ","ZZ"]:
            if op=="SM":
                path = os.path.join(base_dir, f"{process}_{decay}/mc16_13TeV.*.MGPy8EG_aQGCFM0_{op}_1_{process}_{decay}.merge.EVNT.*/DOCUT_YES/Tables/Tables_file_second/hists.root")
            else:
                path = os.path.join(base_dir, f"{process}_{decay}/mc16_13TeV.*.MGPy8EG_aQGC{op}_QUAD_1_{process}_{decay}.merge.EVNT.*/DOCUT_YES/Tables/Tables_file_second/hists.root")
            matches = glob.glob(path)
            #print(matches)
            if not matches:
                print(f"No match found for operator {op} and process {process}")
                continue
            Root_paths[f"{process}_{decay}_{op}"] = matches[0]
    #print(Cutflow_paths)



    # Call the function with the desired parameters
    desired_num_bins = int(opts.bins)  # Replace with your desired number of bins

    output_plot = f"{base_dir_plot}/2Lepton/{op}/"

    # Check if the directory exists, if not, create it
    if not os.path.exists(output_plot):
        os.makedirs(output_plot)

    plot_canvas = plot_histograms(desired_num_bins, Root_paths,output_plot)



    