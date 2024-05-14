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
parser.add_option("--SM", default = "YES")
parser.add_option("--QUAD1", default = "FS1")
parser.add_option("--QUAD2", default = "FT1")
parser.add_option("--bins", default = 25)
opts, _ = parser.parse_args()

all_ops_SM = ["SM","FM0","FM1","FM2","FM3","FM4","FM5","FM7",
        "FS02","FS1",
        "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]

base_dir_plot = "/exp/atlas/salin/ATLAS/VBS_mc/EFT_files_AMI/Plots_op/"

categories = {
    "SM": [ROOT.kViolet, ROOT.kMagenta], 
    "FM": [ROOT.kRed, ROOT.kOrange], 
    "FS": [ROOT.kBlue, ROOT.kCyan],  # Changed from kCyan to kAzure+7 for more contrast
    "FT": [ROOT.kGreen+3, ROOT.kSpring-9]  # Changed from kTeal to kSpring-9 for more contrast
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


def plot_histograms(desired_num_bins, root_files_info, output_plot):
    # Open the first ROOT file to retrieve the list of parameters (branches)
    first_root_file_path = next(iter(root_files_info.values()))
    first_root_file = ROOT.TFile(first_root_file_path, "READ")
    keys = [key.GetName() for key in first_root_file.GetListOfKeys()]
    print(keys)
    first_root_file.Close()

    # Create a dictionary to map each operator to a color
# Create a dictionary to map each operator to a color
    operator_colors = {}
    for operator in all_ops_SM:
        category = operator[:2]  # Get the category from the operator name
        operators_in_category = [op for op in all_ops_SM if op.startswith(category)]  # Get all operators in this category
        color_index = operators_in_category.index(operator)  # Get the index of the operator within its category
        color1 = categories[category][0]  # Get the first color for this category
        color2 = categories[category][1]  # Get the second color for this category
        color = get_gradient_color(color1, color2, len(operators_in_category), color_index)  # Create a gradient color
        operator_colors[operator] = color  # Store the color


    for parameter_to_plot in keys:    
        # Create a new TCanvas
        canvas_name = "canvas_" + parameter_to_plot
        canvas = ROOT.TCanvas(canvas_name, "Stacked Histograms", 2000, 2000)

        # Create a THStack
        hs = ROOT.THStack("hs", "Distribution of " + parameter_to_plot+ ' \nnb of bins:'+str(desired_num_bins))
        
        # Create a legend
        legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
        
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
                
            op=legend_name.split("_")[-1]
            color = operator_colors[op]  # Get the color for this operator
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
        "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]



for process in ["WmZ"]:
    Root_paths = {}

    for decay in ["llqq"]:
        for op in all_ops_SM:
            if op=="SM":
                path = os.path.join(base_dir, f"{process}_{decay}/mc16_13TeV.*.MGPy8EG_aQGCFM0_{op}_1_{process}_{decay}.merge.EVNT.*/DOCUT_YES/Tables/Tables_file_first/hists.root")
            else:
                path = os.path.join(base_dir, f"{process}_{decay}/mc16_13TeV.*.MGPy8EG_aQGC{op}_QUAD_1_{process}_{decay}.merge.EVNT.*/DOCUT_YES/Tables/Tables_file_first/hists.root")
            matches = glob.glob(path)
            #print(matches)
            if not matches:
                print(f"No match found for operator {op} and process {process}")
                continue
            Root_paths[f"{process}_{decay}_{op}"] = matches[0]
    #print(Cutflow_paths)



    # Call the function with the desired parameters
    desired_num_bins = int(opts.bins)  # Replace with your desired number of bins

    output_plot = f"{base_dir_plot}/{process}/"

    # Check if the directory exists, if not, create it
    if not os.path.exists(output_plot):
        os.makedirs(output_plot)

    plot_canvas = plot_histograms(desired_num_bins, Root_paths,output_plot)



    