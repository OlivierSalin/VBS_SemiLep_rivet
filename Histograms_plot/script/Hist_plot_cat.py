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

base_dir_plot = "/exp/atlas/salin/ATLAS/VBS_mc/EFT_files_AMI/Plots_op2/"

categories = {
    "FM0_FM1_FM7": ["FM0", "FM1", "FM7", ROOT.kRed-2, ROOT.kPink+6],
    "FM2_FM3_FM4_FM5": ["FM2", "FM3", "FM4", "FM5", ROOT.kOrange+10, ROOT.kYellow-4],
    "FT0_FT1_FT2": ["FT0", "FT1", "FT2", ROOT.kBlue-9, ROOT.kViolet-6],
    "FT5_FT6_FT7": ["FT5", "FT6", "FT7", ROOT.kAzure, ROOT.kCyan],
    "cat5": ["FT8", "FT9", ROOT.kCyan-10, ROOT.kTeal-6],
    "FS1_FS02": ["FS1","FS02", ROOT.kGreen+3, ROOT.kSpring-9],
    "sm": ["SM", ROOT.kBlack, ROOT.kBlack]
}

categories_bis = {
    "FM0": [ROOT.kMagenta-4,1],
    "FM1": [ROOT.kPink+1,9],
    "FM7": [ROOT.kPink+6,7],
    "FM2": [ROOT.kRed+1,1],
    "FM3": [ROOT.kOrange+10,9],
    "FM4": [ROOT.kOrange-3,7],
    "FM5": [ROOT.kYellow-9,4],
    "FT0": [ROOT.kBlue-9,1],
    "FT1": [ROOT.kViolet-9,9],
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


def plot_histograms(desired_num_bins, root_files_info, output_plot):
    # Open the first ROOT file to retrieve the list of parameters (branches)
    first_root_file_path = next(iter(root_files_info.values()))
    first_root_file = ROOT.TFile(first_root_file_path, "READ")
    keys = [key.GetName() for key in first_root_file.GetListOfKeys()]
    print(keys)
    first_root_file.Close()

    # Create a dictionary to map each operator to a color
    operator_colors = {}
    for category, operators in categories.items():
        for operator in operators[:-2]:  # Exclude the last two items (colors)
            color_index = operators.index(operator)
            color1 = operators[-2]  # Second to last item is the first color
            color2 = operators[-1]  # Last item is the second color
            color = get_gradient_color(color1, color2, len(operators)-2, color_index)  # Create a gradient color
            operator_colors[operator] = color
            
    operator_colors_bis = {}
    operator_style_bis = {}
    for operator, values in categories_bis.items():
        color1 = values[0]  # First item is the color
        style1 = values[1]  # Second item is the style
        operator_colors_bis[operator] = color1
        operator_style_bis[operator] = style1

    for parameter_to_plot in keys:    
        for category, operators in categories.items():
            # Create a new TCanvas
            canvas_name = "canvas_" + parameter_to_plot + "_" + category
            canvas = ROOT.TCanvas(canvas_name, "Stacked Histograms", 2000, 2000)
            
            ROOT.gPad.SetRightMargin(0.2)
            # Create a THStack
            hs = ROOT.THStack("hs", "Distribution of " + parameter_to_plot+ ' \n\nnb of bins:'+str(desired_num_bins))
            
            # Create a legend
            legend = ROOT.TLegend(0.82, 0.75, 1.0, 0.9)
            output_dir = os.path.join(output_plot, category)
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            output_plot2 = os.path.join(output_dir, "")
            # Loop over the files to retrieve and stack the histograms
            for legend_name, file_path in root_files_info.items():
                op = legend_name.split("_")[-1]
                if op not in operators[:-2]:  # Only plot histograms for operators in the current category
                    continue
                
                # Open the ROOT file
                root_file = ROOT.TFile(file_path, "READ")
            
                # Retrieve the histogram
                histogram = root_file.Get(parameter_to_plot)
                histogram.SetDirectory(0)  
                current_num_bins = histogram.GetNbinsX()
                
                rebin_factor = int(current_num_bins / desired_num_bins)
                if rebin_factor > 0:
                    histogram.Rebin(rebin_factor)
                    
                color = operator_colors_bis[op]  # Get the color for this operator
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
            if hs.GetHists():
                hs.Draw("nostack")
                hs.GetXaxis().SetTitle(parameter_to_plot )
                hs.GetYaxis().SetTitle("Events / bin")
                
                # Set log scale for y-axis
                #canvas.SetLogy()
                #ROOT.gPad.SetLogy()
                
                # Draw the legend
                legend.Draw()
                #canvas.SetCanvasSize(1000, 1000)
                #ROOT.gPad.SetLogy()
                # Update the canvas
                canvas.Update()
                
                # Save the canvas to a file
                canvas.SaveAs(output_plot2+parameter_to_plot  + "_hist.png")
            #canvas.SaveAs(output_plot+parameter_to_plot + "_hist.pdf")

# Call the function with Root_paths as argument

base_dir = "/exp/atlas/salin/ATLAS/VBS_mc/EFT_files_AMI/"

all_ops_SM = ["SM","FM0","FM1","FM2","FM3","FM4","FM5","FM7",
        "FS02","FS1",
        "FT0","FT1","FT2","FT5","FT6","FT7"]

all_ops_cat = ["SM","FM0","FM2","FS1","FT1","FT5"]


for process in ["WpZ"]:
    Root_paths = {}

    for decay in ["llqq"]:
        for op in all_ops_SM:
            if op=="SM":
                path = os.path.join(base_dir, f"{process}_{decay}/mc16_13TeV.*.MGPy8EG_aQGCFM0_{op}_1_{process}_{decay}.merge.EVNT.*/DOCUT_YES/Tables/Tables_file_polar/hists.root")
            else:
                path = os.path.join(base_dir, f"{process}_{decay}/mc16_13TeV.*.MGPy8EG_aQGC{op}_QUAD_1_{process}_{decay}.merge.EVNT.*/DOCUT_YES/Tables/Tables_file_polar/hists.root")
            matches = glob.glob(path)
            #print(matches)
            if not matches:
                print(f"No match found for operator {op} and process {process}")
                continue
            Root_paths[f"{process}_{decay}_{op}"] = matches[0]
    #print(Cutflow_paths)



    # Call the function with the desired parameters
    desired_num_bins = int(opts.bins)  # Replace with your desired number of bins

    output_plot = f"{base_dir_plot}/Categorie/Linear/{process}_{decay}/"

    # Check if the directory exists, if not, create it
    if not os.path.exists(output_plot):
        os.makedirs(output_plot)

    plot_canvas = plot_histograms(desired_num_bins, Root_paths,output_plot)



    