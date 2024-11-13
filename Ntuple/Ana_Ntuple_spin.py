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
parser.add_option("--Full_op", default =False)
parser.add_option("--linear", default =True)

parser.add_option("--bins", default = 25)
opts, _ = parser.parse_args()

all_ops_SM = ["SM","FM0","FM1","FM2","FM3","FM4","FM5","FM7",
        "FS02","FS1",
        "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]

Complement_path="/Spin/NoSpin/"

if not opts.linear:
    base_dir_plot = f"/exp/atlas/salin/ATLAS/VBS_mc/plotting/Ntuple/Plots/Plot_stat_spin3/{Complement_path}/bins_{opts.bins}/"
else:
    base_dir_plot = f"/exp/atlas/salin/ATLAS/VBS_mc/plotting/Ntuple/Plots/Plot_stat_spin3/{Complement_path}/linear/bins_{opts.bins}/"

categories = {
    "cat1": ["FM0", "FM1", "FM7", ROOT.kRed-2, ROOT.kPink+6],
    "cat2": ["FM2", "FM3", "FM4", "FM5", ROOT.kOrange+10, ROOT.kYellow-4],
    "cat3": ["FT0", "FT1", "FT2", ROOT.kBlue-9, ROOT.kViolet-6],
    "cat4": ["FT5", "FT6", "FT7", ROOT.kAzure, ROOT.kCyan],
    "cat5": ["FT8", "FT9", ROOT.kCyan-10, ROOT.kTeal-6],
    "cat6": ["FS1","FS02", ROOT.kGreen+3, ROOT.kSpring-9],
    "cat7": ["SM", ROOT.kBlack, ROOT.kBlack]
}

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



def plot_histograms2(desired_num_bins, root_files_info, output_plot,tree_name="Merged"):
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
        # Create a new TCanvas
        canvas_name = "canvas_" + parameter_to_plot
        canvas = ROOT.TCanvas(canvas_name, "Stacked Histograms", 2000, 2000)

        ROOT.gPad.SetRightMargin(0.2)
        # Create a THStack
        hs = ROOT.THStack("hs", "Distribution of " + parameter_to_plot)

        legend = ROOT.TLegend(0.82, 0.6, 1.0, 0.9)
        min_hist=0
        max_hist=float('-inf')
        # Loop over the files to retrieve and stack the histograms
        for legend_name, file_path in root_files_info.items():
            # Open the ROOT file
            root_file = ROOT.TFile(file_path, "READ")
            tree = root_file.Get(tree_name)
            
            if not tree:
                print(f"Error: No TTree named 'Merged' found in file {file_path}. Skipping this file.")
                root_file.Close()
                continue
            
            min_hist=min(min_hist,tree.GetMinimum(parameter_to_plot))
            max_hist=max(max_hist,tree.GetMaximum(parameter_to_plot))

        for legend_name, file_path in root_files_info.items():
            # Open the ROOT file
            root_file = ROOT.TFile(file_path, "READ")
            tree = root_file.Get(tree_name)
            
            if not tree:
                print(f"Error: No TTree named 'Merged' found in file {file_path}. Skipping this file.")
                root_file.Close()
                continue
            
            tolerance = 0.1
            min_hist_=min_hist - tolerance*abs(min_hist)

            max_hist_=max_hist + tolerance*abs(max_hist)
            
            if parameter_to_plot == "merged_ll_mass":
                min_hist_=88
                max_hist_=94
            histogram = ROOT.TH1D("histogram", "title", desired_num_bins, min_hist_, max_hist_)

            tree.Draw(parameter_to_plot + ">>histogram", "", "norm")
            histogram.SetDirectory(0)  
            
            op=legend_name.split("_")[-1]
            color = operator_colors_bis[op]  # Get the color for this operator
            style = operator_style_bis[op]  # Get the color for this operator
            histogram.SetLineColor(color)
            histogram.SetLineWidth(3)
            #histogram.GetXaxis().SetTitle(legend_name)
            hs.Add(histogram)
            #histogram.SetLineStyle(1)

            # Add entry to the legend
            legend.AddEntry(histogram, legend_name, "lep")

            # Close the ROOT file
            root_file.Close()
            
                           
        hs.Draw("nostack HIST")
        X_axis_param = parameter_to_plot.lstrip("merged_")  
        hs.GetXaxis().SetTitle(X_axis_param)           
        # Set log scale for y-axis
        if not opts.linear:
            canvas.SetLogy()
            ROOT.gPad.SetLogy()

        # Draw the legend
        legend.Draw()
        canvas.SetCanvasSize(1000, 1000)
        #ROOT.gPad.SetLogy()
        # Update the canvas
        canvas.Update()
        
        # Save the canvas to a file
        canvas.SaveAs(output_plot+parameter_to_plot + "_hist.png")

# Call the function with Root_paths as argument

base_dir = "/exp/atlas/salin/ATLAS/VBS_mc/eft_files/"
base_dir = "/exp/atlas/salin/ATLAS/VBS_mc/EFT_files_AMI/"

all_ops_SM = ["SM","FM0","FM1","FM2","FM3","FM4","FM5","FM7",
        "FS02","FS1",
        "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]

all_ops_cat = ["FM0","FM2","FS1","FT1","FT5"]
all_ops_cat = ["FM0","FS1","FT1","FT5"]
processes=["WpZ","ZZ","WmZ"]
decays=["llqq"]
#all_ops_cat = ["SM"]

Complement_path="/Spin/Polarization/"
Complement_path="/Spin/NoSpin/"

if not opts.linear:
    base_dir_plot = f"/exp/atlas/salin/ATLAS/VBS_mc/plotting/Ntuple/Plots/Plot_stat_spin3/{Complement_path}/bins_{opts.bins}/"
else:
    base_dir_plot = f"/exp/atlas/salin/ATLAS/VBS_mc/plotting/Ntuple/Plots/Plot_stat_spin3/{Complement_path}/linear/bins_{opts.bins}/"


for process in processes:
    Root_paths = {}

    for decay in decays:
        if opts.Full_op:
            all_op_plot = all_ops_SM
        else:
            all_op_plot = all_ops_cat
        for op in all_op_plot:
            if op=="SM":
                path = os.path.join(base_dir, f"{process}_{decay}/mc16_13TeV.*.MGPy8EG_aQGCFM0_{op}_1_{process}_{decay}.merge.EVNT.*/DOCUT_YES/Tables/All_channel/Stats/combined/ntuple_rivet.root")
                #path = os.path.join(base_dir, f"{process}_{decay}/user.osalin.MadGraph_{process}_{decay}_FM0_SM_EXT0/DOCUT_YES/Tables/{Complement_path}/ntuple_rivet.root")
            else:
                path = os.path.join(base_dir, f"{process}_{decay}/mc16_13TeV.*.MGPy8EG_aQGC{op}_QUAD_1_{process}_{decay}.merge.EVNT.*/DOCUT_YES/Tables/All_channel/Stats/ntuple_rivet.root")
                #path = os.path.join(base_dir, f"{process}_{decay}/user.osalin.MadGraph_{process}_{decay}_{op}_QUAD_EXT0/DOCUT_YES/Tables/{Complement_path}/ntuple_rivet.root")
            matches = glob.glob(path)
            #print(matches)
            if not matches:
                print(f"No match found for operator {op} and process {process}")
                continue
            Root_paths[f"{process}_{decay}_{op}"] = matches[0]
    #print(Cutflow_paths)



    # Call the function with the desired parameters
    desired_num_bins = int(opts.bins)  # Replace with your desired number of bins
    
    if opts.Full_op:
        output_plot = f"{base_dir_plot}/full_op/{process}/"
    else:   
        output_plot = f"{base_dir_plot}/fewer_few_only_op/{process}/"
        #output_plot = f"{base_dir_plot}/only_op_FT_FS_FM/{process}/"
        #output_plot = f"{base_dir_plot}/fewer_op_less/{process}/"
        #output_plot = f"{base_dir_plot}/SM_only/{process}/"

    # Check if the directory exists, if not, create it
    if not os.path.exists(output_plot):
        os.makedirs(output_plot)

    plot_canvas = plot_histograms2(desired_num_bins, Root_paths,output_plot,tree_name="Angle")
    plot_canvas = plot_histograms2(desired_num_bins, Root_paths,output_plot,tree_name="Merged")