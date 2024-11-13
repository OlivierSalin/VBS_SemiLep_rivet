import ROOT
import os
from itertools import product
import array


from optparse import OptionParser
parser = OptionParser()
parser.add_option("--Ana", default = "WpZ_llqq")
parser.add_option("--DOCUT", default = "YES")
parser.add_option("--Full_op", default =False)
parser.add_option("--linear", default =True)

parser.add_option("--bins", default = 50)
opts, _ = parser.parse_args()
# Open the ROOT file
file_path = "Truth_eta/Truth_aQGC/WmZ_llqq_ntuple_truth.root"


    
def plot_histograms1d_quartil(root_files_path, output_plot, tree_name, desired_num_bins, legend_name="WmZ_llqq", algo="", linear=False):
    root_file = ROOT.TFile(root_files_path, "READ")
    tree = root_file.Get(tree_name)
    ROOT.gROOT.SetBatch(True)
    bin_num=desired_num_bins
    if not tree:
        print(f"Error: No TTree named 'Truth' found in file {root_files_path}.")
        return
    keys = [branch.GetName() for branch in tree.GetListOfBranches() if branch.GetName() != "Event Number"]
    print(keys)
    if not os.path.exists(output_plot):
        os.makedirs(output_plot)

    for parameter_to_plot in keys:
        if tree.GetBranch(parameter_to_plot) is None:
            print(f"Error: No branch named '{parameter_to_plot}' found in the tree.")
            continue

        canvas = ROOT.TCanvas(parameter_to_plot, parameter_to_plot, 1000, 1000)
        
        min_hist = tree.GetMinimum(parameter_to_plot) - 0.1 * abs(tree.GetMinimum(parameter_to_plot))
        max_value = 1.1 * tree.GetMaximum(parameter_to_plot)
        # Draw the parameter into a temporary histogram to calculate the 95th percentile
        temp_hist_name = "tempHist"
        tree.Draw(f"{parameter_to_plot}>>{temp_hist_name}(1000,{min_hist},{max_value})")
        temp_hist = ROOT.gDirectory.Get(temp_hist_name)
        quantiles = array.array('d', [0.95])
        qvalues = array.array('d', [0.])
        temp_hist.GetQuantiles(1, qvalues, quantiles)
        max_hist = qvalues[0]  # Use the 95th percentile as max_hist

        histogram_title = f"{tree_name} for Algo {algo} Distribution of {parameter_to_plot}"
        
        if "matched" in parameter_to_plot:
            desired_num_bins = 6
            
        elif "Label" in parameter_to_plot:
            desired_num_bins = 10
        else:
            desired_num_bins = bin_num
            

        histogram = ROOT.TH1F(histogram_title, histogram_title, desired_num_bins, min_hist, max_hist)
        tree.Draw(f"{parameter_to_plot}>>{histogram_title}")

        # Normalize the histogram
        integral = histogram.Integral(0, -1)
        if integral > 0:  # Avoid division by zero
            histogram.Scale(1.0 / integral)

        histogram.SetLineWidth(3)
        histogram.Draw("HIST")
        canvas.Update()

        # Adjust legend position
        legend = ROOT.TLegend(0.9, 0.75, 1.0, 0.9)  # Move legend to the right of the plot
        legend.AddEntry(histogram, legend_name, "l")
        legend.Draw()

        if not linear:
            ROOT.gPad.SetLogy()

        canvas.SaveAs(os.path.join(output_plot, algo + "_" + parameter_to_plot + ".png"))

    root_file.Close()
    
def plot_histograms1d(root_files_path, output_plot,tree_name,desired_num_bins,legend_name="WmZ_llqq",algo="",linear=False):
    root_file = ROOT.TFile(root_files_path, "READ")
    tree = root_file.Get(tree_name)  # Get the TTree named "Merged"
    ROOT.gROOT.SetBatch(True)
    bin_num=desired_num_bins
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

        canvas = ROOT.TCanvas(parameter_to_plot, parameter_to_plot, 1000, 1000)
        #canvas.SetBatch(True)
        #ROOT.gStyle.SetOptStat(0)
        
        
        min_hist= tree.GetMinimum(parameter_to_plot) - 0.1*abs(tree.GetMinimum(parameter_to_plot))
        max_hist= tree.GetMaximum(parameter_to_plot) + 0.1*abs(tree.GetMaximum(parameter_to_plot))
        
        if algo == "Vhad_first_bis":
            algo = "Vhad_first"
        
        histogram_title = f"{tree_name}: Distribution of {parameter_to_plot}"
        #histogram = ROOT.TH1F(parameter_to_plot, parameter_to_plot, desired_num_bins, min_hist, max_hist)
        
        if "matched" in parameter_to_plot:
            desired_num_bins = 6
            max_hist=3
        elif "label" in parameter_to_plot:
            desired_num_bins = 10
        if "n_EligibleJets" in parameter_to_plot:
            min_hist = 0
            max_hist = 8
            desired_num_bins = 8
        if "Delta_pT_Jets_q_candidate_NotMatched_q" in parameter_to_plot:
            min_hist = 0
            max_hist = 250
            desired_num_bins = 50

        else:
            desired_num_bins = bin_num
            
        histogram = ROOT.TH1F(histogram_title, histogram_title, desired_num_bins, min_hist, max_hist)
        histogram.SetXTitle(parameter_to_plot)  
        tree.Draw(f"{parameter_to_plot}>>{histogram_title}")
        
        # Normalize the histogram
        integral = histogram.Integral(0,-1)
        if integral > 0:  # Avoid division by zero
            histogram.Scale(1.0 / integral)
            
        histogram.SetLineWidth(3)  

        histogram.Draw("HIST")
        canvas.Update()

        # Adjust legend position
        legend = ROOT.TLegend(0.9, 0.75, 1.0, 0.9)  # Move legend to the right of the plot        
       # Adjust these coordinates as needed
        legend.AddEntry(histogram, legend_name, "l")
        #legend.Draw()

        if not linear:
            ROOT.gPad.SetLogy()

        canvas.SaveAs(os.path.join(output_plot, algo + "_"+ parameter_to_plot + ".png"))

    root_file.Close()


def plot_histograms1d_op(root_files_info, output_plot,tree_name,desired_num_bins,legend_name="WmZ_llqq",algo="",linear=False):
    first_root_file_path = next(iter(root_files_info.values()))
    first_root_file = ROOT.TFile(first_root_file_path, "READ")
    ROOT.gROOT.SetBatch(True)
    bin_num=desired_num_bins
    tree = first_root_file.Get(tree_name)  # Get the TTree named "Merged"
    if not tree:
        print(f"Error: No TTree named 'Merged' found in file {first_root_file_path}.")
        return
    keys = [branch.GetName() for branch in tree.GetListOfBranches() if branch.GetName() != "Event Number"]  # Get the list of branches
    print(keys)
    first_root_file.Close()
    if not os.path.exists(output_plot):
        os.makedirs(output_plot)
    colors = [ROOT.kBlue,ROOT.kRed,ROOT.kGreen+4,ROOT.kBlack,ROOT.kOrange,ROOT.kMagenta,ROOT.kCyan,ROOT.kYellow,ROOT.kViolet,ROOT.kAzure,ROOT.kTeal,ROOT.kSpring,ROOT.kPink,ROOT.kGray,ROOT.kWhite]
    for parameter_to_plot in keys:    
        # Create a new TCanvas
        canvas_name = "canvas_" + parameter_to_plot
        canvas = ROOT.TCanvas(canvas_name, "Stacked Histograms", 2000, 2000)
        #canvas.SetBatch(True)

        ROOT.gPad.SetRightMargin(0.2)
        # Create a THStack
        hs = ROOT.THStack("hs", f"Distribution of " + parameter_to_plot)

        legend = ROOT.TLegend(0.82, 0.8, 1.0, 0.9)
        min_hist=0
        max_hist=float('-inf')
        
        quantiles = 0.95  # 95th percentile
        q = array.array('d', [0.])  # Array to store the quantile value
        probSum = array.array('d', [quantiles])  # Quantile position
        
        # Loop over the files to retrieve and stack the histograms
        for legend_name, file_path in root_files_info.items():
            # Open the ROOT file
            root_file = ROOT.TFile(file_path, "READ")
            tree = root_file.Get(tree_name)
            
            if not tree:
                print(f"Error: No TTree named 'Truth' found in file {file_path}. Skipping this file.")
                root_file.Close()
                continue
            
            tree.Draw(parameter_to_plot + ">>hTemp", "", "goff")  # Draw with option "goff" to not display the histogram
            hTemp = ROOT.gDirectory.Get("hTemp")
            hTemp.GetQuantiles(1, q, probSum)
            
            min_hist = min(min_hist, hTemp.GetMinimum())
            max_hist = max(max_hist, q[0])  # Use the quantile value as the relevant maximum
            root_file.Close()

        i=0
        tolerance = 0.1
        min_hist_=min_hist - tolerance*abs(min_hist)

        max_hist_=max_hist + tolerance*abs(max_hist)
        
        if "matched" in parameter_to_plot:
            desired_num_bins = 6
            max_hist_=3
        elif "Label" in parameter_to_plot:
            desired_num_bins = 10
        elif "Tagging_jets_VBS_q_truth" in parameter_to_plot:
            min_hist_ = -max_hist_
        elif "absDelta" in parameter_to_plot:
            min_hist_ = 0
        else:
            desired_num_bins = bin_num


        for legend_name, file_path in root_files_info.items():
            # Open the ROOT file
            root_file = ROOT.TFile(file_path, "READ")
            tree = root_file.Get(tree_name)
            
            if not tree:
                print(f"Error: No TTree named 'Truth' found in file {file_path}. Skipping this file.")
                root_file.Close()
                continue
            


            #print(min_hist_,max_hist_)

            histogram = ROOT.TH1D("histogram", "title", desired_num_bins, min_hist_, max_hist_)

            tree.Draw(parameter_to_plot + ">>histogram", "", "norm")
            
            integral = histogram.Integral(0,-1)
            #print(f"Integral of {legend_name} for {parameter_to_plot}: {integral}")
                
            histogram.SetDirectory(0)  
            histogram.SetLineColor(colors[i])
            

            histogram.SetLineWidth(3)
            #histogram.GetXaxis().SetTitle(legend_name)
            hs.Add(histogram)
            #histogram.SetLineStyle(1)

            # Add entry to the legend
            legend.AddEntry(histogram, legend_name, "lep")
            i=i+1

            # Close the ROOT file
            root_file.Close()
            
                           
        hs.Draw("nostack HIST")
        X_axis_param = parameter_to_plot.lstrip("merged_")  
        hs.GetXaxis().SetTitle(X_axis_param)           
        # Set log scale for y-axis
        if not linear:
            canvas.SetLogy()
            ROOT.gPad.SetLogy()

        # Draw the legend
        legend.Draw()
        canvas.SetCanvasSize(1000, 1000)

        canvas.SaveAs(os.path.join(output_plot, algo + "_"+ parameter_to_plot + ".png"))

    root_file.Close()
    
def plot_histograms1d_tree(root_file_path, output_plot, tree_names, desired_num_bins, algo="", linear=False):
    # Open the ROOT file
    root_file = ROOT.TFile(root_file_path, "READ")
    ROOT.gROOT.SetBatch(True)
    print(output_plot)
    bin_num=desired_num_bins
    if not os.path.exists(output_plot):
        os.makedirs(output_plot)
    colors = [ROOT.kRed, ROOT.kGreen+2,ROOT.kBlue,  ROOT.kBlack, ROOT.kOrange, ROOT.kMagenta, ROOT.kCyan, ROOT.kYellow, ROOT.kViolet, ROOT.kAzure, ROOT.kTeal, ROOT.kSpring, ROOT.kPink, ROOT.kGray, ROOT.kWhite]

    # Assuming all trees have similar branches, get branches from the first tree to setup histograms
    first_tree = root_file.Get(tree_names[0])
    if not first_tree:
        print(f"Error: No TTree named '{tree_names[0]}' found in file {root_file_path}.")
        return
    keys = [branch.GetName() for branch in first_tree.GetListOfBranches() if branch.GetName() != "Event Number"]

    for parameter_to_plot in keys:
        global_min = float('inf')
        global_max = -float('inf')
        

        # Calculate global min and max using mean and standard deviation
        for tree_name in tree_names:
            tree = root_file.Get(tree_name)
            if not tree:
                continue

            tree.Draw(f"{parameter_to_plot}>>hist_{tree_name}", "", "goff")


            min_val= tree.GetMinimum(parameter_to_plot) - 0.1*abs(tree.GetMinimum(parameter_to_plot))
            #print("min val:",min_val)
            max_val= tree.GetMaximum(parameter_to_plot) + 0.1*abs(tree.GetMaximum(parameter_to_plot))
            global_min = min(global_min, min_val)
            global_max = max(global_max, max_val)


        # Adjust for a small margin
        margin = (global_max - global_min) * 0.05
        min_hist = global_min - margin
        max_hist = global_max + margin
        
        
        
        if "tag_jets_matched" in parameter_to_plot:
            desired_num_bins = 6
            max_hist=3
            min_hist=0
        if "label" in parameter_to_plot:
            desired_num_bins = 10
        elif "Tagging_jets_VBS_q_truth" in parameter_to_plot:
            min_hist = -max_hist

            


        canvas_name = "canvas_" + parameter_to_plot
        canvas = ROOT.TCanvas(canvas_name, "Comparison of Trees", 800, 800)
        ROOT.gPad.SetRightMargin(0.2)
        hs = ROOT.THStack("hs", f"Distribution of " + parameter_to_plot)
        legend = ROOT.TLegend(0.82, 0.8, 1.0, 0.9)

        for i, tree_name in enumerate(tree_names):
            tree = root_file.Get(tree_name)
            if not tree:
                print(f"Error: No TTree named '{tree_name}' found in file {root_file_path}. Skipping this tree.")
                continue

            histogram_name = f"hist_{tree_name}_{parameter_to_plot}"
            histogram = ROOT.TH1D(histogram_name, parameter_to_plot, desired_num_bins, min_hist, max_hist)
            tree.Draw(f"{parameter_to_plot}>>{histogram_name}", "", "goff")
            
            if histogram.Integral() != 0:  # Avoid division by zero
                histogram.Scale(1.0 / histogram.Integral())
    
            
            histogram.SetLineColor(colors[i % len(colors)])
            histogram.SetLineWidth(3)
            hs.Add(histogram)
            legend.AddEntry(histogram, tree_name, "l")

        hs.Draw("nostack HIST")
        hs.GetXaxis().SetTitle(parameter_to_plot)
        if not linear:
            canvas.SetLogy()
        legend.Draw()
        canvas.SaveAs(os.path.join(output_plot, f"{algo}_{parameter_to_plot}.png"))

    root_file.Close()
    
def plot_histograms2d(root_files_path, output_plot,tree_name, desired_num_bins, pairs, legend_name="WmZ_llqq",algo="", linear=False):
    root_file = ROOT.TFile(root_files_path, "READ")
    ROOT.gStyle.SetOptStat(0)
    tree = root_file.Get(tree_name)
    ROOT.gROOT.SetBatch(True)
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

        canvas = ROOT.TCanvas(f"{x_parameter}_{y_parameter}", f"{legend_name} {algo} for {tree_name} \n {x_parameter} vs {y_parameter}", 1200, 800)
         
        x_min = tree.GetMinimum(x_parameter)
        x_max = tree.GetMaximum(x_parameter)
        x_range = x_max - x_min
        x_tolerance = x_range * 0.05  # 5% tolerance, adjust as needed

        # Calculate the range and apply a tolerance for y_parameter
        y_min = tree.GetMinimum(y_parameter)
        y_max = tree.GetMaximum(y_parameter)
        y_range = y_max - y_min
        y_tolerance = y_range * 0.05  # 5% tolerance, adjust as needed

        # Create the histogram with the adjusted min and max
        histogram = ROOT.TH2F(f"{x_parameter}_{y_parameter}", f"{legend_name} Tagging jet {algo} for {tree_name}", 
                            desired_num_bins, x_min - x_tolerance, x_max + x_tolerance,
                            desired_num_bins, y_min - y_tolerance, y_max + y_tolerance)



        histogram.SetXTitle(x_parameter)
        histogram.SetYTitle(y_parameter)
        
        # Adjust title font size (and possibly style)
        histogram.GetXaxis().SetTitleSize(0.04)  # Adjust X-axis title size
        histogram.GetYaxis().SetTitleSize(0.04)  # Adjust Y-axis title size
        histogram.GetZaxis().SetTitleSize(0.04)  # Adjust Z-axis title size, if needed
        histogram.SetTitleSize(0.04, "t")  # Adjust the main title size
        
        tree.Draw(f"{y_parameter}:{x_parameter}>>{x_parameter}_{y_parameter}", "", "COLZ")

        # Normalize the histogram
        integral = histogram.Integral()
        if integral > 0:  # Avoid division by zero
            histogram.Scale(1.0 / integral)
            
        canvas.Update()  # Force canvas update to ensure the stats box is created

        if not linear:
            ROOT.gPad.SetLogz()

        canvas.SaveAs(os.path.join(output_plot, f"{algo}_{x_parameter}_{y_parameter}.png"))

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

processes = ["WmZ"]
decays = ["llqq"]

all_ops_cat = ["SM","FM0","FM2","FS1","FT1","FT5"]

detail_path=f"/Tables/Mismatch_VBS/Highest_eta/2Lepton/"
run_dir=f"/exp/atlas/salin/ATLAS/VBS_mc/eft_files/{detail_path}/Combined/"
ouptut_dir=f"/exp/atlas/salin/ATLAS/VBS_mc/eft_files/{detail_path}/Plot/"

Algo_name=["Highest_mass","Highest_eta","Vhad_first_"]
Algo_name=["Highest_eta"]
#Algo_name=["Eligible_"]
#Algo_name=["Highest_mass_"]
Tree_name=["Truth","Truth_mismatch","Truth_matched","Bef_cut","Lep_cut"
           ,"Truth_mismatch_noVhad_VBS_found","Truth_mismatch_Vhad_VBS_found","Truth_Merged","Truth_merged_mismatch"]

Tree_name=["Truth","Truth_mismatch","Truth_matched","Lep_cut"
           ,"Truth_mismatch_noVhad_VBS_found","Truth_mismatch_Vhad_VBS_found","Truth_Merged","Truth_merged_mismatch"]
Tree_name= ["Truth_mismatch_noVhad_VBS_found_DR1"]
Tree_name_=["Truth_mismatch","Truth_matched","Truth_mismatch_noVhad"]
Tree_name_=["Truth_mismatch_Vhad_VBS_found","Truth_mismatch_noVhad_VBS_found"]
Tree_name_=["Truth_mismatch_noVhad_VBS_found","Truth_merged_mismatch"]
Tree_name_=["Truth","Truth_mismatch",'Truth_matched']

option_SM_aQGC=["SM","aQGC"]
File_path_={}

detail_path_=f"/Tables/Mismatch_VBS/2Lepton/"
desired_num_bins=int(opts.bins)
for desired_num_bins in [6,50]:
    for tree_name in Tree_name:
        for algo in Algo_name:
            
            detail_path_algo=f"/Tables/Mismatch_VBS/{algo}/2Lepton/"
            
            run_dir=f"/exp/atlas/salin/ATLAS/VBS_mc/eft_files/{detail_path_algo}/Combined/"
            ouptut_dir=f"/exp/atlas/salin/ATLAS/VBS_mc/eft_files/{detail_path_algo}/Plot/"
            Folder= run_dir
            
            file_path = Folder+ "/All_aQGC_ntuple_rivet.root"
            File_path_["aQGC_"+algo]=file_path
            output_plot_1D=f"{ouptut_dir}/Plot_1D/{tree_name}/bins_{desired_num_bins}/"
            output_plot_1D_SM_aQGC=f"{ouptut_dir}/Plot_1D_Algo/"


            process_op="aQGC"  
            plot_histograms1d(file_path, output_plot_1D,tree_name,desired_num_bins,legend_name=process_op,algo=algo,linear=True)

            output_plot_2D=f"{ouptut_dir}/Plot_2D/{tree_name}/"
            desired_num_bins_2D=20



            targets = ['Delta_Eta_Tagging_jets_VBS_q_truth', 'Delta_Mass_Tagging_jets_VBS_q_truth', "tag_jets_matched","Tagging_Jets_mass", "Tagging_Jets_delta_Eta"]
            parameters = ['VBS_Quarks_truth_mass', 'VBS_Quarks_truth_Delta_Eta', 'Tagging_Jets_mass', 'W_pT']
            
            pairs = [list(pair) for pair in product(parameters, targets) if pair[0] != pair[1]]
            
            #plot_histograms2d(file_path, output_plot_2D,tree_name, desired_num_bins_2D, pairs, legend_name=process_op,algo=algo, linear=False)
        
        output_plot_1D_SM_aQGC=f"/exp/atlas/salin/ATLAS/VBS_mc/eft_files/{detail_path_}/Plot/Plot_1D_Algo_/{tree_name}/bins_{desired_num_bins}/"
        #plot_histograms1d_op(File_path_, output_plot_1D_SM_aQGC,tree_name, desired_num_bins,legend_name=process_op,algo=algo,linear=True) 


for algo in Algo_name:
    for bins_number in [4,6,50]:
        detail_path_algo=f"/Tables/Mismatch_VBS/{algo}/2Lepton/"
        
        run_dir=f"/exp/atlas/salin/ATLAS/VBS_mc/eft_files/{detail_path_algo}/Combined/"
        ouptut_dir=f"/exp/atlas/salin/ATLAS/VBS_mc/eft_files/{detail_path_algo}/Plot/"
        file_path = run_dir+ "/All_aQGC_ntuple_rivet.root"
        
        file_path = Folder+ "/All_aQGC_ntuple_rivet.root"
        process_op="aQGC"
        output_dir_tree=f"{ouptut_dir}/Plot_1D_tree_merged_noVhad_mis/bins_{bins_number}/"
        #plot_histograms1d_tree(file_path, output_dir_tree,Tree_name_, bins_number,algo=algo,linear=True)
        