import glob
import os
from array import array
import ROOT
import matplotlib.pyplot as plt
import shutil
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

output_plot_cutflow = f"{base_dir}/{opts.Ana}/Plots/{opts.QUAD1}_{opts.QUAD2}/DOCUT_{opts.DOCUT}/Cutflow_comp/"
# Check if the directory exists, if not, create it
if not os.path.exists(output_plot_cutflow):
    os.makedirs(output_plot_cutflow)
    
Cutflow_paths = {
    "SM": f"{base_dir}/{opts.Ana}/user.osalin.MadGraph_{opts.Ana}_FM0_SM_EXT0/DOCUT_{opts.DOCUT}/",
    opts.QUAD1: f"{base_dir}/{opts.Ana}/user.osalin.MadGraph_{opts.Ana}_{opts.QUAD1}_QUAD_EXT0/DOCUT_{opts.DOCUT}/",
    opts.QUAD2: f"{base_dir}/{opts.Ana}/user.osalin.MadGraph_{opts.Ana}_{opts.QUAD2}_QUAD_EXT0/DOCUT_{opts.DOCUT}/",
}



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