# run rivet together with rivet plotting afterwards, if files were already downloded
import subprocess
import lib_utils as lu
import utils_func as uf
import os
import shutil
import run_chain
import yoda
import glob
import ROOT
from datetime import datetime
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(ROOT.kTRUE)
from array import array
import math
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--conf", default = "user.okurdysh.MadGraph_WmWm_lvlv_FT0_FULL")
parser.add_option("--DOCUT", default = "NO")
parser.add_option("--evtMax", default = 100)
parser.add_option("--redoRivet", default = "no")
parser.add_option("--redoPlots", default = "no")
parser.add_option("--runAgain", default = "yes")
parser.add_option("--name", default = "")
parser.add_option("--keep", default = False)
parser.add_option("--part", default = "")
parser.add_option("--type_MC", default = "")
opts, _ = parser.parse_args()

prod_dec, base_dir = lu.find_prod_dec_and_dir_tres(opts.conf,opts.type_MC)

print(f'prod_dec: {prod_dec}, base_dir: {base_dir}')

conf_dir, _, _ = lu.find_evnt_dir_and_file_bis(base_dir,opts.conf)
print(f'conf_dir: {conf_dir}')
conf_cut_dir = lu.get_conf_cut_dir(conf_dir, opts.DOCUT)
print(f'conf_cut_dir: {conf_cut_dir}')

if opts.part != "": 
    print("will run chain for part", opts.part)
    conf_dir = conf_dir + f"/{opts.part}/"
    


rivet_out_name = conf_cut_dir + f'/MyOutput.yoda.gz'

EFT_op, EFT_type, proc, decay = uf.extract_EFT_op_proces_dec_bis(opts.conf)
print(f'EFT_op: {EFT_op}, EFT_type: {EFT_type}, proc: {proc}, decay: {decay}')


keyy = f"{EFT_op}_{EFT_type}_{proc}_{decay}"
print(f'keyy: {keyy}')    
#xsection_fb = uf.cross_section_fb(EFT_op,EFT_type, proc, decay)
if opts.type_MC == "model" or opts.type_MC == "aqgc":
    print("Will take cross section from txt")
    VBS_txt='VBS_cross_section_aqgc.txt'
    xsection_fb = uf.take_xsec_fb_aqgc(VBS_txt,EFT_op,EFT_type, proc, decay)
else:
    VBS_txt='VBS_xsection_test.txt'
    xsection_fb = uf.take_xsec_fb(VBS_txt,EFT_op,EFT_type, proc, decay)


if xsection_fb is None:
    print(f'No cross section in txt found for {keyy}')
    
    xsection_fb=uf.cross_section_fb(EFT_op,EFT_type, proc, decay)
    if xsection_fb is None:
        print(f'No cross section found in AMI for {keyy}')
        print("Cross section set at 0 to avoid errors but ISSUE")
        
print(f'Cross section in fb for {keyy}: {xsection_fb}')

mydir=conf_dir + f"/DOCUT_{opts.DOCUT}/"

print(f'mydir: {mydir}')

plots_backup_dir = mydir + "Plots_backup/"
os.makedirs(plots_backup_dir, exist_ok=True)

label_plot= f"{proc}_{decay}_{EFT_op}"

dir_plot_ = mydir + "/plots/"
dir_plot_2 = mydir + "/plots/ntuple/"

# Check if the directory exists
if os.path.exists(dir_plot_):
    # Delete the directory
    shutil.rmtree(dir_plot_)
    
#uf.plot_histograms(output_plot=dir_plot_ , desired_num_bins=200, file_path=mydir + "/hists.root", label=label_plot)

#uf.plot_histograms(output_plot=dir_plot_2 , desired_num_bins=200, file_path=mydir + "/ntuple_rivet.root", label=label_plot)


if (opts.keep):
    timestamp_day = datetime.now().strftime("%d_%m")
    day_dir = plots_backup_dir+timestamp_day+ "/"
    os.makedirs(plots_backup_dir+timestamp_day+ "/", exist_ok=True)
    timestamp_hour = datetime.now().strftime("_%Hh%M")
    run_dir = day_dir + opts.name + timestamp_hour + "/"
    os.makedirs(run_dir, exist_ok=True)
    print(f'run_dir: {run_dir}')


    # Copy the repository mydir + "/plots/"
    #shutil.copytree(mydir + "/plots/", run_dir + "plots/")
#print(f'Cross section in fb: {xsec_fb}')


