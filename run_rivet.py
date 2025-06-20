# run rivet together with rivet plotting afterwards, if files were already downloded
import subprocess
import lib_utils as lu
import utils_func as uf
import os
import shutil
import run_chain
import yoda
import rivet
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

if opts.part != "": 
    print("will run chain for part", opts.part)
    conf_dir = conf_dir + f"/{opts.part}/"
    
conf_cut_dir = lu.get_conf_cut_dir(conf_dir, opts.DOCUT)

rivet_out_name = conf_cut_dir + f'/MyOutput.yoda.gz'
if not os.path.exists(rivet_out_name): do_rivet = 1
elif os.path.exists(rivet_out_name) and opts.redoRivet=="yes": do_rivet = 1
else: do_rivet = 0
if do_rivet:
    run_com = "athena rivet_job.py -c 'conf=" + f'"{opts.conf}";DOCUT=' + f'"{opts.DOCUT}";Part=' + f'"{opts.part}";type_MC=' + f'"{opts.type_MC}"' + f"' --evtMax {opts.evtMax} "

    print("#### will run rivet with", run_com)
    subprocess.call(run_com, shell=True)
else:
    print("dont do evnt conversion since redoRivet=", opts.redoRivet, f"and file {rivet_out_name} exists ", os.path.exists(rivet_out_name))

plots_dir = conf_cut_dir + "/rivet-plots/"
if opts.redoPlots=="no": do_plots = 0
elif opts.redoPlots=="yes": do_plots = 1
else: do_plots = 0
if do_plots:
    plot_com = f"rivet-mkhtml  MyOutput.yoda.gz:'Title={prod_dec}' --no-ratio"
    print("#### will run mkhtml in dir", conf_cut_dir, "with com", plot_com)
    subprocess.call(plot_com, shell=True, cwd = conf_cut_dir)
else:
    print("dont do plot from yoda since redoPlots=", opts.redoPlots, f"and dir {plots_dir} exists ", os.path.exists(plots_dir))


def write_to_f(product_file,product):
    f = open(product_file, "w") #since opening yoda a bit slow for 50 configs save when have it
    f.write(str(product))
    f.close()
    
def save_job_infos(DOCUT_str, mydir, prod_dec,xsec_fb):
    yoda_f_str = mydir + "MyOutput.yoda.gz"

    if not os.path.exists(yoda_f_str): 
        print("dont see yoda file in dir ", mydir, ",return")
        return
    
    print("reading from yoda file ", yoda_f_str)
    yoda_f = yoda.read(yoda_f_str)
    print("Weight name:")
    print(set(name for name in yoda_f if 'merged_VlepVhad_mass' in name))
    weight_names= set([rivet.extractWeightName(name) for name in yoda_f])
    #print({rivet.extractWeightName(name) for name in yoda_f if 'quad' in rivet.extractWeightName(name) or 'cross' in rivet.extractWeightName(name)})
    #print({rivet.extractWeightName(name) for name in yoda_f })
    #print("reading from yoda file ", yoda_f_str)


    yoda_f = yoda.read(yoda_f_str)
    #print("reading from yoda file ", yoda_f_str)
    #print("hist yoda keys", [name for name in yoda_f.keys() if "RAW" not in name] )
    all_hists_in_yoda = [iname  for iname in yoda_f.keys() if "RAW" not in iname]
    hists_1h_in_yoda = []
    for i_name in all_hists_in_yoda:
        if yoda_f[i_name].type()=="Histo1D": hists_1h_in_yoda.append(i_name)  
    #print("have 1d hists to be saved in root:", hists_1h_in_yoda, "in yoda file", yoda_f_str)



    root_file = mydir + "/hists.root"
    if os.path.exists(root_file):
        os.remove(root_file)
    
    # save fid xsec
    rivet_dir_name = f"/{prod_dec}:OUTDIR=/{mydir}".replace("//","/")
    #rivet_dir_name = f"/{prod_dec}".replace("//","/")
    print("looking for prefix in counter",rivet_dir_name)
    pos_n_in = yoda_f[f"{rivet_dir_name}/pos_w_initial"].numEntries()
    neg_n_in = yoda_f[f"{rivet_dir_name}/neg_w_initial"].numEntries()
    pos_n_f = yoda_f[f"{rivet_dir_name}/pos_w_final"].numEntries()
    neg_n_f = yoda_f[f"{rivet_dir_name}/neg_w_final"].numEntries()
    
    pos_n_f_merged = yoda_f[f"{rivet_dir_name}/pos_w_final_merged"].numEntries()
    neg_n_f_merged = yoda_f[f"{rivet_dir_name}/neg_w_final_merged"].numEntries()
    
    pos_n_f_resolved = yoda_f[f"{rivet_dir_name}/pos_w_final_resolved"].numEntries()
    neg_n_f_resolved = yoda_f[f"{rivet_dir_name}/neg_w_final_resolved"].numEntries()
    #
    print(f'pos_n_in: {pos_n_in}, neg_n_in: {neg_n_in}, pos_n_f: {pos_n_f}, neg_n_f: {neg_n_f}', 
          f'pos_n_f_merged: {pos_n_f_merged}, neg_n_f_merged: {neg_n_f_merged}',
          f'pos_n_f_resolved: {pos_n_f_resolved}, neg_n_f_resolved: {neg_n_f_resolved}')
    
    frac_cut = (pos_n_f+neg_n_f) / (pos_n_in+neg_n_in) if (pos_n_in+neg_n_in) != 0 else 0
    frac_pos = pos_n_f / pos_n_in if pos_n_in != 0 else 0
    frac_neg = neg_n_f / neg_n_in if neg_n_in != 0 else 0
    frac_cut_er_bar = 1/(pos_n_in+neg_n_in) * math.sqrt(pos_n_in*frac_pos*(1-frac_pos) + neg_n_in*frac_neg*(1-frac_neg)) if (pos_n_in+neg_n_in) != 0 else 0
    
    frac_cut_merged = (pos_n_f_merged+neg_n_f_merged) / (pos_n_in+neg_n_in) if (pos_n_in+neg_n_in) != 0 else 0
    frac_pos_merged = pos_n_f_merged / pos_n_in if pos_n_in != 0 else 0
    frac_neg_merged = neg_n_f_merged / neg_n_in if neg_n_in != 0 else 0
    frac_cut_er_bar_merged = 1/(pos_n_in+neg_n_in) * math.sqrt(pos_n_in*frac_pos_merged*(1-frac_pos_merged) + neg_n_in*frac_neg_merged*(1-frac_neg_merged)) if (pos_n_in+neg_n_in) != 0 else 0
    
    frac_cut_resolved = (pos_n_f_resolved+neg_n_f_resolved) / (pos_n_in+neg_n_in) if (pos_n_in+neg_n_in) != 0 else 0
    frac_pos_resolved = pos_n_f_resolved / pos_n_in if pos_n_in != 0 else 0
    frac_neg_resolved = neg_n_f_resolved / neg_n_in if neg_n_in != 0 else 0
    frac_cut_er_bar_resolved = 1/(pos_n_in+neg_n_in) * math.sqrt(pos_n_in*frac_pos_resolved*(1-frac_pos_resolved) + neg_n_in*frac_neg_resolved*(1-frac_neg_resolved)) if (pos_n_in+neg_n_in) != 0 else 0
  
    # frac_cut_er_bar = frac_cut_unc / 2 
    #
    pos_w_in = yoda_f[f"{rivet_dir_name}/pos_w_initial"].sumW()
    neg_w_in = yoda_f[f"{rivet_dir_name}/neg_w_initial"].sumW()
    pos_w_f = yoda_f[f"{rivet_dir_name}/pos_w_final"].sumW()
    neg_w_f = yoda_f[f"{rivet_dir_name}/neg_w_final"].sumW()
    
    pos_w_f_resolved = yoda_f[f"{rivet_dir_name}/pos_w_final_resolved"].sumW()
    neg_w_f_resolved = yoda_f[f"{rivet_dir_name}/neg_w_final_resolved"].sumW()
    
    pos_w_f_merged = yoda_f[f"{rivet_dir_name}/pos_w_final_merged"].sumW()
    neg_w_f_merged = yoda_f[f"{rivet_dir_name}/neg_w_final_merged"].sumW()

    Dir_info = mydir + "/Info/"
    
    if not os.path.exists(Dir_info):
        # If not, create it
        os.makedirs(Dir_info)
    
    write_to_f(Dir_info + "frac_cuts_merged.txt", frac_cut_merged)
    write_to_f(Dir_info + "frac_after_cuts_error_bar_merged.txt", frac_cut_er_bar_merged)
    
    write_to_f(Dir_info + "frac_cuts_resolved.txt", frac_cut_resolved)
    write_to_f(Dir_info + "frac_after_cuts_error_bar_resolved.txt", frac_cut_er_bar_resolved)
    write_to_f(Dir_info + "Cross_section_fb.txt", xsec_fb)
    
    print("Write to files")


    # save hists in root for further plotting
    root_file = ROOT.TFile(root_file,"UPDATE")
    for i_hist in hists_1h_in_yoda: # they are in format '/WpWm_lvlv:DOCUT=YES/leptons_pids'
        h_yoda =  yoda_f[i_hist]
        h_root = lu.yoda_to_root_1d(h_yoda, i_hist.split("/")[-1])
        h_root.Write("", ROOT.TObject.kOverwrite)
    root_file.Close()

    os.makedirs(mydir + "/plots_bis/", exist_ok=True)

    uf.plot_histograms(output_plot=mydir + "/plots_bis/" , desired_num_bins=25, file_path=mydir + "/hists.root", label='M rwg')

    
    lu.save_xsec_frac_prod(Dir_info,xsec_fb,
                        frac_cut_merged,frac_pos_merged,frac_neg_merged, 
                        pos_w_in, neg_w_in, pos_w_f_resolved, neg_w_f_resolved,pos_w_f_merged, neg_w_f_merged,
                        pos_n_in, neg_n_in, pos_n_f_resolved, neg_n_f_resolved, pos_n_f_merged, neg_n_f_merged)
    
    print("Save frac files")
    
    #plot_root_histograms(mydir + "/hists.root")
    

    ############## draw event and cutflow

    print("saving cutflow resolved as img") 
    
    print("cutflow merged")
        
    cutflow_merged_file = mydir + "cutflow_merged.txt"
    if os.path.exists(cutflow_merged_file):
        cut_merged_names, cut_merged_cumu, cut_merged_incr = lu.get_cutflow_arrays(cutflow_merged_file)
        #lu.draw_cutflows(cut_merged_names, [cut_merged_incr,cut_merged_cumu], ["incremental","cumulative"],mydir+"/cutflow_merged_img.png", prod_dec) 
    print("saving cutflow merged ")    
    cutflow_resolved_file = mydir + "cutflow_resolved.txt"
    if os.path.exists(cutflow_resolved_file):
        cut_resolved_names, cut_resolved_cumu, cut_resolved_incr = lu.get_cutflow_arrays(cutflow_resolved_file)
        #lu.draw_cutflows(cut_resolved_names, [cut_resolved_incr,cut_resolved_cumu], ["incremental","cumulative"],mydir+"/cutflow_resolved_img.png", prod_dec)   


# prod_dec, _ = lu.find_prod_dec_and_dir(job_name)

#save_job_infos("DOCUT='" + f'{opts.DOCUT}', conf_dir + "/DOCUT='" + f'"{opts.DOCUT}"' + "/", prod_dec)

EFT_op, EFT_type, proc, decay = uf.extract_EFT_op_proces_dec_bis(opts.conf)

keyy = f"{EFT_op}_{EFT_type}_{proc}_{decay}" 
print(f'Process studied: {keyy}')   
#xsection_fb = uf.cross_section_fb(EFT_op,EFT_type, proc, decay)
if opts.type_MC == "model" or opts.type_MC == "aqgc":
    print("Will take cross section from aqgc model file")
    VBS_txt='VBS_cross_section_aqgc.txt'
    xsection_fb = uf.take_xsec_fb_aqgc(VBS_txt,EFT_op,EFT_type, proc, decay)
elif opts.type_MC == "run3" or opts.type_MC == "Run3":
    print("Will take cross section from Run3 file")
    VBS_txt='VBS_cross_section_Eboli_run3.txt'
    xsection_fb = uf.take_xsec_fb_aqgc(VBS_txt,EFT_op,EFT_type, proc, decay)
else:
    VBS_txt='VBS_xsection_test.txt'
    xsection_fb = uf.take_xsec_fb(VBS_txt,EFT_op,EFT_type, proc, decay)


if xsection_fb is None:
    print(f'No cross section in txt found for {keyy}')
    xsection_fb = None
    #xsection_fb=uf.cross_section_fb(EFT_op,EFT_type, proc, decay)
    if xsection_fb is None:
        print(f'No cross section found in AMI for {keyy}')
        print("Cross section set at 0 to avoid errors but ISSUE")
        xsection_fb = 0
        
print(f'Cross section in fb for {keyy}: {xsection_fb}')

mydir=conf_dir + f"/DOCUT_{opts.DOCUT}/"

plots_backup_dir = mydir + "Plots_backup/"
os.makedirs(plots_backup_dir, exist_ok=True)

save_job_infos(f"DOCUT={opts.DOCUT}", conf_dir + f"/DOCUT_{opts.DOCUT}/", prod_dec,xsec_fb=xsection_fb)
label_plot= f"{proc}_{decay}_{EFT_op}"

dir_plot_ = mydir + "/plots/"
dir_plot_2 = mydir + "/plots/ntuple/"

# Check if the directory exists
if os.path.exists(dir_plot_):
    # Delete the directory
    shutil.rmtree(dir_plot_)
    


if (opts.keep):
    timestamp_day = datetime.now().strftime("%d_%m")
    day_dir = plots_backup_dir+timestamp_day+ "/"
    os.makedirs(plots_backup_dir+timestamp_day+ "/", exist_ok=True)
    timestamp_hour = datetime.now().strftime("_%Hh%M")
    run_dir = day_dir + opts.name + timestamp_hour + "/"
    os.makedirs(run_dir, exist_ok=True)

    files_to_copy = [
        "ntuple_rivet.root",
        "/Info/frac_cuts_merged.txt",
        "/Info/frac_after_cuts_error_bar_merged.txt",
        "/Info/frac_cuts_resolved.txt",
        "/Info/frac_after_cuts_error_bar_resolved.txt",
        "cutflow_merged.txt",
        "cutflow_resolved.txt",
        "/Info/Cross_section_fb.txt",
    ]

    for file in files_to_copy:
        shutil.copy2(mydir  + file, run_dir)

    # Copy the repository mydir + "/plots/"
    #shutil.copytree(mydir + "/plots/", run_dir + "plots/")
print(f'Cross section in fb: {xsection_fb}')

