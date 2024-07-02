# combine panda downloading, running rivet, making plots for many combinations
from pandaclient import panda_api
c = panda_api.get_api()
import subprocess
import ROOT
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(ROOT.kTRUE)
import os
import glob
from array import array
import matplotlib
import yoda
import lib_utils as lu
from optparse import OptionParser
import math

def prepare_grid_files(i_job_name):
    print("will download+untar files evnt and log for", i_job_name)
    evnt_did, evnt_dir, log_did, log_dir = lu.get_envt_log_names_dirs(base_dir,i_job_name)
    if not os.path.exists(evnt_dir):
        print ("will download", evnt_did, "since dir doesn't exist", evnt_dir)
        subprocess.call(f"rucio download {evnt_did}", shell=True, cwd=base_dir)
    else:
        print("have ", evnt_dir ,"already so dont download")
    if not os.path.exists(log_dir):
        print ("will download", log_did, "since dir doesn't exist", log_dir)
        subprocess.call(f"rucio download {log_did}", shell=True, cwd=evnt_dir)
    else:
        print("have ", evnt_dir ,"already so dont download")
    # if job were reran tehre will be several logs and last one comes from succesfull
    untared_dir_cand = glob.glob(f"{evnt_dir}/*.log/tarball_PandaJob*")
    if len(untared_dir_cand)==0:
        tar_file_candidates = sorted(glob.glob(f"{log_dir}/*log.tgz"))
        tar_file = os.path.basename(tar_file_candidates[-1])
        print("have candidates for tar files", tar_file_candidates)
        print("from them select the file to untar", tar_file)
        subprocess.call(f"tar -xvf {tar_file}", shell=True, cwd=log_dir)
    else:
        print("dont untar since did it before res in ", untared_dir_cand[0])

def save_job_infos(DOCUT_str, mydir, xsec_fb, prod_dec):
    yoda_f_str = mydir + "MyOutput.yoda.gz"

    if not os.path.exists(yoda_f_str): 
        print("dont see yoda file in dir ", mydir, ",return")
        return
    
    yoda_f = yoda.read(yoda_f_str)
    print("reading from yoda file ", yoda_f_str)
    all_hists_in_yoda = [iname  for iname in yoda_f.keys() if "[" not in iname and "RAW" not in iname]
    hists_1h_in_yoda = []
    for i_name in all_hists_in_yoda:
        if yoda_f[i_name].type()=="Histo1D": hists_1h_in_yoda.append(i_name)  
    print("have 1d hists to be saved in root:", hists_1h_in_yoda, "in yoda file", yoda_f_str)

    root_file = mydir + "/hists.root"
    if not os.path.exists(root_file) or not os.path.exists(mydir + "xsec_fb.txt"): proceed = 1
    elif (os.path.exists(root_file) or os.path.exists(mydir + "xsec_fb.txt")) and opts.runAgain=="yes": proceed = 1
    else: proceed = 0

    if not proceed:
        print("dont do xsec and hisst to root")
        return  

    # save fid xsec
    rivet_dir_name = f"/{prod_dec}:OUTDIR=/{mydir}".replace("//","/")
    print("looking for prefix in counter",rivet_dir_name)
    pos_n_in = yoda_f[f"{rivet_dir_name}/pos_w_initial"].numEntries()
    neg_n_in = yoda_f[f"{rivet_dir_name}/neg_w_initial"].numEntries()
    pos_n_f = yoda_f[f"{rivet_dir_name}/pos_w_final"].numEntries()
    neg_n_f = yoda_f[f"{rivet_dir_name}/neg_w_final"].numEntries()
    #
    frac_cut = (pos_n_f+neg_n_f) / (pos_n_in+neg_n_in) 
    frac_pos = pos_n_f / pos_n_in 
    frac_neg = neg_n_f / neg_n_in
    frac_cut_er_bar = 1/(pos_n_in+neg_n_in) * math.sqrt(pos_n_in*frac_pos*(1-frac_pos) + neg_n_in*frac_neg*(1-frac_neg))
    # frac_cut_er_bar = frac_cut_unc / 2 
    #
    pos_w_in = yoda_f[f"{rivet_dir_name}/pos_w_initial"].sumW()
    neg_w_in = yoda_f[f"{rivet_dir_name}/neg_w_initial"].sumW()
    pos_w_f = yoda_f[f"{rivet_dir_name}/pos_w_final"].sumW()
    neg_w_f = yoda_f[f"{rivet_dir_name}/neg_w_final"].sumW()
    #
    lu.save_xsec_frac_prod(mydir,xsec_fb,
                            frac_cut, frac_pos, frac_neg, frac_cut_er_bar,
                            pos_w_in, neg_w_in, pos_w_f, neg_w_f,
                            pos_n_in, neg_n_in, pos_n_f, neg_n_f)
    # save hists in root for further plotting
    root_file = ROOT.TFile(root_file,"UPDATE")
    for i_hist in hists_1h_in_yoda: # they are in format '/WpWm_lvlv:DOCUT=YES/leptons_pids'
        h_yoda =  yoda_f[i_hist]
        h_root = lu.yoda_to_root_1d(h_yoda, i_hist.split("/")[-1])
        h_root.Write("", ROOT.TObject.kOverwrite)
    root_file.Close()

    ############## draw event and cutflow
    if opts.runWithCuts=="yes":
        print("drawing average and distubtioni image of event")
        lu.draw_average_event(mydir, average_im=True)
        lu.draw_average_event(mydir, average_im=False)
        #
        print("saving cutflow as img")
        cutflow_file = mydir + "cutflow.txt"
        if os.path.exists(cutflow_file):
            cut_names, cut_cumu, cut_incr = lu.get_cutflow_arrays(cutflow_file)
            lu.draw_cutflows(cut_names, [cut_incr,cut_cumu], ["incremental","cumulative"],
                            mydir+"/cutflow_img.png", prod_dec)
            
def main():
    parser = OptionParser()
    parser.add_option("--runWithCuts", default = "yes")
    parser.add_option("--runAgain", default = "no")
    parser.add_option("--jobName", default = "")
    parser.add_option("--evtMax", default = 20000)
    global opts
    opts, _ = parser.parse_args()

    global base_dir
    # global evtMax
    evtMax = int(opts.evtMax)

    if len(opts.jobName)==0: return

    print("##################### \n ############# will work on job", opts.jobName)
    prod_dec, base_dir = lu.find_prod_dec_and_dir(opts.jobName) # dir where all files are stored
    # save_hists_log_get_xsec_after_cuts(opts.jobName)
    job_name = opts.jobName
    prepare_grid_files(job_name)
    evnt_file, log_file = lu.get_evnt_log_files(base_dir,job_name)

    if evnt_file==-1 or log_file==-1: return 
    
    if opts.runWithCuts=="yes":
        com_run_rivet = lu.get_rivet_com(job_name, evtMax = evtMax, DOCUT = "YES", redoRivet=opts.runAgain, redoPlots=opts.runAgain)
        print("run rivet+untar in with com \n", com_run_rivet)
        subprocess.call(com_run_rivet, shell=True)
    else:
        com_run_rivet = lu.get_rivet_com(job_name, evtMax = evtMax, DOCUT = "NO", redoRivet=opts.runAgain, redoPlots=opts.runAgain)
        print("run rivet+untar in with com \n", com_run_rivet)
        subprocess.call(com_run_rivet, shell=True)
    # get xsec after cuts
    evnt_dir = os.path.dirname(evnt_file)
    xsec_fb = lu.get_xsec(log_file)
    # prod_dec, _ = lu.find_prod_dec_and_dir(job_name)
    if opts.runWithCuts=="yes": 
        save_job_infos("DOCUT=YES",evnt_dir + "/DOCUT_YES/", xsec_fb, prod_dec)
    else: 
        save_job_infos("DOCUT=NO",evnt_dir + "/DOCUT_NO/", xsec_fb, prod_dec)

    return 0

if __name__ == "__main__":
    main()