import subprocess
from multiprocessing import Pool
import itertools
import glob
import os
import shutil
import time

from optparse import OptionParser
parser = OptionParser()

parser.add_option("--nb_part", default = 0)
parser.add_option("--nb_lep", default = 2)
parser.add_option("--nb_events", default = 1000)
parser.add_option("--All_channel", default= True)
parser.add_option("--aQGC", default= True)
parser.add_option("--Name", default= "All_ops")
parser.add_option("--name_copy", default= "aqgc_prod")
parser.add_option("--EFT_order", default = "QUAD")
parser.add_option('--type_MC', default = "")
parser.add_option("--Channel", default = "")
parser.add_option("--decays", default = "llqq")
parser.add_option("--Conf", default = "user.osalin")
opts, _ = parser.parse_args()

nb_events = int(opts.nb_events)
nb_lepton= int(opts.nb_lep)
#nb_lepton=1

order=opts.EFT_order


all_ops_cat=["FM0","FM"]



Conf_list = {}






all_ops_cat = ["SM","FM0","FM2","FS1","FT1","FT5"]
all_ops_cat = ["FM0","FS0","FT0"]
#all_ops_cat = ["FM0"]
all_ops_cat = ["FM0","FM1","FM2","FM3","FM4","FM5","FM7","FM8","FM9",
            "FS0","FS1","FS2",
            "FT0","FT1","FT2","FT3","FT4","FT5","FT6","FT7"]

#all_ops_cat = ["FM0","FM1"]
all_ops_cat = ["FM","FS","FT"]

#all_ops_cat = ["FS0","FS1","FS2"]
#all_ops_cat = ["FM0","FM1"]
cross_terms_ = [f"{op1}vs{op2}" for op1, op2 in itertools.combinations(all_ops_cat, 2)]
print("Defined cross terms:", cross_terms_)

if opts.EFT_order == "CROSS":
    all_ops_cat = cross_terms_


processes = ["WpZ"]
decays = ['llqq']
proc_decays_tuple = list(itertools.product(processes, decays))
#decays=[opts.decays]
name_run = opts.Name


def find_prod_dec_and_dir_tres(conf, type_MC=None):
    base_path = "/exp/atlas/salin/ATLAS/VBS_mc/eft_files"
    #print("New models")
    
    def extract_prod_dec(conf):
        prod_temp = conf[conf.find("user.osalin.MadGraph_") + len("user.osalin.MadGraph_"):]
        #prod_temp = conf[conf.find("user.osalin.Madgraph_") + len("user.osalin.MadGraph_"):]
        print("start from string", prod_temp)
        prod_dec = prod_temp[:prod_temp.find("qq_") + 2]
        #print("from conf found production dec", prod_dec)
        return prod_dec

    if conf.startswith("user."):
        prod_dec = extract_prod_dec(conf)
        
        if "Run3" in type_MC or "run3" in type_MC:
            base_dir = f"{base_path}/Run3/{prod_dec}/"
        elif "aqgc" in type_MC or "model" in type_MC:
            base_dir = f"{base_path}/aqgc_model/{prod_dec}/"
        #elif "Reweighting" in type_MC or "reweight" in type_MC or "rwg" in type_MC:
            #base_dir = f"{base_path}/Reweighting/{prod_dec}/"
        elif "13p0" in type_MC or "13TeV" in type_MC:
            base_dir = f"{base_path}/13p0/{prod_dec}/"
        elif "MCprod_QGC_R3" in type_MC or "MCprod_aqgc_R3" in type_MC:
            base_dir = f"{base_path}/MCprod/aqgc/13p6/{prod_dec}/"
        elif "MCprod_QGC_R2" in type_MC or "MCprod_aqgc_R2" in type_MC:
            base_dir = f"{base_path}/MCprod/aqgc/13p0/{prod_dec}/"
            
        elif "Reweighting_Madspin" in type_MC:
            base_dir = f"{base_path}/Reweighting/Madspin/{prod_dec}/"
        elif "Reweighting_NoSpin" in type_MC:
            base_dir = f"{base_path}/Reweighting/NoSpin/{prod_dec}/"
        elif "Reweighting_Decay_chain" in type_MC:
            base_dir = f"{base_path}/Reweighting/Decay_chain/{prod_dec}/"
        elif "Reweighthel_ignore" in type_MC or "Reweighting_hel_ignore" in type_MC:
            base_dir = f"{base_path}/Reweighting/Madspin/hel_ignore/{prod_dec}/"
        elif "Reweighthel_aware" in type_MC or "Reweighting_hel_aware" in type_MC:
            base_dir = f"{base_path}/Reweighting/Madspin/hel_aware/{prod_dec}/"
        elif "Reweight_Polarisation" in type_MC or "Reweighting_Polarisation" in type_MC:
            base_dir = f"{base_path}/Reweighting/Polarisation/{prod_dec}/"
            
        elif "EFTDec_Madspin" in type_MC:
            base_dir = f"{base_path}/EFTDec/Madspin//{prod_dec}/"
        elif "EFTDec_NoSpin" in type_MC:
            base_dir = f"{base_path}/EFTDec/NoSpin//{prod_dec}/"
        elif "EFTDec_Decay_chain" in type_MC:
            base_dir = f"{base_path}/EFTDec/Decay_chain//{prod_dec}/"

        else:
            base_dir = f"{base_path}/{prod_dec}/"
        
        print("base_dir: ", base_dir)
        search_com= base_dir + f"/*{conf}*EXT0"
        print("searching for dir with pattern", search_com)
        conf_dir_arr = glob.glob(search_com)
        print("found possibilities for dir", conf_dir_arr)
        conf_dir = conf_dir_arr[0] if len(conf_dir_arr)>=1 else -1  
        if conf_dir == -1: raise ValueError("did not find folder for this config ",search_com)
    
    return prod_dec, conf_dir


path_build_command = "-I/exp/atlas/salin/ATLAS/VBS_mc/vcpkg/installed/x64-linux/include/ `root-config --cflags --libs`"
for process in processes:
    for decay in decays:
        build_command = ["rivet-build", f"Rivet{process}_{decay}.so", f"{process}_{decay}.cc", "EventWeights.cc", f"{path_build_command}"]
        print(f"Build command: {build_command}")
        subprocess.run(build_command)
        
def run_command(proc_decay_op_tuple):
    process, decay, op = proc_decay_op_tuple
    if op == "SM":
        command = ["python", "run_rivet.py", "--evtMax", f"{nb_events}", "--conf", f"user.osalin.MadGraph_{process}_{decay}_FM0_SM", "--DOCUT", "YES","--type_MC",f"{opts.type_MC}", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", f"{name_run}"]
    else:
        if "Run2" in opts.type_MC or "run2" in opts.type_MC:
            command = ["python", "run_rivet.py", "--evtMax", f"{nb_events}", "--conf", f"MGPy8EG_aQGC{op}_{order}_1_{process}_{decay}", "--DOCUT", "YES","--type_MC",f"{opts.type_MC}", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", f"{name_run}"]
        else:
            command = ["python", "run_rivet.py", "--evtMax", f"{nb_events}", "--conf", f"user.osalin.MadGraph_{process}_{decay}_{op}_{order}", "--DOCUT", "YES","--type_MC",f"{opts.type_MC}", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", f"{name_run}"]
    print(f"Running command: {command}")
    subprocess.run(command)

combined_processes_ops = list(itertools.product(processes, decays))
print(combined_processes_ops)

if __name__ == "__main__":
    # Prepare all (process, decay, op) combinations
    proc_decay_op_tuples = [
        (process, decay, op)
        for (process, decay) in proc_decays_tuple
        for op in all_ops_cat
    ]
    with Pool() as p:
        print("test")
        p.map(run_command, proc_decay_op_tuples)

    for process,decay in proc_decays_tuple:
        for op in all_ops_cat:

            Conf= f"user.osalin.MadGraph_{process}_{decay}_{op}_{order}"
            print("Conf: ", Conf)
            try:
                proc, Base_dir = find_prod_dec_and_dir_tres(Conf, opts.type_MC)
            except ValueError as e:
                print(f"Error: {e}")
                continue
            print("Base_dir: ", Base_dir + f"/DOCUT_YES/")
            Dir_ntuple= Base_dir + f"/DOCUT_YES/"
            path_hist_base =f"/exp/atlas/salin/ATLAS/VBS_mc/eft_files/Histograms/Reweighting/Rwg_test/{opts.type_MC}/"
            path_hist= path_hist_base + f"/{process}_{decay}/{op}_{order}//"
            if not os.path.exists(path_hist):
                os.makedirs(path_hist)
            else:
                print("Warning: directory already exists, skipping.")
                timestamp = time.strftime("%m%d-%H%M")
                path_hist = path_hist + f"/Already_exist/{opts.name_copy}_{timestamp}//{process}_{decay}/{op}_{order}/{opts.type_MC}/"
                os.makedirs(path_hist)
            
            # Copy the content of Dir_ntuple into path_hist
            print(f"Copying {Dir_ntuple} to {path_hist}")
            shutil.copytree(Dir_ntuple, path_hist, dirs_exist_ok=True)














