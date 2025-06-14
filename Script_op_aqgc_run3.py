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


if opts.All_channel:
    processes = ["WmZ","WpZ","ZZ","WmWm","WpWm","WpWp"]
    decays = ['lvqq','llqq','vvqq']
else:
    if nb_lepton == 1:
        print("1 lepton")
        processes = ["WmZ","WpZ","WmWm","WpWm","WpWp"]
        decays = ['lvqq']
    elif nb_lepton == 2:
        print("2 lepton")
        processes=["WpZ","WmZ","ZZ"]
        decays = ['llqq']
    elif nb_lepton == 0:
        print("0 lepton")
        processes=["WpZ","WmZ","ZZ"]
        decays = ['vvqq']

proc_decays_tuple = list(itertools.product(processes, decays))

# Filter out invalid combinations
valid_combinations = {
    1: [("WmZ", "lvqq"), ("WpZ", "lvqq"), ("WmWm", "lvqq"), ("WpWm", "lvqq"), ("WpWp", "lvqq")],
    2: [("WpZ", "llqq"), ("WmZ", "llqq"), ("ZZ", "llqq")],
    0: [("WpZ", "vvqq"), ("WmZ", "vvqq"), ("ZZ", "vvqq")]
}

if opts.All_channel:
    proc_decays_tuple = [comb for comb in proc_decays_tuple if comb in valid_combinations[1] + valid_combinations[2] + valid_combinations[0]]
else:
    proc_decays_tuple = [comb for comb in proc_decays_tuple if comb in valid_combinations[nb_lepton]]



all_ops_cat = ["SM","FM0","FM2","FS1","FT1","FT5"]
all_ops_cat = ["FM0","FS0","FT0"]
#all_ops_cat = ["FM0"]


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
        elif "Reweighting" in type_MC or "reweight" in type_MC or "rwg" in type_MC:
            base_dir = f"{base_path}/Reweighting/{prod_dec}/"
        elif "13p0" in type_MC or "13TeV" in type_MC:
            base_dir = f"{base_path}/13p0/{prod_dec}/"
        elif "MCprod_QGC_R3" in type_MC or "MCprod_aqgc_R3" in type_MC:
            base_dir = f"{base_path}/MCprod/aqgc/13p6/{prod_dec}/"
        elif "MCprod_QGC_R2" in type_MC or "MCprod_aqgc_R2" in type_MC:
            base_dir = f"{base_path}/MCprod/aqgc/13p0/{prod_dec}/"
        else:
            base_dir = f"{base_path}/{prod_dec}/"
        
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
        build_command = ["rivet-build", f"Rivet{process}_{decay}.so", f"{process}_{decay}.cc", f"{path_build_command}"]
        subprocess.run(build_command)
        
def run_command(proc_gen_tuple):
    #process, Gen= proc_gen_tuple
    process, decay = proc_gen_tuple
    for op in all_ops_cat:
        if op == "SM":
            command = ["python", "run_rivet.py", "--evtMax", f"{nb_events}", "--conf", f"user.osalin.MadGraph_{process}_{decay}_FM0_SM", "--DOCUT", "YES","--type_MC",f"{opts.type_MC}", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", f"{name_run}"]
        else:
            if "Run2" in opts.type_MC or "run2" in opts.type_MC:
                command = ["python", "run_rivet.py", "--evtMax", f"{nb_events}", "--conf", f"MGPy8EG_aQGC{op}_{order}_1_{process}_{decay}", "--DOCUT", "YES","--type_MC",f"{opts.type_MC}", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", f"{name_run}"]
            else:
                #command = ["python", "run_rivet.py", "--evtMax", f"{nb_events}", "--conf", f"user.osalin.MadGraph_{process}_{decay}_{op}_{order}", "--DOCUT", "YES","--type_MC",f"{opts.type_MC}", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", f"{name_run}"]
                command = ["python", "run_rivet.py", "--evtMax", f"{nb_events}", "--conf", f"user.osalin.MadGraph_{process}_{decay}_{op}_{order}", "--DOCUT", "YES","--type_MC",f"{opts.type_MC}", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", f"{name_run}"]
            
        print(f"Running command: {command}")

        subprocess.run(command)

def run_command_run2(process_decay):

    process, decay = process_decay
            #command = ["python", "run_rivet.py", "--evtMax", f"{nb_events}", "--conf", f"user.osalin.MadGraph_{process}_{decay}_{op}_{order}", "--DOCUT", "YES","--type_MC",f"{opts.type_MC}", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", f"{name_run}"]
    command = ["python", "run_rivet.py", "--evtMax", f"{nb_events}", "--conf", f"user.osalin.MadGraph_{process}_{decay}_FM0_{order}", "--DOCUT", "YES","--type_MC",f"", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", f"{name_run}"]
            
    print(f"Running command Run2 : {command}")

    subprocess.run(command)


combined_processes_ops = list(itertools.product(processes, decays))
print(combined_processes_ops)

if __name__ == "__main__":
    with Pool() as p:
        #proc_gen_tuple = list(itertools.product(processes, Generator_))
        #proc_decays_tuple = list(itertools.product(processes, decays))
        
        p.map(run_command, proc_decays_tuple)
        #p.map(run_command_run2, proc_decays_tuple)

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
            path_hist_base ="/exp/atlas/salin/ATLAS/VBS_mc/eft_files/Histograms/aQGC_Run3/MC_prod/"
            path_hist= path_hist_base + f"/Run3/{process}_{decay}/{op}_{order}/{opts.name_copy}/"
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














