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
parser.add_option("--nb_events", default = 1000)
parser.add_option("--aQGC", default= True)
parser.add_option("--Name", default= "All_ops")
parser.add_option("--name_copy", default= "Py_H7")
parser.add_option("--EFT_order", default = "SM")
parser.add_option('--type_MC', default = "")
parser.add_option("--Channel", default = "")
parser.add_option("--decays", default = "llqq")
opts, _ = parser.parse_args()

nb_events = int(opts.nb_events)

order=opts.EFT_order


all_ops_cat=["FM0","FM"]
Generator_=["Py8","H7"]


Conf_list = {}
processes = ["WZjj","ZZ"]
processes = ["WZjj"]
decays = ["llqq","vvqq","lvqq"]
decays=[opts.decays]
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
            conf_dir = f"{base_path}/13p0/{prod_dec}/"
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
    process, Gen= proc_gen_tuple
    for decay in decays:
        if Gen == "SM":
            command = ["python", "run_rivet.py", "--evtMax", f"{nb_events}", "--conf", f"user.osalin.MadGraph_{process}_{decay}_FM0_SM", "--DOCUT", "YES","--type_MC",f"{opts.type_MC}", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", f"{name_run}"]
        else:
            #command = ["python", "run_rivet.py", "--evtMax", f"{nb_events}", "--conf", f"user.osalin.MadGraph_{process}_{decay}_{op}_{order}", "--DOCUT", "YES","--type_MC",f"{opts.type_MC}", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", f"{name_run}"]
            command = ["python", "run_rivet.py", "--evtMax", f"{nb_events}", "--conf", f"user.osalin.MadGraph_{process}_{decay}_{order}_{Gen}", "--DOCUT", "YES","--type_MC",f"{opts.type_MC}", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", f"{name_run}"]
            
        print(f"Running command: {command}")

        subprocess.run(command)

def run_command_run2(process_decay):

    process, decay = process_decay
            #command = ["python", "run_rivet.py", "--evtMax", f"{nb_events}", "--conf", f"user.osalin.MadGraph_{process}_{decay}_{op}_{order}", "--DOCUT", "YES","--type_MC",f"{opts.type_MC}", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", f"{name_run}"]
    command = ["python", "run_rivet.py", "--evtMax", f"{nb_events}", "--conf", f"user.osalin.MadGraph_{process}_{decay}_FM0_{order}", "--DOCUT", "YES","--type_MC",f"", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", f"{name_run}"]
            
    print(f"Running command Run2 : {command}")

    subprocess.run(command)


combined_processes_ops = list(itertools.product(processes, Generator_))
print(combined_processes_ops)

if __name__ == "__main__":
    with Pool() as p:
        proc_gen_tuple = list(itertools.product(processes, Generator_))
        proc_decays_tuple = list(itertools.product(processes, decays))
        p.map(run_command, proc_gen_tuple)
        #p.map(run_command_run2, proc_decays_tuple)

    for proc_gen in proc_gen_tuple:
        for dec in decays:
            process, Gen = proc_gen
            Conf= f"user.osalin.MadGraph_{process}_{dec}_{order}_{Gen}"
            print("Conf: ", Conf)
            proc,Base_dir=find_prod_dec_and_dir_tres(Conf, opts.type_MC)
            print("Base_dir: ", Base_dir + f"/DOCUT_YES/")
            print("Proc: ", proc)
            Dir_ntuple= Base_dir + f"/DOCUT_YES/"
            path_hist_base ="/exp/atlas/salin/ATLAS/VBS_mc/eft_files/Histograms/SM_Run3/"
            path_hist= path_hist_base + f"/{opts.name_copy}/{process}_{dec}/{Gen}/"
            if not os.path.exists(path_hist):
                os.makedirs(path_hist)
            else:
                print("Warning: directory already exists, skipping.")
                timestamp = time.strftime("%m%d-%H%M")
                path_hist = path_hist_base + f"/Already_exist/{opts.name_copy}_{timestamp}/{process}_{dec}/{Gen}/"
                os.makedirs(path_hist)
            
            # Copy the content of Dir_ntuple into path_hist
            shutil.copytree(Dir_ntuple, path_hist, dirs_exist_ok=True)














