import subprocess
from multiprocessing import Pool
import itertools

from optparse import OptionParser
parser = OptionParser()

parser.add_option("--nb_part", default = 0)
parser.add_option("--nb_events", default = 50000)
parser.add_option("--aQGC", default= True)
parser.add_option("--Name", default= "All_ops")
parser.add_option("--EFT_order", default = "QUAD")
parser.add_option('--type_MC', default = "")
parser.add_option("--Channel", default = "")
opts, _ = parser.parse_args()

nb_events = int(opts.nb_events)

order=opts.EFT_order

all_ops = ["SM","FM0", "FM1", "FM2", "FM3", "FM4", "FM5", "FM7", "FS2", "FS1", "FT0", "FT1", "FT2", "FT5", "FT6", "FT7"]
all_ops_cat=["FM0","FM1","FM2","FS1",'FT0',"FT1","FT5"]

if "aqgc" in opts.type_MC or "model" in opts.type_MC:
    all_ops_cat = ["FM0","FM1","FM2",'FM3','FM4','FM5','FM7',"FM8","FM9",
            "FS0","FS1","FS2",
            "FT0","FT1",'FT2','FT3','FT4',"FT5",'FT6','FT7'
            "FM1odd","FM2odd","FT2odd","FT6odd"]
    
    all_ops_cat = ["FM0","FM1","FM2",'FM3','FM4','FM5','FM7',"FM8","FM9",
            "FS0","FS1","FS2",
            "FT0","FT1",'FT2','FT3','FT4',"FT5",'FT6','FT7']
    #all_ops_cat=["FM0","FM1","FM2","FS1",'FT0',"FT1","FT5","FM8","FM9","FM2odd","FT6odd"]


processes = ["WpZ"]
decays = ["llqq"]
name_run = opts.Name


path_build_command = "-I/exp/atlas/salin/ATLAS/VBS_mc/vcpkg/installed/x64-linux/include/ `root-config --cflags --libs`"
for process in processes:
    for decay in decays:
        build_command = ["rivet-build", f"Rivet{process}_{decay}.so", f"{process}_{decay}.cc", f"{path_build_command}"]
        subprocess.run(build_command)
        
def run_command(proc_op_tuple):
    process, op = proc_op_tuple
    for decay in decays:
        if op == "SM":
            command = ["python", "run_rivet.py", "--evtMax", f"{nb_events}", "--conf", f"user.osalin.MadGraph_{process}_{decay}_FM0_SM", "--DOCUT", "YES","--type_MC",f"{opts.type_MC}", "--redoRivet", "yes", "--redoPlots", "yes", "--keep", "True", "--name", f"{name_run}"]
        else:
            command = ["python", "run_rivet.py", "--evtMax", f"{nb_events}", "--conf", f"user.osalin.MadGraph_{process}_{decay}_{op}_{order}", "--DOCUT", "YES","--type_MC",f"{opts.type_MC}", "--redoRivet", "yes", "--redoPlots", "yes", "--keep", "True", "--name", f"{name_run}"]
        subprocess.run(command)

""" all_ops_SM=["FM0"]
processes_SM = ["WmZ", "WpZ","ZZ"]        
def run_command_SM(process):
    command = ["python", "run_rivet.py", "--evtMax", "10000", "--conf", f"MGPy8EG_aQGCFM0_SM_1_{process}_llqq", "--DOCUT", "YES", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", "All_op_plot"]
    subprocess.run(command) """

combined_processes_ops = list(itertools.product(processes, all_ops_cat))

if __name__ == "__main__":
    with Pool() as p:
        process_op_tuples = list(itertools.product(processes, all_ops_cat))
        p.map(run_command, process_op_tuples)
        #p.map(run_command_SM, processes_SM)
        
        
