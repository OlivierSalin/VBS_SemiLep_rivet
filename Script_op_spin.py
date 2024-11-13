import subprocess
from multiprocessing import Pool
import itertools

all_ops = ["SM","FM0", "FM1", "FM2", "FM3", "FM4", "FM5", "FM7", "FS02", "FS1", "FT0", "FT1", "FT2", "FT5", "FT6", "FT7"]
all_ops_cat=["SM","FM0","FM2","FS1","FT1","FT5"]
processes = ["WpZ","WmZ"]
decays = ["llqq"]
name_run = "VBS_mismatch"



path_build_command = "-I/exp/atlas/salin/ATLAS/VBS_mc/vcpkg/installed/x64-linux/include/ `root-config --cflags --libs`"
for process in processes:
    for decay in decays:
        build_command = ["rivet-build", f"Rivet{process}_{decay}.so", f"{process}_{decay}.cc", f"{path_build_command}"]
        subprocess.run(build_command)
        
def run_command(proc_op_tuple):
    process, op = proc_op_tuple
    for decay in decays:
        if op == "SM":
            command = ["python", "run_rivet.py", "--evtMax", "5000", "--conf", f"user.osalin.MadGraph_{process}_{decay}_FM0_SM", "--DOCUT", "YES", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", f"{name_run}"]
        else:
            command = ["python", "run_rivet.py", "--evtMax", "5000", "--conf", f"user.osalin.MadGraph_{process}_{decay}_{op}_QUAD", "--DOCUT", "YES", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", f"{name_run}"]
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
        
        
