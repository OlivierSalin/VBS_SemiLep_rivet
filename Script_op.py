import subprocess
from multiprocessing import Pool

all_ops = ["SM","FM0", "FM1", "FM2", "FM3", "FM4", "FM5", "FM7", "FS02", "FS1", "FT0", "FT1", "FT2", "FT5", "FT6", "FT7"]
#all_ops=["FM0","FM2"]
processes = ["WmZ", "WpZ","ZZ"]

def run_command(op):
    for process in processes:
        if op == "SM":
            command = ["python", "run_rivet.py", "--evtMax", "10000", "--conf", f"MGPy8EG_aQGCFM0_SM_1_{process}_llqq", "--DOCUT", "YES", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", "All_op_plot"]
        else:
            command = ["python", "run_rivet.py", "--evtMax", "10000", "--conf", f"MGPy8EG_aQGC{op}_QUAD_1_{process}_llqq", "--DOCUT", "YES", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", "All_op_plot"]
        subprocess.run(command)

""" all_ops_SM=["FM0"]
processes_SM = ["WmZ", "WpZ","ZZ"]        
def run_command_SM(process):
    command = ["python", "run_rivet.py", "--evtMax", "10000", "--conf", f"MGPy8EG_aQGCFM0_SM_1_{process}_llqq", "--DOCUT", "YES", "--redoRivet", "yes", "--redoPlots", "no", "--keep", "True", "--name", "All_op_plot"]
    subprocess.run(command) """

if __name__ == "__main__":
    with Pool() as p:
        p.map(run_command, all_ops)
        #p.map(run_command_SM, processes_SM)
        
        