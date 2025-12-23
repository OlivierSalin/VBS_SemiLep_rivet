#!/usr/bin/env python3
import subprocess
import multiprocessing
import sys
from optparse import OptionParser
import os
import shutil
import glob	
import time
parser = OptionParser()
parser.add_option("--nb_lep", default = 2)
parser.add_option("--nb_events", default = 100)
parser.add_option("--EFT_order", default = "QUAD")
parser.add_option("--name", default = "")
parser.add_option("--Mode", default = "Lep")

parser.add_option("--type_mc", default = "")

(opts, args) = parser.parse_args()
nb_lepton= int(opts.nb_lep)
nb_events = int(opts.nb_events)
name_outdir= opts.name

order_EFT= opts.EFT_order
type_mc_rivet = opts.type_mc
# --- Configuration: edit these lists/values directly ---
valid_combinations = {
	4: [("ZZ", "llll")],
	3: [("WmZ", "lllv"), ("WpZ", "lllv")],
	2: [("WmWm", "lvlv"), ("WpWm", "lvlv"), ("WpWp", "lvlv")],
	1: [("Wpy", "lvy")],
}

if opts.Mode == "Lep":
	print("Running in Leptonic mode")
	valid_combinations = {
		4: [("ZZ", "llll")],
		3: [("WmZ", "lllv"), ("WpZ", "lllv")],
		2: [("WmWm", "lvlv"), ("WpWm", "lvlv"), ("WpWp", "lvlv"),("Zy", "lly")],
		1: [("Wpy", "lvy")],
	}


elif opts.Mode == "SemiLep":
	print("Running in Semi-Leptonic mode")
	valid_combinations = {
		1: [("WmZ", "lvqq"), ("WpZ", "lvqq"), ("WmWm", "lvqq"), ("WpWm", "lvqq"), ("WpWp", "lvqq")],
		2: [("WpZ", "llqq"), ("WmZ", "llqq"), ("ZZ", "llqq")],
		0: [("WpZ", "vvqq"), ("WmZ", "vvqq"), ("ZZ", "vvqq")]
	}

valid_combi_aQGC = [f"{proc}_{dec}" for pairs in valid_combinations.values() for proc, dec in pairs]
valid_combi_aQGC_lep = [f"{proc}_{dec}" for proc, dec in valid_combinations[nb_lepton]]
print("Selected processes for nb_lep =", nb_lepton, ":\n", valid_combi_aQGC_lep)

base_path_outdir = f"/exp/atlas/salin/ATLAS/VBS_mc/eft_files/Histograms/Polarisation/Leptonic/Rwg/"
outdir = f"{base_path_outdir}//GOOD_output//{name_outdir}/{opts.Mode}/"

PROCESSES = valid_combi_aQGC_lep         # list of processes (example: ["ZZ_llll"])
#PROCESSES = ["WpWp_lvlv"]          # list of processes (example: ["ZZ_llll"])
OPERATORS = ["FM","FS","FT","FModd","FTodd"] 
POLARISATION= [None]
POLARISATION= ["LL","LT","TL","TT"]      
EVTMAX = nb_events
DOCUT = "YES"
REDO_RIVET = "yes"
REDO_PLOTS = "no"
TYPE_MC = type_mc_rivet
JOBS = 0                         # 0 => use number of jobs = number of operators (per process)
PYTHON_CMD = sys.executable      # use same python interpreter that runs this script
# --- end configuration ---

def find_prod_dec_and_dir_tres(conf, type_MC=None):
    base_path = "/data/atlas/salin/VBS_mc/VBS/eft_Files/"
    #print("New models")
	
    
    def extract_prod_dec(conf):
        prod_temp = conf[conf.find("user.osalin.MadGraph_") + len("user.osalin.MadGraph_"):]
        #print("start from string", prod_temp)
        # take the first two underscore-separated fields as production+decay, e.g. "ZZ_llll", "WpZ_lvlv", "WmZ_lllv"
        parts = prod_temp.split("_")
        if len(parts) >= 2:
            prod_dec = f"{parts[0]}_{parts[1]}"
        else:
            prod_dec = parts[0]
        #print("from conf found production dec", prod_dec)
        return prod_dec

    if conf.startswith("user."):
        prod_dec = extract_prod_dec(conf)
        
        if "Reweight_VBSLep" in type_MC:
            base_dir = f"{base_path}/Reweighting/VBSAll/Pol/Gen/VBS/Leptonic/{prod_dec}/"

        elif "Reweight_VBSSemiLep" in type_MC:
            base_dir = f"{base_path}/Reweighting/VBSAll/Pol/Gen/VBS/SemiLep/{prod_dec}/"

        elif "EFTDec_VBSLeptonic" in type_MC:
            base_dir = f"{base_path}/EFTDec/VBSAll/Pol/Leptonic/{prod_dec}/"   
            
        elif "EFTDec_Madspin" in type_MC:
            base_dir = f"{base_path}/EFTDec/Madspin//{prod_dec}/"

        elif "EFTDec_Polarisation" in type_MC:
            base_dir = f"{base_path}/EFTDec/Polarisation/{prod_dec}/"

        else:
            base_dir = f"{base_path}/test/{prod_dec}/"
        
        #print("base_dir: ", base_dir)
        search_com= base_dir + f"/*{conf}*EXT0"
        #print("searching for dir with pattern", search_com)
        conf_dir_arr = glob.glob(search_com)
        #print("found possibilities for dir", conf_dir_arr)
        conf_dir = conf_dir_arr[0] if len(conf_dir_arr)>=1 else -1  
        if conf_dir == -1: raise ValueError("did not find folder for this config ",search_com)
    
    return prod_dec, conf_dir

def run_command(cmd, shell=False):
	"""Run command and return (returncode, stdout, stderr)."""
	try:
		res = subprocess.run(cmd, check=True, capture_output=True, text=True, shell=shell)
		return 0, res.stdout, res.stderr
	except subprocess.CalledProcessError as e:
		return e.returncode, e.stdout or "", e.stderr or ""



def worker(process, op, pol):
	# ensure operator has _QUAD suffix
	op_suff = op if op.endswith(f"_{order_EFT}") else f"{op}_{order_EFT}"
	if pol is None:
		conf = f"user.osalin.MadGraph_{process}_{op_suff}"
	else:
		conf = f"user.osalin.MadGraph_{process}_{op_suff}_{pol}"
	cmd = [
		PYTHON_CMD,
		"run_rivet.py",
		"--evtMax", str(EVTMAX),
		"--conf", conf,
		"--DOCUT", DOCUT,
		"--redoRivet", REDO_RIVET,
		"--redoPlots", REDO_PLOTS,
		"--type_MC", TYPE_MC,
	]
	code, out, err = run_command(cmd)
	return (process, op, pol, conf, code, out, err)

def main():
	import glob
	# 1) Run rivet-build for each process (keeps backticks for `root-config`)
	for proc in PROCESSES:
		rivet_build_cmd = (
			f"rivet-build Rivet{proc}.so {proc}.cc EventWeights.cc "
			"-I/exp/atlas/salin/ATLAS/VBS_mc/vcpkg/installed/x64-linux/include/ "
			"`root-config --cflags --libs`"
		)
		print(f"Running rivet-build for process '{proc}'...")
		ret, out, err = run_command(rivet_build_cmd, shell=True)
		if ret != 0:
			print(f"rivet-build failed for {proc} (exit {ret}). STDOUT/STDERR:")
			if out: print(out)
			if err: print(err)
			return

	# 2) Prepare job list (process, operator, polarization) triplets
	jobs = [(proc, op, pol) for proc in PROCESSES for op in OPERATORS for pol in POLARISATION]
	if not jobs:
		print("No jobs to run (empty PROCESSES or OPERATORS or POLARISATION). Exiting.")
		return

	# 3) Run jobs in parallel
	# Make sure we don't create more workers than jobs so operators can run concurrently
	n_jobs = JOBS if JOBS > 0 else len(jobs)
	n_jobs = min(n_jobs, len(jobs))
	print(f"Launching {len(jobs)} job(s) with up to {n_jobs} parallel worker(s).")
	pool = multiprocessing.Pool(processes=n_jobs)
	results = []
	def _collect(res):
		# callback executed in parent when a job completes
		results.append(res)

	# submit all jobs asynchronously so multiple operators/polarisations run at the same time
	async_results = []
	for proc, op, pol in jobs:
		print(f"Submitting job: process={proc}, operator={op}, pol={pol}")
		async_results.append(pool.apply_async(worker, args=(proc, op, pol), callback=_collect))

	# close pool to new submissions and wait for all to finish
	pool.close()
	for ar in async_results:
		ar.wait()
	pool.join()

	# 4) Report results
	timestamp = time.strftime("%Y%m%d-%H%M%S")
	pid = os.getpid()
	cwd = os.getcwd()
	log_dir = os.path.join(cwd, "logs")
	os.makedirs(log_dir, exist_ok=True)
	log_path = os.path.join(log_dir, f"rivet_jobs_{type_mc_rivet}_{timestamp}_{pid}.txt")

	with open(log_path, "w") as logf:
		logf.write(f"Run timestamp: {timestamp}\n")
		logf.write(f"PID: {pid}\n")
		logf.write(f"Total jobs: {len(results)}\n\n")
		for proc, op, pol, conf, code, out, err in results:
			logf.write(f"=== {proc} | operator: {op} | pol: {pol} | conf: {conf} | exit: {code} ===\n")
			if out:
				logf.write("STDOUT:\n")
				logf.write(out if out.endswith("\n") else out + "\n")
			if err:
				logf.write("STDERR:\n")
				logf.write(err if err.endswith("\n") else err + "\n")
			logf.write("\n")



	# derive variables from this script's configuration
	path_hist_base = outdir
	order = order_EFT
	name_copy = getattr(opts, "name_copy", "copy")

	# build process/decay tuples from PROCESSES (expecting strings like "WpWp_lvlv")
	proc_decays_tuple = []
	for p in PROCESSES:
		if "_" in p:
			_process, _decay = p.split("_", 1)
		else:
			_process, _decay = p, ""
		proc_decays_tuple.append((_process, _decay))

	# classify operators that require polarisation (FM*) and those that don't (FS/FT)

	polarisations = [pol for pol in POLARISATION if pol is not None]

	# File operations for FM (with polarisation)
	for process, decay in proc_decays_tuple:
		for op in OPERATORS:
			for polarisation in polarisations:
				Conf = f"user.osalin.MadGraph_{process}_{decay}_{op}_{order}_{polarisation}"
				print("Conf: ", Conf)
				try:
					proc_name, Base_dir = find_prod_dec_and_dir_tres(Conf, type_mc_rivet)
				except ValueError as e:
					print(f"Error: {e}")
					continue
				Dir_ntuple = os.path.join(Base_dir, "DOCUT_YES")
				File_cross_section = os.path.join(Base_dir, "cross_section_fb.txt")
				File_log_file = os.path.join(Base_dir, "log.generate")
				path_hist = os.path.join(path_hist_base, f"{process}_{decay}", f"{op}_{order}_{polarisation}")

				if not os.path.exists(path_hist):
					os.makedirs(path_hist, exist_ok=True)
				else:
					print("Warning: directory already exists, creating timestamped fallback.")
					timestamp = time.strftime("%m%d-%H%M")
					path_hist = os.path.join(path_hist, "Already_exist", f"{name_copy}_{timestamp}", f"{process}_{decay}", f"{op}_{order}", type_mc_rivet)
					os.makedirs(path_hist, exist_ok=True)

				# Copy the content of Dir_ntuple into path_hist
				if os.path.exists(Dir_ntuple):
					#print(f"Copying {Dir_ntuple} to {path_hist}")
					shutil.copytree(Dir_ntuple, path_hist, dirs_exist_ok=True)
				else:
					print(f"Warning: ntuple directory not found: {Dir_ntuple}")
				if os.path.exists(File_cross_section):
					shutil.copy(File_cross_section, path_hist)
					#print(f"Copied cross section to {path_hist}")
				else:
					print(f"Warning: cross section file not found: {File_cross_section}")
				if os.path.exists(File_log_file):
					shutil.copy(File_log_file, path_hist)
					#print(f"Copied log file to {path_hist}")
				else:
					print(f"Warning: log file not found: {File_log_file}")
	print(f"Copy was succesful in {path_hist_base}")

	

if __name__ == "__main__":
	main()