#!/usr/bin/env python3
import subprocess
import multiprocessing
import sys
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--nb_lep", default = 2)
parser.add_option("--EFT_order", default = "QUAD")
parser.add_option("--type_mc", default = "")

(opts, args) = parser.parse_args()
nb_lepton= int(opts.nb_lep)
order_EFT= opts.EFT_order
type_mc_rivet = opts.type_mc
# --- Configuration: edit these lists/values directly ---

valid_combinations = {
    4: [("ZZ", "llll")],
    3: [("WmZ", "lllv"), ("WpZ", "lllv")],
    2: [("WmWm", "lvlv"), ("WpWm", "lvlv"), ("WpWp", "lvlv")],
    1: [("Wpy", "lvy")],
}
valid_combi_aQGC = [f"{proc}_{dec}" for pairs in valid_combinations.values() for proc, dec in pairs]
valid_combi_aQGC_lep = [f"{proc}_{dec}" for proc, dec in valid_combinations[nb_lepton]]
print("Selected processes for nb_lep =", nb_lepton, ":\n", valid_combi_aQGC_lep)

PROCESSES = valid_combi_aQGC_lep          # list of processes (example: ["ZZ_llll"])
#PROCESSES = ["WpWp_lvlv"]          # list of processes (example: ["ZZ_llll"])
OPERATORS = ["FS"] 
POLARISATION= [None]
POLARISATION= ["LL","LT","TL","TT"]      
EVTMAX = 1000
DOCUT = "YES"
REDO_RIVET = "yes"
REDO_PLOTS = "no"
TYPE_MC = type_mc_rivet
JOBS = 0                         # 0 => use number of jobs = number of operators (per process)
PYTHON_CMD = sys.executable      # use same python interpreter that runs this script
# --- end configuration ---

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
	elif "FM" in op:
		conf = f"user.osalin.MadGraph_{process}_{op_suff}_{pol}"
	else:
		conf = f"user.osalin.MadGraph_{process}_{op_suff}"
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
	for proc, op, pol, conf, code, out, err in results:
		print(f"=== {proc} | operator: {op} | pol: {pol} | conf: {conf} | exit: {code} ===")
		if out:
			print("STDOUT:")
			print(out)
		if err:
			print("STDERR:")
			print(err)

if __name__ == "__main__":
	main()