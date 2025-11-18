import os
import lib_utils as lu
import matplotlib.pyplot as plt
from statistics import mean
import math
import numpy as np
import pandas as pd
from pandas.plotting import table
import subprocess
import glob
import re
import itertools
from array import array
from pandaclient import panda_api
from optparse import OptionParser

plt.rcParams.update({'text.usetex': True})
c = panda_api.get_api()
parser = OptionParser()

parser.add_option("--All_channel", default=True)
parser.add_option("--EFT_order", default="QUAD")
parser.add_option("--Rwg", default= False)
parser.add_option("--Channel", default="")
parser.add_option("--nb_lep", default=2)
parser.add_option('--type_MC', default="")
parser.add_option('--Add_stats', default=False)
opts, _ = parser.parse_args()


Processes=["WpZ"]
#Processes=["WpZ", "WmZ","ZZ"]
Decay = ['llqq','lvqq',"vvqq"]
Decay = ['llqq']
valid_combinaison_SM=["WZjj_llqq", "ZZjj_llqq", "WZjj_vvqq", "ZZjj_vvqq", "WZjj_lvqq", "WWjj_lvqq"]
# Replace previous hard-coded valid_combi_aQGC lists with a leptonic-only mapping
valid_combinations = {
    4: [("ZZ", "llll")],
    3: [("WmZ", "lllv"), ("WpZ", "lllv")],
    2: [("WmWm", "lvlv"), ("WpWm", "lvlv"), ("WpWp", "lvlv"), ("Zy", "lly")],
    1: [("Wpy", "lvy")],
}

# Flatten mapping to the same string format used elsewhere in the script (e.g. "ZZ_llll")
valid_combi_aQGC = [f"{proc}_{dec}" for pairs in valid_combinations.values() for proc, dec in pairs]

nb_lepton = int(opts.nb_lep)

order = opts.EFT_order
type_MC = opts.type_MC
name_panda= "Olivier Salin"
if "model" in type_MC or "aqgc" in type_MC:
    name_spe_task="aqgcModel_new_test_Produc"
elif "Run3" in type_MC or "run3" in type_MC:
    name_spe_task="Run13p6_Eboli_Prod_Eb"
    print("Run3")

elif "Reweighting_Polarisation" in type_MC:
    name_spe_task="Reweighting_RwgPolar_test01_50k"


elif "Reweighting_CPodd" in type_MC:
    name_spe_task="Rwg_CPodd_Pol_20k"

elif "Reweighting_VBSAll" in type_MC:
    name_spe_task="Rwg_Alt_missingOp_test_10k"

elif "Reweighting_VBSLeptonic" in type_MC:
    name_spe_task="Rwg_VBSLep_all_Pol_30k_30k"


elif "Reweighting_VBSSemiLeptonic" in type_MC:
    name_spe_task="Rwg_SemiLep_Pol_all_test_20k"

elif "Reweighting_INT" in type_MC:
    name_spe_task="Rwg_Rwg_RewgingINT"


elif "EFTDec_Madspin" in type_MC:
    name_spe_task="aqgcModel_EFTDec_Madspin_Core01"
    name_spe_task="aqgcModel_EFTDec_INT_better"


elif "EFTDec_VBSLeptonic" in type_MC:
    name_spe_task="Leptonic_testbis"
    
elif "EFTDec_Polarisation" in type_MC:
    name_spe_task="aqgcModel_EFTDec_Pol_test02"
    name_spe_task="aqgcModel_EFTDec_Pol_Cross_bis"

elif "EFTDec_CPodd" in type_MC:
    name_spe_task="testCross_Valid_ModdTodd"

elif "EFTDec_INT" in type_MC:
    name_spe_task="INT_intCPodd_50k"

else:
    name_spe_task=""
print("name_spe_task", name_spe_task)

#all_ops_cat = ["SM", "FM0", "FM2", "FS1", "FT1", "FT5"]

all_ops_cat = ["FM0","FM1","FM2","FM3","FM4","FM5","FM7","FM8","FM9",
            "FS0","FS1","FS2",
            "FT0","FT1","FT2","FT3","FT4","FT5","FT6","FT7","FT8","FT9"]
all_ops_cat = ["FM0","FM1","FM2","FM3","FM4","FM5","FM7","FM8","FM9"]
#all_ops_cat = ["FT0","FT1","FT2","FT3","FT4","FT5","FT6"]
#all_ops_cat = ["FS0","FS1","FS2"]

all_ops_cat = ["FM","FT","FS"]
all_ops_cat = ["FM0","FM1","FM2","FM3","FM4","FM5","FM7","FM8","FM9",
            "FS0","FS1","FS2",
            "FT0","FT1","FT2","FT3","FT4","FT5","FT6"]

all_ops_cat= ["FM0","FM1","FM2","FM3","FM4","FM5","FM7","FM8","FM9",
            "FS0","FS1","FS2",
            "FT0","FT1","FT2","FT3","FT4","FT5","FT6","FT7","FT8","FT9",
            "FM1odd","FM2odd","FM3odd","FM4odd","FM5odd","FM6odd","FT2odd","FT3odd","FT4odd","FT5odd","FT6odd"]

all_ops_cat = ["FM0","FM1","FM2","FM3","FM4","FM5","FM7"]
all_ops_cat = ["FM1odd","FM2odd","FM3odd","FM4odd","FM5odd","FM6odd","FT2odd","FT3odd","FT4odd","FT5odd","FT6odd"]

all_ops_cat= ["FM0","FM1","FM2","FM3","FM4","FM5","FM7","FM8","FM9",
            "FS0","FS1","FS2",
            "FT0","FT1","FT2","FT3","FT4","FT5","FT6","FT7","FT8","FT9",
            "FM1odd","FM2odd","FM3odd","FM4odd","FM5odd","FM6odd","FT2odd","FT3odd","FT4odd","FT5odd","FT6odd"]

#all_ops_cat = ["FM0","FT0"]
if opts.Rwg:
    #all_ops_cat = ["FS"]
    all_ops_cat = ["FM","FS","FT","FModd","FTodd"]


cross_terms_ = [f"{op1}vs{op2}" for op1, op2 in itertools.combinations(all_ops_cat, 2)]


print("Defined cross terms:", cross_terms_)

if opts.EFT_order == "CROSS":
    all_ops_cat = cross_terms_

#all_ops_cat=cross_terms_ 
#all_ops_cat =["FM0"]
#all_ops_cat= ["FM0","FM2","FS1","FT0","FT1","FT5","FT8"]

if "aqgc" in type_MC or "model" in type_MC:
        all_ops_cat = ["FM0","FM1","FM2","FM3","FM4","FM5","FM7","FM8","FM9",
            "FS0","FS1","FS2",
            "FT0","FT1","FT2","FT3","FT4","FT5","FT6","FT8",
            "FM1odd","FM2odd","FM3odd","FM6odd","FT2odd","FT4odd","FT6odd"]
        #all_ops_cat = ["FM0","FS0","FT0"]


        
base_path = "/exp/atlas/salin/ATLAS/VBS_mc/eft_files/"
base_path = "/data/atlas/salin/VBS_mc/VBS/eft_Files/"
base_dir = f"{base_path}/{'13p0/' if '13p0' in type_MC or '13p0' in type_MC else 'aqgc_model/' if 'aqgc' in type_MC or 'model' in type_MC else ''}"

if "Reweighting_Polarisation" in type_MC:
    base_dir = f"{base_path}/Reweighting/Polarisation/"
    
elif "Reweighting_CPodd" in type_MC:
    base_dir = f"{base_path}/Reweighting/CPodd/Validation/Polarisation/"

elif "Reweighting_VBSAll" in type_MC:
    base_dir = f"{base_path}/Reweighting/VBSAll/Polarisation/test/Gen/"

elif "Reweighting_VBSLeptonic" in type_MC:
    base_dir = f"{base_path}/Reweighting/VBSAll/Pol/Gen/VBS/Leptonic/"

elif "Reweighting_VBSSemiLeptonic" in type_MC:
    base_dir = f"{base_path}/Reweighting/VBSAll/Pol/Gen/VBS/SemiLep/"

elif "Reweighting_INT" in type_MC:
    base_dir = f"{base_path}/Reweighting/INT/Validation/Stats/"

elif "EFTDec_Madspin" in type_MC:
    base_dir = f"{base_path}/EFTDec/Madspin/"

elif "EFTDec_VBSLeptonic" in type_MC:
    base_dir = f"{base_path}/EFTDec/VBSAll/Pol/Leptonic/"

elif "EFTDec_Polarisation" in type_MC:
    base_dir = f"{base_path}/EFTDec/Polarisation/"
elif "EFTDec_CPodd" in type_MC:
    base_dir = f"{base_path}/EFTDec/CPodd/Validation/"

elif "EFTDec_INT" in type_MC:
    base_dir = f"{base_path}/EFTDec/INT/Validation/"

else:
    base_dir = f"{base_path}/BAD_attempt/"
    exception_message = f"Invalid type_MC: {type_MC}. Please check the input."
    raise ValueError(exception_message)
os.makedirs(base_dir, exist_ok=True)
base_dir_SM = f"{base_path}/" if not base_dir else ""
#base_path= "/exp/atlas/salin/ATLAS/VBS_mc/eft_files/Run3/test_WpZ_llqq"

name_run = "Type_MC_test"

def get_evnt_log_files(base_dir, i_job_name):
    evnt_did, evnt_dir, log_did, log_dir_before_untar = lu.get_envt_log_names_dirs(base_dir, i_job_name)
    evnt_file = -1
    log_file = -1
    if os.path.exists(evnt_dir) and os.path.exists(log_dir_before_untar):
        print("directories for evnt and log exist")
        evnt_candidates = glob.glob(evnt_dir + "/*EVNT.root")
        log_candidates = glob.glob(log_dir_before_untar + "/tarball_PandaJob*/log.generate")
        print("evnt candidates of len", len(evnt_candidates), evnt_candidates)
        print("log candidates of len", len(log_candidates), log_candidates)
        if len(evnt_candidates) == 1 and len(log_candidates) == 1:
            evnt_file = evnt_candidates[0]
            log_file = log_candidates[0]
    else:
        print("directories for evnt and log DOESN'T exist")

    print("returning evnt file", evnt_file)
    print("returning log file", log_file)
    return evnt_file, log_file

def extract_jobName_info(job_template):
    parts = job_template.split('.')
    if len(parts) < 3:
        print(f"Invalid job template format: {job_template}")
        return None, None, None, None, None

    main_parts = parts[2].split('_')
    
    if len(main_parts) < 5:
        print(f"Invalid job template format: {job_template}")
        return None, None, None, None, None

    process = main_parts[1]
    decay = main_parts[2]
    operator = main_parts[3]
    order = main_parts[4]
    nickname = '_'.join(main_parts[5:]) if len(main_parts) > 5 else ""

    return process, decay, operator, order, nickname

def prepare_grid_files(Job_name, base_dir=base_dir):
    Event_files = {}
    Log_files = {}
    Cross_sections = {}
    Cross_sections_bef_decay = {}
    X_sections = {}
    X_sections_bef_decay = {}
    
    for i_job_name in Job_name:
        try:
            print("will download+untar files evnt and log for", i_job_name)
            proc, dec, op, order, nickname = extract_jobName_info(i_job_name)
            print(f"Process: {proc}, Decay: {dec}, Operator: {op}, Order: {order}, Nickname: {nickname}") 
            print(f"Base directory: {base_dir}")
            proc_dec = f"{proc}_{dec}" 
            proc_dec_op_order = f"{proc}_{dec}_{op}_{order}"     

            evnt_did, evnt_dir, log_did, log_dir = lu.get_envt_log_names_dirs(base_dir, i_job_name)

            evnt_dir_pattern = re.sub(rf'_{order}.*?_EXT0', f'_{order}*_EXT0', evnt_dir)
            matching_evnt_dirs = glob.glob(evnt_dir_pattern)
            print(f"Event dir: {evnt_dir}\n, pattern: {evnt_dir_pattern},\n matching dirs: {matching_evnt_dirs}")
            print("will download", evnt_did, "since no matching dir found for pattern", evnt_dir_pattern)
            subprocess.call(f"rucio download {evnt_did}", shell=True, cwd=base_dir)
            if not matching_evnt_dirs:
                print("will download", evnt_did, "since no matching dir found for pattern", evnt_dir_pattern)
                #subprocess.call(f"rucio download {evnt_did}", shell=True, cwd=base_dir)
            else:
                evnt_dir = matching_evnt_dirs[0]
                log_dir = matching_evnt_dirs[0] + f"/{log_did}/"
                
                if opts.Add_stats:
                    print("adding stats")
                    os.makedirs(evnt_dir + "/temp/", exist_ok=True)
                    print(f"Creating Temporary directory")
                    subprocess.call(f"rucio download {evnt_did}", shell=True, cwd=evnt_dir + "/temp/")
                    temp_files = glob.glob(f"{evnt_dir}/temp/{evnt_did}/*.root")
                    existing_files = glob.glob(f"{evnt_dir}/*.root")
                    existing_file_names = [os.path.basename(f) for f in existing_files]
                    for temp_file in temp_files:
                        if os.path.basename(temp_file) not in existing_file_names:
                            print("Will copy", temp_file)
                            subprocess.call(f"cp {temp_file} {evnt_dir}", shell=True)
                    subprocess.call(f"rm -r {evnt_dir}/temp/", shell=True)
                    print(f"Removing Temporary directory")
                else:
                    print("have", evnt_dir, "already so don't download")
                    
            event_files = glob.glob(f"{evnt_dir}/*.root")
            Event_files[proc_dec_op_order] = event_files
            evnt_did, evnt_dir, log_did, log_dir = lu.get_envt_log_names_dirs(base_dir, i_job_name)
            log_dir_pattern = re.sub(rf'EXT0/.*?_{order}.*?\.log', f'EXT0/*_{order}_*.log', log_dir)
            matching_log_dirs = glob.glob(log_dir_pattern)
            print(f"Log dir: {log_dir}\n, pattern: {log_dir_pattern},\n matching dirs: {matching_log_dirs}")
            print("LOG ill download", log_did, "since dir doesn't exist", log_dir)
            subprocess.call(f"rucio download {log_did}", shell=True, cwd=evnt_dir)
            if not matching_log_dirs:
                print("LOG ill download", log_did, "since dir doesn't exist", log_dir)
                #subprocess.call(f"rucio download {log_did}", shell=True, cwd=evnt_dir)
            else:
                log_dir = matching_log_dirs[0]
                os.makedirs(evnt_dir + "/Log_other/", exist_ok=True)
                temp_log = glob.glob(evnt_dir + "/Log_other/" + f"{log_did}/")
                print("temp_log", temp_log)

                subprocess.call(f"rucio download {log_did}", shell=True, cwd=evnt_dir + "/Log_other/")
                print("have LOG", "already so download in other folder", log_dir + "/Log_other/")

            
            print(f'Log dir: {log_dir}')
            untared_dir_cand = glob.glob(f"{log_dir}/tarball_PandaJob*")
            print(f"Untared dir candidates: {untared_dir_cand}")


            if len(untared_dir_cand) == 0:
                tar_file_candidates = sorted(glob.glob(f"{log_dir}/*log.tgz"))
                print("tar_file_candidates", tar_file_candidates)
                tar_file = os.path.basename(tar_file_candidates[-1])
                subprocess.call(f"tar -xvf {tar_file}", shell=True, cwd=log_dir)
            else:
                print("dont untar since did it before res in ", untared_dir_cand[0])
                
            log_files = glob.glob(f"{log_dir}/tarball_PandaJob*/log.generate")
            print(f"Log files: {log_files}")
                
            Log_files[proc_dec_op_order] = log_files
            xsec = lu.get_xsec(log_files[0])
            subprocess.call(f"cp {log_files[0]} {evnt_dir}", shell=True)
            ## Cross section before decay
            xsec_value,x_sec_unc = lu.get_xsec_bef_decay(log_files[0])
            xsec_bef_decay= f'{xsec_value} +- {x_sec_unc}'
            
            
            Cross_sections[proc_dec_op_order] = xsec
            X_sections[f"{op}_{order}_{proc_dec}"] = xsec
            
            Cross_sections_bef_decay[proc_dec_op_order] = xsec_bef_decay
            X_sections_bef_decay[f"{op}_{order}_{proc_dec}"] = xsec_bef_decay
            
            with open(os.path.join(evnt_dir, "cross_section_fb.txt"), "w") as f:
                f.write(str(xsec))
        
            print(f"Event file: {event_files}, Log file: {log_files}, xsec file: {xsec}")
        except Exception as e:
            print(f"Error processing job {i_job_name}: {e}")
            continue
    
    dir_xsec = base_dir + "/Cross-section/"
    os.makedirs(dir_xsec, exist_ok=True)
    try:
        with open(os.path.join(dir_xsec, "VBS_xsection_fb.txt"), "w") as f:
            for key, value in X_sections.items():
                f.write(f"{key}: {value}\n")
        with open(os.path.join(dir_xsec, "VBS_xsection_bef_decay_pb.txt"), "w") as f:
            for key, value in X_sections_bef_decay.items():
                f.write(f"{key}: {value}\n")
        
        subprocess.call(f"cp {os.path.join(dir_xsec, 'VBS_xsection_fb.txt')} {os.path.join(dir_xsec, f'VBS_xsection_fb_{nickname}.txt')}", shell=True)
    except Exception as e:
        print(f"Error writing cross section file: {e}")

    # Write the Cross_section dictionary to a text file in the same format as VBS_xsection.txt
    try:
        with open(os.path.join(dir_xsec, f"VBS_cross_section_{nickname}.txt"), "w") as f:
            for prod_dec, sections in Cross_sections.items():
                for key, value in sections.items():
                    f.write(f"{key}_{prod_dec}: {value}\n")
                    
        with open(os.path.join(dir_xsec, f"VBS_bef_decay_cross_section_{nickname}.txt"), "w") as f:
            for prod_dec, sections in Cross_sections_bef_decay.items():
                for key, value in sections.items():
                    f.write(f"{key}_{prod_dec}: {value}\n")
    except Exception as e:
        print(f"Error writing cross section file: {e}")

    return Event_files, Log_files, X_sections, X_sections_bef_decay

valid_combi_aQGC = [f"{proc}_{dec}" for pairs in valid_combinations.values() for proc, dec in pairs]
print("valid_combi_aQGC", valid_combi_aQGC)

TaskName = {}    
Event_files={}
Log_files={}   
Cross_section={}
Cross_section_bef_decay={}

for prod_dec in valid_combi_aQGC:
    process, decay = prod_dec.split("_")
    print(f"Process: {process}, Decay: {decay}")
    if prod_dec in valid_combinaison_SM or prod_dec in valid_combi_aQGC: 
        print(f"{prod_dec} is a good combinaison")           
        tasks = c.get_tasks(limit=100000000, days=13000, username=name_panda, status="done") # get already last try since only retry if it failed
        print(f"Number of tasks: {len(tasks)}")
        print(f"Name special for task: {name_spe_task}")

        task_names = [i_task['taskname'].replace("/","") for i_task in tasks if name_spe_task in i_task['taskname'] and prod_dec in i_task['taskname']]
        print(f"task_names: {task_names}")
        # Filter task_names to keep only those with an operator from all_ops_cat
        print(f"Filtering task names for operators: {all_ops_cat}, eft order: {order}")
        op__= [f"_{op}_" for op in all_ops_cat if op in task_names[0]]
        filtered_task_names = [task for task in task_names if any(op in task for op in all_ops_cat) & any(order_eft in task for order_eft in [order] )]
        print(f"Filtered task names: {filtered_task_names}")
        base_dir_prod = base_dir + f"{prod_dec}/"
        os.makedirs(base_dir_prod, exist_ok=True)
        if filtered_task_names:
            found_ops = [op for op in all_ops_cat if any(op in task for task in filtered_task_names)]
            print(f"Found valid task names for {prod_dec} with operators from {all_ops_cat}: {found_ops}")
            TaskName[f"{prod_dec}"] = filtered_task_names
            print(f"task_names: {filtered_task_names}")
            print("Jobname", extract_jobName_info(filtered_task_names[0]))


            Evnt_f,Log_f,X_sec,X_sec_bef_decay=prepare_grid_files(TaskName[f"{prod_dec}"], base_dir=base_dir_prod)
        else:
            missing_ops = [op for op in all_ops_cat if not any(op in task for task in task_names)]
            print(f"No valid task names found for {prod_dec} with operators from {all_ops_cat}. Missing operators: {missing_ops}")

        Event_files[f"{prod_dec}"] = Evnt_f
        Log_files[f"{prod_dec}"] = Log_f
        Cross_section[f"{prod_dec}"] = X_sec
        Cross_section_bef_decay[f"{prod_dec}"] = X_sec_bef_decay
        

# Write the Cross_section dictionary to a text file in the same format as VBS_xsection.txt
dir_cross_section = base_dir + "/Cross-section/"
os.makedirs(dir_cross_section, exist_ok=True)
try:
    with open(os.path.join(dir_cross_section, f"VBS_cross_section_{name_spe_task}.txt"), "w") as f:
        for prod_dec, sections in Cross_section.items():
            for key, value in sections.items():
                f.write(f"{key}: {value}\n")
                
    with open(os.path.join(dir_cross_section, f"bef_decay_VBS_cross_section_{name_spe_task}.txt"), "w") as f:
        for prod_dec, sections in Cross_section_bef_decay.items():
            for key, value in sections.items():
                f.write(f"{key}: {value}\n")
except Exception as e:
    print(f"Error writing cross section file: {e}")











