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
parser.add_option("--Channel", default="")
parser.add_option("--nb_lep", default=2)
parser.add_option('--type_MC', default="")
parser.add_option('--Add_stats', default=False)
opts, _ = parser.parse_args()

Processes = ["WmZ", "WpZ", "ZZ", "WmWm", "WpWm", "WpWp","WZjj","ZZjj","WWjj"]
Processes=["WmZ", "WpZ", "ZZ", "WmWm", "WpWm", "WpWp"]
Processes=["WmZ", "WpZ", "ZZ"]
Decay = ['llqq','lvqq', 'vvqq']
valid_combinaison_SM=["WZjj_llqq", "ZZjj_llqq", "WZjj_vvqq", "ZZjj_vvqq", "WZjj_lvqq", "WWjj_lvqq"]
valid_combi_aQGC=["WpZ_llqq","WpZ_vvqq","WpZ_lvqq","WmZ_llqq","WmZ_vvqq","WmZ_lvqq","ZZ_llqq","ZZ_vvqq"
                  ,"WmWm_lvqq","WpWm_lvqq","WpWp_lvqq"]
#Decay = ['llqq']

Processes=["WWjj"]
Decay = ["lvqq"]
#Processes=["WZjj"]
#Decay = ["llqq"]

nb_lepton = int(opts.nb_lep)

order = opts.EFT_order
type_MC = opts.type_MC
name_panda= "Olivier Salin"
if "model" in type_MC or "aqgc" in type_MC:
    name_spe_task="aqgcModel"
    name_spe_task="aqgcModel_new_1"
    name_spe_task="aqgcModel_new_test_Produc"
elif "Run3" in type_MC or "run3" in type_MC:
    name_spe_task="Run3_13p6_"
    name_spe_task="Run13p6_prod_100k"
    name_spe_task="Run13p6_Eboli_Prod_Eb"
    print("Run3")
#elif "13p0" in type_MC or "13p0TeV" in type_MC:
    #name_spe_task="Run2_13p0_ProdTest_Batch_Nev50k"
elif "DynScale" in type_MC or "dynScale" in type_MC:
    name_spe_task="Run2_13p0_DynSca_Customsrv_50k"
elif "dyn_nocut" in type_MC:
    name_spe_task="Run2_13p0_ProdAlter_dyn_noCut_Batch_Nev50k"
    name_spe_task="Run2_13p0_dyn2_noCut_Batch5k_Nev50k"
elif "MCprod_13p0" in type_MC:
    name_spe_task="BFilter"
elif "MCprod_13p6" in type_MC:
    name_spe_task="Run3_13p6_Prod_DynScaAll_NEventsBIS_60k"
else:
    name_spe_task="heugebe"
print("name_spe_task", name_spe_task)

#all_ops_cat = ["SM", "FM0", "FM2", "FS1", "FT1", "FT5"]


all_ops_cat = [""]
cross_terms_ = [f"{op1}vs{op2}" for op1, op2 in itertools.combinations(all_ops_cat, 2)]
print("Defined cross terms:", cross_terms_)

#all_ops_cat=cross_terms_ 
#all_ops_cat =["FM0"]
#all_ops_cat= ["FM0","FM2","FS1","FT0","FT1","FT5","FT8"]




        
base_path = "/exp/atlas/salin/ATLAS/VBS_mc/eft_files/"
base_dir = f"{base_path}/{'13p0/' if '13p0' in type_MC or '13p0' in type_MC else 'aqgc_model/' if 'aqgc' in type_MC or 'model' in type_MC else ''}"
if "Run3" in type_MC or "run3" in type_MC:
    base_dir = f"{base_path}/Run3/"
elif "Run2" in type_MC or "run2" in type_MC:
    base_dir = f"{base_path}/Run2/"
elif "DynScale" in type_MC or "dynScale" in type_MC:
    base_dir = f"{base_path}/13p0/DynScale/"
#elif "13p0" in type_MC or "13p0TeV" in type_MC:
    #base_dir = f"{base_path}/13p0/"
elif "aqgc" in type_MC or "model" in type_MC:
    base_dir = f"{base_path}/aqgc_model/"
elif "dyn_nocut" in type_MC or "Dyn_nocut" in type_MC:
    base_dir = f"{base_path}/13p0Alt/dyn_nocut/"
elif "MCprod_13p0" in type_MC:
    base_dir = f"{base_path}/MCprod/SM/13p0/"
elif "MCprod_13p6" in type_MC:
    base_dir = f"{base_path}/MCprod/SM/13p6/"

else:
    base_dir = f"{base_path}/Fails/"
os.makedirs(base_dir, exist_ok=True)
#base_dir_SM = f"{base_path}/" if not base_dir else ""
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
            
            # Check for BFilter or BVeto in the job name
            if "BFilter" in i_job_name:
                proc_dec = f"{proc}_{dec}_BFilter"
            elif "BVeto" in i_job_name:
                proc_dec = f"{proc}_{dec}_BVeto"
            else:
                proc_dec = f"{proc}_{dec}"
            
            proc_dec_op_order = f"{proc}_{dec}_SM"

            evnt_did, evnt_dir, log_did, log_dir = lu.get_envt_log_names_dirs(base_dir, i_job_name)

            evnt_dir_pattern = re.sub(rf'_{order}.*?_EXT0', f'_{order}*_EXT0', evnt_dir)
            matching_evnt_dirs = glob.glob(evnt_dir_pattern)
            print(f"Event dir: {evnt_dir}\n, pattern: {evnt_dir_pattern},\n matching dirs: {matching_evnt_dirs}")
            if not matching_evnt_dirs:
                print("will download", evnt_did, "since no matching dir found for pattern", evnt_dir_pattern)
                subprocess.call(f"rucio download {evnt_did}", shell=True, cwd=base_dir)
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
            
            log_dir_pattern = re.sub(rf'EXT0/.*?_{order}.*?\.log', f'EXT0/*_{order}_*.log', log_dir)
            matching_log_dirs = glob.glob(log_dir_pattern)
            print(f"Log dir: {log_dir}\n, pattern: {log_dir_pattern},\n matching dirs: {matching_log_dirs}")
            if not matching_log_dirs:
                print("LOG ill download", log_did, "since dir doesn't exist", log_dir)
                subprocess.call(f"rucio download {log_did}", shell=True, cwd=evnt_dir)
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
            subprocess.call(f"cp {log_files[0]} {evnt_dir}", shell=True)
            xsec = lu.get_xsec(log_files[0])
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


    return Event_files, Log_files, X_sections, X_sections_bef_decay



processes, decays = Processes, Decay

TaskName = {}    
Event_files={}
Log_files={}   
Cross_section={}
Cross_section_bef_decay={}
for process in processes:
    for decay in decays:
        
        
        prod_dec = f"{process}_{decay}"
        print(f"Process: {process}, Decay: {decay}")
        if prod_dec in valid_combinaison_SM: 
            print(f"{prod_dec} is a good combinaison")           
            tasks = c.get_tasks(limit=100000000, days=13000, username=name_panda, status="done") # get already last try since only retry if it failed
            print(f"Number of tasks: {len(tasks)}")
            print(f"Name special for task: {name_spe_task}")

            task_names = [i_task['taskname'].replace("/","") for i_task in tasks if prod_dec in i_task['taskname']]
            print(f"task_names: {len(task_names)}")
            # Filter task_names to keep only those with an operator from all_ops_cat
            print(f"Filtering task names for proc decay: {prod_dec}")
            filtered_task_names = [task for task in task_names if name_spe_task in task and prod_dec in task]
            print(f"Filtered task names: {filtered_task_names}, len: {len(filtered_task_names)}")
            if len(filtered_task_names) != 0:
                if any("BVeto" in task for task in filtered_task_names):
                    print('BVETO if')
                    base_dir_prod = base_dir + f"{prod_dec}_BVeto/"
                    os.makedirs(base_dir_prod, exist_ok=True)
                    filtered_task_names_BVeto = [task for task in filtered_task_names if "BVeto" in task]
                    TaskName[f"{prod_dec}_BVeto"] = filtered_task_names_BVeto
                    Evnt_f,Log_f,X_sec,X_sec_bef_decay=prepare_grid_files(TaskName[f"{prod_dec}_BVeto"], base_dir=base_dir_prod)            
                # Adjust base directory for BFilter and BVeto
                if "BFilter" in filtered_task_names[0]:
                    print('BFilter if')
                    base_dir_prod = base_dir + f"{prod_dec}_BFilter/"
                    os.makedirs(base_dir_prod, exist_ok=True)
                    filtered_task_names_BFilter = [task for task in filtered_task_names if "BFilter" in task]
                    TaskName[f"{prod_dec}_BFilter"] = filtered_task_names_BFilter
                    print(f"task_names: {filtered_task_names}")
                    print("Jobname", extract_jobName_info(filtered_task_names[0]))
                    #Evnt_f,Log_f,X_sec,X_sec_bef_decay=prepare_grid_files(TaskName[f"{prod_dec}_BFilter"], base_dir=base_dir_prod)


                else:
                    base_dir_prod = base_dir + f"{prod_dec}/"
                
                    os.makedirs(base_dir_prod, exist_ok=True)
                    if filtered_task_names:
                        print(f"Found samples for {prod_dec} and with spe name: {name_spe_task}")
                        TaskName[f"{prod_dec}"] = filtered_task_names
                        print(f"task_names: {filtered_task_names}")
                        print("Jobname", extract_jobName_info(filtered_task_names[0]))


                        Evnt_f,Log_f,X_sec,X_sec_bef_decay=prepare_grid_files(TaskName[f"{prod_dec}"], base_dir=base_dir_prod)
                    else:
                        print(f"Not Found samples for {prod_dec} and with spe name: {name_spe_task}")

                Event_files[f"{prod_dec}"] = Evnt_f
                Log_files[f"{prod_dec}"] = Log_f
                Cross_section[f"{prod_dec}"] = X_sec
                Cross_section_bef_decay[f"{prod_dec}"] = X_sec_bef_decay
            else:
                print(f"No task names found for {prod_dec} with special name: {name_spe_task}")
        

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











