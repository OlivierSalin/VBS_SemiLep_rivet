import glob
import os
from array import array

import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import pandas as pd
import numpy as np
import json
import re
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--conf", default = "user.okurdysh.MadGraph_WmWm_lvlv_FT0_FULL")
parser.add_option("--DOCUT", default = "NO")
opts, _ = parser.parse_args()



Sample_name_AMI_SM="mc16_13TeV.364852.MGPy8EvtGen_NNPDF30NLO_A14NNPDF23LO_WZjj_llqq_EW6.merge.EVNT.e7827_e7400_tid22682503_00"
Sample_name_="user.osalin.MadGraph_WpZ_llqq_FM1_QUAD_EXT0"
Sample_name_AMI_QUAD="MGPy8EG_aQGCFM0_QUAD_1_WmZ_llqq.merge.EVNT"

pattern = r"mc16_13TeV\.(\d+)\.MGPy8EG_aQGC(.*)_(.*)_1_(.*)_(.*)\.merge"

# Use the re.search function to find the pattern in the sample name
match = re.search(pattern, Sample_name_AMI_QUAD)

if match:
    dataset_number = match.group(1)
    operator,operator_type = match.group(2), match.group(3)
    process, decay = match.group(4), match.group(5)
    process_decay=process+"_"+decay


    
    
def extract_EFT_op_proces_dec_bis(conf):
    """
    This function takes conf as a parameter,
    and returns the EFT operator, process, and decay.
    """
    pattern = r"MGPy8EG_aQGC(.*)_(.*)_1_(.*)_(.*)"
    match = re.search(pattern, conf)
    if match:
        EFT_op = match.group(1)
        EFT_type = match.group(2)
        proc = match.group(3)
        decay = match.group(4)
        
    elif conf.startswith("user."):
        # If the conf matches the Sample_name_ pattern
        parts = conf.split('_')
        EFT_op,EFT_type, proc, decay = parts[3], parts[4], parts[1], parts[2]
    else:
        raise ValueError("Invalid conf format")

    return EFT_op, EFT_type, proc, decay

EFT_op, EFT_type, proc, decay = extract_EFT_op_proces_dec_bis(opts.conf)

keyy = f"{EFT_op}_{EFT_type}_{proc}_{decay}"

print("Process studied:", keyy)  


def get_conf_cut_dir(evnt_dir, docut):
    mydir = evnt_dir + f"/DOCUT_{docut}/"
    if not os.path.exists(mydir): os.makedirs(mydir) 
    return mydir

def get_job_name_template(prod,dec,op,EFTmode):
    return f"user.osalin.MadGraph_{prod}_{dec}_{op}_{EFTmode}"



def find_prod_dec_and_dir(conf):

    prod_temp = conf[conf.find("user.osalin.MadGraph_")+len("user.osalin.MadGraph_"):]
    print("start from string", prod_temp)
    prod_dec = prod_temp[:prod_temp.find("_F")]
    print("from conf found production dec", prod_dec)
    conf_dir = f"/exp/atlas/salin/ATLAS/VBS_mc/eft_files/{prod_dec}/"
    print("dir would be", conf_dir)
    return prod_dec, conf_dir

print("\nFind prod and decay dir function:")
def find_prod_dec_and_dir_bis(conf):
    if conf.startswith("user."):
        prod_temp = conf[conf.find("user.osalin.MadGraph_")+len("user.osalin.MadGraph_"):]
        print("start from string", prod_temp)
        prod_dec = prod_temp[:prod_temp.find("_F")]
        print("from conf found production dec", prod_dec)
        conf_dir = f"/exp/atlas/salin/ATLAS/VBS_mc/eft_files/{prod_dec}/"
        print("dir would be", conf_dir)
    else:
        pattern = r"MGPy8EG_aQGC(.*)_(.*)_1_(.*)_(.*)"
        match = re.search(pattern, conf)
        if match:
            prod_dec = match.group(3) + "_" + match.group(4)
            print("From conf found production dec", prod_dec)
            conf_dir = f"/exp/atlas/salin/ATLAS/VBS_mc/EFT_files_AMI/{prod_dec}/"
            print("Dir would be", conf_dir)
        else:
            raise ValueError("Invalid conf format")
    return prod_dec, conf_dir

prod_dec, base_dir =find_prod_dec_and_dir_bis(opts.conf)
print("\nFind evnt dir and file function:")
def find_evnt_dir_and_file_bis(base_dir, conf):
    if conf.startswith("user."):
        search_com= base_dir + f"/*{opts.conf}*EXT0"
        print("searching for dir with pattern", search_com)
        conf_dir_arr = glob.glob(search_com)
        print("found possibilities for dir", conf_dir_arr)
        conf_dir = conf_dir_arr[0] if len(conf_dir_arr)>=1 else -1  
        if conf_dir == -1: raise ValueError("did not find folder for this config ",search_com)

        evnt_file_candidates = glob.glob(conf_dir + "/*EVNT.root")
        print("found possibilities for evnt file", evnt_file_candidates)
        evnt_file = evnt_file_candidates[0] if len(evnt_file_candidates)>0 else -1
        if evnt_file == -1: raise ValueError("did not find EVNT file for this config ",conf_dir)
    
    else:
        pattern = r"MGPy8EG_aQGC(.*)_(.*)_1_(.*)_(.*)"
        match = re.search(pattern, conf)
        if(match):
            search_com = base_dir + f"/*{conf}*"
            print("searching for dir with pattern", search_com)
            conf_dir_arr = glob.glob(search_com)
            print("found possibilities for dir", conf_dir_arr)
            conf_dir = conf_dir_arr[0] if len(conf_dir_arr)>=1 else -1  
            if conf_dir == -1: raise ValueError("did not find folder for this config ",search_com)
            
            evnt_file_candidates = [file for file in glob.glob(conf_dir + "/*EVNT*.pool.root.1") if not file.endswith('.part')]
            print("found possibilities for evnt file", evnt_file_candidates)
            evnt_file = evnt_file_candidates[0] if len(evnt_file_candidates)>0 else -1
        

    return conf_dir, evnt_file, evnt_file_candidates

_,_,_=find_evnt_dir_and_file_bis(base_dir,opts.conf)