import pyAMI.client
import pyAMI.atlas.api as AtlasAPI

import subprocess
import os
import re


from optparse import OptionParser
parser = OptionParser()
opts, _ = parser.parse_args()

client = pyAMI.client.Client('atlas')

AtlasAPI.init()

lifetime= "2592000"
rse="spacetoken=ATLASSCRATCHDISK&cloud=FR"


all_ops = ["FM0","FM1","FM2","FM3","FM4","FM5","FM7",
           "FS02","FS1",
           "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]
all_ops_ = [op + "_QUAD" for op in all_ops]

Processes = ["WmZ","WpZ","ZZ"]
Processes = ["WmZ","WpZ","WmWm","WpWm","WmWp","WpWp"]

Decay = ["vvqq"]
Decay = ["lvqq"]

Processes = ["WmZ","WpZ","ZZ","WmWm","WpWm","WpWp"]
Decay = ["llqq",'lvqq','vvqq']

def possible_process(proc, decay):
    
    valid_combinations = {
        "WmZ": ["llqq", "lvqq", "vvqq"],
        "WpZ": ["llqq", "lvqq", "vvqq"],
        "ZZ": ["llqq", "vvqq"],  # Only ll and vv decays are valid for ZZ
        "WmWm": ["lvqq"],
        "WpWm": ["lvqq"],
        "WpWp": ["lvqq"],
    } 

    return decay in valid_combinations.get(proc, [])

def get_datasets(process, decay):
    pattern = f"%%MGPy8EG_aQGCF%%_{process}_{decay}%%"
    #print(pattern)
    return AtlasAPI.list_datasets(client, patterns=[pattern], 
                                  fields=['ldn', 'cross_section', 'dataset_number'], 
                                  limit=[1, 600], type='EVNT', project='mc16%')

def organize_datasets(datasets):
    organized_info = {}
    for dataset in datasets:
        operator = next((op for op in all_ops_ if op in dataset['ldn']), None)
        if operator:
            for process in Processes:
                for decay in Decay:
                    if possible_process(process, decay):
                        if process in dataset['ldn'] and decay in dataset['ldn']:
                            key = f"{operator}_{process}_{decay}"
                            if key not in organized_info:
                                organized_info[key] = {'name': [], 'cross_section': [], 'dataset_number': []}
                            organized_info[key]['name'].append(dataset['ldn'])
                            organized_info[key]['cross_section'].append(dataset['cross_section'])
                            organized_info[key]['dataset_number'].append(dataset['dataset_number'])
    return organized_info

def get_first_dataset_numbers_and_names(organized_info):
    dataset_numbers = []
    names = []
    for op in all_ops_:
        for process in Processes:
            for decay in Decay:
                if possible_process(process, decay):
                    key = f'{op}_{process}_{decay}'
                    if key in organized_info:
                        dataset_numbers.append(organized_info[key]['dataset_number'][0])
                        names.append(organized_info[key]['name'][0])
                    else:
                        print(f"Warning: No datasets found for {key}")
    return dataset_numbers, names



List_info = {f'{process}_{decay}': get_datasets(process, decay) 
             for process in Processes for decay in Decay 
             if get_datasets(process, decay)}

organized_info = organize_datasets([dataset for datasets in List_info.values() for dataset in datasets])

dataset_numbers, names = get_first_dataset_numbers_and_names(organized_info)

#print(dataset_numbers)
#print(names)
name_test = names[6]

#example command : rucio add-rule --lifetime 2592000 mc16_13TeV.364852.MGPy8EvtGen_NNPDF30NLO_A14NNPDF23LO_WZjj_llqq_EW6.merge.EVNT.e7827_e7400_tid22682503_00 1 spacetoken=ATLASSCRATCHDISK&cloud=FR --skip-duplicates
def run_rucio_command(dataset_name,lifetime,rse):
    command = ["rucio", "add-rule", "--skip-duplicates", "--lifetime", lifetime, dataset_name, "1", rse]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    if process.returncode != 0:
        print(f"Error : {stderr.decode('utf-8')}")
    else:
        print(f"Command output: {stdout.decode('utf-8')}")  

def download_rucio_command(dataset_name, directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    command = ["rucio", "download", dataset_name, "--dir", directory]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    if process.returncode != 0:
        print(f"Error : {stderr.decode('utf-8')}")
    else:
        print(f"Command output: {stdout.decode('utf-8')}") 
        
        
for process in Processes:
    for decay in Decay:
        if possible_process(process,decay):
            key=f'{process}_{decay}'
            print(key)
            output_plot = f"/exp/atlas/salin/ATLAS/VBS_mc/EFT_files_AMI/{process}_{decay}/"
            print(output_plot)
            if not os.path.exists(output_plot):
                os.makedirs(output_plot)
            for name in names:
                if key in name:
                    print(name)
                    #run_rucio_command(name,lifetime,rse)
                    download_rucio_command(name, output_plot)
