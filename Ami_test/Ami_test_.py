import pyAMI.client
import pyAMI.atlas.api as AtlasAPI

import rucio.client


client = pyAMI.client.Client('atlas')

AtlasAPI.init()

all_ops = ["FM0","FM1","FM2","FM3","FM4","FM5","FM7",
           "FS02","FS1",
           "FT0","FT1","FT2","FT5","FT6","FT7"]
all_ops_ = [op + "_QUAD" for op in all_ops]
Processes = ["WmZ","WpZ","ZZ"]
Decay = ["llqq"]

def get_datasets(op, process, decay):
    pattern = f"%%MGPy8EG_aQGCF%%_{process}_{decay}%%"
    print(pattern)
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
                dataset_numbers.append(organized_info[f'{op}_{process}_{decay}']['dataset_number'][0])
                names.append(organized_info[f'{op}_{process}_{decay}']['name'][0])
    return dataset_numbers, names

List_info = {f'{op}_{process}_{decay}': get_datasets(op, process, decay) 
             for op in ["FM"] for process in Processes for decay in Decay 
             if get_datasets(op, process, decay)}

organized_info = organize_datasets([dataset for datasets in List_info.values() for dataset in datasets])

dataset_numbers, names = get_first_dataset_numbers_and_names(organized_info)

print(names)
    