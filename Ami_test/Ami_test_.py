import pyAMI.client
import pyAMI.atlas.api as AtlasAPI

client = pyAMI.client.Client('atlas')

AtlasAPI.init()

all_ops = ["FM0","FM1","FM2","FM3","FM4","FM5","FM7",
            "FS02","FS1",
            "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]
all_ops2 = ["FM","FS","FT"]

Processes = ["WmZ"]
Decay = ["llqq"]

List_info={}
for op in all_ops2:
    for process in Processes:
        for decay in Decay:
            List_ami = []
            #print(f'\n%MGPy8EG_aQGC{op}_QUAD_1_{process}_{decay}%')
            pattern = f'%MGPy8EG_aQGC{op}%_QUAD_1_{process}_{decay}%'
            List_ami = AtlasAPI.list_datasets(client, patterns=[pattern], 
                                              fields=['ldn', 'cross_section', 'dataset_number'], 
                                              limit=[1, 100], type='EVNT', project='mc16%')

            if List_ami:
                #print(List_ami)
                List_info[f'{op}_{process}_{decay}'] = List_ami
            else:
                print(f'No datasets found for {op}_{process}_{decay}\n')
                


organized_info = {}


for key, datasets in List_info.items():
    for dataset in datasets:
        # Extract operator from the dataset name
        operator = next((op for op in all_ops if op in dataset['ldn']), None)
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

print(organized_info['FM1_WmZ_llqq']['dataset_number'][0])
print(organized_info['FM1_WmZ_llqq']['cross_section'][0])