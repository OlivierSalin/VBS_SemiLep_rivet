import pyAMI.client
import pyAMI.atlas.api as AtlasAPI
from lib_utils import find_prod_dec_and_dir
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--conf", default = "user.osalin.MadGraph_WmZ_llqq_FM1_QUAD")
opts, _ = parser.parse_args()

client = pyAMI.client.Client('atlas')

AtlasAPI.init()

prod_dec, base_dir = find_prod_dec_and_dir(opts.conf)
print(f'prod_dec: {prod_dec}')
def extract_EFT_op_proces_dec(conf):
    """
    This function takes conf as a parameter,
    and returns the EFT operator, process, and decay.
    """
    # Split the conf name with the underscore 
    parts = conf.split('_')
    EFT_op,EFT_type, proc, decay = parts[3], parts[4], parts[1], parts[2]
    
    return EFT_op, EFT_type, proc, decay

EFT_op, EFT_type, proc, decay = extract_EFT_op_proces_dec(opts.conf)
print(f'EFT_op: {EFT_op}, EFT type:{EFT_type}, proc: {proc}, decay: {decay}')

def cross_section_fb(EFT_op, EFT_type,  proces, dec):
    """
    This function takes an EFT operator, a process, and a decay as parameters,
    and returns the cross section for the corresponding dataset.
    """
    all_ops = ["FM0","FM1","FM2","FM3","FM4","FM5","FM7",
               "FS02","FS1",
               "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]
    all_ops2 = ["FM","FS","FT"]
    Processes = [proces]
    Decay = [dec]

    organized_info = {}

    if EFT_type == 'QUAD':
        for op in all_ops2:
            for process in Processes:
                for decay in Decay:
                    pattern = f'%MGPy8EG_aQGC{op}%_QUAD_1_{process}_{decay}%'
                    List_ami = AtlasAPI.list_datasets(client, patterns=[pattern], 
                                                    fields=['ldn', 'cross_section', 'dataset_number'], 
                                                    limit=[1, 100], type='EVNT', project='mc16%')

                    if List_ami:
                        for dataset in List_ami:
                            operator = next((op for op in all_ops if op in dataset['ldn']), None)
                            if operator:
                                key = f"{operator}_{EFT_type}_{process}_{decay}"
                                if key not in organized_info:
                                    organized_info[key] = {'name': [], 'cross_section': [], 'dataset_number': []}

                                organized_info[key]['name'].append(dataset['ldn'])
                                organized_info[key]['cross_section'].append(dataset['cross_section']) # Cross section in nb
                                organized_info[key]['dataset_number'].append(dataset['dataset_number'])
                    else:
                        print(f'No datasets found for {op}_{process}_{decay}\n')
                        
    elif (EFT_type == 'SM'):
        for process in Processes:
            for decay in Decay:
                if process == 'WmZ' or process == 'WpZ':
                    process_ = 'WZ'
                pattern = f'%WZjj_llqq_EW6%'
                List_ami = AtlasAPI.list_datasets(client, patterns=[pattern], 
                                                fields=['ldn', 'cross_section', 'dataset_number'], 
                                                limit=[1, 100], type='EVNT', project='mc16%')
                print(List_ami)
                if List_ami:
                    for dataset in List_ami:
                        key = f"{EFT_op}_{EFT_type}_{process}_{decay}"
                        if key not in organized_info:
                            organized_info[key] = {'name': [], 'cross_section': [], 'dataset_number': []}

                        organized_info[key]['name'].append(dataset['ldn'])
                        organized_info[key]['cross_section'].append(dataset['cross_section'])

    # Construct the key from the EFT operator, process, and decay
    key_ = f"{EFT_op}_{EFT_type}_{proces}_{dec}"

    # Check if the key exists in the organized_info dictionary
    if key_ in organized_info:
        # If the key exists, return the first cross section
        xsection_nb = organized_info[key_]['cross_section'][0] # Return the cross section in nb
        xsection_fb = float(xsection_nb) * 1e6 # Convert the cross section to fb
        return xsection_fb
    else:
        # If the key does not exist, print an error message and return None
        print(organized_info[key_]['cross_section'][0])
        print(f'Error: No datasets found for {key_}')
        return None

keyy = f"{EFT_op}_{EFT_type}_{proc}_{decay}"    
xsection_fb = cross_section_fb(EFT_op,EFT_type, proc, decay)
print(f'Cross section in fb for {keyy}: {xsection_fb}')