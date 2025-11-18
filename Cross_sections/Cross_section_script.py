import pyAMI.client
import pyAMI.atlas.api as AtlasAPI
import lib_utils as lu
from optparse import OptionParser
import utils_func as uf
parser = OptionParser()
parser.add_option("--conf", default = "user.osalin.MadGraph_WmZ_llqq_FM5_QUAD")
opts, _ = parser.parse_args()

client = pyAMI.client.Client('atlas')

AtlasAPI.init()

prod_dec, base_dir = lu.find_prod_dec_and_dir(opts.conf)
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


all_ops = ["FM0","FM1","FM2","FM3","FM4","FM5","FM7",
           "FS02","FS1",
           "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]
all_ops_ = [op + "_QUAD" for op in all_ops]
Processes = ["WmZ","WpZ","ZZ","WmWm","WpWm","WmWp","WpWp"]
Decay = ["llqq",'lvqq','vvqq'] 


def cross_section_fb(EFT_op, EFT_type,  proces, dec):
    """
    This function takes an EFT operator, a process, and a decay as parameters,
    and returns the cross section for the corresponding dataset.
    """
    all_ops = ["FM0","FM1","FM2","FM3","FM4","FM5","FM7",
               "FS02","FS1",
               "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]
    all_ops2 = ["FM","FS","FT"]
    all_ops_ = [op + "_QUAD" for op in all_ops]
    Processes = [proces]
    Decay = [dec]

    organized_info = {}

    if EFT_type == 'QUAD':
        for process in Processes:
            for decay in Decay:  
                pattern = f"%%MGPy8EG_aQGCF%%_{process}_{decay}%%"
                List_ami = AtlasAPI.list_datasets(client, patterns=[pattern], 
                                                fields=['ldn', 'cross_section', 'dataset_number'], 
                                                limit=[1, 100], type='EVNT', project='mc16%')

                if List_ami:
                    for dataset in List_ami:
                        operator = next((op for op in all_ops_ if op in dataset['ldn']), None)
                        if operator:
                            key = f"{operator}_{process}_{decay}"
                            if key not in organized_info:
                                organized_info[key] = {'name': [], 'cross_section': [], 'dataset_number': []}

                            organized_info[key]['name'].append(dataset['ldn'])
                            organized_info[key]['cross_section'].append(dataset['cross_section']) # Cross section in nb
                            organized_info[key]['dataset_number'].append(dataset['dataset_number'])
                else:
                    print(f'No datasets found for {process}_{decay}\n')
                        
        #print(organized_info)                
    elif (EFT_type == 'SM'):
        for process in Processes:
            for decay in Decay:
                if process == 'WmZ' or process == 'WpZ':
                    process_ = 'WZ'
                if process == "ZZ":
                    process_= "ZZ"
                if process == "WmWm" or process == "WpWm" or process == "WmWp" or process == "WpWp":
                    process_ = "WW"
                pattern = f'%{process_}jj_{decay}_EW6%'
                List_ami = AtlasAPI.list_datasets(client, patterns=[pattern], 
                                                fields=['ldn', 'cross_section', 'dataset_number'], 
                                                limit=[1, 100], type='EVNT', project='mc16%')
                #print(List_ami)
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
        #print(organized_info[key_]['cross_section'][0])
        print(f'Error: No datasets found for {key_}')
        return None


all_ops = ["FM0","FM1","FM2","FM3","FM4","FM5","FM7",
           "FS02","FS1",
           "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]
all_ops_ = [op + "_QUAD" for op in all_ops]
all_ops_.append("FM0_SM")
Processes = ["WmZ","WpZ","ZZ","WmWm","WpWm","WpWp"]
Decay = ["llqq",'lvqq','vvqq']
 



VBS_xsection = {}

path_log="/exp/atlas/salin/ATLAS/VBS_mc/eft_files/Test/aqgc_model/WpZ_llqq/user.osalin.MadGraph_WpZ_llqq_FM2_QUAD_aqgcModel_new_all_EXT0/user.osalin.MadGraph_WpZ_llqq_FM2_QUAD_aqgcModel_new_1.log/tarball_PandaJob_6437928208_FZK-LCG2/log.generate"


x_sec,xsec_unc= lu.get_xsec_bef_decay(path_log)
print(f'Xsec before decay: {x_sec} +- {xsec_unc}')

# Loop over all combinations of operators, processes, and decays
""" for EFT_op_ in all_ops_:
    EFT_op, EFT_type = EFT_op_.split('_')
    for proc in Processes:
        for decay in Decay:
            # Create the key
            if uf.possible_process(proc, decay):
                keyy = f"{proc}_{decay}_{EFT_op}_{EFT_type}"    
                xsection_fb = cross_section_fb(EFT_op,EFT_type, proc, decay)
                if xsection_fb is not None:
                    # Store the cross section in the dictionary only if it's not None
                    VBS_xsection[keyy] = xsection_fb
                    print(f'Cross section in fb for {keyy}: {xsection_fb}')
                
# Write the dictionary to a text file
with open('VBS_xsection_test.txt', 'w') as f:
    for key, value in VBS_xsection.items():
        f.write(f'{key}: {value}\n')
        
with open('VBS_valid_ops.txt', 'w') as f:
    for key in VBS_xsection.keys():
        f.write(f'{key}\n') """