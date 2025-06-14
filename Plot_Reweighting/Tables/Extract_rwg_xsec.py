import re
import pandas as pd
import os
import csv  
import numpy as np

# Filepath base
filepath_base = '/eos/home-o/osalin/VBS_MC/mcjoboptions/Output/Model_restricted2/'

# Dictionary with nicknames and log paths
Reweight_specific = {
    #"Param_1e12": 'Reweight_M_xsec30k/LOG/Reweight_M_5k_1e12.txt',
    #"No_process": 'Rwg_nochange_process/LOG/Reweight_M_10k_nochangePro_noheli.txt',
    #"Param_sqrte12": 'Reweight_M_xsec3/LOG/Reweight_M_10ksqrt.txt',
    #"Param_invsqrte12": 'Reweight_M_xsec30k/LOG/Reweight_M_5k_parem1e12.txt',
    #"Param_invsqrt": 'Reweight_M_xsec30k/LOG/Reweight_auto_5k.txt',
    #"Helicity_1e12": 'Reweight_M_xsec3/LOG/Reweight_M_10k_1par_helicity.txt',
    #"Helicity_invsqrte12": 'Rwg_M_xsec_param/LOG/Reweight_M_10k_invsqrt.txt',
    #"Helicity_sqrte12": 'Rwg_M_xsec_param/LOG/Reweight_M_10k_helicty.txt',
    #"Nb_evnts_1k": 'Rwg_M_xsec_param/LOG/Reweight_M_1k_param_helicity.txt',
    #"Nb_evnts_10k": 'Rwg_M_xsec_param/LOG/Reweight_M_10k_invsqrt.txt',
    #"Nb_evnts_100k": 'Rwg_M_xsec_param/LOG/Reweight_M_helicity_100k.txt',
    #"TM_helicity_param1": 'Rwg_T_M_xsec/LOG/Reweight_TM_10k_hel_parm1_works.txt',
    #"TM_helicity_InvXsec": 'Rwg_TM_xsec/LOG/Reweight_TM_hel_10k_invxsec_works.txt',
    #"S_M_hel_param1": 'Rwg_S_M_param1/LOG/Reweight_SM_10k_.txt',
    #"S_M_hel_InvXsec": 'Rwg_S_M_xsec/LOG/Reweight_SM_10k_invsec2.txt',
    "S_T_hel_InvXsec": "Rwg_S_T_xsec/LOG/Reweight_ST_10k_invsec.txt",
}

# Filepaths for cross-section data
xsec_fb_filepath = os.path.join(filepath_base, 'Tables/Cross_section/VBS_xsection_fb_aqgcModel_new_op.txt')
xsec_pb_filepath = os.path.join(filepath_base, 'Tables/Cross_section/VBS_xsection_bef_decay_pb.txt')

# Function to read cross-section data from a file
def read_xsec_file(filepath, unit_conversion=1e-3):
    xsec_data = {}
    with open(filepath, 'r') as file:
        for line in file:
            if ':' in line:
                operator, value = line.split(':')
                value = value.strip()
                operator = operator.split('_')[0]
                if '+-' in value:
                    xsec_value, xsec_uncertainty = value.split(' +- ')
                    xsec_value = float(xsec_value) * unit_conversion
                    xsec_uncertainty = float(xsec_uncertainty) * unit_conversion
                    xsec_data[operator.strip()] = (xsec_value, xsec_uncertainty)
                else:
                    xsec_value = float(value) * unit_conversion
                    xsec_data[operator.strip()] = xsec_value
    return xsec_data

# Read cross-section data
xsec_fb_data = read_xsec_file(xsec_fb_filepath)
xsec_pb_data = read_xsec_file(xsec_pb_filepath, unit_conversion=1)
print(xsec_pb_data)

for nickname, log_path in Reweight_specific.items():
    # Read the text file
    filepath = os.path.join(filepath_base, log_path)
    with open(filepath, 'r') as file:
        content = file.read()

    # Regular expression to extract the cross-section information
    pattern_xsec = re.compile(r'INFO: (\w+_\w+) : ([\d\.\-e]+) \+- ([\d\.\-e]+) pb')
    # Regular expression to extract the parameter values
    pattern_param = re.compile(r'INFO\s+(\d+)\s+([\d\.\-e]+)\s+#\s+(\w+)')

    # Find all matches
    matches_xsec = pattern_xsec.findall(content)
    matches_param = pattern_param.findall(content)

    # Create a dictionary for parameter values
    param_dict = {param: value for _, value, param in matches_param}

    Op_list=[]
    for match in matches_xsec:
        operator, _, _ = match
        if '_CROSS' in operator:
            operator = operator.replace('_CROSS', '').replace('_', 'vs')
        else:
            operator = operator.strip('_QUAD')
        Op_list.append(operator)
    # Process and format the values
    
    
    data = []
    for match in matches_xsec:
        operator, xsec_value, xsec_uncertainty = match
        print(f'Operator: {operator}, Xsec Value: {xsec_value}, Xsec Uncertainty: {xsec_uncertainty}')
        if '_CROSS' in operator:
            operator = operator.replace('_CROSS', '').replace('_', 'vs').replace('set', '')
            param_value = ''
        else:
            operator = operator.strip('_QUAD')
            param_value = "{:.2e}".format(float(param_dict.get(operator, '')))
        xsec_value = "{:.2e}".format(float(xsec_value))
        xsec_uncertainty = "{:.2e}".format(float(xsec_uncertainty))
        
        # Get cross-section values from the additional files
        xsec_fb = xsec_fb_data.get(operator)
        xsec_pb = xsec_pb_data.get(operator)
        #print(f'Operator: {operator}, Xsec FB: {xsec_fb}, Xsec PB: {xsec_pb}')

        data.append({
            'Operator': operator,
            'Param original': param_value,
            #'σ nominal': "{:.2e}".format(xsec_fb) if xsec_fb is not None else '',
            'σ bef decay': "{:.2e}".format(xsec_pb[0]) if xsec_pb[0] is not None else '',
            'Unc. σ bef decay': "{:.2e}".format(xsec_pb[1]) if xsec_pb[1] is not None else '',
            'σ Rwg ': xsec_value,
            'Unc. σ Rwg': xsec_uncertainty,
            'Rel Rwg Unc. %': "{:.2f}%".format(np.abs(float(xsec_uncertainty) / float(xsec_value)) * 100),
            'Rel diff in % σ rwg': "{:.2f}%".format(np.abs((float(xsec_value) - xsec_pb[0]) / xsec_pb[0]) * 100) if xsec_pb[0] is not None else '',
        })

    # Convert the list of dictionaries to a pandas DataFrame
    df = pd.DataFrame(data)

    path_tables = os.path.join(filepath_base, 'Tables')
    path_tables_spe = os.path.join(filepath_base, 'Tables/Excel/', nickname)
    os.makedirs(path_tables_spe, exist_ok=True)

    # Save the DataFrame to an Excel file
    output_excel_filepath = os.path.join(path_tables_spe, f'{nickname}_output.xlsx')
    df.to_excel(output_excel_filepath, index=False)

    # Save the DataFrame to a CSV file
    output_csv_filepath = os.path.join(path_tables_spe, f'{nickname}_output.csv')
    df.to_csv(output_csv_filepath, index=False)

    print(f'Data successfully extracted and saved to {output_excel_filepath} and {output_csv_filepath}')