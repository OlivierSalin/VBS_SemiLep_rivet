import re
import pandas as pd
import os
import csv  
import numpy as np
import itertools
import glob
from pathlib import Path


from optparse import OptionParser
parser = OptionParser()

parser.add_option("--nb_lep", default = 2)
parser.add_option("--name", default = "")
parser.add_option("--All_channel", default= True)
parser.add_option("--aQGC", default= True)
parser.add_option("--EFT_order", default = "QUAD")
parser.add_option('--type_MC', default = "")
parser.add_option("--Channel", default = "")
parser.add_option("--decays", default = "llqq")
parser.add_option("--Conf", default = "user.osalin")
parser.add_option("--Rwg", default = False)
parser.add_option("--Pol", default="")
opts, _ = parser.parse_args()

order="QUAD"

all_ops_cat = ["FM0","FM1","FM2","FM3","FM4","FM5","FM7","FM8","FM9",
            "FS0","FS1","FS2",
            "FT0","FT1","FT2","FT3","FT4","FT5","FT6"]

#all_ops_cat = ["FM0","FM1"]  # Reduced set for testing
all_ops_cat_rwg=["FM","FS","FT"]
if opts.Rwg:
    all_ops_cat = ["FM","FS","FT"]

processes = ["WpZ"]
decays = ['llqq']
proc_decays_tuple = list(itertools.product(processes, decays))

cross_terms_ = [f"{op1}vs{op2}" for op1, op2 in itertools.combinations(all_ops_cat, 2)]

if opts.EFT_order == "CROSS":
    all_ops_cat = cross_terms_

Type_MC_all= ["EFTDec_Madspin", "EFTDec_polarisation","Rwg_hel_ignore_polarisation", "Rwg_hel_aware_polarisation"]
Type_MC_all= [ "EFTDec_Madspin","Rwg_hel_ignore_polarisation", "Rwg_hel_aware_polarisation"]

# Filepath base
filepath_base = '/exp/atlas/salin/ATLAS/VBS_mc/plotting/Plot_Reweighting/Tables/'

# Dictionary with nicknames and log paths
Reweight_specific_path = {
    "M_operator": 'Polarisation_xsec/LOG/Reweight_TM_10k_hel_parm1_works.txt',
    "S_operator": 'Rwg_TM_xsec/LOG/Reweight_TM_hel_10k_invxsec_works.txt',
    'T operator': 'Rwg_T_M_xsec/LOG/Reweight_TM_10k_hel_parm1_works.txt',
}
# Filepaths for cross-section data


def get_xsec(log_file):
    with open(log_file) as textf:
        xsec_val, xsec_unit = -999.0 , "fb" # here pb but later for plots will convert to fb
        for line in textf:
            if 'MetaData: cross-section' in line:
                xsec_val = float(line[line.find('=')+1:])
                xsec_unit = line[line.find('(')+1:line.find(')')]
    conv_fact_to_pb = {"mb":1e9, "um":1e6, "nb":1e3, "pb":1, "fb":1e-3}
    xsec_fb = xsec_val * conv_fact_to_pb[xsec_unit] * 1000
    print("found xsec value ",xsec_val,"with unit",xsec_unit,"converting to fb get in fb",xsec_fb)
    return xsec_fb

def pb_to_fb(xsec_pb):
    """
    Convert cross-section from picobarns (pb) to femtobarns (fb).
    1 pb = 1000 fb
    """
    return float(xsec_pb) * 1000

def logfile_find(conf, type_MC=None):
    base_path = "/exp/atlas/salin/ATLAS/VBS_mc/eft_files"
    #print("New models")
    
    def extract_prod_dec(conf):
        prod_temp = conf[conf.find("user.osalin.MadGraph_") + len("user.osalin.MadGraph_"):]
        #prod_temp = conf[conf.find("user.osalin.Madgraph_") + len("user.osalin.MadGraph_"):]
        #print("start from string", prod_temp)
        prod_dec = prod_temp[:prod_temp.find("qq_") + 2]
        #print("from conf found production dec", prod_dec)
        return prod_dec

    if conf.startswith("user."):
        prod_dec = extract_prod_dec(conf)
            
        if "Reweighting_Madspin" in type_MC:
            base_dir = f"{base_path}/Reweighting/Madspin/{prod_dec}/"
        elif "Reweighting_Polarisation" in type_MC or "Rwg_polarisation" in type_MC:
            base_dir = f"{base_path}/Reweighting/Polarisation/{prod_dec}/"
        elif "Reweighthel_ignore" in type_MC or "Rwg_hel_ignore_polarisation" in type_MC:
            base_dir = f"{base_path}/Reweighting/Polarisation/hel_ignore/{prod_dec}/"
        elif "Reweighthel_aware" in type_MC or "Rwg_hel_aware_polarisation" in type_MC:
            base_dir = f"{base_path}/Reweighting/Polarisation/hel_aware/{prod_dec}/"
            
            
        elif "EFTDec_Madspin" in type_MC:
            base_dir = f"{base_path}/EFTDec/Madspin//{prod_dec}/"
        elif "EFTDec_Polarisation" in type_MC or "EFTDec_polarisation" in type_MC:
            base_dir = f"{base_path}/EFTDec/Polarisation/{prod_dec}/"

        else:
            raise ValueError("Unknown type_MC: ", type_MC)
        
        #print("base_dir: ", base_dir)
        search_com= base_dir + f"/*{conf}*EXT0"
        #print("searching for dir with pattern", search_com)
        conf_dir_arr = glob.glob(search_com)
        #print("found possibilities for dir", conf_dir_arr)
        conf_dir = conf_dir_arr[0] if len(conf_dir_arr)>=1 else -1  
        if conf_dir == -1: raise ValueError("did not find folder for this config ",search_com)
    
    return prod_dec, conf_dir


LOG_files= {}
polarisations = ["LL", "LT", "TL", "TT"]
#polarisations = ["LL"]




# --- BEGIN: Automated polarisation scan and Excel export ---
results = []
pattern_xsec_rwg = re.compile(r'INFO: (\w+_\w+) : ([\d.\-e]+) \+- ([\d.\-e]+) pb')
pattern_xsec_eftdec = re.compile(r'Current estimate of cross-section: ([\d.eE\-]+) \+- ([\d.eE\-]+)')

for type_MC in Type_MC_all:
    for process, decay in proc_decays_tuple:
        for op in all_ops_cat:
            op_cross = op
            if "madspin" in type_MC.lower():
                if opts.EFT_order == "CROSS":
                    conf = f"user.osalin.MadGraph_{process}_{decay}_{op}_CROSS"
                else:
                    conf = f"user.osalin.MadGraph_{process}_{decay}_{op}_{order}"
                try:
                    _, conf_dir = logfile_find(conf, type_MC)
                except Exception as e:
                    print(f"[WARN] Could not find dir for {conf} / {type_MC}: {e}")
                    continue
                log_path = os.path.join(conf_dir, "log.generate")
                if not os.path.exists(log_path):
                    print(f"[WARN] Log file does not exist: {log_path}")
                    continue
                with open(log_path, "r") as f:
                    content = f.read()
                match = pattern_xsec_eftdec.search(content)
                #print(f"{op} Searching for {conf} in {type_MC}")
                #print(f"[INFO] Found match for {conf} in {type_MC}")
                #print("Match found: ", match)
                if match:
                    xsec_val, xsec_unc = match.groups()
                    xsec_val_fb = pb_to_fb(xsec_val)
                    xsec_unc_fb = pb_to_fb(xsec_unc)
                    results.append({
                        "Type_MC": type_MC,
                        "Polarisation": None,
                        "Operator": op,
                        "xsec_fb": xsec_val_fb,
                        "uncertainty_fb": xsec_unc_fb,
                    })
            elif "rwg" in type_MC.lower() or "reweight" in type_MC.lower():
                op_rwg = op[:2]  # e.g., FM for FM0, FM1, etc.
                for pol in polarisations:
                    conf = f"user.osalin.MadGraph_{process}_{decay}_{op_rwg}_{order}_{pol}"
                    try:
                        _, conf_dir = logfile_find(conf, type_MC)
                    except Exception as e:
                        print(f"[WARN] Could not find dir for {conf} / {type_MC}: {e}")
                        continue
                    log_path = os.path.join(conf_dir, "log.generate")
                    if not os.path.exists(log_path):
                        print(f"[WARN] Log file does not exist: {log_path}")
                        continue
                    with open(log_path, "r") as f:
                        content = f.read()
                    matches = pattern_xsec_rwg.findall(content)
                    for match in matches:
                        operator, xsec_val, xsec_unc = match
                        # Only keep the operator that matches the current op (e.g., FM0)
                        #if operator.startswith(op) and "CROSS" not in operator and "QUAD_" not in operator:
                        op_= op.replace("vs", "_")  # Replace 'vs' with '_'
                        if operator.startswith(op_)  and "QUAD_" not in operator and "CROSS_" not in operator:
                            operator = op
                            results.append({
                                "Type_MC": type_MC,
                                "Polarisation": pol,
                                "Operator": operator,
                                "xsec_fb": pb_to_fb(xsec_val),
                                "uncertainty_fb": pb_to_fb(xsec_unc),
                                "log_path": log_path
                            })
            elif "eftdec" in type_MC.lower():
                for pol in polarisations:
                    conf = f"user.osalin.MadGraph_{process}_{decay}_{op}_{order}_{pol}"
                    try:
                        _, conf_dir = logfile_find(conf, type_MC)
                    except Exception as e:
                        print(f"[WARN] Could not find dir for {conf} / {type_MC}: {e}")
                        continue
                    log_path = os.path.join(conf_dir, "log.generate")
                    if not os.path.exists(log_path):
                        print(f"[WARN] Log file does not exist: {log_path}")
                        continue
                    with open(log_path, "r") as f:
                        content = f.read()
                    match = pattern_xsec_eftdec.search(content)

                    if match:
                        xsec_val, xsec_unc = match.groups()
                        results.append({
                            "Type_MC": type_MC,
                            "Polarisation": pol,
                            "Operator": op,
                            "xsec_fb": pb_to_fb(xsec_val),
                            "uncertainty_fb": pb_to_fb(xsec_unc),
                            "log_path": log_path
                        })

if results:
    df = pd.DataFrame(results)
    # Remove log_path column if present
    if 'log_path' in df.columns:
        df = df.drop(columns=['log_path'])
    # Add 'Mixte' polarisation as the sum of LT and TL for each group (by Type_MC, Operator, etc.)
# After df has been built in the script
    df["Operator"] = df["Operator"].str.replace(r"(_QUAD|_CROSS).*", "", regex=True)
    df["xsec_fb"] = pd.to_numeric(df["xsec_fb"], errors="coerce")
    df["uncertainty_fb"] = pd.to_numeric(df["uncertainty_fb"], errors="coerce")

    # Création du tableau structuré
    polar_states = ["LL", "LT", "TL", "TT"]
    type_mc_order = ["EFTDec_polarisation", "Rwg_hel_ignore_polarisation", "Rwg_hel_aware_polarisation"]
    df_all = df[(df["Type_MC"] == "EFTDec_Madspin") & (df["Polarisation"].isna())].copy()
    df_all = df_all.rename(columns={"xsec_fb": "xsec_all", "uncertainty_fb": "unc_all"})[["Operator", "xsec_all", "unc_all"]]

    rows = []
    for op in sorted(df["Operator"].unique()):
        all_xsec_val = df_all[df_all["Operator"] == op]["xsec_all"].values
        unc_all_val = df_all[df_all["Operator"] == op]["unc_all"].values
        all_xsec = float(all_xsec_val[0]) if len(all_xsec_val) else np.nan
        unc_all = float(unc_all_val[0]) if len(unc_all_val) else np.nan
        for mc in type_mc_order:
            # For type_MC with 'polarisation' in name, sum over all pol states
            if "rwg" in mc.lower():
                pol_matches = df[(df["Operator"] == op) & (df["Type_MC"] == mc) & (df["Polarisation"].isin(polar_states))]
                pol_xsecs = pol_matches["xsec_fb"].values
                pol_uncs = pol_matches["uncertainty_fb"].values
                all_xsec_pol = np.nansum(pol_xsecs) if len(pol_xsecs) else ""
                unc_all_pol = np.sqrt(np.nansum(np.square(pol_uncs))) if len(pol_uncs) else ""
                row = [op, mc, all_xsec_pol, unc_all_pol]
            else:
                row = [op, mc, all_xsec if mc == "EFTDec_polarisation" else "", unc_all if mc == "EFTDec_polarisation" else ""]
            for pol in polar_states:
                match = df[(df["Operator"] == op) & (df["Type_MC"] == mc) & (df["Polarisation"] == pol)]
                if not match.empty:
                    xsec = float(match["xsec_fb"].values[0])
                    unc = float(match["uncertainty_fb"].values[0])
                    # Use the correct all_xsec for ratio
                    if "polarisation" in mc.lower():
                        ratio = abs(round(xsec / all_xsec_pol, 4)) if all_xsec_pol else ""
                    else:
                        ratio = abs(round(xsec / all_xsec, 4)) if all_xsec else ""
                else:
                    xsec, unc, ratio = "", "", ""
                row.extend([xsec, unc, ratio])
            rows.append(row)

    columns = ["Operator", "Type_MC", "All_xsec", "All_unc"]
    for pol in polar_states:
        columns += [f"{pol}_xsec", f"{pol}_unc", f"{pol}_RatioPolAll"]



    df = pd.DataFrame(rows, columns=columns)  # <-- replace raw df with structured df

    float_cols = [col for col in df.columns if "_xsec" in col or "_unc" in col]
    for col in float_cols:
        df[col] = df[col].apply(lambda x: f"{x:.2e}" if isinstance(x, float) else x)

    # Export as final result
    df.to_excel(f"Tables/Excel/MC_Polarisation_Comparison_Structured_{opts.name}.xlsx", index=False)
    print(f"[INFO] Structured table saved: MC_Polarisation_Comparison_Structured_{opts.name}.xlsx")
# --- END: Automated polarisation scan and Excel export ---