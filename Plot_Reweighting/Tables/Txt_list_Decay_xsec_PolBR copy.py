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

#all_ops_cat = ["FM0","FM1","FM2"]

#all_ops_cat = ["FM0","FM1"]  # Reduced set for testing
all_ops_cat_rwg=["FM","FS","FT"]
if opts.Rwg:
    all_ops_cat = ["FM","FS","FT"]

processes = ["WpZ"]
decays = ['llqq']

Branching_ratios = {
    "WpZ_llqq": 0.069,  
}

EFT_ORDER = opts.EFT_order
print(f"[INFO] Using EFT order: {EFT_ORDER}")
proc_decays_tuple = list(itertools.product(processes, decays))

cross_terms_ = [f"{op1}vs{op2}" for op1, op2 in itertools.combinations(all_ops_cat, 2) if op1[:2] == op2[:2]]

if EFT_ORDER == "CROSS":
    all_ops_cat = cross_terms_
if "both" in EFT_ORDER.lower():
    all_ops_cat = all_ops_cat + cross_terms_

Type_MC_all= ["EFTDec_Madspin", "EFTDec_polarisation","Reweighthel_ignore_polarisation", "Reweighthel_aware_polarisation"]
Type_MC_all= ["EFTDec_Madspin", "EFTDec_polarisation","Rwg_pol_50k","Rwg_pol_100k", "Rwg_InvSqrtXsec_50k", "Rwg_InvXsec_50k"]

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

def xsec_decay(xsec, branching_ratio):
    xsec_ = pb_to_fb(xsec)  
    if branching_ratio <= 0:
        raise ValueError("Branching ratio must be greater than zero.")
    return xsec_ * branching_ratio

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
        elif "Reweighthel_ignore" in type_MC or "Rwg_pol_50k" in type_MC:
            base_dir = f"{base_path}/Reweighting/Polarisation/hel_ignore/{prod_dec}/"
        elif "Reweighthel_aware" in type_MC or "Rwg_hel_aware_polarisation" in type_MC:
            base_dir = f"{base_path}/Reweighting/Polarisation/hel_aware/{prod_dec}/"
            
        elif "Reweighting_100k" in type_MC or "Rwg_pol_100k" in type_MC:
            base_dir = f"{base_path}/Reweighting/Polarisation/hel_ign_100k/{prod_dec}/"
        elif "Reweighting_InvSqrtXsec_50k" in type_MC or "Rwg_InvSqrtXsec_50k" in type_MC:
            base_dir = f"{base_path}/Reweighting/Polarisation/InvSqrtXsec/{prod_dec}/"
        elif "Reweighting_InvXsec_50k" in type_MC or "Rwg_InvXsec_50k" in type_MC:
            base_dir = f"{base_path}/Reweighting/Polarisation/InvXsec/{prod_dec}/"
            
            
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
        br = Branching_ratios.get(f"{process}_{decay}", 1.0)  # Default to 1.0 if not found
        for op in all_ops_cat:
            op_cross = op
            if "madspin" in type_MC.lower():
                if "vs" in op:
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
                    xsec_val_fb = xsec_decay(xsec_val,br)
                    xsec_unc_fb = xsec_decay(xsec_unc,br)
                    results.append({
                        "Type_MC": type_MC,
                        "Polarisation": None,
                        "Operator": op,
                        "xsec_fb": xsec_val_fb,
                        "uncertainty_fb": xsec_unc_fb,
                    })
            elif "rwg" in type_MC.lower() or "reweight" in type_MC.lower():
                
                op_rwg = op[:2]  # e.g., FM for FM0, FM1, etc.
                if "FS" in op_rwg or "FT" in op_rwg and "inv" in type_MC.lower() and "xsec" in type_MC.lower():
                    conf = f"user.osalin.MadGraph_{process}_{decay}_{op_rwg}_{order}"
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
                        op_= op.replace("vs", "_")  # Replace 'vs' with '_'
                        if "vs" in op:
                            Avoid_string_cross = "CROSS_"
                        else:
                            Avoid_string_cross = "CROSS"
                        if operator.startswith(op_)  and "QUAD_" not in operator and Avoid_string_cross not in operator:
                            
                            operator = op
                            results.append({
                                "Type_MC": type_MC,
                                "Polarisation": pol,
                                "Operator": operator,
                                "xsec_fb": xsec_decay(xsec_val,br),
                                "uncertainty_fb": xsec_decay(xsec_unc,br),
                                "log_path": log_path
                            })
                    
                else:
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
                            op_= op.replace("vs", "_")  # Replace 'vs' with '_'
                            if "vs" in op:
                                Avoid_string_cross = "CROSS_"
                            else:
                                Avoid_string_cross = "CROSS"
                            if operator.startswith(op_)  and "QUAD_" not in operator and Avoid_string_cross not in operator:
                                
                                operator = op
                                results.append({
                                    "Type_MC": type_MC,
                                    "Polarisation": pol,
                                    "Operator": operator,
                                    "xsec_fb": xsec_decay(xsec_val,br),
                                    "uncertainty_fb": xsec_decay(xsec_unc,br),
                                    "log_path": log_path
                                })

            elif "eftdec" in type_MC.lower():
                for pol in polarisations:
                    if "vs" in op:
                        conf = f"user.osalin.MadGraph_{process}_{decay}_{op}_CROSS_{pol}"
                    else:
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
                            "xsec_fb": xsec_decay(xsec_val,br),
                            "uncertainty_fb": xsec_decay(xsec_unc,br),
                            "log_path": log_path
                        })

if results:
    df = pd.DataFrame(results)
    if 'log_path' in df.columns:
        df = df.drop(columns=['log_path'])
    df["Operator"] = df["Operator"].str.replace(r"(_QUAD|_CROSS).*", "", regex=True)
    df["xsec_fb"] = pd.to_numeric(df["xsec_fb"], errors="coerce")
    pol_states = ["LL", "LT", "TL", "TT"]
    os.makedirs("Tables/Txt", exist_ok=True)
    # Determine which EFT orders to process
    if "both" in opts.EFT_order.lower():
        print("[INFO] Processing both EFT orders: QUAD and CROSS")
        eft_orders = ["QUAD", "CROSS"]
    else:
        eft_orders = [opts.EFT_order]
    for type_mc in df["Type_MC"].unique():
        txt_lines = []
        for order_eft in eft_orders:
            df_order = df.copy()
            if order_eft == "CROSS":
                df_order = df_order[df_order["Operator"].str.contains("vs")]
            else:
                df_order = df_order[~df_order["Operator"].str.contains("vs")]
            df_mc = df_order[df_order["Type_MC"] == type_mc]
            for op in sorted(df_mc["Operator"].unique()):
                df_op = df_mc[df_mc["Operator"] == op]
                xsec_alone = float(df_op["xsec_fb"].sum())
                df_pol = df_op[df_op["Polarisation"].isin(pol_states)]
                xsec_sum = float(df_pol["xsec_fb"].sum())
                process = "WpZ"
                decay = "llqq"
                if 'polarisation' in type_mc.lower():
                    line_all = f"{process}_{decay}_{op}_{order_eft}: {xsec_sum}"
                else:
                    line_all = f"{process}_{decay}_{op}_{order_eft}: {xsec_alone}"
                txt_lines.append(line_all)
                for _, row in df_pol.iterrows():
                    pol = row.get("Polarisation", None)
                    xsec = float(row.get("xsec_fb", 0))
                    pol_str = f"_{pol}" if pd.notnull(pol) else ""
                    line = f"{process}_{decay}_{op}_{order_eft}{pol_str}: {xsec}"
                    txt_lines.append(line)
        outdir = f"./Tables/Uncertainty/Xsec/test/"
        os.makedirs(outdir, exist_ok=True)
        lower_type_mc = type_mc.lower()
        out_path = f"{outdir}/VBS_cross_section_run2_{lower_type_mc}.txt"
        with open(out_path, "w") as f:
            f.write("\n".join(txt_lines) + "\n")
        print(f"[INFO] TXT table saved: {out_path}")
# --- END: Automated polarisation scan and Excel export ---