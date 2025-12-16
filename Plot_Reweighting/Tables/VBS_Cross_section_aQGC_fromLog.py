#!/usr/bin/env python3
import re
import pandas as pd
import os
import numpy as np
import itertools
import glob

from optparse import OptionParser

# --- CLI (preserve original options) ---
parser = OptionParser()
parser.add_option("--nb_lep", default=2)
parser.add_option("--name", default="")
parser.add_option("--All_channel", default=True)
parser.add_option("--aQGC", default=True)
parser.add_option("--EFT_order", default="QUAD")
parser.add_option("--Mode", default="")
parser.add_option('--type_MC', default="")
parser.add_option("--Channel", default="")
parser.add_option("--decays", default="llqq")
parser.add_option("--Conf", default="user.osalin")
parser.add_option("--Rwg", default=False)
parser.add_option("--Pol", default="")
# NEW: add CLI option to control merged-charge aggregation (string "True"/"False")
parser.add_option("--Merged_charge", default="False")
opts, _ = parser.parse_args()

# --- Constants / configuration ---
order = "QUAD"
# Full operator list (kept as original)

type_MC_all_bis = [
    "Reweighting_Madspin",
    "Reweighting_Polarisation",
    "Reweighting_CPodd",
    "Reweighting_PolCPodd",
    "EFTDec_Madspin",
    "EFTDec_Polarisation",
    "EFTDec_CPodd"
]

if opts.Mode == "Lep":
    type_MC_all_bis= ["Reweighting_VBSLep_pol","eftdec_vbslep"]
elif opts.Mode == "SemiLep":
    type_MC_all_bis= ["reweighting_vbssemilep_pol"]
else:
    type_MC_all_bis= ["Reweighting_VBSLep_pol"]
outdir_path= f"./Cross_section/Validation//{opts.name}/"

all_ops_cat = [
    "FM0","FM1","FM2","FM3","FM4","FM5","FM7","FM8","FM9",
    "FS0","FS1","FS2",
    "FT0","FT1","FT2","FT3","FT4","FT5","FT6",
    "FM1odd","FM2odd","FM3odd","FM4odd","FM5odd","FM6odd",
    "FT1odd","FT2odd","FT3odd","FT4odd","FT5odd","FT6odd"
]

all_ops_cat = [
    "FM1odd","FM2odd","FM3odd","FM4odd","FM5odd","FM6odd"
]

all_ops_cat = [
    "FM0","FM1","FM2","FM3","FM4","FM5","FM7","FM8","FM9",
    "FS0","FS1","FS2",
    "FT0","FT1","FT2","FT3","FT4","FT5","FT6","FT8","FT9",
    "FM1odd","FM2odd","FM3odd","FM4odd","FM5odd","FM6odd",
    "FT1odd","FT2odd","FT3odd","FT4odd","FT5odd","FT6odd"
]
#all_ops_cat = ["FM0"]
all_ops_cat_rwg = ["FM","FS","FT","FModd","FTodd"]
all_ops_cat_rwg = ["FModd"]
#all_ops_cat_rwg = ["FM","FS","FT"]
if opts.Rwg:
    all_ops_cat = all_ops_cat_rwg[:]  # preserve behavior when Rwg flag is set

processes = ["WpZ"]
decays = ['llqq']

if opts.Mode == "Lep":
	print("Running in Leptonic mode")
	valid_combinations = {
		4: [("ZZ", "llll")],
		3: [("WmZ", "lllv"), ("WpZ", "lllv")],
		2: [("WmWm", "lvlv"), ("WpWm", "lvlv"), ("WpWp", "lvlv"),("Zy", "lly")],
		1: [("Wpy", "lvy")],
	}
	valid_combinations = {
		#4: [("ZZ", "llll")],
		#3: [("WmZ", "lllv"), ("WpZ", "lllv")],
		2: [("WpWm", "lvlv")],
		#1: [("Wpy", "lvy")],
	}



elif opts.Mode == "SemiLep":
	print("Running in Semi-Leptonic mode")
	valid_combinations = {
		1: [("WmZ", "lvqq"), ("WpZ", "lvqq"), ("WmWm", "lvqq"), ("WpWm", "lvqq"), ("WpWp", "lvqq")],
		2: [("WpZ", "llqq"), ("WmZ", "llqq"), ("ZZ", "llqq")],
		0: [("WpZ", "vvqq"), ("WmZ", "vvqq"), ("ZZ", "vvqq")]
	}

	valid_combinations = {
		1: [("WmZ", "lvqq"), ("WpZ", "lvqq"), ("WmWm", "lvqq"), ("WpWm", "lvqq"), ("WpWp", "lvqq")],
		2: [("WpZ", "llqq"), ("WmZ", "llqq"), ("ZZ", "llqq")],
		0: [("WpZ", "vvqq"), ("WmZ", "vvqq"), ("ZZ", "vvqq")]
	}

else:
	valid_combinations = {
		2: [("WpZ", "vvqq",), ("WmZ", "vvqq")],
	}


def diboson_br(process: str, decay: str):
    """
    process: 'WpZ','WmZ','WpWp','WpWm','WmWm','ZZ','Zy','Wpy' (charge doesn’t matter for BR)
    decay:   'llll','lllv','lvlv','lvqq','llqq','vvqq','lly','lvy'
    lepton_set: 'emu' or 'all' (all = e, mu, tau)
    returns: branching ratio as float
    """
    # PDG 2024/25 values
    W_e, W_mu, W_tau = 0.1071, 0.1063, 0.1138      # W → eν, μν, τν
    W_had = 1.0 - 3.0*0.1086                       # from universal fit → 0.6742
    Z_e, Z_mu, Z_tau = 0.03363, 0.03366, 0.03369   # Z → ee, μμ, ττ
    Z_had, Z_inv = 0.6991, 0.2000                  # Z → hadrons, Z → νν (sum)


    W_l = W_e + W_mu + W_tau
    Z_ll = Z_e + Z_mu + Z_tau


    # Map final-state tags to products of single-boson BRs
    if decay == "llll":   # ZZ
        return Z_ll * Z_ll
    if decay == "lllv":   # WZ
        return Z_ll * W_l
    if decay == "lvlv":   # WW
        return W_l * W_l
    if decay == "lvqq":   # WZ or WW semileptonic
        if "Z" in process:
            return W_l * Z_had
        else:
            return W_l * W_had
    if decay == "llqq":   # WZ or ZZ semileptonic with dilepton
        if "Z" in process and process.count("Z") == 1:
            return Z_ll * W_had
        else:
            return Z_ll * Z_had
    if decay == "vvqq":   # ZZ → νν qq
        return Z_inv * Z_had
    if decay == "lly":    # Zγ with Z → ℓℓ
        return Z_ll
    if decay == "lvy":    # Wγ with W → ℓν
        return W_l
    raise ValueError("Unknown decay tag")


valid_combi = [f"{proc}_{dec}" for pairs in valid_combinations.values() for proc, dec in pairs]
print(f"[INFO] Valid process+decay combinations: {valid_combi}")

Branching_ratios = {}
for key in valid_combinations:
    for proc, dec in valid_combinations[key]:
        br = diboson_br(proc, dec)
        Branching_ratios[f"{proc}_{dec}"] = br
print(f"[INFO] Branching ratios: {Branching_ratios}")

EFT_ORDER = opts.EFT_order
print(f"[INFO] Using EFT order: {EFT_ORDER}")
proc_decays_tuple = valid_combi

# cross terms: same-family pairs as original behavior
cross_terms_ = [f"{op1}vs{op2}" for op1, op2 in itertools.combinations(all_ops_cat, 2) if op1[:2] == op2[:2] or ("odd" in op1.lower() and "odd" in op2.lower())]
if EFT_ORDER == "CROSS":
    all_ops_cat = cross_terms_
if "both" in EFT_ORDER.lower():
    all_ops_cat = all_ops_cat + cross_terms_

# base path for EFT files (kept)
BASE_EFT_PATH = "/data/atlas/salin/VBS_mc/VBS/eft_Files/"

# regex patterns (kept logic similar to original)
PATTERN_XSEC_RWG = re.compile(r'INFO: (\w+_\w+) : ([\d.\-e]+) \+- ([\d.\-e]+) pb')
PATTERN_XSEC_EFTDEC = re.compile(r'Current estimate of cross-section: ([\d.eE\-]+) \+- ([\d.eE\-]+)')

# polarisation states (original code set polarisations = [] — preserve)
polarisations = ["LL", "LT", "TL", "TT"]
pol_states = ["LL", "LT", "TL", "TT"]
pol_states = ["LL", "LT", "TL", "TT"]

# --- Helper functions ---
def safe_float(val):
    """Convert numeric string to float safely; return NaN on failure."""
    try:
        return float(val)
    except Exception:
        return float("nan")

def pb_to_fb(xsec_pb):
    """Convert picobarns to femtobarns (1 pb = 1000 fb). Accepts numeric or numeric string."""
    return safe_float(xsec_pb) * 1000.0

def xsec_decay(xsec_pb_like, branching_ratio):
    """Convert cross-section in pb-like value to fb and multiply by branching ratio.
    Accepts numeric or string values for xsec_pb_like.
    """
    xsec_fb = pb_to_fb(xsec_pb_like)
    if branching_ratio <= 0:
        raise ValueError("Branching ratio must be greater than zero.")
    return xsec_fb * branching_ratio

def extract_prod_dec(conf):
    """Extract the production+decay substring from configuration name.
    Mirrors previous logic: finds substring after 'user.osalin.MadGraph_' and takes up to 'qq_'.
    """
    marker = "user.osalin.MadGraph_"
    if marker not in conf:
        raise ValueError("Unexpected conf format, cannot extract prod_dec")
    prod_temp = conf[conf.find(marker) + len(marker):]
    parts = prod_temp.split("_")
    if len(parts) >= 2:
        prod_dec = f"{parts[0]}_{parts[1]}"
    else:
        prod_dec = parts[0]
    return prod_dec

def map_typeMC_to_base_dir(type_MC, prod_dec):
    """Map the provided type_MC string to the corresponding directory under BASE_EFT_PATH.
    This centralizes mapping rules, preserving the original branching logic.
    """
    t = type_MC or ""
    type_MC_low = t.lower()
    if "reweighting_madspin" in type_MC_low:
        return os.path.join(BASE_EFT_PATH, "Reweighting", "Madspin", "hel_ignore", prod_dec)

    if "reweighting_cpodd" in type_MC_low:
        return os.path.join(BASE_EFT_PATH, "Reweighting", "CPodd", "Validation", prod_dec)
    if "reweighting_polcpodd" in type_MC_low or "rwg_OddPol" in type_MC_low:
        return os.path.join(BASE_EFT_PATH, "Reweighting", "CPodd", "Validation","Polarisation", prod_dec) 

    if "reweighting_vbslep_pol" in type_MC_low or "reweight_vbs_lep_pol" in type_MC_low:
        return os.path.join(BASE_EFT_PATH, "Reweighting", "VBSAll", "Pol","Gen","VBS","Leptonic", prod_dec)

    if "reweighting_vbssemilep_pol" in type_MC_low or "reweight_vbs_semilep_pol" in type_MC_low:
        return os.path.join(BASE_EFT_PATH, "Reweighting", "VBSAll", "Pol","Gen","VBS","SemiLep", prod_dec)

    if "reweighting_vbssemilep_validation" in type_MC_low or "reweight_vbs_semilep_validation" in type_MC_low:
        return os.path.join(BASE_EFT_PATH, "Reweighting", "Validation", "VBS_SemiLep","20k", prod_dec)

    if "eftdec_madspin" in type_MC_low:
        return os.path.join(BASE_EFT_PATH, "EFTDec", "Madspin", prod_dec)
    if "eftdec_polarisation" in type_MC_low or "eftdec_polarisation" in type_MC_low:
        return os.path.join(BASE_EFT_PATH, "EFTDec", "Polarisation", prod_dec)
    if "eftdec_cpodd" in type_MC_low:
        return os.path.join(BASE_EFT_PATH, "EFTDec", "CPodd", "Validation", prod_dec)

    if "eftdec_vbslep" in type_MC_low:
        return os.path.join(BASE_EFT_PATH, "EFTDec", "VBSAll", "Pol", "Leptonic", prod_dec)
    # fallback: unknown
    raise ValueError("Unknown type_MC: " + str(type_MC))

def logfile_find(conf, type_MC=None):
    """Find the directory for a given configuration name and type_MC.
    Returns (prod_dec, conf_dir) or raises ValueError if not found.
    """
    if not conf.startswith("user."):
        raise ValueError("Expected conf to start with 'user.'")

    prod_dec = extract_prod_dec(conf)
    base_dir = map_typeMC_to_base_dir(type_MC, prod_dec)
    # search the EXT0 pattern as in original code
    search_pattern = os.path.join(base_dir, f"*{conf}*EXT0")
    matches = glob.glob(search_pattern)
    if not matches:
        raise ValueError(f"did not find folder for this config {search_pattern}")
    conf_dir = matches[0]
    return prod_dec, conf_dir

def make_reweight_opname(op):
    """Map operator string to reweight operator name used in filenames/configs.
    This function encapsulates the logic previously spread through the code.
    Rules preserved from original script:
      - For a single operator: use first two chars plus 'odd' if original contains 'odd'
      - For a 'vs' pair: consider prefixes and 'odd' flags and build combined token
    """
    if "vs" not in op:
        return f"{op[:2]}{'odd' if 'odd' in op.lower() else ''}"
    left, right = op.split("vs", 1)
    left_pref, right_pref = left[:2], right[:2]
    left_odd, right_odd = ("odd" in left.lower()), ("odd" in right.lower())
    if left_pref == right_pref:
        # same family: if no odd -> single prefix, otherwise keep both prefixes and add single 'odd'
        if not left_odd and not right_odd:
            return left_pref
        if left_odd and right_odd:
            return f"{left_pref}odd"
        return f"{left_pref}{right_pref}{'odd' if (left_odd or right_odd) else ''}"
    # different families: preserve prefix and odd flags for each
    return f"{left_pref}{'odd' if left_odd else ''}{right_pref}{'odd' if right_odd else ''}"

# --- Main processing ---
def process_all():
    results = []
    # final MC types (preserve last assignment from the original script)
    Type_MC_all = type_MC_all_bis
    #Type_MC_all = ["Reweighting_Madspin"]

    for type_MC in Type_MC_all:
        for process_decay in proc_decays_tuple:
            process, decay = process_decay.split("_", 1)
            br = Branching_ratios.get(f"{process}_{decay}", 1.0)
            for op in all_ops_cat:
                # EFTDec-like processing (original checks: if "madspin" in type_MC.lower() or "eftdec_cpodd" in type_MC.lower())
                type_MC_low = type_MC.lower()
                try:
                    if ("eftdec" in type_MC_low and "pol" not in type_MC_low) or "eftdec_madspin" in type_MC_low or "eftdec_cpodd" in type_MC_low:
                        # construct conf as original
                        if "vs" in op:
                            conf = f"user.osalin.MadGraph_{process}_{decay}_{op}_CROSS"
                        else:
                            conf = f"user.osalin.MadGraph_{process}_{decay}_{op}_{order}"
                        _, conf_dir = logfile_find(conf, type_MC)
                        log_path = os.path.join(conf_dir, "log.generate")
                        if not os.path.exists(log_path):
                            print(f"[WARN] Log file does not exist: {log_path}")
                            continue
                        with open(log_path, "r") as fh:
                            content = fh.read()
                        match = PATTERN_XSEC_EFTDEC.search(content)
                        if match:
                            xsec_val, xsec_unc = match.groups()
                            results.append({
                                "Type_MC": type_MC,
                                "Process_decay": process_decay,
                                "Polarisation": None,
                                "Operator": op,
                                "xsec_fb": xsec_decay(xsec_val, br),
                                "uncertainty_fb": xsec_decay(xsec_unc, br),
                            })
                    # Reweighting block (original code looked for reweighting or rwg text)
                    elif "reweight" in type_MC_low or "rwg" in type_MC_low:
                        op_rwg = make_reweight_opname(op)
                        # debug line (preserved)
                        print(f"[DEBUG] op='{op}' -> op_rwg='{op_rwg}'")
                        # When type_MC mentions CPodd reweighting (preserve original branch)
                        if "reweighting_cpodd" in type_MC_low or "reweighting_madspin" in type_MC_low:
                            conf = f"user.osalin.MadGraph_{process}_{decay}_{op_rwg}_{order}"
                            _, conf_dir = logfile_find(conf, type_MC)
                            log_path = os.path.join(conf_dir, "log.generate")
                            if not os.path.exists(log_path):
                                print(f"[WARN] Log file does not exist: {log_path}")
                                continue
                            with open(log_path, "r") as fh:
                                content = fh.read()
                            matches = PATTERN_XSEC_RWG.findall(content)
                            for operator, xsec_val, xsec_unc in matches:
                                op_ = op.replace("vs", "_")
                                Avoid_string_cross = "CROSS_" if "vs" in op else "CROSS"
                                # apply same filtering as original
                                if operator.startswith(op_) and "QUAD_" not in operator and Avoid_string_cross not in operator:
                                    results.append({
                                        "Type_MC": type_MC,
                                        "Process_decay": process_decay,
                                        "Polarisation": None,
                                        "Operator": op,
                                        "xsec_fb": xsec_decay(xsec_val, br),
                                        "uncertainty_fb": xsec_decay(xsec_unc, br),
                                        "log_path": log_path
                                    })
                        # the original script also had a branch for polarisation-containing reweight types;
                        # preserve that logic, although polarisations list is empty by default (keeps behavior)
                        elif "polarisation" in type_MC_low or "pol" in type_MC_low:
                            for pol in polarisations:
                                print(f"[DEBUG] Processing pol='{pol}' for op='{op}'")
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
                                with open(log_path, "r") as fh:
                                    content = fh.read()
                                matches = PATTERN_XSEC_RWG.findall(content)
                                for operator, xsec_val, xsec_unc in matches:
                                    op_ = op.replace("vs", "_")
                                    Avoid_string_cross = "CROSS_" if "vs" in op else "CROSS"
                                    if operator.startswith(op_) and "QUAD_" not in operator and Avoid_string_cross not in operator and "INT" not in operator:
                                        results.append({
                                            "Type_MC": type_MC,
                                            "Process_decay": process_decay,
                                            "Polarisation": pol,
                                            "Operator": op,
                                            "xsec_fb": xsec_decay(xsec_val, br),
                                            "uncertainty_fb": xsec_decay(xsec_unc, br),
                                            "log_path": log_path
                                        })
                    # EFTDec with polarisation (preserve original branch)
                    elif "eftdec" in type_MC_low and "pol" in type_MC_low:
                        for pol in polarisations:
                            if "vs" in op:
                                conf = f"user.osalin.MadGraph_{process}_{decay}_{op}_CROSS_{pol}"
                            else:
                                conf = f"user.osalin.MadGraph_{process}_{decay}_{op}_{order}_{pol}"
                            _, conf_dir = logfile_find(conf, type_MC)
                            log_path = os.path.join(conf_dir, "log.generate")
                            if not os.path.exists(log_path):
                                print(f"[WARN] Log file does not exist: {log_path}")
                                continue
                            with open(log_path, "r") as fh:
                                content = fh.read()
                            match = PATTERN_XSEC_EFTDEC.search(content)
                            if match:
                                xsec_val, xsec_unc = match.groups()
                                results.append({
                                    "Type_MC": type_MC,
                                    "Process_decay": process_decay,
                                    "Polarisation": pol,
                                    "Operator": op,
                                    "xsec_fb": xsec_decay(xsec_val, br),
                                    "uncertainty_fb": xsec_decay(xsec_unc, br),
                                    "log_path": log_path
                                })
                except Exception as e:
                    # preserve original behavior of warning and continue on failures
                    print(f"[WARN] Could not process op={op} / {type_MC} / {process}_{decay}: {e}")
                    continue

    # --- Post-process results and write outputs (same lists & names as original) ---
    if not results:
        print("[INFO] No results collected.")
        return

    df = pd.DataFrame(results)
    if 'log_path' in df.columns:
        df = df.drop(columns=['log_path'])
    # remove QUAD/CROSS suffixes if they ended up in operator names
    df["Operator"] = df["Operator"].str.replace(r"(_QUAD|_CROSS).*", "", regex=True)
    print(f"[INFO] Collected {len(df)} cross-section entries.")
    #print(df.head())
    df["xsec_fb"] = pd.to_numeric(df["xsec_fb"], errors="coerce")

    # --- NEW: merged-charge aggregation ---
    merge_charge_flag = str(opts.Merged_charge).lower() in ("true", "1", "yes")
    if merge_charge_flag and not df.empty:
        # robustly split Process_decay into proc and decay (use keyword n and expand)
        temp = df['Process_decay'].str.split('_', n=1, expand=True)
        # ensure we always have two columns
        if temp.shape[1] == 1:
            temp[1] = ""
        temp.columns = ['__proc_tmp', '__decay_tmp']
        df = pd.concat([df, temp], axis=1)
        # mapping of merged name -> list of charged process prefixes to sum
        merge_map = {
            "WZ": ["WpZ", "WmZ"],
            "ssWW": ["WmWm", "WpWp"],
            "osWW": ["WpWm"]
        }
        new_rows = []
        # iterate per Type_MC, decay, operator, polarisation
        for type_mc in df['Type_MC'].unique():
            df_tm = df[df['Type_MC'] == type_mc]
            for decay in df_tm['__decay_tmp'].unique():
                for merged_name, parts in merge_map.items():
                    df_match = df_tm[(df_tm['__proc_tmp'].isin(parts)) & (df_tm['__decay_tmp'] == decay)]
                    if df_match.empty:
                        continue
                    # group by Operator and Polarisation (treat NaN/None polarisation as group)
                    group_cols = ['Operator', 'Polarisation']
                    grouped = df_match.groupby(group_cols, dropna=False)
                    for (op, pol), grp in grouped:
                        xsec_sum = grp['xsec_fb'].sum()
                        # uncertainties combined in quadrature
                        unc_vals = grp['uncertainty_fb'].astype(float).values
                        unc_quad = float(np.sqrt(np.nansum(np.square(unc_vals)))) if unc_vals.size else 0.0
                        new_rows.append({
                            "Type_MC": type_mc,
                            "Process_decay": f"{merged_name}_{decay}",
                            "Polarisation": pol if pd.notnull(pol) else None,
                            "Operator": op,
                            "xsec_fb": xsec_sum,
                            "uncertainty_fb": unc_quad
                        })
        if new_rows:
            df = pd.concat([df, pd.DataFrame(new_rows)], ignore_index=True, sort=False)
            # update valid_combi so outputs include merged names
            for merged_name in merge_map.keys():
                for decay in df['__decay_tmp'].unique():
                    merged_entry = f"{merged_name}_{decay}"
                    if merged_entry not in valid_combi:
                        valid_combi.append(merged_entry)
        # cleanup temp columns if present
        for c in ('__proc_tmp', '__decay_tmp'):
            if c in df.columns:
                df = df.drop(columns=[c])
    # --- END NEW ---

    # create output directories and file names exactly as original code
    outdir = outdir_path
    os.makedirs(outdir, exist_ok=True)

    # determine EFT orders to write (preserve original both/both-like behavior)
    if "both" in opts.EFT_order.lower():
        eft_orders = ["QUAD", "CROSS"]
    else:
        eft_orders = [opts.EFT_order]

    for type_mc in df["Type_MC"].unique():
        txt_lines = []
        for process_decay in valid_combi:
            process, decay = process_decay.split("_", 1)

            for order_eft in eft_orders:
                df_order = df.copy()
                if order_eft == "CROSS":
                    df_order = df_order[df_order["Operator"].str.contains("vs")]
                else:
                    df_order = df_order[~df_order["Operator"].str.contains("vs")]
                df_mc = df_order[(df_order["Type_MC"] == type_mc) & (df_order["Process_decay"] == process_decay)]
                for op in sorted(df_mc["Operator"].unique()):
                    df_op = df_mc[df_mc["Operator"] == op]
                    xsec_alone = float(df_op["xsec_fb"].sum())
                    df_pol = df_op[df_op["Polarisation"].isin(pol_states)]
                    xsec_sum = float(df_pol["xsec_fb"].sum()) if not df_pol.empty else 0.0
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

        lower_type_mc = type_mc.lower()
        out_path = f"{outdir}/VBS_cross_section_run2_{lower_type_mc}_{opts.Mode}.txt"
        with open(out_path, "w") as fh:
            fh.write("\n".join(txt_lines) + "\n")
        print(f"[INFO] TXT table saved: {out_path}")

    for type_mc in df["Type_MC"].unique():
        for process_decay in valid_combi:           # uncertainties
            unc_lines = []
            for order_eft in eft_orders:
                df_order = df.copy()
                if order_eft == "CROSS":
                    df_order = df_order[df_order["Operator"].str.contains("vs")]
                else:
                    df_order = df_order[~df_order["Operator"].str.contains("vs")]
                df_mc = df_order[(df_order["Type_MC"] == type_mc) & (df_order["Process_decay"] == process_decay)]
                for op in sorted(df_mc["Operator"].unique()):
                    df_op = df_mc[df_mc["Operator"] == op]
                    unc_alone = float(df_op["uncertainty_fb"].sum())
                    df_pol = df_op[df_op["Polarisation"].isin(pol_states)]
                    pol_uncs = df_pol["uncertainty_fb"].values.astype(float) if not df_pol.empty else np.array([])
                    unc_all_pol = np.sqrt(np.nansum(np.square(pol_uncs))) if pol_uncs.size else ""
                    if 'polarisation' in type_mc.lower():
                        line_all = f"{process}_{decay}_{op}_{order_eft}: {unc_all_pol}"
                    else:
                        line_all = f"{process}_{decay}_{op}_{order_eft}: {unc_alone}"
                    unc_lines.append(line_all)
                    for _, row in df_pol.iterrows():
                        pol = row.get("Polarisation", None)
                        unc = float(row.get("uncertainty_fb", 0))
                        pol_str = f"_{pol}" if pd.notnull(pol) else ""
                        line = f"{process}_{decay}_{op}_{order_eft}{pol_str}: {unc}"
                        unc_lines.append(line)

        unc_out_path = f"{outdir}/VBS_Uncertainty_cross_section_run2_{lower_type_mc}_{opts.Mode}.txt"
        with open(unc_out_path, "w") as fh:
            fh.write("\n".join(unc_lines) + "\n")
        print(f"[INFO] Uncertainty TXT table saved: {unc_out_path}")

if __name__ == "__main__":
    process_all()