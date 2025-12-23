import glob
import os
import re
import ROOT
import pandas as pd
import numpy as np

ROOT.gROOT.SetBatch(True)

def parse_options():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--Mode", default="Lep")          # Lep | SemiLep
    parser.add_option("--EFT_order", default="QUAD")    # QUAD
    parser.add_option("--outdir", default="./acceptance/")
    parser.add_option(
        "--xsec_file",
        default="/exp/atlas/salin/ATLAS/VBS_mc/VBS_Pol_Rivet/VBS_rivet/Plot_Reweighting/Tables/Cross_section/Validation/Validation/VBS_cross_section_run2_reweighting_vbs_ALL.txt",
        help="Path to VBS_cross_section_*.txt"
    )
    parser.add_option("--verbose", default=True, help="Verbose prints (True/False)")
    parser.add_option(
        "--pol_base_dir",
        default="/exp/atlas/salin/ATLAS/VBS_mc/eft_files/Histograms/Polarisation/Leptonic/Rwg/GOOD/Pol/Leptonic",
        help="Base dir that contains folders like <process>_<decay>/<OPGROUP>_QUAD[_POL]/**/ntuple_rivet.root"
    )
    parser.add_option("--write_xlsx", default=True, help="Also write an .xlsx next to the .csv (True/False)")
    opts, _ = parser.parse_args()
    return opts

def _as_bool(x):
    if isinstance(x, bool):
        return x
    return str(x).strip().lower() in ("1", "true", "yes", "y", "on")

def _log(verbose, msg):
    if verbose:
        print(msg)

def filter_operators_by_group(all_ops, op_group):
    g = (op_group or "").strip()
    if g == "FM":
        return [op for op in all_ops if op.startswith("FM") and not op.endswith("odd")]
    if g == "FS":
        return [op for op in all_ops if op.startswith("FS")]
    if g == "FT":
        return [op for op in all_ops if op.startswith("FT") and not op.endswith("odd")]
    if g == "FModd":
        return [op for op in all_ops if op.startswith("FM") and op.endswith("odd")]
    if g == "FTodd":
        return [op for op in all_ops if op.startswith("FT") and op.endswith("odd")]
    return list(all_ops)

def take_xsec_fb(VBS_txt, op, order_EFT, process, decay, pol=None):
    # your exact key-matching function
    if pol in ["LL", "LT", "TL", "TT"]:
        key = f"{process}_{decay}_{op}_{order_EFT}_{pol}"
    else:
        key = f"{process}_{decay}_{op}_{order_EFT}"
    try:
        with open(VBS_txt, "r") as f:
            for line in f:
                info_xsec = line.strip().split(": ")
                if len(info_xsec) > 1:
                    key_file, xsection_fb = info_xsec[0], info_xsec[1]
                    if key_file == key:
                        return float(xsection_fb)
    except OSError:
        return None
    return None

def _sum_branch_rdf(root_path, tree_name, branch_name):
    """
    Fast sum using RDataFrame (C++ side). Returns 0.0 on missing tree/branch.
    """
    if not os.path.exists(root_path):
        return 0.0
    f = ROOT.TFile.Open(root_path, "READ")
    if not f or f.IsZombie():
        return 0.0
    t = f.Get(tree_name)
    if not t:
        f.Close()
        return 0.0
    if not t.GetBranch(branch_name):
        f.Close()
        return 0.0
    f.Close()
    try:
        df = ROOT.RDataFrame(tree_name, root_path)
        return float(df.Sum(branch_name).GetValue())
    except Exception:
        return 0.0

def _folder_for(op_rwg, order_eft, pol):
    # structure you described: <op_rwg>_QUAD and <op_rwg>_QUAD_LL, etc.
    if pol in ("LL", "LT", "TL", "TT"):
        return f"{op_rwg}_{order_eft}_{pol}"
    return f"{op_rwg}_{order_eft}"

def _after_tree_for_mode(mode):
    return "SR" if mode == "Lep" else "Merged"


def _acceptance(root_path, before_tree, after_tree, weight_branch):
    sum_before = _sum_branch_rdf(root_path, before_tree, weight_branch)
    sum_after = _sum_branch_rdf(root_path, after_tree, weight_branch)
    if sum_before == 0.0:
        return None, sum_before, sum_after
    return (sum_after / sum_before), sum_before, sum_after


def _format_fiducial_cell(xsec_fid_fb, acc):
    if xsec_fid_fb is None:
        return ""
    if acc is None:
        return f"{xsec_fid_fb:.6g}"
    return f"{xsec_fid_fb:.6g} ({acc*100:.2f}%)"

def main():
    opts = parse_options()
    verbose = _as_bool(getattr(opts, "verbose", False))
    write_xlsx = _as_bool(getattr(opts, "write_xlsx", False))

    _log(verbose, f"[CFG] Mode={opts.Mode} EFT_order={opts.EFT_order}")
    _log(verbose, f"[CFG] pol_base_dir={opts.pol_base_dir} exists={os.path.isdir(opts.pol_base_dir)}")
    _log(verbose, f"[CFG] xsec_file={opts.xsec_file} exists={os.path.exists(opts.xsec_file)}")
    _log(verbose, f"[CFG] outdir={opts.outdir}")

    if opts.Mode == "Lep":
        valid_combinations = {
            4: [("ZZ", "llll")],
            3: [("WZ", "lllv")],
            2: [("ssWW", "lvlv"), ("osWW", "lvlv"), ("Zy", "lly")],
        }

    elif opts.Mode == "SemiLep":
        valid_combinations = {
            1: [("WZ", "lvqq"), ("ssWW", "lvqq"), ("osWW", "lvqq")],
            2: [("WZ", "llqq"), ("ZZ", "llqq")],
            0: [("WZ", "vvqq"), ("ZZ", "vvqq")],
        }

    else:
        valid_combinations = {2: [("WpZ", "vvqq"), ("WmZ", "vvqq")]}

    proc_decay_pairs = [pd for pairs in valid_combinations.values() for pd in pairs]
    _log(verbose, f"[CFG] proc_decay_pairs={proc_decay_pairs}")

    all_operators_all = [
        "FM0","FM1","FM2","FM3","FM4","FM5","FM7","FM8","FM9",
        "FS0","FS1","FS2",
        "FT0","FT1","FT2","FT3","FT4","FT5","FT6","FT7",
        "FM1odd","FM2odd","FM3odd","FM4odd","FM5odd","FM6odd",
        "FT1odd","FT2odd","FT3odd","FT4odd","FT5odd","FT6odd"
    ]

    #all_operators_all = ["FM0","FS0","FT0","FM3"]
    #all_operators_all = ["FM1odd","FM2odd","FM3odd","FM4odd","FM5odd","FM6odd",]
    

    operators_rwg = ["FM", "FS", "FT", "FModd", "FTodd"]
    pol_states = ["LL", "LT", "TL", "TT", "unpolarised"]
    before_tree = "Bef_cut"
    after_tree = _after_tree_for_mode(opts.Mode)

    if opts.Mode == "Lep":
        pol_directory = "/exp/atlas/salin/ATLAS/VBS_mc/eft_files/Histograms/Polarisation/Leptonic/Rwg/GOOD/Pol/Leptonic"
    elif opts.Mode == "SemiLep":
        pol_directory = "/exp/atlas/salin/ATLAS/VBS_mc/eft_files/Histograms/Polarisation/Leptonic/Rwg/GOOD/Pol/SemiLep/"

    rows = []
    for proc_decay in proc_decay_pairs:
        process, decay = proc_decay
        for op_rwg in operators_rwg:
            ops_this = filter_operators_by_group(all_operators_all, op_rwg)
            _log(verbose, f"\n[GRP] op_rwg={op_rwg} -> n_ops={len(ops_this)} ops={ops_this}")

            for pol in pol_states:
                op_folder = _folder_for(op_rwg, opts.EFT_order, pol)
                
                proc_dir = os.path.join(pol_directory, f"{process}_{decay}")
                base_dir = os.path.join(proc_dir, op_folder)
                _log(verbose, f"[POL] pol={pol} folder={op_folder}")
                _log(verbose, f"[POL] base_dir={base_dir} exists={os.path.isdir(base_dir)}")
                if not os.path.isdir(base_dir):
                    continue


                # accept any depth under process_decay that contains the ntuple
                ntuples = glob.glob(os.path.join(base_dir, "**", "ntuple_rivet.root"), recursive=True)
                _log(verbose, f"[PROC] found ntuples={len(ntuples)}")
                if not ntuples:
                    continue

                for ntuple_path in ntuples:
                    # compute per-operator acceptances + xsecs
                    for op in ops_this:
                        wbr = f"EventWeight_{op.lower()}_{opts.EFT_order.lower()}"
                        acc, sum_before, sum_after = _acceptance(ntuple_path, before_tree, after_tree, wbr)

                        pol_key = pol if pol in ("LL", "LT", "TL", "TT") else None
                        print(f"Taking xsec for process={process} decay={decay} op={op} pol={pol_key}")
                        print(f"Acceptance: {acc} sum_before={sum_before} sum_after={sum_after} wbr={wbr}")
                        xsec_pol = take_xsec_fb(opts.xsec_file, op, opts.EFT_order, process, decay, pol=pol_key)
                        xsec_comb = take_xsec_fb(opts.xsec_file, op, opts.EFT_order, process, decay, pol=None)

                        if acc is None:
                            _log(verbose, f"[WARN] missing acceptance op={op} pol={pol} proc={process} decay={decay} branch={wbr}")
                            # still keep xsec-only rows

                        xsec_fid_pol = (xsec_pol * acc) if (xsec_pol is not None and acc is not None) else None
                        xsec_fid_comb = (xsec_comb * acc) if (xsec_comb is not None and acc is not None) else None

                        rows.append({
                            "mode": opts.Mode,
                            "process_decay": f"{process}_{decay}",
                            "process": process,
                            "decay": decay,
                            "operator": op,
                            "order": opts.EFT_order,
                            "op_rwg": op_rwg,
                            "pol_state": pol,
                            "ntuple_path": ntuple_path,
                            "weight_branch": wbr,
                            "sumW_before": sum_before,
                            "sumW_after": sum_after,
                            "acceptance": acc,
                            "xsec_pol_fb": xsec_pol,
                            "xsec_comb_fb": xsec_comb,
                            "xsec_fid_pol_fb": xsec_fid_pol,
                            "xsec_fid_comb_fb": xsec_fid_comb,
                        })

    df = pd.DataFrame(rows)
    if df.empty:
        _log(True, "[WARN] no rows found (no ntuples or no matches)")
        return

    df = df.sort_values(["process_decay", "operator", "pol_state", "ntuple_path"]).reset_index(drop=True)

    # If multiple ntuples are found per (process_decay, operator, pol_state), sum weights are per-ntuple.
    # For acceptances, take a sum-weighted average across ntuples: acc = sum_after/sum_before.
    grp_cols = ["mode", "process_decay", "process", "decay", "operator", "order", "pol_state"]
    agg = df.groupby(grp_cols, dropna=False).agg(
        sumW_before=("sumW_before", "sum"),
        sumW_after=("sumW_after", "sum"),
        xsec_pol_fb=("xsec_pol_fb", "first"),
        xsec_comb_fb=("xsec_comb_fb", "first"),
    ).reset_index()
    agg["acceptance"] = agg.apply(
        lambda r: (r["sumW_after"] / r["sumW_before"]) if (r["sumW_before"] and r["sumW_before"] != 0.0) else None,
        axis=1,
    )
    agg["xsec_fid_pol_fb"] = agg.apply(
        lambda r: (r["xsec_pol_fb"] * r["acceptance"]) if (r["xsec_pol_fb"] is not None and r["acceptance"] is not None) else None,
        axis=1,
    )
    agg["xsec_fid_comb_fb"] = agg.apply(
        lambda r: (r["xsec_comb_fb"] * r["acceptance"]) if (r["xsec_comb_fb"] is not None and r["acceptance"] is not None) else None,
        axis=1,
    )

    pol_order = ["unpolarised", "LL", "LT", "TL", "TT"]
    agg["pol_state"] = pd.Categorical(agg["pol_state"], categories=pol_order, ordered=True)

    # -------------------------
    # Table A: inclusive xsec + pol fraction (in %)
    # Output columns (flat):
    #   unpolarised_xsec_fb, LL_xsec_fb, ..., and LL_pol_frac, ...
    # where pol_frac is relative to unpolarised in %.
    # -------------------------
    xsec_piv = agg.pivot_table(
        index=["process_decay", "operator"],
        columns="pol_state",
        values="xsec_pol_fb",
        aggfunc="first",
    )

    # fractions vs unpolarised (exclude unpolarised itself)
    if "unpolarised" in xsec_piv.columns:
        base = xsec_piv["unpolarised"]
        frac_piv = (xsec_piv.div(base, axis=0) * 100.0).round(3)
    else:
        frac_piv = xsec_piv.copy() * 0.0

    # build final flat table
    table_xsec = pd.DataFrame(index=xsec_piv.index)
    table_xsec["unpolarised_xsec_fb"] = xsec_piv.get("unpolarised")
    for pol in ["LL", "LT", "TL", "TT"]:
        table_xsec[f"{pol}_pol_frac"] = frac_piv.get(pol)
    table_xsec = table_xsec.reset_index()

    # -------------------------
    # Table B: fiducial (ATLAS fiducial volume)
    # Requirements:
    #  - remove fiducial fraction for unpolarised
    #  - provide unpolarised fiducial xsec and acc_unpolarised
    #  - for each polarisation state: include fiducial xsec (with acc%), pol_frac in SR, and acc_<pol>
    # where pol_frac in SR is (fid_xsec_pol / fid_xsec_unpolarised) * 100.
    # -------------------------
    fid_num = agg.pivot_table(
        index=["process_decay", "operator"],
        columns="pol_state",
        values="xsec_fid_pol_fb",
        aggfunc="first",
    )
    # unpolarised total cross section BEFORE cuts (from txt file)
    xsec_tot = agg.pivot_table(
        index=["process_decay", "operator"],
        columns="pol_state",
        values="xsec_pol_fb",
        aggfunc="first",
    )
    acc_piv = agg.pivot_table(
        index=["process_decay", "operator"],
        columns="pol_state",
        values="acceptance",
        aggfunc="first",
    )
    fid_str = fid_num.copy()
    for c in fid_str.columns:
        fid_str[c] = [
            _format_fiducial_cell(fid_num.loc[idx, c], acc_piv.loc[idx, c])
            for idx in fid_num.index
        ]

    if "unpolarised" in fid_num.columns:
        comb_fid = fid_num["unpolarised"]
        fid_frac = (fid_num.div(comb_fid, axis=0) * 100.0).round(3)
    else:
        fid_frac = fid_num * 0.0

    table_fid = pd.DataFrame(index=fid_num.index)

    # unpolarised first
    table_fid["xsec_fb"] = xsec_tot.get("unpolarised")
    
    # each polarisation: pol fraction (SR), acc
    for pol in ["LL", "LT", "TL", "TT"]:
        table_fid[f"{pol}_pol_frac"] = fid_frac.get(pol)

    table_fid["acc_unpolarised"] = (acc_piv.get("unpolarised") * 100.0).round(3) if "unpolarised" in acc_piv.columns else None
    for pol in ["LL", "LT", "TL", "TT"]:
        table_fid[f"acc_{pol}"] = (acc_piv.get(pol) * 100.0).round(3) if pol in acc_piv.columns else None

    table_fid = table_fid.reset_index()

    def _flatten_columns_for_excel(frame: pd.DataFrame) -> pd.DataFrame:
        out = frame.copy()
        if isinstance(out.columns, pd.MultiIndex):
            out.columns = ["__".join([str(x) for x in col if str(x) != ""]) for col in out.columns]
        return out

    os.makedirs(opts.outdir, exist_ok=True)
    tag = f"vbs_{opts.Mode.lower()}_{opts.EFT_order.lower()}"

    out_xsec_csv = os.path.join(opts.outdir, f"{tag}_pol_table.csv")
    out_fid_csv = os.path.join(opts.outdir, f"{tag}_fiducial_pol_table.csv")
    table_xsec.to_csv(out_xsec_csv, index=False)
    table_fid.to_csv(out_fid_csv, index=False)
    _log(True, f"[OK] wrote -> {out_xsec_csv}")
    _log(True, f"[OK] wrote -> {out_fid_csv}")


    if write_xlsx:
        out_xlsx = os.path.join(opts.outdir, f"{tag}_pol_tables.xlsx")
        out_xlsx_fidu = os.path.join(opts.outdir, f"{tag}_fiducial_pol_tables.xlsx")
        with pd.ExcelWriter(out_xlsx) as writer:
            _flatten_columns_for_excel(table_xsec).to_excel(writer, sheet_name="xsec", index=False)
        with pd.ExcelWriter(out_xlsx_fidu) as writer:
            _flatten_columns_for_excel(table_fid).to_excel(writer, sheet_name="fiducial", index=False)
        _log(True, f"[OK] wrote -> {out_xlsx}")
        _log(True, f"[OK] wrote -> {out_xlsx_fidu}")

if __name__ == "__main__":
    main()