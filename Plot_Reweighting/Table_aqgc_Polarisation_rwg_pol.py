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
    parser.add_option("--Mode", default="Lep")
    parser.add_option("--EFT_order", default="QUAD")
    parser.add_option("--outdir", default="/exp/atlas/salin/ATLAS/VBS_mc/Plots/Plot_Reweighting/MCRequest/Validation/Tables/acceptance/")
    parser.add_option("--invert", default=False)
    parser.add_option(
        "--xsec_file",
        default="/exp/atlas/salin/ATLAS/VBS_mc/VBS_Pol_Rivet/VBS_rivet/Plot_Reweighting/Tables/Cross_section/Validation/Validation/Validation/VBS_cross_section_run2_reweighting_vbslep_pol_Lep.txt",
        help="Path to VBS_cross_section_*.txt"
    )
    parser.add_option("--verbose", default=True, help="Verbose prints (True/False)")
    parser.add_option(
        "--pol_base_dir",
        default="/exp/atlas/salin/ATLAS/VBS_mc/eft_files/Histograms/Polarisation",
        help="Base dir holding Polarisation/<Mode>/Rwg/GOOD/Pol/<Mode>/<op_rwg>/..."
    )
    opts, _ = parser.parse_args()
    return opts

def _as_bool(x):
    if isinstance(x, bool):
        return x
    return str(x).strip().lower() in ("1", "true", "yes", "y", "on")

def _log(verbose, msg):
    if verbose:
        print(msg)

# -------------------------
# Acceptance computation
# -------------------------

def _weight_branch_for_op(op, order_eft="QUAD"):
    op_order = f"{op}_{order_eft}".lower().replace("vs", "_")
    return f"EventWeight_{op_order}"

def _tree_sum_of_branch(tree, branch_name):
    if not tree or not tree.GetBranch(branch_name):
        return 0.0
    s = np.float64(0.0)
    for ev in tree:
        try:
            s += np.float64(getattr(ev, branch_name))
        except Exception:
            continue
    return float(s)

def _infer_process_decay_from_path(path):
    m = re.search(r"/([A-Za-z0-9]+)_([A-Za-z0-9]+)/", path)
    if not m:
        return None, None
    return m.group(1), m.group(2)

def _infer_pol_state_from_path(path):
    # detect explicit polarisation tag; otherwise treat as combined
    for pol in ("LL", "LT", "TL", "TT"):
        if re.search(rf"(^|[_/]){pol}([_/]|$)", path):
            return pol
    return "COMB"

def find_pol_ntuple_paths(pol_base_dir, mode):
    """
    Minimal discovery:
      <pol_base_dir>/<process>_<decay>/**/ntuple_rivet.root
    """
    rows = []
    pattern = os.path.join(pol_base_dir, "*_*", "**", "ntuple_rivet.root")
    for nt in glob.glob(pattern, recursive=True):
        process, decay = _infer_process_decay_from_path(nt)
        if not process or not decay:
            continue
        rows.append({
            "mode": mode,
            "process": process,
            "decay": decay,
            "pol_state": _infer_pol_state_from_path(nt),
            "ntuple_path": nt,
        })
    # de-dup
    seen = set()
    uniq = []
    for r in rows:
        if r["ntuple_path"] in seen:
            continue
        seen.add(r["ntuple_path"])
        uniq.append(r)
    return uniq

def find_pol_ntuple_paths_for_process_decay(pol_base_dir, mode, process, decay):
    """
    Discovery restricted to one process+decay:
      <pol_base_dir>/<process>_<decay>/**/ntuple_rivet.root
    """
    rows = []
    pattern = os.path.join(pol_base_dir, f"{process}_{decay}", "**", "ntuple_rivet.root")
    for nt in glob.glob(pattern, recursive=True):
        rows.append({
            "mode": mode,
            "process": process,
            "decay": decay,
            "pol_state": _infer_pol_state_from_path(nt),
            "ntuple_path": nt,
        })
    # de-dup
    seen = set()
    uniq = []
    for r in rows:
        if r["ntuple_path"] in seen:
            continue
        seen.add(r["ntuple_path"])
        uniq.append(r)
    return uniq

def _parse_xsec_key(key_str):
    """
    Key formats in file:
      <proc>_<decay>_<OP>_<ORDER>
      <proc>_<decay>_<OP>_<ORDER>_<POL>
    Example:
      ZZ_llll_FM0_QUAD
      ZZ_llll_FM0_QUAD_LL
    """
    parts = key_str.strip().split("_")
    if len(parts) < 4:
        return None
    proc = parts[0]
    decay = parts[1]
    op = parts[2]
    order = parts[3]
    pol = parts[4] if len(parts) >= 5 else "COMB"
    return proc, decay, op, order, pol

def read_cross_sections_txt(path):
    """
    Returns dict: (proc, decay, op, order, pol_state) -> xsec_fb (float)
    """
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or ":" not in line:
                continue
            k, v = line.split(":", 1)
            parsed = _parse_xsec_key(k.strip())
            if not parsed:
                continue
            try:
                out[parsed] = float(v.strip())
            except ValueError:
                continue
    return out

def _lookup_xsec_fb(xsec_map, process, decay, operator, order, pol_state):
    # exact match first
    v = xsec_map.get((process, decay, operator, order, pol_state))
    if v is not None:
        return v
    # common alias: COMB vs missing pol suffix (already encoded as COMB)
    if pol_state != "COMB":
        return None
    return xsec_map.get((process, decay, operator, order, "COMB"))

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
    if pol in ["LL", "LT", "TL", "TT"]:
        key = f"{process}_{decay}_{op}_{order_EFT}_{pol}"
    else:
        key = f"{process}_{decay}_{op}_{order_EFT}"
    with open(VBS_txt, "r") as f:
        for line in f:
            info_xsec = line.strip().split(": ")
            if len(info_xsec) > 1:
                key_file, xsection_fb = info_xsec[0], info_xsec[1]
                if key_file == key:
                    return float(xsection_fb)
    return None

def compute_acceptance_table(
    ntuple_infos,
    all_operators,
    mode,
    order_eft="QUAD",
    invert_ratio=False,
    xsec_file=None,
    verbose=False,
    op_rwg=None,
):
    before_tree = "Bef_cut"
    after_tree = "SR" if mode == "Lep" else "Merged"
    xsec_map = xsec_map or {}
    _log(verbose, f"[ACC] mode={mode} before_tree={before_tree} after_tree={after_tree} order={order_eft} invert={invert_ratio}")

    _log(verbose, f"[ACC] processing {len(ntuple_infos)} ntuples")

    out = []
    for info in ntuple_infos:
        _log(verbose, f"[ACC] ntuple={info['ntuple_path']}")
        f = ROOT.TFile.Open(info["ntuple_path"], "READ")
        if not f or f.IsZombie():
            _log(verbose, f"[WARN] cannot open ROOT file: {info['ntuple_path']}")
            continue

        t_before = f.Get(before_tree)
        t_after = f.Get(after_tree)
        if not t_before or not t_after:
            _log(verbose, f"[WARN] missing trees in {info['ntuple_path']} (have Bef_cut? {bool(t_before)} ; have {after_tree}? {bool(t_after)})")
            f.Close()
            continuepol_base_dir=

        # quick branch listing once (helps debugging naming mismatches)
        if verbose:
            bnames_before = set([b.GetName() for b in t_before.GetListOfBranches()])
            bnames_after = set([b.GetName() for b in t_after.GetListOfBranches()])
            _log(verbose, f"[ACC] branches: Bef_cut={len(bnames_before)}  after={len(bnames_after)}  pol_state={info['pol_state']}")

        for op in all_operators:
            br = _weight_branch_for_op(op, order_eft=order_eft)

            # prints only for first few ops to keep stdout readable
            do_op_print = verbose and (j_op < 3)

            if do_op_print:
                _log(verbose, f"  [OP] {op} branch={br}")

            # branch existence hints
            if verbose:
                if not t_before.GetBranch(br):
                    _log(verbose, f"  [WARN] missing in Bef_cut: {br}")
                if not t_after.GetBranch(br):
                    _log(verbose, f"  [WARN] missing in {after_tree}: {br}")

            sum_before = _tree_sum_of_branch(t_before, br)
            sum_after = _tree_sum_of_branch(t_after, br)


            acc = None
            if sum_before != 0.0 and sum_after != 0.0:
                acc = (sum_before / sum_after) if invert_ratio else (sum_after / sum_before)

            xsec_fb = None
            if xsec_file and os.path.exists(xsec_file):
                pol = info["pol_state"] if info["pol_state"] in ("LL", "LT", "TL", "TT") else None
                xsec_fb = take_xsec_fb(xsec_file, op, order_eft, info["process"], info["decay"], pol=pol)
                if xsec_fb is None and pol is not None:
                    xsec_fb = take_xsec_fb(xsec_file, op, order_eft, info["process"], info["decay"], pol=None)

            print(f"xsec: {xsec_fb} fb")
            print(f"    sum_before={sum_before:.6g} sum_after={sum_after:.6g} acc={acc if acc is None else f'{acc:.6g}'} xsec_fb={xsec_fb}")
            if do_op_print:
                _log(verbose, f"    sum_before={sum_before:.6g} sum_after={sum_after:.6g} acc={acc if acc is None else f'{acc:.6g}'} xsec_fb={xsec_fb}")

            if verbose and xsec_fb is None:
                _log(verbose, f"  [WARN] missing xsec for key=({info['process']},{info['decay']},{op},{order_eft},{info['pol_state']})")

            out.append({
                "mode": info["mode"],
                "op_rwg": op_rwg,
                "process": info["process"],
                "decay": info["decay"],
                "pol_state": info["pol_state"],
                "operator": op,
                "order": order_eft,
                "xsec_fb": xsec_fb,
                "weight_branch": br,
                "sumW_before": sum_before,
                "sumW_after": sum_after,
                "acceptance": acc,
                "ntuple_path": info["ntuple_path"],
            })

        f.Close()

    df = pd.DataFrame(out)
    # stable ordering
    if not df.empty:
        df = df.sort_values(["process", "decay", "pol_state", "operator"]).reset_index(drop=True)
    _log(verbose, f"[ACC] produced {len(df)} rows")
    return df

def main():
    opts = parse_options()
    verbose = _as_bool(getattr(opts, "verbose", False))

    _log(verbose, f"[CFG] Mode={opts.Mode} EFT_order={opts.EFT_order} invert={opts.invert}")
    _log(verbose, f"[CFG] outdir={opts.outdir}")
    _log(verbose, f"[CFG] xsec_file={opts.xsec_file}")
    _log(verbose, f"[CFG] pol_base_dir={opts.pol_base_dir}")

    if opts.Mode == "Lep":
        print("Running in Leptonic mode")
        valid_combinations = {
            4: [("ZZ", "llll")],
            3: [("WZ", "lllv")],
            2: [("ssWW", "lvlv"), ("osWW", "lvlv"), ("Zy", "lly")],
        }
    elif opts.Mode == "SemiLep":
        print("Running in Semi-Leptonic mode")
        valid_combinations = {
            1: [("WZ", "lvqq"), ("ssWW", "lvqq"), ("osWW", "lvqq")],
            2: [("WZ", "llqq"), ("ZZ", "llqq")],
            0: [("WZ", "vvqq"), ("ZZ", "vvqq")],
        }
    else:
        valid_combinations = {2: [("WpZ", "vvqq"), ("WmZ", "vvqq")]}

    proc_decay_pairs = []
    for pairs in valid_combinations.values():
        proc_decay_pairs.extend(pairs)
    _log(verbose, f"[CFG] will scan {len(proc_decay_pairs)} (process,decay) pairs: {proc_decay_pairs}")

    all_operators_all = [
        "FM0","FM1","FM2","FM3","FM4","FM5","FM7","FM8","FM9",
        "FS0","FS1","FS2",
        "FT0","FT1","FT2","FT3","FT4","FT5","FT6","FT7",
        "FM1odd","FM2odd","FM3odd","FM4odd","FM5odd","FM6odd",
        "FT2odd","FT3odd","FT4odd","FT5odd","FT6odd"
    ]

    operators_rwg = ["FM", "FS", "FT", "FModd", "FTodd"]
    Pol_state= ["LL", "LT", "TL", "TT", "COMB"]
    mode_folder = "Leptonic" if opts.Mode == "Lep" else "SemiLep"

    dfs = []
    for op_rwg in operators_rwg:
        for pol in Pol_state:
            _log(verbose, f"[GRP] op_rwg={op_rwg} pol={pol}")
            if pol == "COMB":
                pol_state = "all"
                end_path_op=f"{op_rwg}_QUAD"
                pol_base_dir = os.path.join(opts.pol_base_dir, mode_folder, "Rwg", "GOOD", "Pol", mode_folder, end_path_op)
            else:
                pol_state = pol
                end_path_op=f"{op_rwg}_QUAD_{pol}"
                pol_base_dir = os.path.join(opts.pol_base_dir, mode_folder, "Rwg", "GOOD", "Pol", mode_folder, end_path_op)
                
            ops_this = filter_operators_by_group(all_operators_all, op_rwg)
            _log(verbose, f"[GRP] op_rwg={op_rwg} -> {len(ops_this)} operators")

            
            _log(verbose, f"[GRP] pol_base_dir={pol_base_dir}")
            if not os.path.isdir(pol_base_dir):
                _log(verbose, f"[WARN] missing directory for op_rwg={op_rwg}: {pol_base_dir}")
                continue

            for (process, decay) in proc_decay_pairs:
                ntuple_infos = find_pol_ntuple_paths_for_process_decay(pol_base_dir, opts.Mode, process, decay)
                if not ntuple_infos:
                    continue
                dfs.append(
                    compute_acceptance_table(
                        ntuple_infos=ntuple_infos,
                        all_operators=ops_this,
                        polarisation=pol_state,
                        mode=opts.Mode,
                        order_eft=opts.EFT_order,
                        invert_ratio=bool(opts.invert),
                        xsec_file=opts.xsec_file,
                        verbose=verbose,
                        op_rwg=op_rwg,
                    )
                )

    df = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()

    # ...existing output writing (keep as-is)...
    tag = f"vbs_{opts.Mode.lower()}_acceptance_{opts.EFT_order.lower()}_ALLGROUPS"
    csv_path = os.path.join(opts.outdir, f"{tag}.csv")
    txt_path = os.path.join(opts.outdir, f"{tag}.txt")
    tex_path = os.path.join(opts.outdir, f"{tag}.tex")

    df.to_csv(csv_path, index=False)
    with open(txt_path, "w") as f:
        f.write(df.to_string(index=False))
        f.write("\n")
    with open(tex_path, "w") as f:
        f.write(df.to_latex(index=False, float_format="%.6g"))

    print(f"[OK] wrote {len(df)} rows")
    print(f"  {csv_path}")
    print(f"  {txt_path}")
    print(f"  {tex_path}")

if __name__ == "__main__":
    main()