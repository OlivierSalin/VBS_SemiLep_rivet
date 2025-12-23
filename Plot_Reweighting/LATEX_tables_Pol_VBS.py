import pandas as pd
import numpy as np
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--Fidu", default="True")
(opts, args) = parser.parse_args()

def fmt_xsec(x):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return ""
    ax = abs(float(x))
    # 2 significant digits; keep scientific for very small/large
    if ax != 0.0 and (ax < 1e-2 or ax >= 1e3):
        return r"\num{" + f"{x:.1e}" + "}"
    if ax >= 10 and ax < 1e3:
        return r"\num{" + f"{x:.3g}" + "}"
    return r"\num{" + f"{x:.2g}" + "}"

def fmt_pct(pct):
    """Format percent value (already in %) with 4 significant digits, never scientific."""
    if pct is None or (isinstance(pct, float) and np.isnan(pct)):
        return ""
    v = float(pct)
    # Use general format, but forbid scientific
    s = f"{v:.4g}"
    if "e" in s or "E" in s:
        # fallback: fixed-point with enough decimals to avoid scientific for typical % ranges
        s = f"{v:.6f}".rstrip("0").rstrip(".")
    return r"\num{" + s + "}"

def get_op_type(op):
    op_u = op.upper()
    if "ODD" in op_u:
        if op_u.startswith("FM"): return "FModd"
        if op_u.startswith("FT"): return "FTodd"
    else:
        if op_u.startswith("FS"): return "FS"
        if op_u.startswith("FM"): return "FM"
        if op_u.startswith("FT"): return "FT"
    return "Other"

def generate_latex(df, caption, label, process_title_map, process_order, is_fiducial=False):
    if df.empty: return ""
    df = df.copy()
    process_order = [p for p in process_order if p in set(df["process_decay"].unique())]
    operators = list(dict.fromkeys(df["operator"].tolist()))

    pol_rows = [
        ("LL", "LL_pol_frac", "acc_LL"),
        ("LT", "LT_pol_frac", "acc_LT"),
        ("TL", "TL_pol_frac", "acc_TL"),
        ("TT", "TT_pol_frac", "acc_TT"),
    ]

    cols_per_proc = 2 if is_fiducial else 1
    colspec = "ll" + ("r" * cols_per_proc) * len(process_order)

    lines = []
    # longtable: allows splitting over pages (vertical fit)
    lines.append(r"\begin{table}[h]")
    lines.append(r"\caption{" + caption + r"}\label{" + label + r"}") 
    # Restrict scientific notation usage to cross sections only; keep % as plain numbers
    lines.append(r"\setlength{\tabcolsep}{5pt}")
    lines.append(r"\renewcommand{\arraystretch}{1.15}")
    # Fit horizontally by resizing the tabular content
    lines.append(r"\resizebox{\textwidth}{!}{%")
    lines.append(r"\begin{tabular}{" + colspec + r"}")

    lines.append(r"\toprule")

    header_top = [r"\textbf{operator}", r"\textbf{pol}"]
    for p in process_order:
        header_top.append(r"\multicolumn{" + str(cols_per_proc) + r"}{c}{" + process_title_map.get(p, p) + r"}")
    lines.append(" & ".join(header_top) + r" \\")
    
    header_sub = ["", ""]
    for _ in process_order:
        header_sub.append(r"\textbf{pol frac [\%]}")
        if is_fiducial:
            header_sub.append(r"\textbf{acc [\%]}")
    lines.append(" & ".join(header_sub) + r" \\")
    lines.append(r"\midrule")

    def get_row(proc, op):
        sub = df[(df["process_decay"] == proc) & (df["operator"] == op)]
        if len(sub) == 0:
            return None
        return sub.iloc[0]

    nrows = 5

    for op in operators:
        lines.append(r"\multirow{" + str(nrows) + r"}{*}{" + op + r"}")

        row_cells = ["", "unpol"]
        for proc in process_order:
            r0 = get_row(proc, op)
            if r0 is None:
                row_cells += [""] * cols_per_proc
                continue
            # Handle different column names for xsec
            unpol = float(r0["unpolarised_xsec_fb"]) if "unpolarised_xsec_fb" in r0 else float(r0.get("xsec_fb", 0.0))
            unpol_str = fmt_xsec(unpol)
            if unpol_str:
                row_cells.append(r"$\sigma = " + unpol_str + r"$ fb")
            else:
                row_cells.append("")
            
            if is_fiducial:
                acc = r0.get("acc_unpolarised", None)
                row_cells.append(fmt_pct(acc))
        lines.append(" & ".join(row_cells) + r" \\")

        for _, (pol_label, pol_col, acc_col) in enumerate(pol_rows):
            row_cells = ["", pol_label]
            for proc in process_order:
                r0 = get_row(proc, op)
                if r0 is None:
                    row_cells += [""] * cols_per_proc
                    continue
                
                frac_pct = float(r0[pol_col])
                row_cells.append(fmt_pct(frac_pct))
                if is_fiducial:
                    acc = r0.get(acc_col, None)
                    row_cells.append(fmt_pct(acc))
            lines.append(" & ".join(row_cells) + r" \\")

        lines.append(r"\midrule")

    lines.pop()
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"}% end resizebox")
    lines.append(r"\end{table}")

    return "\n".join(lines)

def make_split_tables(df, out_base, caption_base, label_base, process_title_map, process_order, is_fiducial=False):
    df = df.copy()
    df["op_type"] = df["operator"].apply(get_op_type)
    
    categories = ["FS", "FM", "FT", "FModd", "FTodd"]
    all_latex = []
    
    for cat in categories:
        sub_df = df[df["op_type"] == cat]
        if sub_df.empty:
            continue
            
        tex_content = generate_latex(
            sub_df, 
            f"{caption_base} ({cat})" if not is_fiducial else f"Fiducial {caption_base} ({cat})", 
            f"{label_base}_{cat}", 
            process_title_map, 
            process_order,
            is_fiducial=is_fiducial
        )
        
        out_file = f"{out_base}_{cat}.tex"
        with open(out_file, "w", encoding="utf8") as f:
            f.write(tex_content)
            print(f"Wrote table to {out_file}")
        
        all_latex.append(tex_content)
        
    if all_latex:
        if not is_fiducial:
            with open(f"{out_base}_all.tex", "w", encoding="utf8") as f:
                f.write("\n\\newpage\n".join(all_latex))
                print(f"Wrote combined table to {out_base}_all.tex")
        elif is_fiducial:
            with open(f"{out_base}_fiducial_all.tex", "w", encoding="utf8") as f:
                f.write("\n\\newpage\n".join(all_latex))
                print(f"Wrote combined fiducial table to {out_base}_fiducial_all.tex")


base_path = "/exp/atlas/salin/ATLAS/VBS_mc/VBS_Pol_Rivet/VBS_rivet/Plot_Reweighting"
path_to_tables_csv =f"{base_path}/acceptance/GoodTables/"

if opts.Fidu.lower() in ["true", "1", "yes"]:
    outdir_path = f"{base_path}/acceptance/LATEX_tables/Fiducial"
    lep_csv_name = "vbs_lep_quad_fiducial_pol_table.csv"
    semi_csv_name = "vbs_semilep_quad_fiducial_pol_table.csv"
else:
    outdir_path = f"{base_path}/acceptance/LATEX_tables/"
    lep_csv_name = "vbs_lep_quad_pol_table.csv"
    semi_csv_name = "vbs_semilep_quad_pol_table.csv"


lep_df  = pd.read_csv(f"{path_to_tables_csv}/{lep_csv_name}")
semi_df = pd.read_csv(f"{path_to_tables_csv}/{semi_csv_name}")

lep_process_title_map = {
    "ZZ_llll":      r"$ZZ \to \ell\ell\ell\ell$",
    "WZ_lllv":      r"$WZ \to \ell\ell\ell\nu$",
    "ssWW_lvlv":    r"$W^{\pm}W^{\pm} \to \ell\nu\ell\nu$",
    "osWW_lvlv":    r"$W^{+}W^{-} \to \ell\nu\ell\nu$",
    "Zy_lly":       r"$Z\gamma \to \ell\ell\gamma$",
}

semi_process_title_map = {
    "ZZ_llqq":      r"$ZZ \to \ell\ell qq$",
    "WZ_llqq":      r"$WZ \to \ell\ell qq$",
    "ssWW_lvqq":    r"$W^{\pm}W^{\pm} \to \ell\nu qq$",
    "osWW_lvqq":    r"$W^{+}W^{-} \to \ell\nu qq$",
    "ZZ_vvqq":      r"$ZZ \to \nu\nu qq$",
    "WZ_vvqq":      r"$WZ \to \nu\nu qq$",
}

lep_order  = ["ZZ_llll", "WZ_lllv", "ssWW_lvlv", "osWW_lvlv", "Zy_lly"]
semi_order = ["ZZ_llqq", "WZ_llqq", "ssWW_lvqq", "osWW_lvqq", "ZZ_vvqq", "WZ_vvqq"]

make_split_tables(
    lep_df,
    out_base=f"{outdir_path}/table_leptonic",
    caption_base=r"Leptonic final states",
    label_base="tab:lep_xsec_pol",
    process_title_map=lep_process_title_map,
    process_order=lep_order,
    is_fiducial=("fiducial" in lep_csv_name or "fidu" in lep_csv_name)
)

make_split_tables(
    semi_df,
    out_base=f"{outdir_path}/table_semileptonic",
    caption_base=r"Semi leptonic final states",
    label_base="tab:semi_xsec_pol",
    process_title_map=semi_process_title_map,
    process_order=semi_order,
    is_fiducial=("fiducial" in semi_csv_name or "fidu" in semi_csv_name)
)
