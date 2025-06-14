import utils_func as uf



def format_model_name(model_name,op_rwg=None):
    if "aqgc" in model_name:
        model_name = model_name.replace("aqgc", "aQGC").replace("Aqgc", "aQGC").replace("new", "new basis").replace("_", " ")
    elif "rwg" in model_name or "reweight" in model_name:
        model_name = model_name.replace("rwg", "Rwg").replace("_", " ").replace("fs", "").replace("fm", "").replace("ft", "")
        model_name = model_name.replace("Nohel", "Hel ignorant").replace("hel_aware", "Hel aware")
    elif "Reweighting" in model_name:
        model_name = model_name.replace("Reweighting",f"Reweighting {op_rwg}").replace("_", " ")
    elif "EFTDec" in model_name:
        model_name = model_name.replace("EFTDec", "EFT Decomposition").replace("_", " ")
    return model_name

def format_title_name(title_name):
    return title_name.replace("WpZ", "W^{+}Z").replace("WmZ", "W^{-}Z").replace("llqq", "#rightarrow llqq")

def format_param_name(param_name):
    param_name = param_name.lstrip("merged_").replace("VlepVhad", "VV").replace("_", " ")
    if param_name == "VV mass":
        return "m_{VV} (GeV)"
    elif param_name == "CS V cos theta":
        return "cos#theta_{CS}"
    elif "phi" in param_name:
        return param_name.replace("phi", "#phi")
    elif param_name == "cos theta star":
        return "cos#theta*"
    elif "pt" in param_name:
        return "p_{T}(" + param_name.replace(" pt", "") + ")"
    elif "DeltaEta" in param_name or "delta eta" in param_name:
        return "#Delta#eta(" + param_name.replace("DeltaEta", "").replace("delta eta", "") + ")"
    elif "DeltaPhi" in param_name:
        return "#Delta#phi(" + param_name.replace("DeltaPhi", "") + ")"
    return param_name

def format_Y_axis(X_param, nb_bins, max_hist):
    per_bins = round(max_hist / nb_bins, 2)
    return f"Normalized Entry / {per_bins} (GeV)" if "pt" in X_param or "mass" in X_param else f"Normalized Entry / {per_bins}"


def get_cross_section(EFT_op, EFT_type, proc, decay, name_model_):
    path = "/exp/atlas/salin/ATLAS/VBS_mc/plotting/"
    if "aqgc_new" in name_model_ or "model" in name_model_ or "new" in name_model_:
        VBS_txt = f'{path}/VBS_cross_section_aqgc.txt'
    elif "Eboli_Run3" in name_model_ or "Run3" in name_model_ or "run3" in name_model_:
        VBS_txt = f'{path}/VBS_cross_section_Eboli_run3.txt'
    elif "run2" in name_model_ or "Run2" in name_model_:
        VBS_txt = f'{path}/VBS_xsection_test.txt'
    else:
        VBS_txt = f'{path}/VBS_cross_section_aqgc.txt'
    return uf.take_xsec_fb_aqgc(VBS_txt, EFT_op, EFT_type, proc, decay)