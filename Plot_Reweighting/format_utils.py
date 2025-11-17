import utils_func as uf
import os
import re



def format_model_name(model_name,op_rwg=None):
    if "aqgc" in model_name:
        model_name = model_name.replace("aqgc", "aQGC").replace("Aqgc", "aQGC").replace("new", "new basis").replace("_", " ")
    elif "rwg" in model_name or "reweight" in model_name:
        model_name = model_name.replace("rwg", "Rwg").replace("_", " ").replace("fs", "").replace("fm", "").replace("ft", "")
        model_name = model_name.replace("Nohel", "Hel ignorant").replace("hel_aware", "Hel aware")
    elif "Reweighting" in model_name or "Reweight" in model_name:
        #model_name = model_name.replace("Reweighting",f"Rwg {op_rwg}").replace("_", " ")
        #model_name = model_name.replace("Reweight",f" {op_rwg} Rwg").replace("_", " ")
        #model_name = model_name.replace("Reweighthel_ignore", f"{op_rwg} Rwg hel ignore: ")
        #model_name = model_name.replace("Reweighthel_aware", f"{op_rwg} Rwg hel aware: ")
        model_name = model_name.replace("Reweighthel_ignore", f"{op_rwg} Rwg hel ignore: ")
        model_name = model_name.replace("Reweighthel_aware", f"{op_rwg} Rwg hel aware: ")
        model_name = model_name.replace("Polarisation", f"Pol")
        

    elif "EFTDec" in model_name:
        model_name = model_name.replace("EFTDec", "EFT Dec").replace("_", " ")
    elif "Polarisation" in model_name:
        model_name = model_name.replace("Polarisation", "Pol").replace("_", " ")
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

def take_xsec_fb(VBS_txt,op, order_EFT, process, decay,pol=None):
    # Create the key
    if pol in ["LL", "LT", "TL", "TT"]:
        key = f"{process}_{decay}_{op}_{order_EFT}_{pol}"
    else:
        key = f"{process}_{decay}_{op}_{order_EFT}" 
    
    with open(VBS_txt, 'r') as f:
        for line in f:
            # Split the line into key and value
            info_xsec= line.strip().split(': ')
            if len(info_xsec)>1:
                key_file, xsection_fb = info_xsec[0], info_xsec[1]
                if key_file == key:
                    return float(xsection_fb)
    
    # If the key was not found in the file, return None
    return None




#def get_cross_section_reweighting(path ,EFT_op, EFT_type, proc, decay, name_model_,pol=None):
#    if "reweighting" in name_model_.lower() or "rwg" in name_model_.lower():
#        log_path= f"{path}/log.generate"
#        if not os.path.exists(log_path):
#            print(f"Log file {log_path} does not exist.")
#            return None
#        PATTERN_XSEC_RWG = re.compile(r'INFO: (\w+_\w+) : ([\d.\-e]+) \+- ([\d.\-e]+) pb')
#        with open(log_path, 'r') as file:
#            content = file.read()
#        matches = PATTERN_XSEC_RWG.findall(content)







def get_cross_section(EFT_op, EFT_type, proc, decay, name_model_,pol=None):
    path = "/exp/atlas/salin/ATLAS/VBS_mc/plotting//Plot_Reweighting/Tables/Cross_section/"
    lower_name_model_ = name_model_.lower()
    if "polarisation" in lower_name_model_ or "pol" in lower_name_model_ or "xsec" in lower_name_model_:
        xsec_text_name_spe=lower_name_model_
        if "reweighting_hel" in lower_name_model_:
            xsec_text_name_spe = xsec_text_name_spe.replace("reweighting_hel", "reweighthel")
        VBS_txt =  f'{path}/VBS_cross_section_run2_{xsec_text_name_spe}.txt'
        
    else: 
        if "aqgc_new" in name_model_ or "model" in name_model_ or "new" in name_model_:
            VBS_txt = f'{path}/VBS_cross_section_aqgc.txt'
        elif "Eboli_Run3" in name_model_ or "Run3" in name_model_ or "run3" in name_model_:
            VBS_txt = f'{path}/VBS_cross_section_Eboli_run3.txt'
        elif "run2" in name_model_ or "Run2" in name_model_:
            VBS_txt = f'{path}/VBS_xsection_test.txt'

        else:
            VBS_txt = f'{path}/VBS_cross_section_aqgc.txt'
    return take_xsec_fb(VBS_txt, EFT_op, EFT_type, proc, decay,pol)

def get_cross_section_cpodd(EFT_op, EFT_type, proc, decay, name_model_,pol=None):
    path = "/exp/atlas/salin/ATLAS/VBS_mc/plotting//Plot_Reweighting/Tables/Cross_section/"
    lower_name_model_ = name_model_.lower()
    if "eftdec" in lower_name_model_:
        if "polarisation"  in lower_name_model_:
            VBS_txt = f'{path}/VBS_cross_section_run2_eftdec_polarisation.txt'
        elif "int" in lower_name_model_.lower():
            VBS_txt = f'{path}/VBS_cross_section_run2_eftdec_int.txt'
        else:
            VBS_txt = f'{path}/VBS_cross_section_run2_eftdec.txt'
    elif "reweight" in lower_name_model_ or "rwg" in lower_name_model_:
        if "polarisation"  in lower_name_model_ or "pol" in lower_name_model_:
            VBS_txt = f'{path}/VBS_cross_section_run2_rwg_pol.txt'
        elif "int" in lower_name_model_.lower():
            VBS_txt = f'{path}/VBS_cross_section_run2_rwg_int.txt'
        else:
            VBS_txt = f'{path}/VBS_cross_section_run2_rwg.txt'
    else:
        VBS_txt = f'{path}/VBS_cross_section_run2_eftdec.txt'
    print(f"VBS_txt: {VBS_txt} ") 

    return take_xsec_fb(VBS_txt, EFT_op, EFT_type, proc, decay,pol)