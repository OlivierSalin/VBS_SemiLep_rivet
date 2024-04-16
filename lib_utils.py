import glob
import os
from array import array
import ROOT
import itertools
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import pandas as pd
import numpy as np
import math
import json
import pyAMI.client
import pyAMI.atlas.api as AtlasAPI
client = pyAMI.client.Client('atlas')
AtlasAPI.init()

def get_im_color(particle_name, for_distribution=False):
    color = "black"
    if "jet" in particle_name: 
        color="blue" if not for_distribution else "deepskyblue"
    elif "lepton" in particle_name:
        color="red" if not for_distribution else "lightcoral"
    elif "photon" in particle_name: 
        color = "orange" if not for_distribution else "yellow"
    elif "MET" in particle_name:
        color = "dimgrey" if not for_distribution else "darkgrey"

    return color

def get_fitted_plot(prod_dec):
    mystr=""
    if prod_dec=="Zy_vvy":mystr="pt_photon"
    return mystr

def get_var_latex(varname): # to be used with ROOT latex
    latexstr=varname
    if varname=="pt_photon": latexstr="p_{T}^{\gamma} [GeV]"
    return latexstr

def latex_ana_str(prod_dec):
    mystr=""
    if "Zy_lly" in prod_dec: mystr = "Z(" + r"$\rightarrow$" + f"ll)y"
    return mystr

def draw_average_event(files_dir, average_im=True):
    op = files_dir[files_dir.find("_F")+1 : files_dir.find("_EXT0")]
    prod_dec = files_dir[files_dir.find("/eft_files/")+11 : files_dir.find("/user.osalin")]
    print("found op and prod_dec", op, prod_dec)

    f_str = files_dir + "info_for_image.csv" 
    # "/exp/atlas/kurdysh/vbs_cross_terms_study/eft_files/Zy_lly/user.okurdysh.MadGraph_Zy_lly_FT7vsFT8_CROSS_EXT0/DOCUT_YES/Rivet.csv"
    df = pd.read_csv(f_str,delimiter=";")
    particles = set([i_col[i_col.find("_")+1:] for i_col in  list(df.head()) if "Unnamed" not in i_col])
    print("found particles", particles, "in file", f_str)

    fig, ax = plt.subplots(figsize=(8, 6))
    plt.clf()
    if average_im:
        for i_part in sorted(list(particles)):
            i_eta = round(np.mean(df["eta_" + i_part]), 2)
            i_phi = round(np.mean(df["phi_" + i_part]), 2)
            i_pt = round(np.mean(df["pt_" + i_part]), 2)
            i_color = get_im_color(i_part)
            print("drawing", i_part, "with eta phi pt", i_eta, i_phi, i_pt, "and color", i_color)
            plt.plot([0, i_pt*math.sinh(i_eta)], [0, i_pt*math.sin(i_phi)], 
                    color = i_color, label=f"{i_part} $<p_T>={i_pt:.1f}$,$<\eta>={i_eta:.1f}$")
        plt.legend(loc="upper left")
    else:
        for _, row in df.iterrows():
            for i_part in sorted(list(particles)):
                i_eta = round(row["eta_" + i_part], 2)
                i_phi = round(row["phi_" + i_part], 2)
                i_pt = round(row["pt_" + i_part], 2)
                i_color = get_im_color(i_part, for_distribution=True)
                plt.plot([0, i_pt*math.sinh(i_eta)], [0, i_pt*math.sin(i_phi)], color = i_color, alpha=0.4)
    plt.xlabel('beam Z')
    plt.ylabel('Y')
    if average_im: my_title = f"{latex_ana_str(prod_dec)} SR av.img for {op}"
    else: my_title = f"{latex_ana_str(prod_dec)} SR distr.img for {op}" 
    plt.title(my_title)
    # take super long to save pdf switch to png
    save_path_no_ext = files_dir+"average_image" if average_im else files_dir+"distribution_image"  
    plt.savefig(save_path_no_ext+".png", bbox_inches='tight') 
    plt.savefig(save_path_no_ext+".svg", bbox_inches='tight') 
    return

def get_cutflow_arrays(cutflow_file):
    names = []
    cumu = []
    incr = []
    with open(cutflow_file) as f:
        for i_line in f.readlines():
            line_arr = [i_thing.strip().replace("%", "") for i_thing in i_line.split(" ") if i_thing != ""]
            print(line_arr)
            if len(line_arr) != 5: continue
            names.append(line_arr[1])
            cumu.append(float(line_arr[3]))
            incr.append(float(line_arr[4]))
    return names, cumu, incr

def draw_cutflows(cut_names, y_arrays, labels_arr, outname, prod_dec=""):
    cycle_colors = itertools.cycle('bgrcmk')
    plt.clf()
    fig, ax = plt.subplots(figsize=(6, 4))
    for y_arr,label in zip(y_arrays, labels_arr):
        plt.plot(range(len(cut_names)), y_arr, marker='o', label=label, color=next(cycle_colors))
    plt.xticks(range(len(cut_names)), cut_names, fontsize=7)
    plt.xticks(rotation=90)
    for xc in cut_names: ax.axvline(x=xc, color='0.3', linestyle='--', linewidth=0.3)
    plt.ylim(0, 110)
    plt.legend()
    plt.ylabel('fraction $[\%]$')
    plt.xlabel('cut')
    op = outname[outname.find("_F")+1 : outname.find("_EXT0")]
    plt.title(latex_ana_str(prod_dec)+" "+op)
    plt.savefig(outname, bbox_inches='tight')

# @TODO maybe dublication with find_evnt_dir_and_file
def get_evnt_log_files(base_dir,i_job_name):
    evnt_did, evnt_dir, log_did, log_dir_before_untar = get_envt_log_names_dirs(base_dir,i_job_name)
    evnt_file = -1
    log_file = -1
    if os.path.exists(evnt_dir) and os.path.exists(log_dir_before_untar):
        print("directories for evnt and log exist")
        evnt_candidates = glob.glob(evnt_dir + "/*EVNT.root")
        log_candidates = glob.glob(log_dir_before_untar + "/tarball_PandaJob*/log.generate")
        print("evnt candidates of len",len(evnt_candidates), evnt_candidates)
        print("log candidates of len",len(log_candidates), log_candidates)
        if len(evnt_candidates)==1 and len(log_candidates)==1:
            evnt_file = evnt_candidates[0]
            log_file =  log_candidates[0]
    else:
        print("directories for evnt and log DOESN:T exist")

    print("returning evnt file", evnt_file)
    print("returning log file", log_file)
    return evnt_file, log_file

def get_plotdir(prod_dec, DOCUT_str):
    my_dir = f"/exp/atlas/salin/ATLAS/VBS_mc/plotting/plots/{prod_dec}/{DOCUT_str}/" 
    if not os.path.exists(my_dir): os.makedirs(my_dir)
    return my_dir

def get_pair_str(op1,op2):
    mypair = sorted([op1,op2])
    return f"{mypair[0]}vs{mypair[1]}"

def get_bookletdir(start_path, normalized=""):
    my_dir = start_path + "/booklets/"
    if len(normalized)>0: my_dir = my_dir[:-1] + "_normalized/"
    if not os.path.exists(my_dir): os.makedirs(my_dir)
    if not os.path.exists(my_dir + "/svg/"): os.makedirs(my_dir + "/svg/")
    return my_dir

def get_ops(include_fs0_2):
    all_ops = ["FM0","FM1","FM2","FM3","FM4","FM5","FM7",
                "FS02","FS1",
                "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]
    if include_fs0_2: 
        all_ops += ["FS0","FS2"]
        all_ops.remove("FS02")
    op_pairs = sorted(list(itertools.combinations(all_ops,2)))
    return all_ops, op_pairs

def get_hists_bounds_cuts(prod_dec):
    with open(f"{prod_dec}_hists.json") as fo: total_h = json.load(fo)
    # always load analysis specific hists and jet ones
    with open("jet_hists.json") as fo: jet_h = json.load(fo)
    with open("photon_hists.json") as fo: photon_h = json.load(fo)
    with open("lepton_hists.json") as fo: lepton_h = json.load(fo)
    total_h.update(jet_h)
    # update others depending on ana
    if prod_dec in ["Zy_lly", "Wmy_lvy"]:
        total_h.update(lepton_h)
        total_h.update(photon_h)
    elif prod_dec=="Zy_vvy":
        total_h.update(photon_h)

    # this loop is unncessecary
    return_dict = {} # key is hist name, value[nbin,min,max,cut,cutdir="+,-"]
    for i_hist, i_h_arr in total_h.items():
        if len(i_h_arr)!=3: 
            print("----were not able to  for",i_hist,"will not be there in kin plots")
            continue
        return_dict[i_hist]=i_h_arr
    # replace original bin params with what want to have [rebin_x, x_low, x_up]
    params={}
    params["pt_tagjet1"]  = [10, 0, 2500]
    params["pt_tagjet2"] = [10, 0, 1000]
    params["eta_tagjets"] = [3, -1, -1]
    params["phi_tagjets"] = [5, 0, 6]
    params["m_tagjets"] = [10, -1, -1]
    params["dy_tagjets"]= [5, 0, 9]
    for i_dphi in ["dphi_tagjets","dphi_MET_photon", "dphi_MET_tagjet"]:
        params[i_dphi]= [-1, 0, 3.5]
    params["pt_lepton"]  = [10, 0, 2500]
    params["eta_lepton"] = [3, -3, 3]
    params["pt_photon"]  = params["pt_MET"] = [10, 0, 4000]
    params["eta_photon"] = [3, -3, 3]
    params["m_ll"] = [10,0,300]
    params["m_lly"] =  params["m_ly"] = [10,-1,1]
    for i_cent in ["centrality_lly","centrality_jjy","centrality_jjly"]:
        params[i_cent] = [-1,0,0.8]
    params["cone_frac_photon"] = [-1, 0, 0.1]
    params["m_W_T"] = [3, 0, 150]
    params["dR_lepton_photon"] = [2, 0, 6]
    params["dR_tagjets"] = [2,0,10]
    default_params = [-1,-1,-1]
    for i_hist, i_h_arr in return_dict.items():
        if i_hist in params.keys(): replace_params=params[i_hist]
        else: replace_params = default_params
        i_h_arr[0]=replace_params[0]
        i_h_arr[1]=replace_params[1]
        i_h_arr[2]=replace_params[2]

    print("returning plotting dict")
    return return_dict 


# def get_root_hist_param(plot_name):
#     params = {}
#     params["photon_iso_pt"] = params["m_lly"] = [-1,-1,10]
#     params["all_lep_pairs_m_ll"] = [0, 2500, 10] 
#     params["m_ll_of_pairs_best_quadruplet"] = [0, 200, 1]


def find_prod_dec_and_dir(conf):
    prod_temp = conf[conf.find("user.osalin.MadGraph_")+len("user.osalin.MadGraph_"):]
    print("start from string", prod_temp)
    prod_dec = prod_temp[:prod_temp.find("_F")]
    print("from conf found production dec", prod_dec)
    conf_dir = f"/exp/atlas/salin/ATLAS/VBS_mc/eft_files/{prod_dec}/"
    print("dir would be", conf_dir)
    return prod_dec, conf_dir

def find_evnt_dir_and_file(search_com):
    conf_dir_arr = glob.glob(search_com)
    print("found possibilities for dir", conf_dir_arr)
    conf_dir = conf_dir_arr[0] if len(conf_dir_arr)>=1 else -1  
    if conf_dir == -1: raise ValueError("did not find folder for this config ",search_com)

    evnt_file_candidates = glob.glob(conf_dir + "/*EVNT.root")
    evnt_file = evnt_file_candidates[0] if len(evnt_file_candidates)>0 else -1

    return conf_dir, evnt_file

def get_conf_cut_dir(evnt_dir, docut):
    mydir = evnt_dir + f"/DOCUT_{docut}/"
    if not os.path.exists(mydir): os.makedirs(mydir) 
    return mydir

def get_job_name_template(prod,dec,op,EFTmode):
    return f"user.osalin.MadGraph_{prod}_{dec}_{op}_{EFTmode}"

def find_last_match_job(names_arr, str_op_EFTmode):
    matches_temp = []
    for i_name in names_arr:
        if str_op_EFTmode in i_name: matches_temp.append(i_name)
    matches = sorted(matches_temp)
    print("found possible jobs", matches)
    last_try_job = matches[-1] if len(matches)>0 else -1 
    print("select the last job ", last_try_job)
    return last_try_job

def get_envt_log_names_dirs(base_dir,i_job_name):
    evnt_did = i_job_name + "_EXT0"
    evnt_dir = base_dir + "/" + evnt_did 
    log_did = i_job_name + ".log"
    log_dir = evnt_dir + "/" + log_did
    print("returnning envt did and dir",evnt_did, evnt_dir)
    print("returnning log did and dir",log_did, log_dir)
    return evnt_did, evnt_dir, log_did, log_dir  

def get_rivet_com(job_name, evtMax=-1, redoRivet=-1, redoPlots=-1, DOCUT=-1):
    mycom = f'python run_rivet.py --conf="{job_name}" '
    if evtMax!=-1: mycom += f' --evtMax {evtMax} '
    if redoRivet!=-1: mycom += f' --redoRivet "{redoRivet}" '
    if redoPlots!=-1: mycom += f' --redoPlots "{redoPlots}" '
    if DOCUT!=-1: mycom += f' --DOCUT "{DOCUT}" '
    return mycom

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


def cross_section_fb(EFT_op, proces, dec):
    """
    This function takes an EFT operator, a process, and a decay as string param,
    and returns the cross section in fb.
    """
    all_ops = ["FM0","FM1","FM2","FM3","FM4","FM5","FM7",
               "FS02","FS1",
               "FT0","FT1","FT2","FT5","FT6","FT7","FT8","FT9"]
    
    all_ops2 = ["FM","FS","FT"]
    Processes = [proces]
    Decay = [dec]

    organized_info = {}

    for op in all_ops2:
        for process in Processes:
            for decay in Decay:
                pattern = f'%MGPy8EG_aQGC{op}%_QUAD_1_{process}_{decay}%'
                List_ami = AtlasAPI.list_datasets(client, patterns=[pattern], 
                                                  fields=['ldn', 'cross_section', 'dataset_number'], 
                                                  limit=[1, 100], type='EVNT', project='mc16%')

                if List_ami:
                    for dataset in List_ami:
                        operator = next((op for op in all_ops if op in dataset['ldn']), None)
                        if operator:
                            key = f"{operator}_{process}_{decay}"
                            if key not in organized_info:
                                organized_info[key] = {'name': [], 'cross_section': [], 'dataset_number': []}

                            organized_info[key]['name'].append(dataset['ldn'])
                            organized_info[key]['cross_section'].append(dataset['cross_section']) # Cross section in nb
                            organized_info[key]['dataset_number'].append(dataset['dataset_number'])
                else:
                    print(f'No datasets found for {op}_{process}_{decay}\n')

    # Construct the key from the EFT operator, process, and decay
    key_ = f"{EFT_op}_{proces}_{dec}"

    # Check if the key exists in the organized_info dictionary
    if key_ in organized_info:
        # If the key exists, return the first cross section
        xsection_nb = organized_info[key_]['cross_section'][0] # Return the cross section in nb
        xsection_fb = float(xsection_nb) * 1e6 # Convert the cross section to fb
        return xsection_fb
    else:
        # If the key does not exist, print an error message and return None
        print(organized_info[key_]['cross_section'][0])
        print(f'Error: No datasets found for {key_}')
        return None

def get_sumw_initial(log_file):
    with open(log_file) as textf:
        sumw_neg, sumw_pos, filt_ef = -999.0 , -999.0, -999.0
        for line in textf:
            if 'MetaData: sumOfPosWeights =' in line:
                sumw_pos = float(line[line.find('=')+1:])
            if "MetaData: sumOfNegWeights =" in line:
                sumw_neg = float(line[line.find('=')+1:])
            if "MetaData: GenFiltEff =" in line:
                filt_ef = float(line[line.find('=')+1:])    
    sumw_in = (sumw_pos + sumw_neg) * filt_ef 
    print(f"found sumw_in: {sumw_in} build from sumw_neg, pos and filt_ef: {sumw_neg}, {sumw_pos}, {filt_ef}")
    return sumw_in

def save_xsec_frac_prod(savedir,xsec_fb,
                        frac,frac_pos,frac_neg, frac_er_bar, 
                        pos_w_in, neg_w_in, pos_w_f, neg_w_f,
                        pos_n_in, neg_n_in, pos_n_f, neg_n_f):
    write_to_f(savedir + "xsec_fb.txt",xsec_fb)
    #
    write_to_f(savedir + "frac_after_cuts.txt",frac)
    write_to_f(savedir + "frac_after_cuts_pos.txt",frac_pos)
    write_to_f(savedir + "frac_after_cuts_neg.txt",frac_neg)
    write_to_f(savedir + "frac_after_cuts_error_bar.txt", frac_er_bar)
    #
    write_to_f(savedir + "xsec_times_frac_fb.txt",xsec_fb*frac)
    #
    write_to_f(savedir + "pos_sumw_in.txt", pos_w_in)
    write_to_f(savedir + "neg_sumw_in.txt", neg_w_in)
    write_to_f(savedir + "pos_sumw_f.txt", pos_w_f)
    write_to_f(savedir + "neg_sumw_f.txt", neg_w_f)
    #
    write_to_f(savedir + "pos_w_n_in.txt", pos_n_in)
    write_to_f(savedir + "neg_w_n_in.txt", neg_n_in)
    write_to_f(savedir + "pos_w_n_f.txt", pos_n_f)
    write_to_f(savedir + "neg_w_n_f.txt", neg_n_f)
    

def write_to_f(product_file,product):
    f = open(product_file, "w") #since opening yoda a bit slow for 50 configs save when have it
    f.write(str(product))
    f.close()

def yoda_to_root_1d(h_yoda, yoda_title):
    print("converting to root hist with name", yoda_title)
    root_title =  yoda_title
    mjj_h_root = ROOT.TH1D(root_title, root_title, h_yoda.numBins(), array('d', h_yoda.xEdges()))
    mjj_h_root.Sumw2()
    mjj_h_root.GetXaxis().SetTitle(yoda_title)
    rtErrs = mjj_h_root.GetSumw2()
    for i in range(mjj_h_root.GetNbinsX()):
        mjj_h_root.SetBinContent(i + 1, h_yoda.bin(i).sumW())
        rtErrs.AddAt(h_yoda.bin(i).sumW2(), i+1)
    mjj_h_root.SetDirectory(0)
    return mjj_h_root.Clone()

def get_op_from_dir(mydir,prod_dec):
    temp1 = mydir[len(f"user.osalin.MadGraph_{prod_dec}_"):]
    temp2 = temp1[:temp1.find("_EXT0")]
    if "try" in temp2: temp3 = temp2[:temp2.find("_try")]
    else: temp3  = temp2
    arr = temp3.split("_")
    print("vec for this str", arr)
    ops = arr[0]
    regime = arr[1]
    ops_arr = []
    if "vs" in ops:
        ops_arr.append(ops[:ops.find("vs")])
        ops_arr.append(ops[ops.find("vs")+2:])
    else:
        ops_arr.append(ops)
    return sorted(ops_arr), regime

def save_plot(plot,path_to_save, draw_option = "text45", log_scale = False, legend=False):
    c=ROOT.TCanvas()
    plot.Draw(draw_option)
    if log_scale: ROOT.gPad.SetLogy()
    if legend: 
        l=ROOT.gPad.BuildLegend()
        l.SetFillColorAlpha(0, 0) # transparent label
    c.Modified()
    c.Update()
    c.Show()
    c.SaveAs(path_to_save)

def read_hists(root_file_name, h_names_arr):
    hists = {}
    if os.path.exists(root_file_name):
        h_file =  ROOT.TFile.Open(root_file_name, "READ")
        hists_in_file = [i_obj.GetName() for i_obj in h_file.GetListOfKeys()]
        # print("list of hists in this file", hists_in_file, "in file ", root_file_name)
        for i_name in h_names_arr:
            if i_name in hists_in_file: 
                i_hist = h_file.Get(i_name)
                i_hist.SetDirectory(0)
                hists[i_name] = i_hist 
            else:
                print("there is no hist", i_name, "in file ", root_file_name)
        h_file.Close()
    return hists

def dress_hist(my_hist, my_title, my_color, my_norm = 1.0):
    my_hist.SetTitle(my_title)
    my_hist.SetName(my_title)
    my_hist.SetLineColor(my_color)
    my_hist.SetMarkerColor(my_color)
    hist_integ = my_hist.Integral()
    # print("normalize", my_hist.GetName(), my_hist.GetTitle(), "to", my_norm, "start from integ", hist_integ)
    if hist_integ!=0:
        my_hist.Scale(my_norm/hist_integ)
        # print("get back integ", my_hist.Integral())
    # else: 
    #     print("dont normalize since integral is 0")
    return my_hist.Clone()

def make_stack(hist_arr, title="",norm=-1):
    my_stack = ROOT.THStack(f"{title}/{norm}", f"{title}/{norm}")
    for i_plot in hist_arr:
        plot_copy = i_plot.Clone() # since will use several time with/without normalzietions  
        if norm!=-1 and plot_copy.Integral()!=0: plot_copy.Scale(norm/plot_copy.Integral())
        my_stack.Add(plot_copy)
    return my_stack.Clone()


### for comparisons
# title = "pt_photon"
# rebin = 10
# def get_hist(op_mode, color):
#     root_file_name =  f"/exp/atlas/kurdysh/vbs_cross_terms_study/eft_files/Zy_vvy/user.okurdysh.MadGraph_Zy_vvy_{op_mode}_EXT0/DOCUT_YES/hists.root"
#     root_file = ROOT.TFile.Open(root_file_name, "READ")
#     hist_in = root_file.Get(title)
#     hist = hist_in.Clone()
#     hist.SetDirectory(0)
#     hist.SetTitle(op_mode)
#     hist.SetName(op_mode)
#     hist.SetLineColor(color)
#     hist.SetMarkerColor(color)
#     hist.RebinX(rebin)
#     hist.Sumw2()
#     return hist
# plot_q1 = get_hist("FT1_QUAD", 2)
# plot_q2 = get_hist("FT2_QUAD", 3)
# plot_c = get_hist("FT1vsFT2_CROSS", 1)


def get_hist_dict(i_h):
    dict_q1 = {}  # key is bin x-axis L edge, value is [bin content, error]
    for ix in range(1, i_h.GetNbinsX() + 1):
        l_edge = i_h.GetBinLowEdge(ix)
        cont = i_h.GetBinContent(ix)
        err = i_h.GetBinError(ix)
        if cont != 0: dict_q1[l_edge] = [cont, err]
    return dict_q1

def get_ratio_plot_tests(hist_1, hist_2): # _2 is the one with to respect to which hist_1 is taken
    dict_1 = get_hist_dict(hist_1)
    dict_2 = get_hist_dict(hist_2)
    bigger_xaxis_dict = dict_1 if len(dict_1) > len(dict_2) else dict_2
    smaller_xaxis_dict = dict_2 if bigger_xaxis_dict == dict_1 else dict_1
    ratio_hist_for_gr_range = ROOT.TH1F("dummy", "dummy", hist_1.GetNbinsX(), 
                                        hist_1.GetXaxis().GetXmin(),hist_1.GetXaxis().GetXmax())
    ratio_plot = ROOT.TGraphErrors(ratio_hist_for_gr_range)
    ratio_test_points = []
    for num_p,ix in enumerate(bigger_xaxis_dict.keys()):
        if ix not in smaller_xaxis_dict.keys(): continue
        i_1 = bigger_xaxis_dict[ix]
        i_2 = smaller_xaxis_dict[ix]
        ratio = i_1[0]/i_2[0]
        ratio_plot.SetPoint(num_p, ix, ratio)
        ratio_plot.SetPointError(num_p, 0, math.sqrt((i_1[1])**2+(i_2[1])**2))
        if ratio>=1: 
            ratio_test_points.append(ratio-1)
        else:
            ratio_test_points.append(1/ratio-1)
    ratio_test = np.sum(ratio_test_points)/len(ratio_test_points)
    rchi2 = hist_1.Chi2Test(hist_2, "CHI2/NDF")
    print("from ratios", ratio_test_points, "get RT", ratio_test, "with chi2/ndf", rchi2)
    ratio_name = f"{hist_1.GetName()}/{hist_2.GetName()} RT={ratio_test:.2f} rChi2={rchi2:.2f}"
    ratio_plot.SetName(ratio_name)
    ratio_plot.SetTitle(ratio_name)
    ratio_color = hist_1.GetLineColor() # notice it's hist_1 where color is taken
    ratio_plot.SetMarkerColor(ratio_color)
    ratio_plot.SetLineColor(ratio_color)
    return ratio_plot.Clone(), ratio_test, rchi2

def draw_stack_with_ratio(my_stack, mg_ratios, xtitle, outname, stack_x_range=[]):
    c = ROOT.TCanvas()
    y_divide = 0.4
    pad_1 = ROOT.TPad("pad1","pad1", 0.0, y_divide, 1.0, 1.0)
    pad_2 = ROOT.TPad("pad1","pad1", 0.0, 0.0, 1.0, y_divide)
    pad_1.Draw()
    pad_2.Draw()

    pad_1.cd()
    my_stack.Draw("nostack")
    if len(stack_x_range)!=0: 
        my_stack.GetXaxis().SetRangeUser(stack_x_range[0], stack_x_range[1])
    ROOT.gPad.SetLogy()
    my_stack.GetXaxis().SetTitle(xtitle)
    l = ROOT.gPad.BuildLegend()
    l.SetFillColorAlpha(0, 0)

    pad_2.cd()
    mg_ratios.Draw("AP")
    # without range can be slight visual different in range wrt to top plot
    if len(stack_x_range)==0:
        mg_ratios.GetXaxis().SetRangeUser(my_stack.GetXaxis().GetXmin(),my_stack.GetXaxis().GetXmax())
    else:
        mg_ratios.GetXaxis().SetRangeUser(stack_x_range[0], stack_x_range[1])
    l = ROOT.gPad.BuildLegend()
    l.SetFillColorAlpha(0, 0)
    c.Modified()
    c.Update()
    c.Show()
    c.SaveAs(outname)