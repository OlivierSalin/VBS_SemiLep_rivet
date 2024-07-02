import ROOT
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gROOT.LoadMacro("~/atlasstyle/AtlasStyle.C")
ROOT.SetAtlasStyle()
ROOT.gROOT.LoadMacro("~/atlasstyle/AtlasLabels.C")
ROOT.gROOT.LoadMacro("~/atlasstyle/AtlasUtils.C")
ROOT.gStyle.SetLegendBorderSize(0)
ROOT.gStyle.SetLegendFillColor(0)
# ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetLegendTextSize(0.04)
import os
import lib_utils as lu
import matplotlib.pyplot as plt
from statistics import mean
plt.rcParams.update({'text.usetex': True})
import math
import numpy as np
import pandas as pd
from pandas.plotting import table
import subprocess
from pandaclient import panda_api
c = panda_api.get_api()
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--tProd", default = "WmWm")
parser.add_option("--tDec", default = "lvlv")
parser.add_option("--runSMAndFULL", default = "no")
parser.add_option("--runQUADAndCROSS", default = "no")
parser.add_option("--numJobsParallel", default = 15)
parser.add_option("--runWithCuts", default = "yes")
parser.add_option("--runAgain", default = "no")
parser.add_option("--sumPlotsOnly", default = "yes")
parser.add_option("--SMOnSumPlots", default = "no")
parser.add_option("--drawCutflow", default = "no")
opts, _ = parser.parse_args()

prod_dec = f"{opts.tProd}_{opts.tDec}"
_, top_files_dir = lu.find_prod_dec_and_dir(f"user.okurdysh.MadGraph_{prod_dec}_FM0_SM")
docut_dir = "DOCUT_YES" if opts.runWithCuts=="yes" else "DOCUT_NO"

plot_dir = lu.get_plotdir(prod_dec, docut_dir)
big_pairs_cutoff = 1.05
plots_no_logy=["n_","ph","dp","et","dy", "dR"]

plots_to_save_info_dict = lu.get_hists_bounds_cuts(prod_dec) 
plots_to_save = plots_to_save_info_dict.keys()

tasks = c.get_tasks(limit=100000000, days=13000, username="Oleksii Kurdysh", status="done") # get already last try since only retry if it failed
task_names = [i_task['taskname'].replace("/","") for i_task in tasks if "MadGraph" in i_task['taskname'] and opts.tProd in i_task['taskname'] and opts.tDec in i_task['taskname']]
all_ops, op_pairs = lu.get_ops(include_fs0_2=True)
splitedSize = int(opts.numJobsParallel)
blocks_of_single_ops = [all_ops[x:x + splitedSize] for x in range(0, len(all_ops), splitedSize)]
blocks_of_op_pairs = [op_pairs[x:x + splitedSize] for x in range(0, len(op_pairs), splitedSize)]
############
# get xsec*frac and save hists to root
##############
def get_com(jobname):
    return f'python run_chain.py --jobName "{jobname}" --runAgain "{opts.runAgain}" --runWithCuts "{opts.runWithCuts}"'

def call_bloc_proc(op_blocks, eft_config):
    for i_num_block,i_ops_block in enumerate(op_blocks,start=1): # loops over blocks of 5
        print ("############## \n ####calling processing of ops", i_ops_block, f"block {i_num_block} out of {len(op_blocks)}")
        i_block_bashcoms = []
        for i_op in i_ops_block:
            if eft_config!="CROSS":
                i_job = lu.find_last_match_job(task_names, f"{i_op}_{eft_config}")
                if i_job!=-1: i_block_bashcoms.append(get_com(i_job))
                else: print("-------------- did not find done job for", i_job)
            else:
                i_job_order_1 = lu.find_last_match_job(task_names, f"{i_op[0]}vs{i_op[1]}_{eft_config}")
                i_job_order_2 = lu.find_last_match_job(task_names, f"{i_op[1]}vs{i_op[0]}_{eft_config}")
                if i_job_order_1!=-1: i_block_bashcoms.append(get_com(i_job_order_1))
                elif i_job_order_2!=-1: i_block_bashcoms.append(get_com(i_job_order_2))
                else: print("-------------- did not find done in both orders for", i_job_order_1, i_job_order_2)
        procs = [subprocess.Popen(i_bash, shell=True) for i_bash in i_block_bashcoms]
        # wait until all processes are finished
        if len(procs)>0: 
            for p in procs: p.wait()
        print (f"############## \n ####finished with block ops num {i_num_block} out of {len(op_blocks)}")

if opts.sumPlotsOnly!="yes":
    if opts.runSMAndFULL=="yes": 
        call_bloc_proc([["FM0"]], "SM")
        call_bloc_proc(blocks_of_single_ops, "FULL")
        
    if opts.runQUADAndCROSS=="yes": 
        call_bloc_proc(blocks_of_single_ops, "QUAD")
        call_bloc_proc(blocks_of_op_pairs, "CROSS")

###########
# build summary plots out of files created above
################


plots_dict = {'SM':{}, 'FULL':{}, 'QUAD':{}, 'CROSS':{}, 'INT':{}}
for i_eft in plots_dict.keys():
    for i_plot in plots_to_save:
        plots_dict[i_eft][i_plot] = {}
xsec_dict = {'SM':{}, 'FULL':{}, 'QUAD':{}, 'CROSS':{}, 'INT':{}}
frac_dict = {'SM':{}, 'FULL':{}, 'QUAD':{}, 'CROSS':{}, 'INT':{}}
frac_er_dict = {'SM':{}, 'FULL':{}, 'QUAD':{}, 'CROSS':{}, 'INT':{}}
cutflow_file_dict = {'SM':{}, 'FULL':{}, 'QUAD':{}, 'CROSS':{}, 'INT':{}}
for op_dir in [i_obj for i_obj in os.listdir(top_files_dir) if os.path.isdir(top_files_dir + "/" + i_obj)]:
    full_op_dir = os.path.join(top_files_dir,op_dir,docut_dir,"")
    print("#### plotting new conf", op_dir, "full path", full_op_dir)
    _, log_file = lu.get_evnt_log_files(top_files_dir,op_dir[:op_dir.find("_EXT0")])
    ops_arr, regime = lu.get_op_from_dir(op_dir, prod_dec)
    print("this is for ops", ops_arr, "in regime", regime)
    product_file = full_op_dir + "xsec_times_frac_fb.txt"
    frac_file = full_op_dir + "frac_after_cuts.txt"
    frac_er_file = full_op_dir + "frac_after_cuts_error_bar.txt"
    hists_file = full_op_dir + "hists.root"
    cutflow_file = full_op_dir + "cutflow.txt"

    if not os.path.exists(product_file) or not os.path.exists(frac_file) or not os.path.exists(hists_file) or not os.path.exists(log_file):
        print("didnt find the txt xsec fid and/or root file or frac file or log file for", op_dir)
        continue 
    if not os.path.exists(cutflow_file) or not os.path.exists(frac_er_file):
        print("didnt find the txt cutflow  and/or frac error file", op_dir)
        continue

    # xsec organization
    with open(product_file, 'r') as f: fid_xsec_fb = float(f.read())
    print("reading back fid_xsec in fb", fid_xsec_fb)
    with open(frac_file, 'r') as f: frac = float(f.read())    
    print("reading back efficiency", frac)
    with open(frac_er_file, 'r') as f: frac_er = float(f.read())    
    print("reading back efficiency error", frac)
    if regime=="CROSS":
        my_op1, my_op2  = ops_arr[0], ops_arr[1]
        if my_op1 not in xsec_dict[regime].keys(): xsec_dict[regime][my_op1] = {}
        xsec_dict[regime][my_op1][my_op2] = fid_xsec_fb
        #
        if my_op1 not in frac_dict[regime].keys(): frac_dict[regime][my_op1] = {}
        frac_dict[regime][my_op1][my_op2] = frac
        #
        if my_op1 not in frac_er_dict[regime].keys(): frac_er_dict[regime][my_op1] = {}
        frac_er_dict[regime][my_op1][my_op2] = frac_er
        #
        if my_op1 not in cutflow_file_dict[regime].keys(): cutflow_file_dict[regime][my_op1] = {}
        cutflow_file_dict[regime][my_op1][my_op2] = cutflow_file
        
    else:
        my_op = ops_arr[0]
        xsec_dict[regime][my_op] = fid_xsec_fb
        frac_dict[regime][my_op] = frac
        frac_er_dict[regime][my_op] = frac_er
        cutflow_file_dict[regime][my_op] = cutflow_file
    # hists organization
    op_hists = lu.read_hists(hists_file, plots_to_save)
    for i_hist_name in plots_to_save:
        if i_hist_name not in op_hists.keys(): continue

        if regime=="CROSS":
            my_op1,my_op2 = ops_arr[0], ops_arr[1]
            if my_op1 not in plots_dict[regime][i_hist_name].keys(): 
                plots_dict[regime][i_hist_name][my_op1] = {}
            plots_dict[regime][i_hist_name][my_op1][my_op2] = op_hists[i_hist_name] 
        else:
            my_op = ops_arr[0]
            plots_dict[regime][i_hist_name][my_op] = op_hists[i_hist_name]
    
print("xsec SM", xsec_dict["SM"])
print("xsec FULL", xsec_dict["FULL"])
print("xsec QUAD", xsec_dict["QUAD"])
for i_key1 in xsec_dict["CROSS"].keys():
    print("xsec CROSS",i_key1,xsec_dict["CROSS"][i_key1])
#
print("saved hists SM", [[i_plot, plots_dict["SM"][i_plot].keys()] for i_plot in plots_to_save])
print("saved hists FULL", [[i_plot, plots_dict["FULL"][i_plot].keys()] for i_plot in plots_to_save])
print("saved hists QUAD", [[i_plot, plots_dict["QUAD"][i_plot].keys()] for i_plot in plots_to_save])
print("saved hists CROSS")
for i_plot in plots_to_save:
    print("for plot", i_plot)
    for i_key1 in plots_dict["CROSS"][i_plot].keys():
        print("have hists", plots_dict["CROSS"][i_plot][i_key1].keys())

all_ops =  sorted(list(xsec_dict["QUAD"].keys()))
print("all ops", all_ops)
nbins = len(all_ops)
if opts.runSMAndFULL=="yes" and len(xsec_dict["SM"])>=1 and len(xsec_dict["FULL"])>=1:
    SM_ref =  xsec_dict["SM"][list(xsec_dict["SM"].keys())[0]] 
    print("SM xsec in fb", SM_ref)
    lu.write_to_f(top_files_dir+"/SM_xsec_times_frac_fb.txt", SM_ref)
    ############# FULL and in comparsion with SM
    FULL_h = ROOT.TH2F("FULL_h","FULL_h", nbins,0,nbins,1,0,1)
    FULL_ratio_SM_h = ROOT.TH2F("FULL_ratio_SM_h","FULL_ratio_SM_h", nbins,0,nbins,1,0,1)
    for num_bin,i_op in enumerate(all_ops, start=1):
        FULL_h.GetXaxis().SetBinLabel(num_bin,i_op)
        FULL_ratio_SM_h.GetXaxis().SetBinLabel(num_bin,i_op)
        if i_op in xsec_dict["FULL"].keys(): 
            FULL_h.SetBinContent(num_bin, 1, xsec_dict["FULL"][i_op])
            FULL_ratio_SM_h.SetBinContent(num_bin, 1, xsec_dict["FULL"][i_op] / SM_ref)
    lu.save_plot(FULL_h,plot_dir + "FULL.pdf")
    lu.save_plot(FULL_ratio_SM_h,plot_dir + "FULL_ratio_SM_h.pdf")

if opts.runQUADAndCROSS=="yes":
    ##########
    # make tables and efficiency plot
    ###########
    QUAD_h = ROOT.TH2F("QUAD_h","QUAD_h", nbins,0,nbins,1,0,1)
    for num_bin,i_op in enumerate(all_ops, start=1):
        QUAD_h.GetXaxis().SetBinLabel(num_bin,i_op)
        if i_op in xsec_dict["QUAD"].keys(): QUAD_h.SetBinContent(num_bin, 1, xsec_dict["QUAD"][i_op])
    lu.save_plot(QUAD_h,plot_dir + "QUAD.pdf")

    CROSS_h = ROOT.TH2F("CROSS_h","CROSS_h", nbins,0,nbins,nbins,0,nbins)
    CROSS_geom_QUAD_h = ROOT.TH2F("CROSS_geom_QUAD_h","CROSS_geom_QUAD_h", nbins,0,nbins,nbins,0,nbins)
    CROSS_el_area_ratio_h = ROOT.TH2F("CROSS_el_area_ratio_h","CROSS_el_area_ratio_h", nbins,0,nbins,nbins,0,nbins)
    pairs_str_avalaible = []
    area_ratios = []
    fracs_quad1 = []
    fracs_quad2 = []
    fracs_cross = []
    fracs_ave = []
    fracs_envelope = []
    fracs_errors_quad1 = []
    fracs_errors_quad2 = []
    fracs_errors_cross = []
    fracs_errors_ave = []
    cutflow_plot_dir = plot_dir + "/cutflows_comp/"
    if not os.path.exists(cutflow_plot_dir): os.makedirs(cutflow_plot_dir)
    for i_op1 in all_ops:
        if i_op1 not in xsec_dict["CROSS"].keys(): continue 

        bin_x = all_ops.index(i_op1) + 1
        CROSS_h.GetXaxis().SetBinLabel(bin_x,i_op1)
        CROSS_geom_QUAD_h.GetXaxis().SetBinLabel(bin_x,i_op1)
        CROSS_el_area_ratio_h.GetXaxis().SetBinLabel(bin_x,i_op1)
        for i_op2 in xsec_dict["CROSS"][i_op1]:
            bin_y = all_ops.index(i_op2) + 1
            CROSS_h.GetYaxis().SetBinLabel(bin_y,i_op2)
            CROSS_geom_QUAD_h.GetYaxis().SetBinLabel(bin_y,i_op2)
            CROSS_el_area_ratio_h.GetYaxis().SetBinLabel(bin_y,i_op2)
            cross = xsec_dict["CROSS"][i_op1][i_op2]
            CROSS_h.SetBinContent(bin_x, bin_y, cross)
            # fill geometric average

            if i_op1 not in xsec_dict["QUAD"].keys() or i_op2 not in xsec_dict["QUAD"].keys(): continue

            quad1 = xsec_dict["QUAD"][i_op1]
            quad2 = xsec_dict["QUAD"][i_op2]
            geom_average = cross / math.sqrt(quad1*quad2) if quad1*quad2>0 else -1
            CROSS_geom_QUAD_h.SetBinContent(bin_x, bin_y, geom_average)
            ##### ellipses
            lumi = 140 # lumi of run2 in 1/fb, expect xsec to be in fb
            el_q1_c = lumi * quad1 / 3 # dont square xsec since already squared from madgraph
            el_q2_c = lumi * quad2 / 3 # divide by 3 since formula for are expect factor to be 1
            el_cross_c = lumi * cross / 3
            den_with_cross2 = 4*el_q1_c*el_q2_c - el_cross_c**2
            if den_with_cross2 > 0:
                area_no_cross = 2 * math.pi / math.sqrt(4*el_q1_c*el_q2_c)
                area_with_cross = 2 * math.pi / math.sqrt(den_with_cross2)
                area_ratio = area_with_cross/area_no_cross
                print(f"for i_op1 i_op2 {i_op1} {i_op2} got areas no cross {area_no_cross:.2f} with cross {area_with_cross:.2f} ratio {area_ratio:.2f} used ellipse q1 q2 cross {el_q1_c:.2f} {el_q2_c:.2f} {el_cross_c:.2f}")
            else:
                print("for i_op1 i_op2", i_op1, i_op2, "cannot do area sqrt in this case",den_with_cross2)
                area_ratio = -99
            CROSS_el_area_ratio_h.SetBinContent(bin_x, bin_y, round(area_ratio,2))
            ##### fractions
            pair_str = lu.get_pair_str(i_op1,i_op2)
            pairs_str_avalaible.append(pair_str)
            area_ratios.append(area_ratio)
            i_frac_quad1 = frac_dict["QUAD"][i_op1]
            i_frac_quad2 = frac_dict["QUAD"][i_op2]
            i_frac_cross = frac_dict["CROSS"][i_op1][i_op2]
            fracs_quad1.append(i_frac_quad1)
            fracs_quad2.append(i_frac_quad2)
            fracs_cross.append(i_frac_cross)
            mean_quads = mean([i_frac_quad1, i_frac_quad2])
            geo_ave = i_frac_cross/mean_quads if mean_quads>0 else 0.0 
            fracs_ave.append(geo_ave)
            fracs_envelope.append(max([abs(i_frac_quad1-i_frac_quad2), 
                                        abs(i_frac_quad1-i_frac_cross), 
                                        abs(i_frac_quad2-i_frac_cross)]))
            print("for pair", pair_str, "founds fracs q1 q2 cross", i_frac_quad1, i_frac_quad2, i_frac_cross)
            print(f"and mean of q1 q2 {mean_quads} then c/ave is {geo_ave}")
            i_frac_quad1_unc = frac_er_dict["QUAD"][i_op1] 
            i_frac_quad2_unc = frac_er_dict["QUAD"][i_op2]
            i_frac_cross_unc = frac_er_dict["CROSS"][i_op1][i_op2] 
            fracs_errors_quad1.append(i_frac_quad1_unc)
            fracs_errors_quad2.append(i_frac_quad2_unc)
            fracs_errors_cross.append(i_frac_cross_unc)
            fracs_errors_ave.append(math.sqrt(i_frac_quad1_unc**2 + i_frac_quad2_unc**2 + i_frac_cross_unc**2))
            print("for pair", pair_str, "founds fracs errors q1 q2 cross", i_frac_quad1_unc, i_frac_quad2_unc, i_frac_cross_unc)
            
            ########### cutflow comparison
            if opts.drawCutflow=="yes":
                i_cutflow_quad1 = cutflow_file_dict["QUAD"][i_op1]
                i_cutflow_quad2 = cutflow_file_dict["QUAD"][i_op2]
                i_cutflow_cross = cutflow_file_dict["CROSS"][i_op1][i_op2]
                #
                cut_names, _, cut_incr_quad1 = lu.get_cutflow_arrays(i_cutflow_quad1)
                _, _, cut_incr_quad2 = lu.get_cutflow_arrays(i_cutflow_quad2)
                _, _, cut_incr_cross = lu.get_cutflow_arrays(i_cutflow_cross)
                
                lu.draw_cutflows(cut_names, [cut_incr_quad1,cut_incr_quad2,cut_incr_cross], 
                                [f"incr {i_op1}",f"incr {i_op2}", f"incr {i_op1}vs{i_op2}"],
                                cutflow_plot_dir+f"/{i_op1}vs{i_op2}.pdf", prod_dec)

    # tables
    lu.save_plot(CROSS_h, plot_dir + "CROSS.pdf")
    lu.save_plot(CROSS_geom_QUAD_h, plot_dir + "CROSS_geom_QUAD.pdf")
    lu.save_plot(CROSS_el_area_ratio_h, plot_dir + "CROSS_el_area_ratio_h.pdf", draw_option = "text colz")

    root_file = ROOT.TFile(plot_dir + "/el_area_ratio.root","UPDATE")  # to divide between ana
    CROSS_el_area_ratio_h.Write("", ROOT.TObject.kOverwrite)
    root_file.Close()
    # fractions 
    for i_pair, i_frac_quad1, i_frac_quad2, i_frac_cross in zip(pairs_str_avalaible,fracs_quad1,fracs_quad2,fracs_cross):
        print(f"for pair {i_pair} q1 q2 c fractions are {i_frac_quad1}, {i_frac_quad2}, {i_frac_cross}")
    
    # def draw_efficiency(area_ratio_min):
    #     new_pairs_str = []
    #     new_fracs_quad1, new_fracs_errors_quad1 = [], []
    #     new_fracs_quad2, new_fracs_errors_quad2 = [], []
    #     new_fracs_cross, new_fracs_errors_cross = [], []
    #     new_fracs_ave, new_fracs_errors_ave = [], []
    #     new_fracs_envelope = []
    #     for i,i_ratio in enumerate(area_ratios):
    #         if i_ratio < area_ratio_min: continue
    #         new_pairs_str.append(pairs_str_avalaible[i])
    #         new_fracs_quad1.append(fracs_quad1[i])
    #         new_fracs_errors_quad1.append(fracs_errors_quad1[i])
    #         new_fracs_quad2.append(fracs_quad2[i])
    #         new_fracs_errors_quad2.append(fracs_errors_quad2[i])
    #         new_fracs_cross.append(fracs_cross[i])
    #         new_fracs_errors_cross.append(fracs_errors_cross[i])
    #         new_fracs_ave.append(fracs_ave[i])
    #         new_fracs_errors_ave.append(fracs_errors_ave[i])
    #         new_fracs_envelope.append(fracs_envelope[i])
    #     new_ops_num = range(len(new_pairs_str))
    #     #
    #     plt.clf()        
    #     fig, (ax1, ax2) = plt.subplots(2, figsize=(16,9), gridspec_kw={'height_ratios': [2, 1]})
    #     ax1.errorbar(new_ops_num, new_fracs_quad1, yerr=new_fracs_errors_quad1, fmt='v', mfc="blue", mec="blue", ecolor="blue", label="QUAD1", alpha=0.7)
    #     ax1.errorbar(new_ops_num, new_fracs_quad2, yerr=new_fracs_errors_quad2, fmt='^', mfc="green", mec="green", ecolor="green", label="QUAD2", alpha=0.7)
    #     ax1.errorbar(new_ops_num, new_fracs_cross, yerr=new_fracs_errors_cross, fmt='*', mfc="black", mec="black", ecolor="black", label="CROSS", alpha=0.7)
    #     ax1.set_ylabel('fraction of weights')
    #     ax1.set_title(r"$A_{cross}/A_{nocross}>$"+str(area_ratio_min))
    #     ax1.legend()
    #     ax1.set_xticks(new_ops_num, minor=False)
    #     ax1.set_xticklabels(new_pairs_str, fontdict=None, minor=False,fontsize = 9)
    #     ax1.tick_params(axis='x', labelrotation=90)
    #     #
    #     color_ave = 'teal'
    #     ax2.set_ylabel(r"$C /<Q1,Q2>$", color=color_ave)
    #     ax2.errorbar(new_ops_num, new_fracs_ave, yerr=new_fracs_errors_ave, mfc=color_ave, mec=color_ave, ecolor=color_ave, fmt='o')
    #     ax2.tick_params(axis='y', labelcolor=color_ave)
    #     ax2.set_xticks(new_ops_num, minor=False)
    #     plt.setp(ax2.get_xticklabels(), visible=False)
    #     ax2.set_xlabel('operator pair')
    #     #
    #     ax3 = ax2.twinx()
    #     color_env = 'purple'
    #     ax3.set_ylabel('envelope', color=color_env)
    #     ax3.scatter(new_ops_num, new_fracs_envelope, color=color_env)
    #     ax3.tick_params(axis='y', labelcolor=color_env)
    #     #
    #     for i_ax in [ax1,ax2]:
    #         for xc in new_ops_num: i_ax.axvline(x=xc, color='0.3', linestyle='--', linewidth=0.3)
    #     fig.tight_layout()
    #     fig.savefig(plot_dir + f"/frac_w_after_cuts_above_{area_ratio_min}.pdf")
    
    # draw_efficiency(0.99)
    # draw_efficiency(1.05)
    
    # ##########
    # # make plots kinematics
    # ###########

    # # for each distrib draw on same plot QUAD1,QUAD2,CROSS - normalied to same area and with each xsec*filt
    # # stacks_arr_per_pair_normalized = {}
    # bookdirbase = lu.get_bookletdir(plot_dir, normalized="yes")
    # for i_plot_name in plots_to_save:
    #     display_params = plots_to_save_info_dict[i_plot_name]
    #     quad_plot_ops = plots_dict["QUAD"][i_plot_name].keys()
    #     for i_key1 in plots_dict["CROSS"][i_plot_name]:
    #         i_key1_keys2 = plots_dict["CROSS"][i_plot_name][i_key1].keys()
    #         if len(i_key1_keys2)==0: continue
    #         for i_key2 in i_key1_keys2:
    #             # for each pairs with this plot - get array (q,q,c,?sm)
    #             print("have hist", i_plot_name, "for CROSS pair", [i_key1,i_key2])
    #             if i_key1 not in quad_plot_ops or i_key2 not in quad_plot_ops: continue
    #             print("also have CROSS plot", i_plot_name, "for both of these ops")
    #             i_plot_cross = plots_dict["CROSS"][i_plot_name][i_key1][i_key2]
    #             i_plot_cross = lu.dress_hist(i_plot_cross, f"CROSS_{i_key1}_{i_key2}", 1, xsec_dict["CROSS"][i_key1][i_key2])
    #             #
    #             i_plot_quad1 = plots_dict["QUAD"][i_plot_name][i_key1]
    #             i_plot_quad1 = lu.dress_hist(i_plot_quad1, f"QUAD_{i_key1}", 2, xsec_dict["QUAD"][i_key1])
    #             #
    #             i_plot_quad2 = plots_dict["QUAD"][i_plot_name][i_key2]
    #             i_plot_quad2 = lu.dress_hist(i_plot_quad2, f"QUAD_{i_key2}", 3, xsec_dict["QUAD"][i_key2])
    #             #
    #             i_op_pair = lu.get_pair_str(i_key1,i_key2)
    #             i_plot_hist_cross_quads_sm = [i_plot_cross.Clone(), i_plot_quad1.Clone(), i_plot_quad2.Clone()]
    #             if "FM0" in plots_dict["SM"][i_plot_name].keys() and opts.SMOnSumPlots=="yes":
    #                 i_plot_sm = plots_dict["SM"][i_plot_name]["FM0"]
    #                 i_plot_sm = lu.dress_hist(i_plot_sm, f"SM", 4, xsec_dict["SM"]["FM0"])
    #                 i_plot_hist_cross_quads_sm.append(i_plot_sm)
    #             # convert this array to stack and draw individual plot
    #             for i_hist in i_plot_hist_cross_quads_sm:
    #                 if display_params[0]!=-1: i_hist.RebinX(display_params[0]) 
    #             i_pair_plot_stack = lu.make_stack(i_plot_hist_cross_quads_sm, title = i_plot_name,norm = 1)
    #             c=ROOT.TCanvas()
    #             i_pair_plot_stack.Draw("nostack")
    #             if display_params[1]!=-1: i_pair_plot_stack.GetXaxis().SetRangeUser(display_params[1], display_params[2]) 
    #             l=ROOT.gPad.BuildLegend()
    #             l.SetFillColorAlpha(0, 0) # transparent background
    #             if i_pair_plot_stack.GetTitle()[:2] not in plots_no_logy: ROOT.gPad.SetLogy()
    #             c.Modified()
    #             c.Update()
    #             c.Show()
    #             i_pair_dir = f"{bookdirbase}/{i_op_pair}/" 
    #             if not os.path.exists(i_pair_dir): os.makedirs(i_pair_dir)
    #             c.SaveAs(f"{i_pair_dir}/{i_plot_name}.pdf")
                
    # # compare for each distribution quads shape
    # quaddir = plot_dir + "/quad_kin/"
    # if not os.path.exists(quaddir): os.makedirs(quaddir)
    # for i_plot_name in plots_to_save:
    #     display_params = plots_to_save_info_dict[i_plot_name]
    #     i_plot_ops_arr =[]
    #     for col, (quad_plot_op, quad_plot_in) in enumerate(plots_dict["QUAD"][i_plot_name].items(),start=1):
    #         quad_plot = lu.dress_hist(quad_plot_in, f"QUAD_{quad_plot_op}", col)
    #         if display_params[0]!=-1: quad_plot.RebinX(display_params[0])
    #         i_plot_ops_arr.append(quad_plot)
    #         i_stack = lu.make_stack(i_plot_ops_arr, title=i_plot_name)
    #         c=ROOT.TCanvas()
    #         i_stack.Draw("nostack")
    #         if display_params[1]!=-1: i_stack.GetXaxis().SetRangeUser(display_params[1], display_params[2])
    #         l=ROOT.gPad.BuildLegend()
    #         l.SetFillColorAlpha(0, 0) # transparent background
    #         if i_stack.GetTitle()[:2] not in plots_no_logy: ROOT.gPad.SetLogy() # for plots like n_jets don't need log scale
    #         c.Modified()
    #         c.Update()
    #         c.Show()
    #         c.SaveAs(f"{quaddir}/{i_plot_name}.pdf")

    # check what can use to replace pair
    fit_plot_str = lu.get_fitted_plot(prod_dec)
    display_params_fitv = plots_to_save_info_dict[fit_plot_str]
    quad_ops_fit = list(plots_dict["QUAD"][fit_plot_str].keys())
    quadpairsdir = plot_dir + "/quad_pair_ratios/"
    if not os.path.exists(quadpairsdir): os.makedirs(quadpairsdir)
    rt_arr = []
    rchi2_arr = []
    pairs_arr = []
    for i_key1 in quad_ops_fit:
        for i_key2 in quad_ops_fit:
            if i_key1==i_key2 or quad_ops_fit.index(i_key1)>quad_ops_fit.index(i_key2): continue
            i_quad1 = lu.dress_hist(plots_dict["QUAD"][fit_plot_str][i_key1], f"QUAD_{i_key1}", 2)
            i_quad2 = lu.dress_hist(plots_dict["QUAD"][fit_plot_str][i_key2], f"QUAD_{i_key2}", 3)
            if display_params_fitv[0]!=-1: 
                i_quad1.RebinX(display_params_fitv[0])
                i_quad2.RebinX(display_params_fitv[0])
            i_stack = lu.make_stack([i_quad1, i_quad2], title=fit_plot_str)
            ratio_plot, rt, rchi2 = lu.get_ratio_plot_tests(i_quad1, i_quad2)
            rt_arr.append(rt)
            rchi2_arr.append(rchi2)
            sorted_keys = sorted([i_key1, i_key2])
            pair_str = f"{sorted_keys[0]}vs{sorted_keys[1]}"
            pairs_arr.append(pair_str)
            i_plot_path = quadpairsdir + f"/{pair_str}.pdf"
            stack_x_range=[] if display_params_fitv[1]==-1 else [display_params_fitv[1],display_params_fitv[2]]
            lu.draw_stack_with_ratio(i_stack, ratio_plot, 
                                    lu.get_var_latex(fit_plot_str), i_plot_path, stack_x_range)
    # save custom RT vs chi2 
    plt.clf()
    plt.plot(rchi2_arr,rt_arr, "o", color="black")
    plt.xlabel('chi2/ndf')
    plt.ylabel('custom ratio test')
    plt.savefig(plot_dir + "/quads_rt_vs_chi2.png")
    # save table for tests
    df = pd.DataFrame(index = quad_ops_fit, columns = quad_ops_fit)
    for i_pair, i_rt, i_rchi2 in zip(pairs_arr, rt_arr, rchi2_arr):
        print("filling table with", i_pair)
        i_op1 = i_pair[:i_pair.find("vs")]
        i_op2 = i_pair[i_pair.find("vs")+2:]
        df.at[i_op1, i_op2] = f"{i_rt:.2f};{i_rchi2:.2f}"
    df = df.fillna('')
    plt.clf()
    ax = plt.subplot(111, frame_on=False)  # no visible frame
    ax.xaxis.set_visible(False)  # hide the x axis
    ax.yaxis.set_visible(False)  # hide the y axis
    table=table(ax, df, loc="center")
    table.set_fontsize(14)
    plt.savefig(plot_dir + "/quads_tests_table.pdf")
    # for latex
    print("got pairs for quads", pairs_arr)


