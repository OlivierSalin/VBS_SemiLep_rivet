import ROOT
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(ROOT.kTRUE)
import lib_utils as lu
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--Ana1", default = "Zy_lly")
parser.add_option("--Ana2", default = "Zy_vvy")
opts, _ = parser.parse_args()

all_ops, op_pairs = lu.get_ops(include_fs0_2=True)
n_bins = len(all_ops)

def get_bin_value_for_pair(my_hist, op_1, op_2):
    found_ok = 0
    for i_bin_x in range(1, my_hist.GetNbinsX() + 1):
        for i_bin_y in range(1, my_hist.GetNbinsY() + 1):
            label_x = my_hist.GetXaxis().GetBinLabel(i_bin_x)
            label_y = my_hist.GetYaxis().GetBinLabel(i_bin_y)
            if (label_x==op_1 and label_y==op_2) or (label_x==op_2 and label_y==op_1):
                x_y_val = my_hist.GetBinContent(i_bin_x, i_bin_y)
                # print("found for pair op_1, op_2 value", x_y_val) 
                return x_y_val
    if not found_ok: return -1

def find_num_bin_x_y(my_hist, op_1, op_2):
    found_ok = 0
    for i_bin_x in range(1, my_hist.GetNbinsX() + 1):
        for i_bin_y in range(1, my_hist.GetNbinsY() + 1):
            label_x = my_hist.GetXaxis().GetBinLabel(i_bin_x)
            label_y = my_hist.GetYaxis().GetBinLabel(i_bin_y)
            if (label_x==op_1 and label_y==op_2) or (label_x==op_2 and label_y==op_1):
                # print("found for pair op_1, op_2 binx biny", i_bin_x,i_bin_y)
                return i_bin_x,i_bin_y
    if not found_ok: return -1, -1

base_dir = "/exp/atlas/kurdysh/vbs_cross_terms_study/plots/"

h_file_1 =  ROOT.TFile.Open(f"{base_dir}/{opts.Ana1}/DOCUT_YES/el_area_ratio.root", "READ")
h_file_2 =  ROOT.TFile.Open(f"{base_dir}/{opts.Ana2}/DOCUT_YES/el_area_ratio.root", "READ")

hist_1 = h_file_1.Get("CROSS_el_area_ratio_h")
hist_2 = h_file_2.Get("CROSS_el_area_ratio_h")

hist_ratio = ROOT.TH2F("ratio","ratio",n_bins,0,n_bins, n_bins,0,n_bins)
for num_bin, i_op in enumerate(all_ops, start=1): hist_ratio.GetXaxis().SetBinLabel(num_bin,i_op)
for num_bin, i_op in enumerate(all_ops, start=1): hist_ratio.GetYaxis().SetBinLabel(num_bin,i_op)

for op_pair in op_pairs:
    # print("doing pair", op_pair)
    hist_1_val = get_bin_value_for_pair(hist_1, op_pair[0], op_pair[1])
    hist_2_val = get_bin_value_for_pair(hist_2, op_pair[0], op_pair[1])
    if hist_1_val<=0 or hist_2_val<=0: continue 
    sum_bin_x,sum_bin_y = find_num_bin_x_y(hist_ratio, op_pair[0], op_pair[1])
    if sum_bin_x<0 or sum_bin_y<0: continue 
    hist_ratio.SetBinContent(sum_bin_x, sum_bin_y, abs(hist_1_val/hist_2_val-1))
hist_ratio.SetTitle(f"|{opts.Ana1}/{opts.Ana2}-1|")

c=ROOT.TCanvas()
hist_ratio.Draw("text45 colz")
c.Modified()
c.Update()
c.Show()
c.SaveAs(base_dir + f"/ratio_{opts.Ana1}_{opts.Ana2}.pdf")


