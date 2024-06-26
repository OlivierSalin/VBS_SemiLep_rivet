import utils_func as uf

processes = ["WpZ","ZZ"]
decays = ["llqq"]

all_ops_cat = ["SM","FM0","FM2","FS1","FT1","FT5"]
#all_ops_cat = ["FM0","FM2","FS1","FT1","FT5"]
#all_ops_cat = ["SM"]
#all_ops_cat = ["SM","FM0","FM2","FS1"]

detail_path="Tables/Truth/Truth_eta_40/"
detail_path="Tables/BDT_sample/BDT_sample_40000/"
output_dir="/exp/atlas/salin/ATLAS/VBS_mc/EFT_files_AMI/BDT_sample/BDT_sample_4000/"

uf.combine_ntuples(detail_path,processes, decays, all_ops_cat,output_dir,truth=False)