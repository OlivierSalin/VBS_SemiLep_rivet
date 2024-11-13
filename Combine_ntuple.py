import utils_func as uf

processes = ["WpZ","ZZ","WmZ"]
decays = ["llqq"]

all_ops_cat = ["SM","FM0","FM2","FS1","FT1","FT5"]

detail_path="Tables/BDT_sample/BDT_sample_40000/"
detail_path="Tables/Spin/Polarization/"
#output_dir="/exp/atlas/salin/ATLAS/VBS_mc/EFT_files_AMI/BDT_sample/Polar/"
output_dir="/exp/atlas/salin/ATLAS/VBS_mc/EFT_files_AMI/BDT_sample/Polar/"

uf.combine_ntuples_bis(detail_path,processes, decays, all_ops_cat,output_dir,truth=False)
uf.write_xsex_merged(detail_path,processes, decays, all_ops_cat,output_dir)