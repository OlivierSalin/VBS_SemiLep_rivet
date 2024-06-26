import utils_function as uf

processes = ["WmZ"]
decays = ["llqq"]

all_ops_cat = ["SM","FM0","FM2","FS1","FT1","FT5"]
#all_ops_cat = ["FM0","FM2","FS1","FT1","FT5"]
#all_ops_cat = ["SM"]
#all_ops_cat = ["SM","FM0","FM2","FS1"]
option_truth=["SM","aQGC","All"]


output_dir_="/exp/atlas/salin/ATLAS/VBS_mc/EFT_files_AMI/Truth/"

detail_path_="Tables/Truth/"

Option_eta = "Truth_bis/Nocut/"
output_dir = output_dir_ + Option_eta  
detail_path = detail_path_ + Option_eta

for option in option_truth:
    if option == "SM":
        all_ops_SM = ["SM"]
        outdir_SM = output_dir + "Truth_SM/"
        uf.combine_ntuples(detail_path,processes, decays, all_ops_SM,outdir_SM,truth=True)
    elif option == "aQGC":
        all_ops_aQGC = ["FM0","FM2","FS1","FT1","FT5"]
        outdir_aQGC = output_dir + "Truth_aQGC/"
        uf.combine_ntuples(detail_path,processes, decays, all_ops_aQGC,outdir_aQGC,truth=True)
    elif option == "All":
        uf.combine_ntuples(detail_path,processes, decays, all_ops_cat,output_dir,truth=True)
