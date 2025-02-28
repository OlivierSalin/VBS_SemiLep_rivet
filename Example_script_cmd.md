## Download file from the grid in a structed manner withoutt overriding current folder with event files and have more stats
python Panda_script_bis.py --Channel "WpZ_llqq" --nb_lep 2 --type_MC "aqgc" --Add_stats True

python Save_run_rivet_type_mc.py --Channel "WpZ_llqq" --nb_lep 2 --type_MC "run3" --Name "aqgc_Eboli_Run3" --EFT_order "QUAD"


## Rivet
rivet-build RivetWpZ_llqq.so WpZ_llqq.cc -I/exp/atlas/salin/ATLAS/VBS_mc/vcpkg/installed/x64-linux/include/ `root-config --cflags --libs` 

python run_rivet.py --evtMax 100 --conf "user.osalin.MadGraph_WpZ_llqq_FM0_QUAD" --DOCUT "YES" --redoRivet "yes" --redoPlots "no" --type_MC ""

## Run Rivet 
python Script_op_model.py --nb_events 20000  --type_MC "aqgc" --Name "aqgc_new" --EFT_order "QUAD"
python Save_run_rivet_type_mc.py --Channel "WpZ_llqq" --nb_lep 2 --type_MC "aqgc" --Name "aqgc_new" --EFT_order "QUAD"

python Script_op_model.py --nb_events 20000  --type_MC "run3" --Name "aqgc_Run3" --EFT_order "QUAD"
python Save_run_rivet_type_mc.py --Channel "WpZ_llqq" --nb_lep 2 --type_MC "run3" --Name "aqgc_Run3" --EFT_order "QUAD"

