python Script_op_model.py --nb_events 20000  --type_MC "aqgc" --Name "aqgc_new" --EFT_order "QUAD"
python Save_run_rivet_type_mc.py --Channel "WpZ_llqq" --nb_lep 2 --type_MC "aqgc" --Name "aqgc_new" --EFT_order "QUAD"

python Script_op_model.py --nb_events 20000  --type_MC "run3" --Name "aqgc_Run3" --EFT_order "QUAD"
python Save_run_rivet_type_mc.py --Channel "WpZ_llqq" --nb_lep 2 --type_MC "run3" --Name "aqgc_Run3" --EFT_order "QUAD"


