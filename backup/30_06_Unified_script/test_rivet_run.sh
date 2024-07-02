rivet-build RivetWpZ_llqq.so WpZ_llqq.cc -I/exp/atlas/salin/ATLAS/VBS_mc/vcpkg/installed/x64-linux/include/
python run_rivet.py --evtMax 10000 --conf "user.osalin.MadGraph_WpZ_llqq_FM0_SM" --DOCUT "YES" --redoRivet "yes" --redoPlots "yes"
python run_rivet.py --evtMax 10000 --conf "user.osalin.MadGraph_WpZ_llqq_FS1_QUAD" --DOCUT "YES" --redoRivet "yes" --redoPlots "yes"
python run_rivet.py --evtMax 10000 --conf "user.osalin.MadGraph_WpZ_llqq_FT1_QUAD" --DOCUT "YES" --redoRivet "yes" --redoPlots "yes"

python run_rivet.py --evtMax 10000 --conf "user.osalin.MadGraph_WpZ_llqq_FM0_SM" --DOCUT "NO" --redoRivet "yes" --redoPlots "yes"
python run_rivet.py --evtMax 10000 --conf "user.osalin.MadGraph_WpZ_llqq_FS1_QUAD" --DOCUT "NO" --redoRivet "yes" --redoPlots "yes"
python run_rivet.py --evtMax 10000 --conf "user.osalin.MadGraph_WpZ_llqq_FT1_QUAD" --DOCUT "NO" --redoRivet "yes" --redoPlots "yes"