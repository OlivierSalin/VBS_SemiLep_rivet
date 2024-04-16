# example of running:   athena rivet_example.py -c 'conf="user.okurdysh.MadGraph_WmWm_lvlv_FT0_FULL";DOCUT="YES"'

theApp.EvtMax = -1
print("#####received conf from cmd through -c 'conf=X;DOCUT=Y':  ",conf, DOCUT)
import lib_utils
prod_dec, base_dir = lib_utils.find_prod_dec_and_dir(conf)
evnt_conf_dir,evnt_file  = lib_utils.find_evnt_dir_and_file(base_dir + f"/*{conf}_EXT0")
conf_cut_dir = lib_utils.get_conf_cut_dir(evnt_conf_dir, DOCUT)

import AthenaPoolCnvSvc.ReadAthenaPool
svcMgr.EventSelector.InputCollections = [ evnt_file ]

from AthenaCommon.AlgSequence import AlgSequence
job = AlgSequence()

from xAODEventInfoCnv.xAODEventInfoCnvConf import xAODMaker__EventInfoCnvAlg
job += xAODMaker__EventInfoCnvAlg()

from Rivet_i.Rivet_iConf import Rivet_i
rivet = Rivet_i()
import os
rivet.AnalysisPath = os.environ['PWD']

# rivet.Analyses += [f'{prod_dec}:DOCUT={DOCUT}']
rivet.Analyses += [f'{prod_dec}:OUTDIR={conf_cut_dir}']
rivet.RunName = ''
rivet.HistoFile = conf_cut_dir + f'/MyOutput.yoda.gz'
rivet.CrossSection = 1.0 #xsec_pb
#rivet.IgnoreBeamCheck = True
job += rivet

