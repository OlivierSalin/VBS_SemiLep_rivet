# example of running:   athena rivet_example.py -c 'conf="user.okurdysh.MadGraph_WmWm_lvlv_FT0_FULL";DOCUT="YES"'

theApp.EvtMax = -1
print(f"##### Received conf from cmd through -c 'conf={conf}; DOCUT={DOCUT}; Part={Part}; type_MC={type_MC}'")
import lib_utils
prod_dec, base_dir = lib_utils.find_prod_dec_and_dir_tres(conf,type_MC)
#evnt_conf_dir,evnt_file  = lib_utils.find_evnt_dir_and_file(base_dir + f"/*{conf}_EXT0")


evnt_conf_dir,evnt_file,evnt_files  = lib_utils.find_evnt_dir_and_file_bis(base_dir,conf)

if Part != "":
    part_number = int(Part.split('_')[1])
    file_per_part=5
    evnt_conf_dir,evnt_file,evnt_files  = lib_utils.find_evnt_dir_and_file_part(base_dir,conf,part_number,file_per_part)
    evnt_conf_dir = evnt_conf_dir + f'/{Part}/'
    
conf_cut_dir = lib_utils.get_conf_cut_dir(evnt_conf_dir, DOCUT)

print("conf_cut_dir: ", conf_cut_dir)

import AthenaPoolCnvSvc.ReadAthenaPool
#svcMgr.EventSelector.InputCollections = [ evnt_file ]
svcMgr.EventSelector.InputCollections = evnt_files
print("Evnts to run over: ", svcMgr.EventSelector.InputCollections)

from PyUtils import AthFile
af = AthFile.fopen(svcMgr.EventSelector.InputCollections[0]) #opens the first file from the InputCollections list
af.fileinfos

metadata = af.fileinfos['metadata']


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
#rivet.CrossSection = 1.0 #xsec_pb
#rivet.IgnoreBeamCheck = True
job += rivet

