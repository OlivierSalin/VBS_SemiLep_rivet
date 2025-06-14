# example of running:   athena rivet_example.py -c 'conf="user.okurdysh.MadGraph_WmWm_lvlv_FT0_FULL";DOCUT="YES"'

theApp.EvtMax = -1
print(f"##### Received conf from cmd through -c 'conf={conf}; DOCUT={DOCUT}; Part={Part}; type_MC={type_MC}'")
import lib_utils
import os
import json

print(f"type_MC: {type_MC}")
prod_dec, base_dir = lib_utils.find_prod_dec_and_dir_tres(conf,type_MC)
#evnt_conf_dir,evnt_file  = lib_utils.find_evnt_dir_and_file(base_dir + f"/*{conf}_EXT0")
print("base_dir: ", base_dir)

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
if '/Generation/Parameters' in metadata:
    genpars=metadata['/Generation/Parameters']
    #print("GenPars: ", genpars)
    if 'HepMCWeightNames' in genpars:
        systWeights=genpars['HepMCWeightNames']
        #print("SystWeights: ", systWeights)

Weight_name_rwg=[]
print('len(systWeights):',len(systWeights))
#print('systWeights:',systWeights)
# Write systWeights to a JSON file
systWeights_dict = {item.strip(): idx for item, idx in systWeights.items()}
with open(os.path.join(conf_cut_dir, 'systWeights.json'), 'w') as f:
    json.dump(systWeights_dict, f, indent=4)

for item in systWeights:
    systName=''
    for n,s in enumerate(item.split()):
        while s.find('.') >= 0:
            s = s.replace('.','p')
        s.strip()
        if n:
            systName=systName+'_'+s
        else:
            systName=s
        #print(f'item: {item}, systName: {systName}')
        if 'DYNSCALE' in systName or 'PDF303000' in systName:
            #print('weight name:',item,', output name',systName)
            Weight_name_rwg.append(systName)
            #print(f'Event_Weight_{systName}')
        if 'QUAD' in systName or 'cross' in systName:
            #print('weight name:',item,', output name',systName)
            Weight_name_rwg.append(systName)
            print(f'{systName}')


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
rivet.RunName = f''
rivet.HistoFile = conf_cut_dir + f'/MyOutput.yoda.gz'
rivet.MatchWeights = '.*PDF303000gugu.*'
#rivet.CrossSection = 1.0 #xsec_pb
#rivet.IgnoreBeamCheck = True
job += rivet

