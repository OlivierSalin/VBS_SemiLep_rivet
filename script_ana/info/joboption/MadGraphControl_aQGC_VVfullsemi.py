from MadGraphControl.MadGraphUtils import *

mode=2
nevents=5000
njobs=8
gridpack_mode=True
gridpack_dir='madevent/'
stringy=""

def aQGC_Generation(
                   aqgcorder="QUAD",
                   aqgctypes=["FT0"],
                   aqgcvalue=1.0,
                   VVgentype="ZZ",
                   VVdecaytype="llqq",
                   rand_seed=1234,
                   nevts=5000,
                   beamEnergy=6500.,
                   pdf='lhapdf',
                   lhaid=260000):

  #---------------------------------------------------------------------------
  # MG5 Proc card
  #---------------------------------------------------------------------------
  procbase=''
  if VVgentype=="WpWp":
    procbase='generate p p > j j w+ w+ QCD=0 '
  elif VVgentype=="WpWm":
    procbase='generate p p > j j w+ w- QCD=0 '
  elif VVgentype=="WmWm":
    procbase='generate p p > j j w- w- QCD=0 '
  elif VVgentype=="WpZ":
    procbase='generate p p > j j w+ z QCD=0 '
  elif VVgentype=="WmZ":
    procbase='generate p p > j j w- z QCD=0 '
  elif VVgentype=="ZZ":
    procbase='generate p p > j j z z QCD=0 '
  else :
    mglog.error('please specify VVgentyle from WW, WZ, or ZZ')
    return -1

  orders={}
  orders["QUAD"]="NP==1"
  orders["INT"]="NP^2==1"
  orders["CROSS"]="NPa^2==1 NPb^2==1"
  orders["FULL"]="NP=1"

  if  "FS02" in aqgctypes : 
    if aqgcorder=="CROSS" : 
      procname = procbase+orders[aqgcorder].replace("NPa","S0").replace("NPb^2==1","S1^2==1")+"\n"
      procname += procbase.replace("generate","add process")+orders[aqgcorder].replace("NPa","S2").replace("NPb^2==1","S1^2==1")
    elif aqgcorder=="INT" : 
      procname = procbase+orders[aqgcorder].replace("NP^2==1","S0^2==1 S2=0")+"\n"
      procname += procbase.replace("generate","add process")+orders[aqgcorder].replace("NP^2==1","S2^2==1 S0=0")
    elif aqgcorder=="QUAD": 
      procname = procbase+orders[aqgcorder].replace("NP","S0")+"\n"
      procname += procbase.replace("generate","add process")+orders[aqgcorder].replace("NP","S2")+"\n"
      procname += procbase.replace("generate","add process")+orders[aqgcorder].replace("NP==1","S0^2==1 S2^2==1")
    elif aqgcorder=="FULL":
      procname = procbase+orders[aqgcorder].replace("NP=1","S0=1 S2=1")+"\n"
  else :
    if aqgcorder=="CROSS" : 
      procname = procbase+orders[aqgcorder].replace("NPa",aqgctypes[0][1:]).replace("NPb",aqgctypes[1][1:])
    else : 
      procname = procbase+orders[aqgcorder].replace("NP",aqgctypes[0][1:])

  process_dir = new_process("""import model SM_Ltotal_Ind5v2020v2_UFO
define p = g u c d s b u~ c~ d~ s~ b~
define j = g u c d s b u~ c~ d~ s~ b~
%s
output -f""" % (procname))


  #---------------------------------------------------------------------------
  # Number of Events
  #---------------------------------------------------------------------------
  safefactor = 1.1
  if nevts>0: nevents = nevts*safefactor
  else: nevents = 20000 * safefactor
  nevents=int(nevents) ## !! must be an integral
  
  #---------------------------------------------------------------------------
  # MG5 Run Card
  #---------------------------------------------------------------------------
  extras = {
    'nevents' : nevents,
    'ebeam1' : beamEnergy,
    'ebeam2' : beamEnergy,
    'pdlabel': pdf,
    'lhaid': lhaid,
    'dynamical_scale_choice': 2,
    #'lhaid':"260000 264000 265000 266000 267000 25100 13100",
    #'reweight_scale' : "True",
    #'reweight_PDF' : "True False False False False False False",
    'asrwgtflavor':"5",
    'lhe_version':"3.0",
    'ptj':"15",
    'ptb':"15",
    'pta':"0",
    'ptl':"4",
    'misset':"0",
    'etaj':"5",
    'etab':"5",
    'etal':"3.0",
    'drjj':"0",
    'drll':"0",
    'draa':"0",
    'draj':"0",
    'drjl':"0",
    'dral':"0",
    'mmjj':"10",
    'mmbb':"10",
    'mmll':"40",
    'maxjetflavor':"5" ,
    'cut_decays'  :'T',
    'use_syst'    :'T',
    'auto_ptj_mjj': 'F',
    'systematics_program': 'systematics',  
    'systematics_arguments': "['--mur=0.5,1,2', '--muf=0.5,1,2', '--dyn=-1,1,2,3,4', '--pdf=errorset,NNPDF30_nlo_as_0119@0,NNPDF30_nlo_as_0117@0,CT14nlo@0,MMHT2014nlo68clas118@0,PDF4LHC15_nlo_30_pdfas']",
    #'ktdurham' : "-1", #Needed for MG2.7
  }
  modify_run_card(process_dir=process_dir,settings=extras)

  #---------------------------------------------------------------------------
  # MG5 param Card
  #---------------------------------------------------------------------------
  ## params is a dictionary of dictionaries (each dictionary is a separate block)
  params={}

  ## default aQGC couplings, set all as 0
  anoinputs={}
  anoinputs['FS0']="0e-12"
  anoinputs['FS1']="0e-12"
  anoinputs['FS2']="0e-12"
  anoinputs['FM0']="0e-12"
  anoinputs['FM1']="0e-12"
  anoinputs['FM2']="0e-12"
  anoinputs['FM3']="0e-12"
  anoinputs['FM4']="0e-12"
  anoinputs['FM5']="0e-12"
  anoinputs['FM6']="0e-12"
  anoinputs['FM7']="0e-12"
  anoinputs['FT0']="0e-12"
  anoinputs['FT1']="0e-12"
  anoinputs['FT2']="0e-12"
  anoinputs['FT3']="0e-12"
  anoinputs['FT4']="0e-12"
  anoinputs['FT5']="0e-12"
  anoinputs['FT6']="0e-12"
  anoinputs['FT7']="0e-12"
  anoinputs['FT8']="0e-12"
  anoinputs['FT9']="0e-12"

  ## update with user defined aQGC couplings
  for aqgc in aqgctypes:
    if aqgc=="FS02" : 
      anoinputs["FS0"]=str(aqgcvalue)+"e-12"   
      anoinputs["FS2"]=str(aqgcvalue)+"e-12"   
    elif aqgc in anoinputs.keys(): anoinputs[aqgc]=str(aqgcvalue)+"e-12"
  params['anoinputs']=anoinputs

  #Update SM parameters to match SM model
  sminputs={'aEWM1':"1.323489e+02"}
  params["SMINPUTS"]=sminputs
  massinputs={'MMU':"0.0","MT":"1.725000e+02"}
  params["MASS"]=massinputs
  yukawainputs={"ymt":"1.725000e+02"}
  params["YUKAWA"]=yukawainputs
  widthinputs={"WT":"1.320000e+00"}
  params["DECAY"]=widthinputs

  modify_param_card(process_dir=process_dir,params=params)

  #--------------------------------------------------------------------------
  # Add reweight Reweight
  #--------------------------------------------------------------------------
  if os.path.exists("reweight_card.dat"):
    copy=subprocess.Popen(["cp","reweight_card.dat","madevent/Cards/reweight_card.dat"])
    copy.wait()
  
  #---------------------------------------------------------------------------
  # MG5 + Pythia8 setup and process (lhe) generation
  #---------------------------------------------------------------------------
  generate(process_dir=process_dir,grid_pack=gridpack_mode,runArgs=runArgs)
  opts.nprocs = 0 #multi-core cleanup after generation

  #--------------------------------------------------------------------------------------------------------------------
  # decay massive particle with MadSpin
  #--------------------------------------------------------------------------------------------------------------------
  if "lvqq"==VVdecaytype : 
    if VVgentype is "WpWp":
      decaysyntax = "decay w+ > l+ vl; decay w+ > j j"
    elif VVgentype is "WmWm":
      decaysyntax = "decay w- > l- vl~; decay w- > j j"
    elif VVgentype is "WpZ":
      decaysyntax = "decay w+ > l+ vl; decay z > j j"
    elif VVgentype is "WmZ":
      decaysyntax = "decay w- > l- vl~; decay z > j j"
    elif VVgentype is "WpWm":
      decaysyntax = "decay w+ > l+ vl; decay w- > j j"
    elif VVgentype is "WmWp":
      decaysyntax = "decay w- > l- vl~; decay w+ > j j"

  elif "llqq"==VVdecaytype : 
    if VVgentype is "WpZ":
      decaysyntax = "decay w+ > j j; decay z > l+ l-"
    elif VVgentype is "WmZ":
      decaysyntax = "decay w- > j j; decay z > l+ l-"
    elif VVgentype is "ZZ":
      decaysyntax = "decay z > j j; decay z > l+ l-"

  elif "vvqq"==VVdecaytype : 
    if VVgentype is "WpZ":
      decaysyntax = "decay w+ > j j; decay z > vl vl~"
    elif VVgentype is "WmZ":
      decaysyntax = "decay w- > j j; decay z > vl vl~"
    elif VVgentype is "ZZ":
      decaysyntax = "decay z > j j; decay z > vl vl~"

  elif "qqqq"==VVdecaytype : 
    if VVgentype is "WpWp":
      decaysyntax = "decay w+ > j j; decay w+ > j j"
    elif VVgentype is "WpWm":
      decaysyntax = "decay w+ > j j; decay w- > j j"
    elif VVgentype is "WmWm":
      decaysyntax = "decay w- > j j; decay w- > j j"
    elif VVgentype is "WpZ":
      decaysyntax = "decay w+ > j j; decay z > j j"
    elif VVgentype is "WmZ":
      decaysyntax = "decay w- > j j; decay z > j j"
    elif VVgentype is "ZZ":
      decaysyntax = "decay z > j j; decay z > j j"

  #Write madspin card, note we revert back to SM
  madspin_card='madspin_card.dat'
  mscard = open(madspin_card,'w')
  mscard.write("""#************************************************************
#*                        MadSpin                           *
#*                                                          *
#*    P. Artoisenet, R. Frederix, R. Rietkerk, O. Mattelaer *
#*                                                          *
#*    Part of the MadGraph5_aMC@NLO Framework:              *
#*    The MadGraph5_aMC@NLO Development Team - Find us at   *
#*    https://server06.fynu.ucl.ac.be/projects/madgraph     *
#*                                                          *
#************************************************************
#Some options (uncomment to apply)
#set Nevents_for_max_weigth 75 # number of events for the estimate of the max. weight
#set BW_cut 15                # cut on how far the particle can be off-shell
#set max_weight_ps_point 400  # number of PS to estimate the maximum for each event
set seed %i

#specify the decay for the final state particles
define l+ = e+ mu+ ta+    
define l- = e- mu- ta-    
define vl = ve vm vt      
define vl~ = ve~ vm~ vt~  
define j = g u c d s b u~ c~ d~ s~ b~
set spinmode none
%s

# running the actual code
launch"""%(rand_seed,decaysyntax))
  mscard.close()

  #Madspin doesn't run with squared couplings, here is hack to trick Madspin by removing the NP coupling order in the process
  gridpack_loc="./madevent/"
  lhetar_path=gridpack_loc+"Events/GridRun_%i/unweighted_events.lhe.gz" %(rand_seed)
  lhe_path=gridpack_loc+"Events/GridRun_%i/unweighted_events.lhe" %(rand_seed)
  proccard_path=gridpack_loc+'/Cards/proc_card_mg5.dat'

  mglog.info('Deleting squred order of couplings at lhe and proc_card.')
  unzip = subprocess.Popen(['gunzip',lhetar_path])
  unzip.wait()
  sed=subprocess.Popen(['sed','-i','-e','s/[SMT][0-9]=\+[0-9eE]*//g',lhe_path])
  sed.wait()
  sed=subprocess.Popen(['sed','-i','-e','s/[SMT][0-9]^2=\+[0-9eE]*//g',lhe_path])
  sed.wait()
  sed=subprocess.Popen(['sed','-i','-e','s/add process.*//g',lhe_path])
  sed.wait()
  zipP = subprocess.Popen(['gzip',lhe_path])
  zipP.wait()
  sed = subprocess.Popen(['sed','-i','-e','s/[SMT][0-9]=\+[0-9eE]*//g',proccard_path])
  sed.wait()
  sed = subprocess.Popen(['sed','-i','-e','s/[SMT][0-9]^2=\+[0-9eE]*//g',proccard_path])
  sed.wait()
  sed=subprocess.Popen(['sed','-i','-e','s/add process.*//g',proccard_path])
  sed.wait()

  #Run madspin
  add_madspin(madspin_card=madspin_card, process_dir=process_dir)

  #--------------------------------------------------------------------------------------------------------------------
  # arrange output file
  #--------------------------------------------------------------------------------------------------------------------
  outputDS=""
  try:
      outputDS=arrange_output(runArgs=runArgs)
  except:
      mglog.error('Error arranging output dataset!')
      return -1

  mglog.info('All done generating events!!')
  return outputDS

#--------------------------------------------------------------------------------------------------------------------
# call the aQGC_Generation
#--------------------------------------------------------------------------------------------------------------------
# extract dataset short name from filename, should be of the form
# "mc.MGPy8EG_aQGCFT0_QUAD_1_ZZjj_qqqq.py" or 
# "mc.MGPy8EG_aQGCFT0vsFT1_INT_1_ZZjj_qqqq.py"

#shortname=runArgs.jobConfig[0].split('/')[-1].split('.')[1]
shortname=jofile.split('/')[-1].split('.')[1]
aqgcorder=""
if "QUAD" in shortname: aqgcorder="QUAD"
elif "INT" in shortname: aqgcorder="INT"
elif "CROSS" in shortname: aqgcorder="CROSS"
elif "FULL" in shortname: aqgcorder="FULL"
else:
  print 'Error=> unrecognized process'
  print '  should be in QUAD, INT, CROSS, or FULL'
  
aqgctypes=[]
if aqgcorder is "QUAD" or aqgcorder is "INT" or aqgcorder is "FULL":  
  aqgctypes.append(shortname.split("_")[1].replace("aQGC",""))
elif aqgcorder is "CROSS":
  aqgctypetmp = shortname.split("_")[1].replace("aQGC","")
  aqgctypes=aqgctypetmp.split("vs")

aqgcvalue=shortname.split("_")[3]
if aqgcvalue[0] == "0":
  aqgcvalue=aqgcvalue.replace("0","0.",1)

VVgentype=""
if "ZZ" in shortname: VVgentype="ZZ"
elif "WpZ" in shortname: VVgentype="WpZ"
elif "WmZ" in shortname: VVgentype="WmZ"
elif "WpWp" in shortname: VVgentype="WpWp"
elif "WpWm" in shortname: VVgentype="WpWm"
elif "WmWm" in shortname: VVgentype="WmWm"
else:
  print 'Error=> unrecognized process'
  print '  should be in ZZ, WpZ, WmZ, WpWp, WpWm, or WmWm'

VVdecaytype=""
if "llqq" in shortname: VVdecaytype="llqq"  
elif "vvqq" in shortname: VVdecaytype="vvqq"  
elif "lvqq" in shortname: VVdecaytype="lvqq"  
elif "qqqq" in shortname: VVdecaytype="qqqq"  
else:
  print 'Error=> unrecognized VV decay'
  print '  should be in llqq, vvqq, lvqq, or qqqq'

# PDF information, in MadGraph's PDF naming scheme.  
# Note that if you change these numbers, you'll probably want to 
# change the "sys_pdf" tag in the run card too.  That's not done
# automatically yet.
pdf='lhapdf'
lhaid=260000 # NNPDF30_nlo_as_0118

#----------------------------------------------------------------------------
# Random Seed
#----------------------------------------------------------------------------
randomSeed = 0
if hasattr(runArgs,'randomSeed'): randomSeed = runArgs.randomSeed

# Run MadGraph!
outputDS=aQGC_Generation( aqgcorder,             # QUAD, INT, CROSS
                          aqgctypes,              # [FT0], [FT1,FT2],,,,
                          aqgcvalue,             # Value if operrator point
                          VVgentype,             # Which VV process
                          VVdecaytype,           # How the W and Z decay
                          randomSeed,            # random seed
                          runArgs.maxEvents,     # number of events for MadGraph to generate
                          runArgs.ecmEnergy/2.,  # Beam energy
                          pdf,                   # PDF information
                          lhaid
                          )


#--------------------------------------------------------------------------------------------------------------------
# multi-core cleanup
#--------------------------------------------------------------------------------------------------------------------
# multi-core running, if allowed!
if 'ATHENA_PROC_NUMBER' in os.environ:
    evgenLog.info('Noticed that you have run with an athena MP-like whole-node setup.  Will re-configure now to make sure that the remainder of the job runs serially.')
    njobs = os.environ.pop('ATHENA_PROC_NUMBER')
    # Try to modify the opts underfoot
    if not hasattr(opts,'nprocs'): mglog.warning('Did not see option!')
    else: opts.nprocs = 0
    print opts


#--------------------------------------------------------------------------------------------------------------------
# Shower
#--------------------------------------------------------------------------------------------------------------------
runArgs.inputGeneratorFile=outputDS
include("Pythia8_i/Pythia8_A14_NNPDF23LO_EvtGen_Common.py")
genSeq.Pythia8.Commands += ["SpaceShower:dipoleRecoil=on"]
include("Pythia8_i/Pythia8_MadGraph.py")

evgenConfig.generators = ["MadGraph", "Pythia8", "EvtGen"]
evgenConfig.description = 'MadGraph aQGC semileptonic or full-hadronic final states'
evgenConfig.keywords+=['SM','ZZ','2jet','VBS']
evgenConfig.inputfilecheck = outputDS
evgenConfig.contact = ['Robert Les <robert.les@cern.ch>','Tatsumi Nitta <tatsumi.nitta@cern.ch>']
