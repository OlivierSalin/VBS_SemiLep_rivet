# for json in c++
lsetup git
cd /exp/atlas/salin/ATLAS/VBS_mc/vcpkg/
./bootstrap-vcpkg.sh
./vcpkg install nlohmann-json
cd /exp/atlas/salin/ATLAS/VBS_mc/plotting/

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh 

asetup 23.6.22,AthGeneration
source setupRivet

lsetup "panda 1.5.68"
lsetup "pyami"
lsetup "rucio -w"
voms-proxy-init -voms atlas --valid 48:0 # 48h grid certif and not 12h default

#source /exp/atlas/kurdysh/vbs_cross_terms_study/python_packages/setup.sh # this is for sympy installed with pipInstall
#source /exp/atlas/kurdysh/vbs_cross_terms_study/python_package2/setup.sh # for plotting modules of sympy 
