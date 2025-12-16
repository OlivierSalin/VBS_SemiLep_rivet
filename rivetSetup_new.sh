# for json in c++
lsetup git
cd /exp/atlas/salin/ATLAS/VBS_mc/vcpkg/
./bootstrap-vcpkg.sh
./vcpkg install nlohmann-json
cd /exp/atlas/salin/ATLAS/VBS_mc/VBS_Pol_Rivet/VBS_rivet

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh 
asetup 23.6.40,AthGeneration
source setupRivet

lsetup "panda"
lsetup "pyami"
lsetup "rucio -w"
lsetup "astyle"
lsetup "git"


voms-proxy-init -voms atlas --valid 48:0 # 48h grid certif and not 12h default

