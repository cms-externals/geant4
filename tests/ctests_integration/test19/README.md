# geant4-dev/tests/ctests_integration/test19
**Geant4 Validation application that is meant for testing mainly simulation of pion production in hadron-nuclear interactions, in the range from several GeV/c to 158 GeV/c of beam momenta.**

<!--Project desription-->
## Project description
The application simulates h+A interactions instanciated by such beam as proton, pi+, or pi-, 
with the beam momenta ranging from 3 GeV/c and up to 158 GeV/c, on a various nuclear targets.
Subsequently it looks at the double differential spectra of secondaries, mainly charges pions,
but in some cases also kaons, Lambda, proton,s antiprotons, and neutrons. 

The hadronic models involved in these tests are Bertini (up to 12 GeV/c of beam momentum),
FTFP and/or FTFP_tuneX (both for full range of beam momenta), and INCLXX (up to 12 GeV/c). 

Similated results are compared vs available experimental data, as follows:

- 3, 5, 8, or 12 GeV/c proton, pi+, or pi- beam on vatious nuclear targets :
Momentum spectra of secondary charged pions or protons in different bins of the polar angle of the secondary. \
As a standard candle, interactions of 8 GeV/c protons, pi-, or pi+ with Carbon (C) or Tantalum (Ta)
are simulated (see later on the interface to geant-val).
But the application can also do 3, 5, or 12 GeV/c beam on such nuclear targets as Al, Be, Cu, or Pb.   
Experimental data from HARP : \
M. Apollonio et al., Nucl. Phys. A821 118, 2009 \
M. Apollonio et al., Phys.Rev.C80 065207, 2009 \
M. Apollonio et al., Phys.Rev.C80 035208, 2009 \
M.G. Catanesi et al., Phys.Rev.C77 055207, 2008 \
M. Appolonio et al., Phys.Rev. C82 (2010) 045208

Please note that we have selected 8 GeV/c hadron (proton, pi+, pi-) interacting with Carbon (C) or Tantalum (Ta) nucleus study cases as standard candles to be featured via geant-val (see details later in this document).
However, we stress that the application can model all cases listed above. 

- 31 GeV/c proton or 60 GeV/c pi+ on Carbon : Momentum spectra of secondary chagred pions (primary focus), 
as well charged or neutral kaons, Lambda, protons, antiprotons. \
Experimantal data from NA61 : \
N. Abgrall et al. ,  Eur.Phys.J.C 76, 2016 (proton beam) \
A. Aduszkiewicz et al. , Phys.Rev.D100 112004, 2019 (pion beam, only data on C used so far)

- 158 GeV/c proton on Carbon :  Spectra such as average multiplicity ot average pT as a function Feynman variable xF,
as well as double differential pT-spectra in different bins of xF, for such secondaries as charged pions, protons, antiprotons, or neutrons. \
Experimental data from NA49 : \
http://spshadrons.web.cern.ch/spshadrons

<!--Authors and contacts-->
## Authors and contacts
- Julia Yarba (Fermilab) - yarba_j@fnal.gov

<!--Geant Val integration-->
## Geant Val integration
Geant Val [![Website shields.io](https://img.shields.io/website-up-down-green-red/http/shields.io.svg)](https://geant-val.cern.ch/) is the Geant4 testing and validation suite. It is a project hosted on gitlab.cern.ch [![Website shields.io](https://img.shields.io/website-up-down-green-red/http/shields.io.svg)](https://gitlab.cern.ch/GeantValidation) used to facilitate the maintenance and validation of Geant4 applications, referred to as <em> tests</em>.\
The following are instructions to use geant4-dev/tests/ctests_integration/test19-<EXP>-<target> within Geant Val, 
from batch submission to website deployment.
These instructions have been tested on lxplus.cern.ch as well as on geant4gpvm-test-al9.fnal.gov

1. First of all, one needs to git clone the Geant Val/geant-config-generator as well as geant4-dev/tests/ctests_integration/test23 (common use software) and geant4-dev/tests/ctests_integration/test19 (hadron-nucleus test). 
In order to check out only tests/test23 and tests/test19 directories without checking out the whole content of the geant4-dev repository one needs to use the sparse-checkout feature which is available starting git 2.25.x series.
On most modern resources proper version of git is availble by default, e.g. on lxplus,cern.ch it is 2.39.3.
   ```sh
   git clone ssh://git@gitlab.cern.ch:7999/GeantValidation/geant-config-generator.git
   git clone ssh://git@gitlab.cern.ch:7999/geant4/geant4-dev.git --no-checkout 
   cd geant4-dev
   git sparse-checkout init --cone
   git sparse-checkout set tests/ctests_integration/test23 tests/ctests_integration/test19 # etc.
   git checkout
   cd ..
   ```
2. Copy geant4-dev/tests/ctests_integration/test19/geantval_scripts into geant-config-generator/tests/geant4/
   ```sh
   cp -r geant4-dev/tests/ctests_integration/test19/geantval_scripts/test19-HARP-Carbon geant-config-generator/tests/geant4/
   cp -r geant4-dev/tests/ctests_integration/test19/geantval_scripts/test19-HARP-Tantalum geant-config-generator/tests/geant4/
   ```
3. Build the validation application for whih one needs to setup compiler, ROOT (see also remarks at the end of this section), and Geant4 install; extend PATH to the area where the excutable is
   ```sh
   source /cvmfs/sft.cern.ch/lcg/contrib/gcc/13/x86_64-el9/setup.sh
   source  /cvmfs/sft.cern.ch/lcg/releases/LCG_104b/ROOT/6.28.08/x86_64-el9-gcc13-opt/ROOT-env.sh
   VERSION="11.2.p01"
   PLATFORM="x86_64-el9-gcc13-optdeb"
   source /cvmfs/geant4.cern.ch/geant4/$VERSION/${PLATFORM}/bin/geant4.sh
   export G4INSTALL=/cvmfs/geant4.cern.ch/geant4/${VERSION}/${PLATFORM}
   cd geant4-dev/tests/ctests_integration
   mkdir test19-build
   cd test19-build
   cmake -DCMAKE_PREFIX_PATH=$G4INSTALL -DCMAKE_CXX_COMPILER=g++ ../test19
   make
   export PATH=$PWD:$PATH
   cd ../..
   ```

   With regards to ROOT, one should bear in mind that not every version/build of ROOT (or other external)
   is going to be compatible with a specific release/build of Geant4.
   In order to figure out, one can look into the LCG_view (stack of dependencies).
   One possible way to determine if is to look at the LD_LIBRARY_PATH environment variable, 
   upon setting up Geant4.
   In this specific case, the LD_LIBRARY_PATH wil contain the following "views" : \
   /cvmfs/sft.cern.ch/lcg/views/LCG_104b_geant4ext20231106/x86_64-el9-gcc13-opt \
   One can subsequently execute the following: \
   source /cvmfs/sft.cern.ch/lcg/views/LCG_104b_geant4ext20231106/x86_64-el9-gcc13-opt/setup.sh \
   then inspect the ROOTSYS environment variable which will point to the right release of ROOT. \
   In principle, settings can be done as described right above. However, we opt for explicitly 
   setting ROOT only.

4. Move to geant-config-generator and prepare environment for the test. Geant4 release 11.2.p01 is used here as an example, and test75 is to be executed against it; thus one should make sure that file ```11.2.p01.sh``` exists into ```configs/geant/```, of the following content :
   ```sh
   #!/bin/bash

   # compiler
   # NOTE: in principle, compiler gets setup upon setiing up ROOT, as shown below,
   #       thus this step is somewhat redundant
   source /cvmfs/sft.cern.ch/lcg/contrib/gcc/13/x86_64-el9/setup.sh

   # ROOT
   source  /cvmfs/sft.cern.ch/lcg/releases/LCG_104b/ROOT/6.28.08/x86_64-el9-gcc13-opt/ROOT-env.sh

   # Geant4
   VERSION="11.2.p01"
   PLATFORM="x86_64-el9-gcc13-optdeb"
   source /cvmfs/geant4.cern.ch/geant4/$VERSION/${PLATFORM}/bin/geant4.sh

   # extend path to the area where the test75 Val executable is located
   export PATH=</path/to>/test19-build:$PATH   
   ```
5. Create macros and metadata for Geant Val execution
   ```sh
   python mc-config-generator.py submit -t test19-HARP-Carbon -d OUTPUT -v 11.2.p01 -q "testmatch" -r
   python mc-config-generator.py submit -t test19-HARP-Tantalum -d OUTPUT -v 11.2.p01 -q "testmatch" -r
   ```
   this command creates the Geant Val files for batch submission using HTCondor 
   under the ```OUTPUT``` folder, using test19-HARP-Carbon and/or test19-HARP-Tantalum, 
   Geant4.11.2.p01 and the ```testmatch``` job flavour.

6. Execute the analysis on the ROOT files in the ```OUTPUT``` folder to create Geant Val JSON output files
    ```sh
    python mc-config-generator.py parse -t test19-HARP-Carbon -d OUTPUT 
    python mc-config-generator.py parse -t test19-HARP-Tantalum -d OUTPUT 
    ```
    the analysis is coded in ```tests/geant4/test19-HARP-Carbon/parser.py``` 
    and ```tests/geant4/test19-HARP-Tantalum/parser.py```. 
    The ```OUTPUTJSON``` folder is created with the corresponding JSON files.

7. Deploy the results on Geant4 Val. The layout of test19-HARP-Carbon or test19-HARP-Tantalum 
on Geant4 Val is defined by ```test19-HARP-Carbon.xml``` or ```test19-HARP-Tantalum.xml``` 
that can be found in the ```gitlab.com/thegriglat/geant-val-layouts``` repository. \
   Deploy JSON files on the Geant Val database
   ```sh
    find . -name '*.json' | while read i; do curl -H "Content-Type: application/json" -H "token: askauthor" --data @$i https://geant-val.cern.ch/upload; echo; done



