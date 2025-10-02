# geant4-dev/tests/ctests_integration/test75
**Geant4 Validation application that is meant for testing simulation of gamma-nuclear interactions.**

<!--Project desription-->
## Project description
The application simulates gamma+A interactions, and subsequently one looks 
at the double differential spectra of secondary proton, pi-, or pi+, that
are compared vs available experimental data.

Currently modeled are the following cases:

- 300 MeV gamma on Cu target : Spectra of secondary protons, at theta=45, 90, or 135 degree, 
as a function of proton's kinetic energy.
Experimental data from: R. Schumacher et al., Phys. Rev. C 25, 2269 (1982)

- 668 MeV gamma on Cu or Pb target : In the case of Cu target, we look at the spectra of secondary 
pi+/-, at theta=28.4 or 44.2 degree, as a function of pion momentum.
In the case of Pb target, we look at the spectra of pi-/+ at 44.2 degree, as a function of pion momentum.
Experimental data from: K. Baba et al., Nucl. Phys. A322, 349 (1979)
   
   
<!--Authors and contacts-->
## Authors and contacts
- Julia Yarba (Fermilab) - yarba_j@fnal.gov

<!--Geant Val integration-->
## Geant Val integration
Geant Val [![Website shields.io](https://img.shields.io/website-up-down-green-red/http/shields.io.svg)](https://geant-val.cern.ch/) is the Geant4 testing and validation suite. It is a project hosted on gitlab.cern.ch [![Website shields.io](https://img.shields.io/website-up-down-green-red/http/shields.io.svg)](https://gitlab.cern.ch/GeantValidation) used to facilitate the maintenance and validation of Geant4 applications, referred to as <em> tests</em>.\
The following are instructions to use geant4-dev/tests/ctests_integration/test75 within Geant Val, from batch submission to website deployment.
These instructions have been tested on lxplus.cern.ch as well as on geant4gpvm-test-al9.fnal.gov

1. First of all, one needs to git clone the Geant Val/geant-config-generator as well as geant4-dev/tests/ctests_integration/test23 (common use software) and geant4-dev/tests/ctests_integration/test75 (gamma-nuclear test). 
In order to check out only tests/test23 and tests/test75 directories without checking out the whole content of the geant4-dev repository one needs to use the sparse-checkout feature which is available starting git 2.25.x series.
On most modern resources proper version of git is availble by default, e.g. on lxplus,cern.ch it is 2.39.3.
   ```sh
   git clone ssh://git@gitlab.cern.ch:7999/GeantValidation/geant-config-generator.git
   git clone ssh://git@gitlab.cern.ch:7999/geant4/geant4-dev.git --no-checkout 
   cd geant4-dev
   git sparse-checkout init --cone
   git sparse-checkout set tests/test23 tests/test75 # etc.
   git checkout
   cd ..
   ```
2. Copy geant4-dev/tests/ctests_integration/test75/geantval_scripts into geant-config-generator/tests/ctests_integration/geant4/
   ```sh
   cp -r geant4-dev/tests/ctests_integration/test75/geantval_scripts/test75 geant-config-generator/tests/ctests_integration/geant4/
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
   mkdir test75-build
   cd test75-build
   cmake3 -DCMAKE_PREFIX_PATH=$G4INSTALL -DCMAKE_CXX_COMPILER=g++ ../test75
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
   then insprct the ROOTSYS environment variable which will point to the right release of ROOT. \
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
   export PATH=</path/to>/test75-build:$PATH   
   ```

5. Create macros and metadata for Geant Val execution
   ```sh
   python mc-config-generator.py submit -t test75 -d OUTPUT -v 11.2.p01 -q "testmatch" -r
   ```
   this command creates the Geant Val files for batch submission using HTCondor under the ```OUTPUT``` folder, using test75, Geant4.11.2.p01 and the ```testmatch``` job flavour.

6. Execute the analysis on the ROOT files in the ```OUTPUT``` folder to create Geant Val JSON output files
    ```sh
    python mc-config-generator.py parse -t test75 -d OUTPUT 
    ```
    the analysis is coded in ```tests/geant4/test75/parser.py```. The ```OUTPUTJSON``` folder is created with the corresponding JSON files.
    
7. Deploy the results on Geant4 Val. The layout of test75 on Geant4 Val is defined by ```test75.xml``` which can be found in the ```gitlab.com/thegriglat/geant-val-layouts``` repository. \
   Deploy JSON files on the Geant Val database
   ```sh
    find . -name '*.json' | while read i; do curl -H "Content-Type: application/json" -H "token: askauthor" --data @$i https://geant-val.cern.ch/upload; echo; done
   ```
 


