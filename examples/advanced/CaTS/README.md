![alt text](https://github.com/hanswenzel/CaTS/blob/master/images/CaTS.png)
# CaTS: Calorimeter and Tracker Simulation

CaTS is a flexible and extend-able framework (based on geant4 and ROOT)
for the simulation of calorimeter and tracking detectors. 
It also serves as an Example that demonstrates how to use opticks
from within Geant4 for the creation and propagation of optical photons.
https://bitbucket.org/simoncblyth/opticks.git.

The components of CaTS are:


Detector Description:      described in gdml input file (e.g. crystalcal.gdml)
(Geometry, Materials,
 optical properties,
 sensitive detector)

Input modules:                

    GPS
    Particle Gun
    HEPMC (Pythia)


Physics Lists:                  choice of all Reference Physics Lists
                                          optical physics processes (Cerenkov, Rayleigh,
                                          Scintillation etc.) are added (talk to)  
                
Sensitive Detectors:        (+ corresponding Hit classes)        

    TrackerSD(Hit)  registers the step points of charged particles.
    CalorimeterSD(Hit) registering energy deposit
    DRCalorimeterSD(Hit) besides registering energy deposit counts produced Cerenkov photons
    MsCSD(Hit) used to study multiple scattering on a thin layer
    PhotonSD(Hit) sensitive detector that registers optical photons.
    lArTPCSD(Hit) sensitive detector that registers ionization and collects Gensteps (Scintillation and Cerenkov) to be processed by Opticks.



    RadiatorSD no Hits just collection of Gensteps  (Scintillation and Cerenkov) to be processed by Opticks.

## Opticks option (default engine)

This requires the opticks environment to be set up properly, see https://github.com/hanswenzel/CaTS/blob/master/Instructions.md

To get started:

```bash
git clone https://github.com/hanswenzel/CaTS.git
cd CaTS/

# to setup the opticks environment source the setup_opticks.sh file described in:
# https://github.com/hanswenzel/CaTS/blob/master/Instructions.md


source (path to opticks WORK_DIR)/setup_opticks.sh 
cd ../
mkdir CaTS-build
cd CaTS-build

cmake -DCMAKE_BUILD_TYPE=Debug \
  -DWITH_ENGINE=opticks \
  -DCMAKE_PREFIX_PATH="${LOCAL_BASE}/opticks/externals;${LOCAL_BASE}/opticks" \
  -DOPTICKS_PREFIX=${LOCAL_BASE}/opticks \
  -DCMAKE_MODULE_PATH=${OPTICKS_HOME}/cmake/Modules \
  -DCMAKE_INSTALL_PREFIX=../CaTS-install \
  ../CaTS

make install
cd ../CaTS-install/bin
source OpticksEnv # Environment Variables for Opticks
time ./CaTS -g  simpleLArTPC_Opticks.gdml -pl 'FTFP_BERT+OPTICAL+STEPLIMIT'  -m time.mac

#Only the -g command line variable is mandatory! If you don't specify the macro file interactive mode is assumed:

./CaTS -g simpleLArTPC_Opticks.gdml -pl 'FTFP_BERT+OPTICAL+STEPLIMIT' 

```

## Celeritas option

Celeritas is a Monte Carlo particle transport code for simulating High Energy Physics (HEP) detectors on general purpose GPUs. See https://celeritas-project.github.io/celeritas/user/index.html for more details.

In this context, Celeritas is an alternative option to the original Opticks + Nvidia's OptiX. It has been integrated into this CaTS example as of v11.5.  

To build and install the Celeritas package, follow the instructions described in https://celeritas-project.github.io/celeritas/user/usage/installation.html

Most Celeritas environment variables can be used here to modify the default behavior.  In particular, `CELER_DEVICE_DISABLE=1` will force Celeritas to run on CPUs only, while `CELER_DISABLE=1` will disable Celeritas, using purely Geant4 optical processes.

Once Celeritas package is installed, to compile CaTS with Celeritas, do:

```bash
cd examples/advanced/CaTS

#.. please ensure that Geant4 and ROOT environments are properly setup

export Celeritas_DIR=/path/to/celeritas/v07.01
mkdir build
cd build

cmake -DCMAKE_BUILD_TYPE=Release \
   -DCMAKE_CUDA_HOST_COMPILER=${CXX} \
   -DCMAKE_CXX_STANDARD=17 \
   -DCMAKE_CXX_FLAGS="-Wall -Wextra -pedantic" \
   -DMPI_C_COMPILER=${MPICC} \
   -DMPI_CXX_COMPILER=${MPICXX} \
   -DGeant4_DIR=${G4INSTALL}/lib64/cmake/Geant4 \
   -DWITH_ENGINE=celeritas \
   -DWITH_ROOT=ON \
   ${CATS_SRC}

#.. run it
export CELER_LOG=info   # or debug, warning, error
export CELER_LOG_LOCAL=warning  # from one of the GPU threads
export CELER_DEVICE_DISABLE=0   # 1: run Celeritas on CPU only
export CELER_DISABLE=0          # 1: equivalent to Geant4 only

time ./CaTS -g ../gdml/simpleLArTPC.gdml \
    -pl 'FTFP_BERT+OPTICAL+STEPLIMIT' \
    -m ../macros/test_celer.mac
```

## For Geant4 only

To compile CaTS with Geant4 only (without Opticks or Celeritas), e.g. for validation purposes, do:

```bash
cmake -DCMAKE_BUILD_TYPE=Debug \
     -DWITH_ENGINE=geant4 \
    -DCMAKE_MODULE_PATH="../CaTS/cmake/Modules" \
    -DCMAKE_INSTALL_PREFIX=../CaTS-install \
    ${CATS_SRC}
```

If the -pl argument is not provided at runtime, the default physics list configuration 'FTFP_BERT+OPTICAL+STEPLIMIT' is used.

The `WITH_ENGINE` CMake option selects the simulation engine. It accepts
`Opticks`, `Celeritas`, or `Geant4`; values are case-insensitive and are
normalized to lowercase. The default is `Geant4`.


![alt text](https://github.com/hanswenzel/CaTS/blob/master/images/display.png)
to look at the hit collection and make a few histograms:

    ./readPhotonHits NewHits_point_Run0.root PhotonHistos.root Det

The 3 arguments here are: name of input file (hits), name of output file for the histograms and the logical Volume that sensitive detector (PhotonSD)is attached to.
One can then use root to look at the plots:

```bash
root histos.root
   ------------------------------------------------------------------
  | Welcome to ROOT 6.22/06                        https://root.cern |
  | (c) 1995-2020, The ROOT Team; conception: R. Brun, F. Rademakers |
  | Built for linuxx8664gcc on Dec 13 2020, 13:28:00                 |
  | From tags/v6-22-06@v6-22-06                                      |
  | Try '.help', '.demo', '.license', '.credits', '.quit'/'.q'       |
   ------------------------------------------------------------------

root [0] 
Attaching file histos.root as _file0...
(TFile *) 0x559578f76f60
root [1] TBrowser b
(TBrowser &) Name: Browser Title: ROOT Object Browser
```

![alt text](https://github.com/hanswenzel/CaTS/blob/master/images/position.png)

For comparison one might want to disable Opticks and use Geant4 to generate and propagate optical photons:
    
    time ./CaTS -g  simpleLArTPC_Opticks.gdml  -pl 'FTFP_BERT+OPTICAL+STEPLIMIT'  -m time_G4.mac
