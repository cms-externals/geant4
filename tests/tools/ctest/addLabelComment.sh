#!/usr/bin/bash

# set NightlyTesting label and add comments to staged Mrs.

which python3 && PYTHON=$(which python3) \
              || PYTHON=python

if echo $MODE | grep -qi release ; then 
  echo "Release build: skip checking labels in gitlab"
else
  export Repository=$(basename ${SourceRepoHttpUrl} .git)   
  export Source_dir="geant4"
  ${PYTHON} scripts/tests/tools/ctest/pre-nightly.py

  export Repository=$(basename ${VerificationRepoHttpUrl} .git)
  export Source_dir="geant4/verification"
  ${PYTHON} scripts/tests/tools/ctest/pre-nightly.py

  export Repository=$(basename ${BenchmarksRepoHttpUrl} .git)
  export Source_dir="geant4/benchmarks"
  ${PYTHON} scripts/tests/tools/ctest/pre-nightly.py
fi
