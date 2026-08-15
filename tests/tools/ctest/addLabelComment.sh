#!/usr/bin/bash

# set NightlyTesting label and add comments to staged Mrs.

which python3 && PYTHON=$(which python3) \
              || PYTHON=python

if echo $MODE | grep -qi release ; then 
  echo "Release build: skip checking labels in gitlab"
else
  sourceBranch=${SourceBranch:-master}
  verificationBranch=${VerificationSourceBranch:-master}
  benchmarksBranch=${BenchmarksSourceBranch:-master}

  export Repository=$(basename ${SourceRepoHttpUrl} .git)   
  export Source_dir="geant4"
  export Branch=$(echo $SourceBranch | awk -F/ 'NF==1 {print $1} NF>1 {print $2}')
  ${PYTHON} scripts/tests/tools/ctest/pre-nightly.py

  export Repository=$(basename ${VerificationRepoHttpUrl} .git)
  export Source_dir="geant4/verification"
 export Branch=$(echo $VerificationSourceBranch | awk -F/ 'NF==1 {print $1} NF>1 {print $2}')
  ${PYTHON} scripts/tests/tools/ctest/pre-nightly.py

  export Repository=$(basename ${BenchmarksRepoHttpUrl} .git)
  export Source_dir="geant4/benchmarks"
  export Branch=$(echo $BenchmarksSourceBranch | awk -F/ 'NF==1 {print $1} NF>1 {print $2}')
  ${PYTHON} scripts/tests/tools/ctest/pre-nightly.py
fi
