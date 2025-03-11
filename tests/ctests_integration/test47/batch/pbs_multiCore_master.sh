#!/usr/bin/env bash
#
# determine number of cores per host
coresPerNode=`cat /proc/cpuinfo | grep -c processor`

# 10M events instead of 1M in most cases, due to low xsec
ntotal=10000000
nevents=$(( ntotal/coresPerNode ))

# NOTE: no need to explicitly pass in those input arguments that are given
# to the master via the -v field of qsub - they'll automatically propagate
# to the parallel jobs 
#
# seq $coresPerNode | /usr/local/bin/parallel -j $coresPerNode "${PBS_O_WORKDIR}/${App} ${G4RELEASE} ${target} ${momentum} ${physlist} ${nevents} {}"
#
seq $coresPerNode | /usr/local/bin/parallel -j $coresPerNode "${PBS_O_WORKDIR}/${App} ${nevents} {}"
