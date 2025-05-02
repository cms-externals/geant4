#!/bin/bash
touch $WORKSPACE/controlfile
echo "Geant4 nightly tag script running...."
date
/bin/pwd

dateFMT='%Y/%m/%d-%H'
#--------- create a tag for nightly
UTC_now=$(date -u +${dateFMT})
UTC_nxt=$(date -u +${dateFMT} --date="next hour")

tokenRepo_source=$(echo ${SourceRepoHttpUrl}             | sed -e"s=https://=https://gitlab:${NightlyTagToken}@=")
tokenRepo_verification=$(echo ${VerificationRepoHttpUrl} | sed -e"s=https://=https://gitlab:${NightlyTagToken}@=")
tokenRepo_benchmarks=$(echo ${BenchmarksRepoHttpUrl}     | sed -e"s=https://=https://gitlab:${NightlyTagToken}@=")

sourceBranch=${SourceBranch:-master}
sourceBranch=$(echo $sourceBranch | awk -F/ 'NF==1 {print $1} NF>1 {print $2}')

verificationBranch=${VerificationSourceBranch:-master}
verificationBranch=$(echo $verificationBranch | awk -F/ 'NF==1 {print $1} NF>1 {print $2}')

benchmarksBranch=${BenchmarksSourceBranch:-master}
benchmarksBranch=$(echo $benchmarksBranch | awk -F/ 'NF==1 {print $1} NF>1 {print $2}')

maxwait=1800 # half an hour is generous. Must beless than an hour, or we will risk to never match UTC_*
wait=30      # wait time between trials

tagRepository()
{
  tokenRepo=$1
  branch=$2
  Repo=$(basename $(dirname $tokenRepo))"/"$(basename $tokenRepo .git)    # ie. geant4/.... 
  echo ""
  echo "Will tag ${Repo} with ${UTC_now}, using branch ${branch}"
  cat << EOF > ghostflow-tag-stage.json
   {
     "$Repo": {
       "$branch": {
         "ref_date_format": "${dateFMT}",
         "reason": "nightly",
         "policy": "keep_topics"
       }
     }
   }
EOF

  curl -XPOST --header "X-Ghostflow-Tag-Stage: $Repo" --data @ghostflow-tag-stage.json https://g4ghostflow-dev.web.cern.ch/geant4-ghostflow-jobs
  sleep 15   # be gentle to ghostflow, and give it a chance to do its work before we check
  rm ghostflow-tag-stage.json
}

tag_avail()
{
  tag=$(git ls-remote $tokenRepo | grep stage | grep ${branch} | grep -e ${UTC_now} -e ${UTC_nxt})
  latest=$(git ls-remote $tokenRepo | grep stage | grep ${branch} | grep -e latest)
  test -n "$latest" && num_latest=$(echo ${latest} | awk '{ print $1 }') || num_latest=1
  test -n "$tag" &&  num_tag=$(echo ${tag} | awk '{ print $1 }') || num_tag=0
  test "$num_latest" == "$num_tag"
  return $rc 
}

checkTagAvail()
{
  tokenRepo=$1
  branch=$2
  maxloop=$((maxwait / wait))
  count=0
  while test $count -lt $maxloop; do
    Repo=$(basename $tokenRepo)
    count=$((count + 1))
    if tag_avail ; then
      echo "tag is available for $(basename $tokenRepo) and branch ${branch}"
      break
    elif test $count -eq $maxloop ; then
      echo "Error: failing to create nightly tag for $(basename $tokenRepo), or taking too long"
      exit 1  
    else
      sleep $wait
    fi     
  done
}

tagRepository ${tokenRepo_source} ${sourceBranch}
echo "sent request to tag $(basename $tokenRepo_source .git)"
checkTagAvail ${tokenRepo_source} ${sourceBranch}

tagRepository ${tokenRepo_verification} ${verificationBranch}
echo "sent request to tag $(basename $tokenRepo_verification .git)"
checkTagAvail ${tokenRepo_verification} ${verificationBranch}

tagRepository ${tokenRepo_benchmarks} ${benchmarksBranch}
echo "sent request to tag $(basename $tokenRepo_benchmarks .git)"
checkTagAvail ${tokenRepo_benchmarks} ${benchmarksBranch}
