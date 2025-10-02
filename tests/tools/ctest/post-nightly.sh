#!/bin/bash

cd geant4

# At this point, we're at the tip of Nightly
# 1. Find commit where the stage was applied
#    Ideally uses merge-base, but need to check
stageBaseCommit=`git merge-base refs/remotes/origin/master HEAD`
stageTipCommit=`git rev-parse HEAD`

# 2. Extract logs of all merge commits between there and the stage tip
#    Some of these may be merges within MRs, but at least limits search
#    Newer gits can extract the trailers directly
stagedTopics=$(git log $stageBaseCommit..HEAD | grep "Topic-id:" | cut -d: -f2)

# 5. Report back to Gitlab
for merge_request_iid in $stagedTopics
do
  noteBody="Jenkins completed the Nightly build using this Topic on Stage:\n\nStage-Base: $stageBaseCommit\n\nStage-Tip: $stageTipCommit\n\nStaged-Topics: $linksToStagedTopics"
  echo "{\"body\": \"$noteBody\"}" > payload.json
  cat payload.json
  curl -f -XPOST \
    --header "Private-Token: $GitLabToken" \
    --header "Content-Type: application/json" \
    -d @payload.json \
    $GITLAB_API/$GITLAB_MR_API/$merge_request_iid/notes
done
