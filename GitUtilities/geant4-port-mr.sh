#!/bin/bash
# Port a Merge Request's commits to a single patch commit on
# a new Topic Branch on a different target branch.
# Currently only supports porting MRs targeting `master` to the
# geant4-10-05_patches branch.
#
# usage:
#   geant4-port-mr.sh <IID>
#
# Where <IID> is the Merge Request ID number on geant4/geant4-dev
# - See https://gitlab.cern.ch/geant4/geant4-git-migration/issues/35
# - TODO: Maybe add as local git alias.

# - The remote where we will fetch MR refs from
#   Maybe set as part of the SetupForDevelopment process, and/or
#   take from .git/config
readonly gitlabRemote="upstream"
readonly gitlabRemoteSSH="ssh://git@gitlab.cern.ch:7999/geant4/geant4-dev.git"

if ! git remote | grep -qE "^$gitlabRemote$" ; then
  echo "No git remote named '$gitlabRemoteSSH'"
  exit 1
fi

if [[ "$(git remote get-url $gitlabRemote)" != "$gitlabRemoteSSH" ]] ; then
  echo "Incorrect URL for '$gitlabRemote' remote" >&2
  exit 1
fi

# - The MR IID and refspec to be fetched and cherry-picked
#   Former is a command-line argument, and must be an integer
readonly mrIID="$1"
readonly mrRefspec="refs/merge-requests/$mrIID/head"

if [[ ! $mrIID =~ [1-9][0-9]* ]] ; then
  echo "Invalid merge request IID '${mrIID}'" >&2
  exit 1
fi

# - Branch which MR originally targeted
#   Eventually as a command line argument, or from GitLab API
#   https://docs.gitlab.com/ee/api/merge_requests.html#get-single-mr
readonly oldTargetBranchName="master"
readonly oldTargetBranchBase="$gitlabRemote/$oldTargetBranchName"

# - Branch to target new cherry-picked MR
#   Eventually as a command line argument
readonly newTargetBranchName="geant4-10-05_patches"
readonly newTargetBranchBase="$gitlabRemote/$newTargetBranchBase"

# 1. To have the MR history available, branch and MR refs need to be fetched 
#    from the authoratative repo
if ! git fetch $gitlabRemote $oldTargetBranchName $newTargetBranchName "+$mrRefspec:$mrRefspec" ; then
  echo "Failed to fetch refspec needed to port MR '${mrRefspec}'" >&2
  exit 1
fi

# 2. General MR commit graph:
#
#                                    refs/merge-requests/IID/head
#                                    |
#            X - Y - ... - ... - J - Z
#           /                   /      \
#... - A - mrMergeBase - ... - I  ... - mrMergeCommit - D - ... - M
#                                                                 |
#                                               oldTargetBranchBase
#
# 3. Only have refs for Z and M, want to obtain list of commits 
#    X-Y-...-Z excluding J (merge commit)
# i. Need to get the merge point, commit "mrMergeCommit"
#    If there is no commit, the topic either isn't merged or targetted
#    a different branch (Use GitLab API earlier to pick this up?)
readonly mrMergeCommit=$(git rev-list --ancestry-path $mrRefspec..$oldTargetBranchBase | tail -n1)

if [[ -z "$mrMergeCommit" ]] ; then
  echo "Cannot cherry-pick MR ${mrIID} as it is not merged into ${oldTargetBranchBase}" >&2
  exit 1
fi

# ii. Find commit "mrMergeBase", the "merge base" of the MR and master:
readonly mrMergeBase=$(git merge-base ${mrMergeCommit}^1 $mrRefspec)

# Should always have a merge-base unless soemthing really weird happens
if [[ -z "$mrMergeBase" ]] ; then
  echo "Cannot cherry-pick MR ${mrIID} as it has no merge base on ${oldTargetBranchBase}" >&2
  exit 1
fi

# iii. Finally, we use rev-list to find all commits reachable from Z,
#      excluding any reachable from mrMergeBase and any merge commits
#      TODO: Assumed that merge-commits are from master and thus unrelated
#            to the MR, but this *might* not be so under all circumstances
#            e.g. merge between same topic branch on two clones
#      TODO: Add a review option to get and display all commits in MR!
#      (again, see the rev-list docs):
readonly mrCommitList="$(git rev-list --reverse --no-merges ${mrRefspec} ^${mrMergeBase})"

# Something is really wrong if this is empty
if [[ -z "$mrCommitList" ]] ; then
  echo "Cannot cherry-pick MR ${mrIID} as it has no non-merge commits" >&2
  exit 1
fi

# TODO: Want a "just print log of each commit oneline" option
#       See also the need to check merge commits above.
# <script> --print-history 42

#-----------------------------------------------------------------------
# Creation of MR patch branch
#-----------------------------------------------------------------------
# 4. Preparation of cherry-pick-branch is as follows:
# i.   Start cherry-pick-branch from MR's merge-base (on the oldTargetBranch)
#      - Error if branch exists already
# ii.  Cherry-pick commit list onto cherry-pick-branch without creating commit(s)
# iii. Commit the changes as a single commit (effective squash)
# iv.  Rebase the cherry-pick-branch onto the new Target branch from the MR merge-base
#      See: https://stackoverflow.com/questions/29914052/i-cant-understand-the-behaviour-of-git-rebase-onto
# v.   If there are conflicts, Git should allow normal interactive fixup etc
# vi.  Clean application may still require fixup (History etc).

readonly cherrypickBranch="cherry-pick-MR$mrIID-to-$newTargetBranchName"
if git rev-parse --verify $cherrypickBranch >/dev/null 2>&1 ; then
  echo "Cannot cherry-pick MR ${mrIID} as target branch ${cherrypickBranch} exists" >&2
  exit 1
fi

git checkout -b $cherrypickBranch $mrMergeBase
git cherry-pick -n $mrCommitList
git commit -m "Port Merge-Request ${mrIID} onto ${newTargetBranchBase}"

if git rebase --onto $newTargetBranchBase $mrMergeBase $cherrypickBranch ; then
  echo "Port of Merge-Request $mrIID onto ${newTargetBranchName} succeeded, using commits:"
  git show --oneline -s $mrCommitList
  echo ""
  echo "Applied on $newTargetBranchName as commit:"
  git show --oneline --stat HEAD
  echo ""
  echo "Fixup History files as appropriate for targeted release ${newTargetBranchName}"
  echo "See full patch/diff with 'git show HEAD'"
  echo ""
  git status
else
  echo "Rebase of $mrIID onto $cherrypickBranch failed..."
  echo "Review git's reporting above and fix problems before running 'git rebase --continue'"
  exit 1
fi

