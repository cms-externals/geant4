import sys
import json
import logging
import requests
import time
import os
import pprint
from git import Repo

Repository = os.environ["Repository"]
# for testing...
# Repository="geant4-devops-playground"

GitLabToken = os.environ["Geant4GitLabToken"]
headers = {"Private-Token": GitLabToken}
URL = (
    "https://gitlab.cern.ch/api/v4/projects/geant4%2F" + Repository + "/merge_requests"
)
# URL=os.environ["GITLAB_API"]+os.environ["GITLAB_MR_API"]

Local_Source = os.environ["WORKSPACE"] + "/" + os.environ["Source_dir"]
Branch       = os.environ["Branch"]

NightlyLabel = "Nightly Testing"
baseCommit = ""
tipCommit = ""

# Default to INFO so all current messages are printed
# TODO: Add as a commline arg along with repo and source.
logging.basicConfig(level=logging.INFO)


def get_staged_MRs():
    logging.info("Local_Source : " + Local_Source)
    repo = Repo(Local_Source)

    global baseCommit, tipCommit
    if Branch == "master":
       baseCommit = str(repo.merge_base("refs/remotes/origin/" + Branch, "HEAD")[0])
    else:
       baseCommit = str(repo.merge_base("refs/stage/" + Branch + "/nightly", "HEAD")[0])
    tipCommit = str(repo.rev_parse("HEAD"))
    topics = repo.git.log(baseCommit + "..HEAD")
    logging.info("baseCommit: " + baseCommit + ", tipCommit: " + tipCommit)
    lines = topics.splitlines()
    mr_iids = []
    for line in lines:
        if line.find("Topic-id") > 0:
            [ignore, mr] = line.split(":")
            mr_iids.append(mr.strip())
    return mr_iids


def addNightlyNote(mr_iids):
    logging.info("Adding Note to staged MRs")
    linksToStagedTopics = ""
    for mr in mr_iids:
        linksToStagedTopics = linksToStagedTopics + "!" + mr + " "

    text = "Jenkins started Nightly testing using this Topic on Stage:"
    text += "\n\nStage-Base: " + baseCommit
    text += "\n\nStage-Tip:  " + tipCommit
    text += "\n\nStaged-Topics: " + linksToStagedTopics
    note = {"body": text}

    for mr in mr_iids:
        logging.info("adding a note to MR" + mr)
        r = requests.post(URL + "/" + mr + "/notes", headers=headers, data=note)
        logging.info("return code for adding Note: " + str(r.status_code))
        time.sleep(1)


def removeNightlyLabel(mr_iids):
    logging.info("Removing Nightly Labels from unstaged/merged MRs")
    params = {"labels": NightlyLabel}
    r = requests.get(URL, headers=headers, params=params)
    res = r.json()

    # remove nightly label on all non staged  MRs
    for req in res:
        if str(req["iid"]) not in mr_iids and str(req["target_branch"]) == Branch:
            newlabels = filter((lambda x: x != NightlyLabel), req["labels"])
            update = {"labels": ",".join(newlabels)}
            logging.info("updating labels for %s to %s", str(req["iid"]), str(update))
            r = requests.put(
                URL + "/" + str(req["iid"]), headers=headers, params=params, data=update
            )
            logging.info(
                "return code for removing nightly label: " + str(r.status_code)
            )
            time.sleep(1)


def addNightlyLabel(mr_iids):
    if not mr_iids:
        logging.info("No staged MRs to label")
        return

    logging.info("Adding Nightly Label to staged MRs")
    pp = pprint.PrettyPrinter(indent=4)
    res = requests.get(URL, headers=headers, params={"iids[]": mr_iids}).json()
    for req in res:
        pp.pprint(req)
        newlabels = req["labels"]
        if NightlyLabel not in req["labels"]:
            newlabels.append(NightlyLabel)
            update = {"labels": ",".join(newlabels)}
            r = requests.put(URL + "/" + str(req["iid"]), headers=headers, data=update)
            logging.info("return code to add nightly label: " + str(r.status_code))
            time.sleep(1)


merge_requests = get_staged_MRs()

removeNightlyLabel(merge_requests)

if merge_requests:
   logging.info("Merge requests: " + str(merge_requests))
   addNightlyNote(merge_requests)
   addNightlyLabel(merge_requests)
else:
   logging.info("No staged merge requests")
