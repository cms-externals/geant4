Managing Development/Testing Shifts
=====

Code is integrated into Geant4 via developers submitting their changes for
testing through GitLab `Merge Requests <https://gitlab.cern.ch/help/user/project/merge_requests/index.md>`_.
Geant4 Working Group Coordinators undertake weekly shifts to manage and monitor
the progress of Merge Requests from submission to eventual integration.
This document provides a guide for working these shifts, which cover the ongoing
mainline of Geant4 development. Shifts during testing and validation of new Releases
and Patches are undertaken by the Release Manager, and are described in the
`Release Management Guide <./RELEASE_MANAGEMENT.rst>`_. In both cases, it is assumed
that you are familiar with the general git workflow for Geant4 described in the
main `Contribution Guide <./CONTRIBUTING.rst>`_.

If you have any questions about, or during, shifts, please `start an issue <https://gitlab.cern.ch/geant4/geant4-dev/issues>`_.
For any problems encountered with specific Merge Requests during a shift, add a
comment to the Merge Request comment thread and "mention" `@bmorgan <https://gitlab.cern.ch/bmorgan>`_ to
notify them and get help, e.g. "*@bmorgan, there's a problem with <describe it here>, could you take look
please?*". Note that throughout Shifts, you should keep all notes and discussions
on the Merge Request comment threads as these act as the Shift Log that
other Shifters and developers will need to consult.


Registering for Shifts
=====

A Google Calendar is available for shift registration. Please contact `@gunter <https://gitlab.cern.ch/gunter>`_
for sign up information. Shifts are taken in 7 day blocks, with each block starting
Monday 00:00 (Geneva) each week.


Needed Tools and Documentation
=====
The majority of tasks during a Shift will only require a web browser with
access to the internet and the pages:

- https://gitlab.cern.ch/geant4/geant4-dev, the authoritative repository for Geant4,
- http://cdash.cern.ch/index.php?project=Geant4, the CDash dashboard where results from
  Continuous and Nightly testing are posted
- *Optional*: https://epsft-jenkins.cern.ch, the Jenkins site which manages testing jobs.
  Only needed if highly detailed tracing of build/test failures is required.

Full documentation and help on GitLab is available at https://gitlab.cern.ch/help,
with additional CERN-specific help at http://information-technology.web.cern.ch/services/git-service.

The CDash pages are mostly intuitive and the guides below will highlight the
major details needed. If further help is needed, a page is available on `the CDash wiki <https://public.kitware.com/Wiki/CDash:Documentation>`_.

Whilst Shift tasks should not require any use of command line Git or
local build/test, you should be familiar with the tools and workflow for Geant4 described
in the `Contribution Guide <./CONTRIBUTING.rst>`_.


Managing and Monitoring Merge Requests in a Shift
=====
Through the Shift period, developers will make changes to Geant4 on
Topic Branches and submit these to the authoritative repository through
Merge Requests, as documented in the `Contribution Guide <./CONTRIBUTING.rst>`_.
All Open (i.e. In Progress) Merge Requests can be viewed on the *Merge Requests*
page for the authoritative repository:

- https://gitlab.cern.ch/geant4/geant4-dev/merge_requests

The Shifter's main responsibilities throughout the week are to:

1. Steer and monitor open Merge Requests through the stages of integration as documented
   below on a daily basis
2. Guide the Merge Request Submitter through any steps needed to fix
   issues identified by the integration testing
3. Check Bugzilla (https://bugzilla-geant4.kek.jp) for bugs that are still in
   the new, unassigned state two weeks after the initial report.

   - This should be done once each week.
   - To get the list of new problems, follow the link to new problems:

     - https://bugzilla-geant4.kek.jp/buglist.cgi?bug_status=NEW&order=component%2Cbug_id

The major stages of integrating a Merge Request are:

- Submission
- *Optional* Code Review
  - Check higher level aspects like style, usefulness, etc
- Continuous Testing
  - Document any issues raised by Continuous or Review, work with Submitter to fix and repeat until resolved
- Nightly Testing
  - Document any issues raised by Nightly, repeat Continuous and then Nightly until resolved
- Acceptance

All stages are handled through the *Discussion* thread that GitLab creates on
the page for each individual Merge Request. Through posting comments on this thread,
you can manage the Merge Request integration in coordination with two "robot" systems:

- Geant4Bot, which runs Git operations on Merge Requests
- Jenkins, which orchestrates Build/Test operations for Merge Requests

The discussion thread should also be used to log testing failures, along with
discussions with the Submitter and others on needed fixes from testing failures
or code review. In other words, the Merge Request discussion thread forms the
*Shift Log*.

Comments to communicate with Geant4Bot and Jenkins are covered in the sections
below. All text in the Discussion thread should be written using *GitLab Flavoured Markdown*.
Discussions and GFM provide many tools to make review and management easier, see:

- https://gitlab.cern.ch/help/user/discussions/index.md
- https://gitlab.cern.ch/help/user/markdown.md
- https://gitlab.cern.ch/help/user/project/quick_actions.md

for help and examples.



Merge Request Submission
-----

On each new Merge Request, you should carry out the following basic checks:

1. Click on the Merge Request to take you to its page
2. Assign yourself by clicking the ``assign yourself`` in the ``Assignee`` block
   on the top right hand of the page
3. Check that the title of the Merge Request follows the recommended format (``category_name-V10-05-XX``,
   where ``category_name`` is name of the top directory the Merge Request has changes for)
   and any additional information is clear.
4. Check that the Description is clear

   - Check that the proposed changes are coherent, e.g. it is not adding/fixing two or
     more orthogonal features/bugs.

5. Click on the ``Changes`` tab at the top of the ``Discussion`` section and
   check that:

   - The ``History`` files for the categories touched by the Merge Request have been updated
   - The number of files changed is reasonable for the scope of the proposed change
     and does not change files outside the categories listed in the Description. There
     is no hard rule here on a reasonable number, other than it should not be more
     than the number of files in the targeted categories!

If there are issues or doubts with any of the above, add a comment to the Discussion asking
the Submitter to fix any identified problems, e.g.

  @bmorgan could you elaborate the description please? It's not quite clear what
  "I'm sorry I can't do that Dave" means.

Geant4Bot will also run some sanity checks to confirm

1. The Topic Branch has no conflicts with the current ``master`` branch
2. Commits on the Topic Branch do not introduce files larger than 2MB

and will leave a comment on the thread reporting any failures. If there are no
failures, then it will *not* comment, but it will post a report as a GitLab "Pipeline"
job at the top of the Merge Request, showing a green tickmark if all checks have passed.
It's expected that Geant4Bot completes these checks in under 5 minutes even for quite
large Topic Branches. If you want to double check this, post the following comment in the
discussion to run the checks:

  Do: check

This is the standard format for Geant4Bot commands, a line *at the end of the comment* with
"Do: " (note the space is important!) followed by the task to run. In this case, Geant4Bot will report that
it's running the check by posting the ``:robot_face:`` emoji on your comment.


Continuous Testing
-----
Provided the Geant4Bot checks pass, the Merge Request will be automatically queued for
*Continuous Testing* by the Jenkins CI service. This process:

1. Confirms the Merge Request has no conflicts with the current `master` branch
2. Temporarily merges the ``master`` and Merge Request branches
3. Configures, Builds, and runs a core set of Tests on a minimal set of platforms

This step is intended to quickly check that the Merge Request has no glaring
commit, build, or runtime issues.

If Geant4Bot or Jenkins report failures due to Merge Conflicts, this means that the Submitter's
Topic Branch has diverged from the current ``master`` branch. In this case, ask
them to fetch/rebase their Topic Branch on ``master`` per the instructions in the
`Contribution Guide <./CONTRIBUTING.rst>`_ to resolve the conflicts.

If build/test failures are reported, work with the Submitter to resolve them. Jenkins
will post a link to the CDash page for the build results presented in a table
with:

- The results are subdivided in **Configure**, **Build**, and **Test**
- A green box means success; a red box means failure; an orange box means warning (only for "Configure" and "Build")

  - The number in each box indicates the number of errors in that category
  - For example in **Tests**, the number indicates the number of failing tests
  - By clicking on the number appearing in a red box you get the list of failures
  - By clicking then on the Failed Status, you can see the log file, including the error message(s)

The Submitter should use this information to fix the issues locally on their Topic Branch
before pushing the new commits to their fork. Continuous Testing will automatically restart
and report results as described above on any new pushes to the Topic Branch. You can repeat
this process as many times as needed to get the Merge Request passing Continuous Testing.

During this phase you can also undertake general Code Review with the Submitter and
other Collaborators if required. To involve other people in the discussion, simply
mention them in a comment in the discussion thread to notify them, e.g.

  @gcosmo, @bmorgan, this Merge Request covers code in global and some CMake scripting.
  Could you quickly  review the changes please? Tests are passing so all clear there.

As above, if you or others request changes, the Submitter only need add commits to their
Topic Branch in response to these and push for Continuous to test the updates.

Should testing or review identify fixes needing more detailed investigation,
you can put the Merge Request into the "Work In Progress" state. This is a useful
to mark the work as not ready for use/integration via Nightly Testing. To do this,
simply edit the Merge Request Title and prepend `WIP:`.


Nightly Testing
-----
Once you and any other reviewers and the Submitter  are satisified with the development state
of the Merge Request and it is passing Continuous Testing, it requires *staging*
for Nightly Testing. If the Merge Request has been marked as `WIP:`, this should be removed
before staging. This can be done by editting the Merge Request Title to remove the
initial `WIP:` text, or by using the "Resolve WIP status" button.

To stage a Merge Request for Nightly Testing, simply post a comment on
it containing:

  Do: stage

and Geant4Bot will pick this up and run the needed operations. It's recommended
to add a note with this to confirm that you and the reviewers are satisfied with
Content/Continuous Testing state, e.g.

  Looks ready for Nightly Testing, @gcosmo and @bmorgan have also reviewed and o.ked
  content

  Do: stage

Geant4Bot understands commands in longer comments, and as above, it will add ``:robot_face:``
to the comment to tell you it picked that task up. Should errors be encountered,
it will report back on the cause. Any errors are mostly likely due to Merge
Conflicts, i.e a Merge Request "A" being staged touches one or more of the
same files as an existing staged Merge Request "B".

In this case, the best course of action is to delay staging of "A" until
"B" is accepted and merged (here is a good place to use `GFM's cross-referencing markup <https://gitlab.cern.ch/help/user/markdown.md#special-gitlab-references>`_)
Once "B" is merged, the conflict should be fixed
by the Submitter *rebasing* their Topic Branch onto the new ``master`` and pushing. Check
that Continuous Testing still passes for the rebased commits before trying
to restage.

Note that if the Submitter pushes new commits to the Topic Branch of a Merge
Request that has been staged, Geant4Bot will **automatically unstage it**. This
is to guarantee that Nightly Tests are known to use code already passing Continuous.

Nightly Testing is launched by Jenkins around 01:00 (Geneva) and generally completes
by early to mid morning on the same day/timezone. Testing jobs:

1. Check out the current Nightly Testing stage branch
2. Configures, Builds, and runs a **full** set of Tests on a **full** set of platforms

Jenkins will post comments on the staged Merge Requests indicating
success or failure, together with the list of staged Merge Requests. Note that
as Nightly Testing involves a set of Merge Requests, it will post the same information on each.
This guarantees full information at a small cost in repetition. The information will
tell you *which* Merge Requests were staged and thus tested. It cannot however
indicate which Merge Request caused any failures. For this you will need
to review the logs on the CDash dashboard via following the link posted by Jenkins to the
Merge Request page. The Nightly results are presented in the same tabular form
as for Continuous builds described earlier, there are simply more rows in the
table due to the wider set of platforms tested. In addition, more tests on each platform.
However, warnings and errors may be found and viewed just as for Continuous Testing.

If there are warnings or errors, you will need to check the logs provided by CDash
to triage them and identify the Merge Request(s) responsible. Due a Nightly test
combining several Merge Requests, this may take some time to track down. Some
general guidelines are:

- Warnings or errors at the *Configure* stage almost certainly relate to either

  - A Merge Request that has modified files under `cmake/` or any other `CMakeLists.txt`
    or `sources.cmake` files
  - A general infrastructure issue on the testing machine, which should be reported
    to `@gunter <https://gitlab.cern.ch/gunter>`_ and `@bmorgan <https://gitlab.cern.ch/bmorgan>`_

- If a *Build* stage error can be traced to a specific file, see if any of the staged
  Merge Requests modified that file

  - This also applies to any Tests that fail at their build stage

- Runtime failures in the *Test* stage may give hints to the cause of the failure if
  there are Exceptions


If you identify a given Merge Request as responsible for the failures, **do not close it**!
Rather, remove it from Nightly Testing by posting a comment on the Merge Request
to ask Geant4Bot to unstage it:

  Do: unstage

As above, it's recommended to add additional information here to record why the
unstage was done, and to formally notify the submitter e.g.

  Failures in Nightly for test00 (see <you could post a link to CDash here>) traced to this
  Merge Request. @<thesubmmiter> could you take a look at the error log above and fix
  please?

  Do: unstage

After this operation, the Merge Request re-enters the Continuous Testing phase,
where the Submitter should use the CDash information to help in fixing the problems identified in Nightly just as in
the Continuous phase. The Submitter is responsible for fixing the issue with new commits pushed to the Topic Branch,
which will be automatically re-tested in Continuous. Once this passes, with iteration if needed,
you can restage the Merge Request as before, e.g.

  Proposed fixes are passing Continuous, so retry the Nightlies

  Do: stage

Like Continuous Testing, a Merge Request can be iterated through Nightly Testing to
get it to pass.

Once you are happy that a Merge Request is not responsible for any
failures in the *Test* stage of Nightly Testing it is ready for `Acceptance and Merging`_.
This decision is at your discretion, and don't hesitate to cross check with the Submitter and
other collaborators here. You should *not* however accept and merge *any* Merge
Request if there are failures during the *Configure* or *Build* stages of Nightly Testing.
Failures here prevent the later stage(s) running, and hence all staged Merge Requests
will not have completed full test coverage.



Acceptance and Merging
-----
After a Merge Request has been staged and passed Nightly Testing, it is
ready for final integration with the main ``master`` branch. All you need to do
here is ask Geant4Bot to run the operation by posting

  Do: merge

as a comment on the Merge Request to be merged to ``master``. As before you
may want to add additional information, especially if some unrelated failures
Nightly Testing occur, e.g.

  Nightly testing now passing, and review complete so good to go on master.
  Heisenbug reappeared in test00, but not caused by the changes here.

  Do: merge

Geant4Bot **must** be used here as it runs all needed merge, book-keeping and
tidying operations for you.

The bot will report on progress and automatically close the Merge Request on
completion. It will also unstage the now merged Merge Request from the
Nightly Testing stage and rebase this on the new ``master`` branch, so
no action is required on your part here to update any still open Merge Requests.
For example, say the previous Nightly had Merge Requests ``X``, ``Y``, and ``Z``
staged. The commit history in ``geant4-dev`` thus looks like::

               stage
                 |
         X - Y - Z
       /
  V - W
      |
    master

Let's say that some Nightly Tests failed, with the culprit traced to `X` and `Z`.
You're happy that ``Y`` is passing testing, so you merge it. Geant4Bot then performs
the merge and *rebases the *``stage``* branch*, leaving::

               stage
                 |
            X' - Z'
           /
  V - W - Y'
          |
        master

We've put a superscript on the commits to show these are different due to the
rebase, but this is a detail. The key thing is that Geant4Bot has rebuilt the `stage`
branch automatically for you, so you never need to cross-manage the staged state of
individual Merge Requests.


Advanced Topics
=====

Most Shifters, Reviewers, and Submitters will only need the workflow described
in this document and `CONTRIBUTING.rst <CONTRIBUTING.rst>`_. Should testing reveal
issues that are difficult to triage whether from complexity or from reproducibility,
it is possible to obtain the exact code *content* in a Local Clone for more
detailed test and evaluation (NB: this does not yet reproduce the exact *build and runtime
environments*, which will come once CERN's testing infrastructure is fully containerized).


Checking out Merge Requests and Stages Locally
-----
GitLab together with Geant4Bot implement a so-called "Git Hosted Workflow" (or "Ghostflow").
All this means is that Git itself is used as an effective database of what to
test at the Continuous and Nightly phases. This uses Git's object store combined
with namespaced refs to refer to these points in development (we defer
to the excellent Git documentation for discussion of the Git internals used here).

These refs are not available by default from the ``geant4-dev`` repository,
but by adding it as a remote of your Local Clone, it is possible to fetch them.
By doing this, you can checkout both the current state of a given Merge Request,
or the specific code used by the Nightly Tests on a given date. This is an
advanced topic and generally only used to confirm issues or triage complex ones.

You should first confirm that you have added the authoritative Geant4 repository
as a remote of your Local Clone, e.g.::

  $ git remote add upstream ssh://git@gitlab.cern.ch:7999/geant4/geant4-dev.git

we will use the name ``upstream`` here, following usual GitLab/Hub convention.

Merge Requests
^^^^^
To fetch the current state of a Merge Request's Topic Branch, first find its ID number. We can then
get the current commit being tested by doing::

  $ git fetch upstream ""+refs/merge-requests/<ID>/head:refs/merge-requests/<ID>/head"

This purely fetches the objects, but we can then checkout the commit by then::

  $ git checkout refs/merge-requests/<ID>/head
  Note: checking out 'refs/merge-requests/<ID>/head'.

  You are in 'detached HEAD' state. You can look around, make experimental
  changes and commit them, and you can discard any commits you make in this
  state without impacting any branches by performing another checkout.

  If you want to create a new branch to retain commits you create, you may
  do so (now or later) by using -b with the checkout command again. Example:

  git checkout -b <new-branch-name>

  HEAD is now at <locally dependent output>
  $

The "detached HEAD" state is important here to maintain a temporary triage/test
checkout, rather than for ongoing development. The ``git-worktree`` `command <https://git-scm.com/docs/git-worktree>`_
available from Git 2.6 onwards can also be used to create a checkout in a separate
directory without interfering with your own ongoing developments.

Nightly Testing Stages
^^^^^
In Ghostflow, Nightly Testing uses the so-called "stage" which is nothing more than
a Git ref constructed by merging a set of open Merge Requests with the current tip of the
targeted branch (``master`` for ongoing development, ``geant4-MAJOR-MINOR_patches`` for
patches to existing ``MAJOR-MINOR`` releases).
Geant4Bot automates this procedure, as documented above, in response to requests
by the Shifter. When Jenkins starts the Nightly Tests, its first action is to
snapshot the current Git ref for the stage with the date on which the Nightly
is being started.

To fetch the stage used by the Nightly Testing for the targeted branch on a
given date, we use a similar syntax to that for Merge Requests. For example, to fetch the stage for the
``master`` branch on the 31st December, 2018, we would do::

  $ git fetch upstream "+refs/stage/master/nightly/2018/12/31:refs/stage/master/nightly/2018/12/31"

and to check it out::

  $ git checkout refs/stage/master/nightly/2018/12/31
  ...

As with Merge Requests, this is a "detached HEAD" checkout. This is especially
important for Nightly stages, as once created they are not modified and so
guarantee to represent the actual content tested on the given date. You can use
this feature to try and reproduced observed failures, and hence report back
to the Merge Request triaged as responsible for the failure.
If you want to fetch the stage for a patch branch, simply replace ``master``
in the above commands with ``geant4-MAJOR-MINOR_patches`` for the ``MAJOR-MINOR``
release of interest, e.g. ``10-04``.

As above, `git-worktree`` `<https://git-scm.com/docs/git-worktree>`_
can also be used to create a checkout in a separate
directory without interfering with your own ongoing developments.
