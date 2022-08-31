Contributing to Geant4
=====
The authoritative Geant4 Git repository is hosted on the `CERN GitLab`_ instance at
https://gitlab.cern.ch/geant4/geant4-dev, accessible to members of the `Geant4 Collaboration`_.
Contributions to Geant4 from outside the Collaboration are welcome, subject to
successful testing and validation, via pull requests submitted thorough our
GitHub mirror of published releases hosted at https://github.com/Geant4/geant4.

If you want to try out the Workflow documented below, the https://gitlab.cern.ch/geant4/geant4-dev
repository provides a "sandbox" which can be used for experimentation.

If you have any issues with or questions on the workflow, please start an issue
on the Issues Board at https://gitlab.cern.ch/geant4/geant4-dev/issues.

Basic Tools
======
Git
-----

`Git <https://git-scm.com>`_ is the version control system used by Geant4. If it
is not installed on your system already, it is usually available through
package managers for your system (apt, yum, Homebrew, vcpkg) or comes
bundled with IDEs such as XCode or Visual Studio. Failing this, the Git homepage provides
a range of `binary bundles for Linux, MacOS, and Windows <https://git-scm.com/downloads>`_.

If you are new to git or need additional information, a wide range of documentation
and tutorials are available, including:

- The `Atlassian Tutorials <https://www.atlassian.com/git/tutorials>`_
  are a good starting point, and also cover more advanced topics.
- The official `Pro Git book <https://git-scm.com/book/en/v2>`_ offers greater
  depth and breadth, Chapters 2 and 3 covering the basics.

Many other resources are available, with a selected few being:

- `Version Control with Git <http://swcarpentry.github.io/git-novice/>`_ tutorial by the
  `Software Carpentry <https://software-carpentry.org>`_ project
- `Help by GitHub <https://help.github.com>`_
- `Learning Resources by GitHub <https://try.github.io>`_. See in particular the
  `Learn Git Branching <https://learngitbranching.js.org>`_ tutorial on the very
  important concept of branches

Whilst the remainder of this guide will focus on the use of ``git`` on the
command line, `several GUIs are available <https://git-scm.com/downloads/guis/>`_
allowing viewing and/or interaction with local and remote repositories.


GitLab
-----

Geant4 uses CERN's `GitLab`_ instance to host the authoritative repository
for Geant4 and developer forks. Here, "authoritative" means the Collaboration's main repository
containing the current development line of Geant4 and from which approved and
validated releases are made. CERN's `help pages <https://gitlab.cern.ch/help>`_
provide a good guide and overall reference to GitLab.

`GitLab Basics <https://gitlab.cern.ch/help/gitlab-basics/README.md>`_
covers most of the major topics, in particular use of `ssh keys for authentication <https://gitlab.cern.ch/help/gitlab-basics/create-your-ssh-keys.md>`_.
We'll refer to pages under this later in the guide of specific areas. `Full
documentation <https://gitlab.cern.ch/help>`_ on GitLab is available, including links to CERN's ServiceNow pages. The latter can be used to check the
status of the service and to report general issues with GitLab functionality.

GitLab has quite an intuitive interface though, so just going to the main repository and clicking around files, branches or commits is quite instructive.

GitLab's Web GUI also provides an intuitive interface to browse the source code,
branches, commits, and in-progress work.


Coding Guidelines
-----
Contributions to Geant4 should follow the set of guidelines linked below for
style, C++ feature usage, Toolkit use, and CMake usage:

- `Geant4 Coding Guidelines <http://geant4.web.cern.ch/Collaboration/coding_guidelines>`_
- `Geant4 C++ Feature Use Guidelines <http://geant4.web.cern.ch/node/82>`_
- `Geant4 Guide for Toolkit Developers <http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/ForToolkitDeveloper/html/index.html>`_
- `CMake Documentation v3.8 <https://cmake.org/cmake/help/v3.8/>`_

C/C++ Compiler and CMake
-----
System and software requirements, plus instructions for building and testing Geant4
are described in the `Geant4 Installation Guide <http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/InstallationGuide/html/index.html>`_.


Geant4 Git Workflow
=====
The following sections walk through the basic Git Workflow for Geant4. 

Basic Git Configuration
-----
If you are using git **for the first time**, it is important to set your git name
and email address which will be used by git to mark authorship of your commits. Failure
to do this properly may cause problems when you try to push to the CERN GitLab
server. You can check that your name and email address are properly set by running

.. code-block:: console

  $ git config --list

and looking for the ``user.name`` and ``user.email`` entries. If these are not
set, they can be set using the following commands

.. code-block:: console

  $ git config --global user.name "Your Name"
  $ git config --global user.email "your.email@domain.xyz"

So that GitLab can track your contributions, ``user.email`` should be an address
that you have `registered in your GitLab profile <https://gitlab.cern.ch/help/user/profile/index.md#profile-settings>`_

We strongly suggest that you `set up SSH keys on GitLab <https://gitlab.cern.ch/help/ssh/README.md#ssh>`_
to authenticate when using command line ``git``.

We also strongly recommend that you also add the following setting (provided
that you have Git version 1.7.11 or newer):

.. code-block:: console

  $ git config --global push.default simple

There are many other configuration settings for git, and the `Atlassian Git Config Tutorial <https://www.atlassian.com/git/tutorials/setting-up-a-repository/git-config>`_
gives a good overview of the most useful ones.

Note that the above are "global" settings, applying to all Git repositories
on your system. If you need different settings for other projects, then
after cloning your fork of the Geant4 repository (see `Clone the Forked Repository Locally`_),
change directory into it and run:

.. code-block:: console

  $ git config --local user.name "Your Name"
  $ git config --local user.email "your.name@cern.ch"

These settings will then only apply to your clone of Geant4.

The above settings are the minimum required for using Geant4's Git/GitLab
Workflow. Additional setup of git aliases and hooks to assist the workflow
is possible through the use of a helper script available after making a local
clone. Use and behaviour of this is covered in `Clone the Forked Repository Locally`_
below.


Fork the Main Repository
-----
The main Geant4 git repository is https://gitlab.cern.ch/geant4/geant4-dev. Follow that link and then
create your own copy of the repository by pressing the ``Fork`` button on the project’s front page.

If you are presented with a number of *group options* when you start to fork, pick your personal
account. There are some good uses for group forks in sub-domains to work on large
developments, but we won’t cover that workflow here. Your fork will be created on GitLab for your
personal account at the URL ``https://gitlab.cern.ch/[YOUR_CERN_ID]/geant4-dev``.

You now have your own personal copy of the Geant4 code on GitLab that you can use to clone
from, and publish to, for your own developments. You will also be unperturbed by other developers
working on their forked copies, and likewise your work will not affect theirs.

The main Geant4 repository is private to Collaboration members. GitLab will
enforce private access on your fork, and in fact by default it will not be
visible to anyone other than you. When you submit a Merge Request, the changes
you have made in this will be visible to the Collaboration as that is necessary
for testing and integration!

To keep informed about ongoing developments and discussions you should set up
GitLab notifications by visiting the main repository at https://gitlab.cern.ch/geant4/geant4-dev
and clicking on the "Bell" icon in the top right hand corner next to the "Star" and
"Fork" buttons. It's recommended that you set the notification level to
at least "Participate" so they you will receive emails for all Issues/Merge Requests
you submit or comment on, as well as when others mention you on Issues/Merge Requests.


Add the Geant4 Robot
^^^^^
When you submit your developments to the main repository for integration, they will be checked
and tested by an automated Continuous Integration system. To enable this and to keep you informed
of progress, the "robot" that oversees the system requires access to your private fork.

To give the robot access, go to ``Project information -> Members`` on the menu that appears on the GitLab
page for your fork (i.e. ``https://gitlab.cern.ch/[YOUR_CERN_ID]/geant4-dev``). In the ``Invite Member`` form add the `geant4bot <https://gitlab.cern.ch/geant4bot>`_ account 
as a Member with ``Maintainer`` access in ``Select a role`` and click ``Invite``.


Clone the Forked Repository Locally
-----
You now have a GitLab fork of your own, but this is used primarily for sharing your changes
with others. To make a change you need a local copy that you can edit yourself. In git
this is done by *cloning* your fork.

The ``git clone`` command requires a URL pointing to the repository you want to clone, in
this case your fork. If you have set up `ssh keys <https://gitlab.cern.ch/help/ssh/README.md#ssh>`_
for your GitLab account, then the following will take a copy of the repository and check out a working
copy for you.

.. code-block:: console

  $ git clone ssh://git@gitlab.cern.ch:7999/[YOUR_CERN_ID]/geant4-dev.git
  ... output ...

This will give you a Local Clone your fork under a ``geant4-dev`` subdirectory of the directory where you run this command.

Note that GitLab offers different ways to authenticate, so there are a few URLs you could
use to clone the repository. Use the ``KRB5/HTTPS/SSH`` option on your fork’s front page
to view the available URLs (but we recommend SSH for security and ease of use). For consistency
through this guide we will always present URLs in their ``ssh`` form.

Once you have created the Local Clone, you should change directory into it
and run the ``GitUtilities/SetupForDevelopment.sh`` program, e.g.

.. code-block:: console

  $ cd geant4-dev
  $ ls
  ...
  $ ./GitUtilities/SetupForDevelopment.sh

It will prompt you to confirm the User and Email you set up earlier are those
you want to use for Geant4 development. Simply follow the instructions to
confirm this or change as you require. It will then install local aliases (shorthands
for useful git commands) and hooks (basic checks on commits) to provide local
versions of the checks the later Merge Request process will enforce. For details
of these, see `Git Hooks and Aliases Installed by SetupForDevelopment.sh`_ below.


Add the Authoritative Repository as a Remote
-----
The ``HEAD`` of  the ``master`` branch of the **authoritative** Geant4 repository, https://gitlab.cern.ch/geant4/geant4-dev,
always has the latest tested and approved code. Your local developments therefore need to be able to
start from, and track the changes on, this (See also `Create and Manage a Topic Branch`_ below).
We enable this by adding the authoritative repository to our Local Clone as a *remote*:

.. code-block:: console

  $ cd geant4-dev
  $ git remote add upstream ssh://git@gitlab.cern.ch:7999/geant4/geant4-dev.git


Your Local Clone now has two remotes, ``origin`` (your fork) and ``upstream`` (the main repository):

.. code-block:: console

  $ git remote -v show
  origin ssh://git@gitlab.cern.ch:7999/[YOUR_CERN_ID]/geant4-dev.git (fetch)
  origin ssh://git@gitlab.cern.ch:7999/[YOUR_CERN_ID]/geant4-dev.git (push)
  upstream ssh://git@gitlab.cern.ch:7999/geant4/geant4-dev.git (fetch)
  upstream ssh://git@gitlab.cern.ch:7999/geant4/geant4-dev.git (push)

We use the name ``upstream`` purely by convention. You are free use another name
if that makes the workflow easier, simply use that name in place of ``upstream``
when running the tasks documented below.


Develop Code
-----
Create and Manage a Topic Branch
^^^^^
In Git, different lines of development are separated by *branches*. The first step of development
is thus to create a new branch to hold your changes:

.. code-block:: console

  $ git checkout master
  ...
  $ git fetch upstream
  ... information on status and updates ...
  $ git rebase upstream/master
  ... information on rebasing master ...
  $ git push origin master
  $ git checkout -b topic-name master

Why the ``checkout``, ``fetch``, ``rebase``, and ``push`` operations first?
As noted above, your fork and local clone are independent copies of the main Geant4 git repository.
GitLab keeps a record of this to help in book-keeping and making Merge Requests, but
it *does not* automatically update your fork or clone with any changes made to the main
repository such as updates to branches or new tags. These commands ensure that
you are starting from the correct source for your Topic Branch, that your local
repository has the newest commits from the main repository, and that the ``master``
branches of your local clone and remote fork are correctly synchronised with the
main repository before creating the branch.

The ``checkout`` command creates and switches your working copy to a branch called
``topic-name``. Whilst branches can have any name, for clarity you should name your 
branch using the code category where the primary changes are being made, the release version base, 
and the development cycle of the category. For example, ``geometry-V11-00-02`` for a branch making the 2nd set of
changes to code under ``source/geometry``, on top of the ``11.0`` release.
The category codes can be found in the ``History`` file inside the directory where
the changes are being made. Read more on ``History`` files in `Local Development Cycle`_.
If the changes will be across more than one category, do not make separate branches,
rather make all changes on a single branch and choose the category part of the name to
reflect the primary category being changed.

You should have one branch per topic you are working on (hence the term
"Topic Branch") to keep the changes focussed and ease the later Merge Request
process. Git branches are cheap in CPU/Storage, so use them to organise your
development work into focussed tasks.

The last element in the ``checkout`` command is the *start point* for your work.
We use ``upstream/master`` as this points to the current stable development line
in Geant4. That is, all code that has passed Nightly Testing, and has
been approved and merged to the release repository's ``master``. You should always start new
Topic Branches from ``upstream/master``, **unless the development addresses a bug fix
for an existing Geant4 release**. In this case you should
start your branch from the upstream branch for the targeted release.
These are named as follows:

- ``master``: Branch targeting "next" release
- ``geant4-MAJOR-MINOR_patches``: Branch targeting next patch release for version ``MAJOR.MINOR``
   - For example, ``geant4-10-05_patches`` is the release branch from which Topic Branches for
     patches to Geant4 10.5 should be made.

The sections below concentrate on "next" development based off of the ``master``
branch. For the changes that are needed to develop from the patch branch(es),
please see `Preparing and Submitting Patch Merge Requests`_.

Because the ``master`` and ``geant4-MAJOR-MINOR_patches`` branches are used for integrating
code for release, you *should not* make commits directly to them, only ever committing to
Topic Branches. This helps both your work and the upstream integration avoid conflicts and
ease testing.

You can see what branches you are currently on plus others available in your working copy
with the ``git branch`` command:

.. code-block:: console

  $ git branch
  ...
    master
  * my-topic
  ...

You can switch between branches at any point using the ``git checkout`` command

.. code-block:: console

  $ git checkout master
  Switched to branch 'master'
  $ git branch
  ...
  * master
    my-topic
  ...

You can also list branches present on the remote repositories using the ``-a`` flag
to ``git branch``:

.. code-block:: console

  $ git branch -a
  ...
  * master
    my-topic
  ...
  remotes/origin/master
  ...
  remotes/upstream/master

When created, your Topic Branch only exists in your Local Clone. To publish it
on your fork on GitLab use ``git push``:

.. code-block:: console

  $ git push -u origin my-topic-branch

the ``-u`` flag is only required the first time you ``push`` your branch and
associates the local Topic Branch with one of the same name on GitLab (or sets up your local branch to
"track" that on your fork).

If you develop on several machines (e.g. Desktop and Laptop), simply repeat the above process
to create Local Clones on them. If you have created and published a
Topic Branch on your fork, it can be checked out later in any of these
clones using:

.. code-block:: console

  $ git fetch origin
  ... updates
  $ git branch -a
  ...
  * whatever-branch-you-are-currently-on
  ...
  remotes/origin/my-topic-branch
  ...
  $ git checkout -t origin/my-topic-branch

Git will checkout and "track" your Topic Branch, so you can work on and
share changes between as many clones as you require.


Local Development Cycle
^^^^^
Once you are on your Topic Branch, you can start making changes to the code
with your preferred editor/IDE. Whilst working on your Topic Branch, you should
regularly synchronize it with the ``master`` branch of the authoritative Geant4 repository (i.e.
``upstream/master`` as we've set the remotes up) to ensure your work builds on the latest
developments. This is done in your local repository using the `fetch` and `rebase` commands:

.. code-block:: console

  $ git checkout my-topic-branch
  ... ensures we are on the correct branch
  $ git status
  On branch my-topic-branch
  $ git fetch upstream
  ...
  $ git rebase upstream/master

You should run this regularly (**at least once per day**), and also **before** submitting your Topic
Branch as a Merge Request. This minimizes any fixes you may have to do as a result
of conflicts between your changes and those in other Merge Requests. If conflicts do occur, then
they are easy to resolve as described in the `GitHub tutorial on resolving conflicts <https://docs.github.com/en/github/getting-started-with-github/resolving-merge-conflicts-after-a-git-rebase>`__.
There is **no need to delete or restart a branch when conflicts occur**, and you **should not** do
this for branches already submitted as Merge Requests. You can easily introduce additional errors
or accidentally revert existing changes by doing this rather than resolving the conflict.
Whilst resolving conflicts can seem intimidating, the `GitHub guide <https://docs.github.com/en/github/getting-started-with-github/resolving-merge-conflicts-after-a-git-rebase>`__
is very clear and helpful.
Be careful with ``rebase`` after you have pushed your Topic Branch to your fork of
Geant4 as ``rebase`` changes the commit history. We'll discuss how to handle this later.

You should also frequently check that your changes compile and pass needed tests (see `Running Tests Locally`_).
Please consult the main `Geant4 Installation Guide <http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/InstallationGuide/html/index.html>`_
for complete details on configuration, build, and testing options.

Do remember that you have the power of a local git repository at your fingertips
during the development process so:

* You can commit code locally anytime, to make it easy to rollback to a known working point
* You can even do more speculative development on a branch from your main Topic Branch - merge it in if it works, throw it away otherwise

We strongly recommend that you *commit early, commit often* for every coherent,
compilable piece of the overall topic you are working on. Git commits code in two phases:

1. You add changed files to a staging area with git add
2. You commit the staging area with git commit

.. code-block:: console

  # ... assuming you are in the git repository
  # ... always good to check where we are, what's been done
  $ git status
  ... content...

  # ... more detail if we need
  $ git diff

  # ... add files
  $ git add <files we want to add to the staging area>
  $ git status
  ... show status output

``git add`` can take wildcards and it’s possible to add all changed files and
new files automatically with ``git add -A``, **but be extremely careful not to add unwanted
files to the repository such as editor or build temporaries**. Always use ``git status`` to
check what will be/had been staged **before you commit**. Note that ``git status`` also helpfully tells
you what to do if you want to rollback changes to a file or unstage files:

- ``git reset`` will unstage changes for you if you did an add of something in error.
- ``git checkout FILE`` will revert ``FILE`` to the original version.
- ``git reset --hard`` will unstage and discard all local changes, use with caution

For more interactive addition of the changes one can use ``git add -p``, which creates an
opportunity to carefully review the differences one by one before adding the modifications.

Once you're happy with the changes that have been staged, we can commit:

.. code-block:: console

  $ git commit

Git will first run a ``pre-commit`` hook to confirm that:

- Hooks and Aliases installed by ``GitUtilities/SetupForDevelopment.sh`` are up-to-date
- Each file being committed is smaller than 2MB (2e6B) in size
- *Formatting checks will be added in 2022*

If any of these checks fail, git will inform you of the problem and not proceed
further. Use the information provided to identify and resolve the problem before trying the ``commit`` again.
If you edit or modify files during this work, remember to use ``git add`` to make
sure they are staged for the commit.

After the ``pre-commit`` checks are passed, ``git`` will open up an editor (determined by ``git config`` or
the ``EDITOR`` environment variable) to allow you to write a message describing the commit.
It is **vital that you write good commit messages**, both for your own records and
to help your colleagues now and in the future. Think of the commit log as the "experimental
logbook" for Geant4!

Commit messages should be informative yet concise. It is suggested that you write
a short one line explanation, followed by a blank line, then some additional longer
paragraphs explaining the changes. Use less than 72 characters per line - commit
logs get read at the terminal.

Further guidelines are in `How to Write a Git Commit Message <https://chris.beams.io/posts/git-commit/>`_
and `5 Useful Tips For A Better Commit Message <https://robots.thoughtbot.com/5-useful-tips-for-a-better-commit-message>`_.

An example of a good commit message would be:

.. code-block:: console

  Fixed uninitialised value and leak in alignment tool

  It was possible to have the m_alignmentOffset variable
  uninitialised in certain situations outlined in Jira.
  This patch corrects for that and changes the behaviour
  to return StatusCode::FAILURE immediately if the
  setup is inconsistent at initialise.

  The memory leak from failing to delete the pointer returned
  from the foo tool was fixed.

The commit message should say *why* something was changed, not *what* was changed
(so it’s not necessary to list changes by filename). You can also cross-reference
GitLab Issues and Merge Requests using `GitLab Markdown References <https://docs.gitlab.com/ee/user/markdown.html#special-gitlab-references>`_,
which can help to keep track of developments and progress.

After you have written the commit message, Git will run a ``commit-msg`` hook
to check that:

- The message is not empty (*some information must be supplied!*)
- There is at least a title of two informative words, e.g.

  - Bad: "*Fix*" is not empty, but not informative (what was fixed?)
  - Good: "*Fix typos*" is simple, but sufficiently informative
  - Bad: "*Fix bug*" is two words, but not sufficiently informative (what bug, how is it fixed?)

Should these checks fail, simply run ``git commit`` again and write a clearer message.

To assist in the Release Management process, you **must** also update the ``History`` file
for each category your Topic Branches touches with an entry comprising a short high level description
of the purpose of the changes made. Each entry is written using Markdown `Markdown <https://en.wikipedia.org/wiki/Markdown>`__
using a `Level 2 heading <https://www.markdownguide.org/basic-syntax/#headings>`__ for date/author/tag,
and `Unordered Lists <https://www.markdownguide.org/basic-syntax/#unordered-lists>`__ for information on the changes.
Each entry **must** be added to the file in reverse chronological order (so newest at the top) and be formatted
as follows:

.. code-block:: markdown

  ## 2021-12-21 Ben Morgan (tagname-V11-00-01)
  - Entry starts with a [Markdown level 2 heading](https://www.markdownguide.org/basic-syntax/#headings) whose
    title must contain, in order
    - The date in `YYYY-MM-DD`
    - Then name of commiter/author
    - Then the tag name in braces, whose format is `(tagname-VXX-YY-PP)` where
      - `tagname` is the identifier for the Category, e.g. `run`.
      - `XX` is the major version of the preceeding Geant4 release
      - `YY` is the minor release of the preceding Geant4 release
      - `PP` is the change to this category starting at `00` since Geant4 release `XX-YY`
  - A [Markdown unordered list](https://www.markdownguide.org/basic-syntax/#unordered-lists) **must**
    be used to structure the description of the change
    - Other than this requirement, you can use other Markdown such as *emphasis*, 
      [links](https://www.markdownguide.org/basic-syntax/#links), even:

      ```cpp
      auto code = blocks();
      ```

  ## 2021-12-20 Ben Morgan (tagname-V11-00-00)
  - Older entries follow, i.e. must be in reverse chronological order
    - If two changes are made on the same date, just increment the `PP` number
      as appropriate
    - Merge Requests touching the same category on the same day are merged in
      the order they are submitted.

**You must not change the number of blank lines between the "preamble" and first tag entry
as changing this can result in Git merging tag entries out of order.** 

**The History files must never be used as a substitute for writing good commit messages**.
These files exist only to assist the Release Manager in preparing final release notes
for release/patch tags when sufficient information has not been supplied in the
commits or Merge Request Description. **It is also critical that you follow the above
formatting rules to ensure easy extraction of information, and that GitLab will render
the information correctly**.

Git is a distributed version control system, so all of the commits made so far
are only present in your local working copy. If you haven't published your Topic Branch
already on your fork on GitLab, as documented earlier, use ``git push``:

.. code-block:: console

  $ git push -u origin my-topic-branch

the ``-u`` flag is only required the first time you ``push`` your branch and
associates the local Topic Branch with that on your GitLab fork (or sets up your local branch to
"track" that on your fork). Subsequent pushes can be done with a simple ``git push``.

If you have pushed to your fork and then rebased your Topic Branch, the next
push will need to be forced as rebasing changes the Git history. In this case,
you should push using the ``--force`` argument:

.. code-block:: console

  $ git push --force origin my-topic-branch

If you work on your Topic Branch across several Local Clones (e.g. Desktop/Laptop), you
should synchronize these using the same procedure of fetch followed by rebase. For example,
using a Local Clone of your fork on your Desktop, you have made some commits to your Topic Branch
and pushed them to your fork. Later on, you move to your Laptop and continue to work on the
Topic Branch, beginning work by doing:

.. code-block:: console

  $ git branch
  ...
  * my-topic-branch
  ...
  $ git fetch origin
  ... notifies you about new commits on my-topic-branch
  $ git rebase origin/my-topic-branch my-topic-branch

Consistent use of fetch/rebase will help to maintain a consistent linear history on
your Topic Branch. It will not prevent true conflicts however. These should be dealt with
as described in the `GitHub tutorial <https://help.github.com/articles/resolving-merge-conflicts-after-a-git-rebase/>`_.

Running Tests Locally
^^^^^
It is possible to run the automated tests locally, for example to troubleshoot failing tests
during the `Merge Request Review and Testing`_ process.

To enable tests, configure the build by running CMake with the following parameters:

* ``GEANT4_ENABLE_TESTING=ON``
* ``GEANT4_USE_GDML=ON`` (requires Xerces-C++, see `Geant4 software prerequisites
  <https://geant4-userdoc.web.cern.ch/UsersGuides/InstallationGuide/html/gettingstarted.html#gdml-xml-geometry-support/>`_).
* ``GEANT4_INSTALL_DATA=ON`` (or ``GEANT4_INSTALL_DATA_DIR`` set to a directory containing valid data).
* Some tests will fail unless `low-energy cross-sections
  <https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/TrackingAndPhysics/physicsProcess.html#cross-sections-for-low-energy-charged-particle-transport/>`_
  are provided and ``G4PARTICLEHPDATA`` is set accordingly.
* Note that some tests are disabled if they require some optional component (eg ROOT)
  which is not enabled.

After building Geant4 as usual, tests are run using ctest. To run all tests on `N` cores, run

.. code-block:: console

  $ ctest -jN

Making a Merge Request
-----
Once you have developed the Topic Branch to your satisfaction, it needs
to be submitted as a GitLab `Merge Request <https://docs.gitlab.com/ee/user/project/merge_requests/>`_
to the authoritative Geant4 repository for full testing and integration.

Submitting your Topic Branch
^^^^^

1. Ensure your Topic Branch is pushed to your fork on GitLab with all commits you want to submit
2. Go to the GitLab interface **for your fork of Geant4**, i.e. ``https://gitlab.cern.ch/[YOUR_CERN_ID]/geant4-dev``
3. Click on the ``Merge Requests`` link in the panel on the left, and then on the green ``New Merge Request`` button
4. In the ``Source Branch`` column, make sure the ``Select source project`` drop down **points to your
   fork**, i.e. ``[YOUR_CERN_ID]/geant4-dev``, and in ``Select source branch`` pick the Topic Branch you want to submit
5. In the ``Target Branch`` column, **make sure** the ``Select target project`` drop down points to ``geant4/geant4-dev``
   and that ``Select target branch`` is ``master``
6. Click on the ``Compare branches and continue`` button

This last step will take you to the main ``New Merge Request`` form that allows you to declare details about
the Merge Request. There are several pieces of information required that **you must provide accurate details for**:

- ``Title`` must be edited to correspond to the Topic Branch name, e.g. ``geometry-V10-05-02``, and optionally
  include a short description, e.g. ``geometry-V10-05-02: Fix Bug #1234``
- ``Description`` should be filled in with an overview of, and motivation for, the proposed changes,
  including cross-references for any GitLab or Bugzilla Issues addressed

  - We recommend that you use the *same* information entered into each ``History``
    file touched by the Topic Branch
  - This leads to some duplication of information, but is used for familiarity
    pending automatic generation of ``History`` files from ``Description``
- If you can, you **must** set ``Assignee`` to the Coordinator for the primary category touched by your changes
  (not every user is allowed to make that change)
  - If you *are* the Coordinator, you can assign your Deputy or other member of your Working Group!
- You do not need set the ``Milestone``, ``Labels``, or ``Approvers`` entries, these will
  be handled by the Assignee and Shifter
- Double check that ``Source branch`` and ``Target Branch`` are correct
- It is **recommended** to tick the ``Remove source branch when Merge Request is accepted`` box as
  this will simplify cleanup once the Merge Request is merged

Finally, click on the ``Submit Merge Request`` button, which will take you to the review and
integration page for your new Merge Request. This page is used to manage and track progress
of the Topic Branch through the review and test process.


Merge Request Review and Testing
^^^^^

After submitting your Topic Branch as a Merge Request it is tested and reviewed, this
process being managed through its page on GitLab under https://gitlab.cern.ch/geant4/geant4-dev/merge_requests.
This page lists all in progress Merge Requests, so simply find yours by the ``Title``
you entered for it. The page for your Merge Request is basically a thread of comments to which the Shifter,
Geant4 Robot, you, and Geant4 collaborators can post to report and communicate on the status of
testing. GitLab should notify you via email of all new posts and updates, so you can
monitor and respond to updates. Please keep all discussions online in the Merge Request
comments to help the Shifter and Working Group Coordinators in the loop on changes,
issues, and progress.

The stages of review and testing are listed below:

1. The Assignee will review the basic information such as ``Title`` and ``Description``,
   and make requests for changes. They may also ask Working Group Coordinators for
   any additional categories your Topic Branch changes touch to be involved in the review.

2. Formatting/file size/conflict checks

   - Geant4Bot runs these checks and posts comments should issues be found
   - You can use the information supplied to resolve issues by making new
     commits on your Topic Branch and pushing them to your fork

3. Continuous Build and Test against latest ``master``

   - These are launched automatically by the Geant4 Jenkins CI system on Merge Request submission
     and any subsequent commits to the Topic Branch
   - It will report back on success/failure, together with a link to the output and results
     on the `Geant4 CDash`_ Dashboard
   - Should there be failures, use the posted link to the CDash Dashboard to track down the cause of the failure
   - If the failure is caused by changes on your Topic Branch, simply add a commit
     to address the issue in your Local Clone and push to your fork. A *single commit must be
     used to prevent overloading the Jenkins CI workers*.
   - If updates have been made to the ``master`` branch, you should
     rebase your Topic Branch against it as documented in `Local Development Cycle`_ above and push.
   - Geant4Bot and Jenkins will test the updated Topic Branch automatically
   - Should any failures be ambiguous, contact the assigned Shifter, or "mention"
     other collaborators for help using their CERN ID, e.g. "@bmorgan, could you look at the failure here please?".
     When you start typing "@.." in GitLab comments, it will pop up a list of matches to
     help you find the right person!

4. Code Review

   - Optionally managed by Assignee and Shifter
   - You are welcome to ask any other collaborators to review the changes
   - Use the Merge Request comments to keep discussion focussed and
     recorded alongside the commits!

5. Nightly Build and Test against latest ``master`` plus all other open Merge Requests

   - Once your Merge Request successfully passes Continuous testing, and any initial Code
     Review issues have been addressed with the Assignee, the Shifter will "stage" your Merge Request for scheduled Nightly testing
   - The Shifter will contact you if staging is not possible due to conflicts, and here you will need
     to rebase your Topic Branch against the latest ``master`` as documented in `Local Development Cycle`_ above,
     fixing the conflicts, and pushing.
   - During this phase, you *must avoid making additional commits to your Topic Branch* other
     than conflict resolution, as changes will result in it being **automatically unstaged** from Nightly testing
   - The Geant4 Jenkins CI launches Nightly tests automatically at around 01:00 Geneva time and will report back
     success/failure, together with a link to the output and results on the `Geant4 CDash`_ Dashboard,
     as a comment on the Merge Request page
   - The Shifter will review any failures with you, and unstage the Topic Branch if traced
     to its changes
   - You can then make a commit on the Topic Branch locally to fix the issues, and force push
     to your GitLab fork to relaunch the *Continuous* tests. As before, you should address the fix in
     a single commit to avoid overload of Jenkins.
   - After successful Continuous testing, the Shifter will restage the Topic Branch for Nightly Testing

6. Acceptance and merge to ``master``

   - Once your Merge Request has passed Nightly Testing, the Shifter will merge the Topic Branch
     into the ``master`` branch if no conflicts occur.
   - If there are conflicts, the Shifter will contact you to request rebasing and push of your
     Topic Branch against the latest ``master`` as documented in `Local Development Cycle`_ above.
   - You do not need to do anything further in GitLab, but you should tidy up your Local Clone(s)
     and fork as documented below in `Finishing Up`_


It's important to note that steps 2-5 are iterative and that you can update your Topic Branch
to address any issues identified during each test phase. You should not close the Merge
Request or delete the associated Topic Branch in response to failures. The workflow
aims to gradually tune and improve the code responsively without requiring a start from scratch.
Generally, Merge Requests are *only* closed without a merge to the ``master`` branch if the work reaches a dead-end or
a better solution to the topic is identified.


Finishing Up
^^^^^

After your Merge Request is accepted and merged to the authoritative ``master`` branch,
you should delete the Topic Branch from your fork and Local Clone, and update them with
the new ``master``.

Whilst deleting the Topic Branch is not required, it avoids clutter and clarifies your active work.
No information is lost, as the commit history is part of ``master`` after merging.
If you ticked the ``Remove source branch when Merge Request is accepted`` box when making
the Merge Request, GitLab will automatically delete the Topic Branch from your fork.
Otherwise, go to the GitLab page for your fork and click on ``Repository -> Branches``
from the left hand pane. This will bring up a list of branches present in your fork,
and you can delete any you no longer require using the red trashbin icon on the right
hand side.

To delete the Topic Branch in your Local Clone, first fetch and rebase the latest
changes from the authoritative ``master``:

.. code-block:: console

  $ git checkout master
  $ git fetch upstream
  ...
  $ git rebase upstream/master master

To syncronize your Local Clone and fork, run a fetch followed by push:

.. code-block:: console

  $ git fetch origin --prune
  ...
  $ git push origin master
  ...
  $ git push origin --tags

Here, the fetch is used with ``--prune`` to remove any local references to
branches deleted on your fork. We do two push operations, one for the ``master``
branch, and one for any new tags (for example, reference, patch or release).
Even if you `Automatically Synchronize your Fork with the Main Repository`_, it
is worthwhile to follow the procedure above to guarantee syncronization in
your Local Clone.

After syncronizing the ``master`` branch, the local Topic Branch can be deleted:

.. code-block:: console

  $ git branch --delete <NAME_OF_TOPIC_BRANCH>

If you need to delete a Topic Branch that has not been merged to ``master``,
you may need to force deletion:

.. code-block:: console

  $ git branch --delete --force <NAME_OF_TOPIC_BRANCH>

Make absolutely sure you want the branch deleted before running either
delete operation!



Draft Merge Requests
-----

The workflow documented above covers the most common use case of submitting
a Topic Branch as a Merge Request at the point development is felt to be complete.
In certain cases, it may be more useful to create a so-called "Draft" (or "WIP") Merge Request,
where an empty or partially complete Topic Branch is submitted. These can be useful for
work such as testing features on platforms you may not have access to (e.g. Windows),
or to identify a regression and then fix it using the testing results.

Other than the point in the development work at which the Topic Branch is submitted
as a Merge Request, the workflow is identical to that described in `Merge Request Review and Testing`_.
However, to ensure that the Shifter and Working Group Coordinators are aware of the
intent of the Merge Request, you **must** prefix the ``Title`` of the Merge Request
with ``Draft:``, e.g. ``Draft: run-V10-05-04: Resolve regression reported in Bug XYZ``.
This helps the Shifter to avoid staging or merging the Topic Branch before work is complete.
During the "Work in Progress" phase, only Continuous testing will be run on the
Merge Request.

Once you feel development is complete, communicate with the Shifter on this, and
either click on the "Mark as Ready" button next to the title, or edit the ``Title`` to remove the ``Draft:`` prefix. 
The Merge Request will then proceed through the additional Nightly test and review and cycle.


Preparing and Submitting Patch Merge Requests
-----

Developing patches for previous releases of Geant4 follows an almost identical workflow
to that for the "next" release, only differing in:

1. Where Topic Branches are started from
2. The Target branch for Merge Requests

As described in `Create and Manage a Topic Branch`_, the ``master`` branch
targets the "next" release, whereas the branches named ``geant4-MAJOR-MINOR_patches``
target patches to the ``MAJOR.MINOR`` release. Thus to start a Topic Branch to
make the sixth patch to the ``run`` category of Geant4 10.4 we would do

.. code-block:: console

  $ git fetch upstream
  ...
  $ git checkout -b run-V10-03-06 upstream/geant4-10-04_patches

Note the change to the Topic Branch naming that reflects the release series targeted
by the patch. Though the patch branches are not as frequently updated
as the ``master``, you must still rebase your Topic Branch regularly:

.. code-block:: console

  $ git checkout run-V10-03-06
  ...
  $ git fetch upstream
  ...
  $ git rebase upstream/geant4-10-04_patches

Be careful at the rebase step to pick the correct branch being targeted!
Development of the patch and publication on your fork proceed exactly as documented
in `Local Development Cycle`_ other than the change in upstream branch name.

Submitting the patch Topic Branch follows the same workflow as described in
`Making a Merge Request`_, except for the change in the branch your Merge Request will target. The only
change you need to make in the procedure outlined in `Submitting your Topic Branch`_
is in Step 5, where the ``Target branch`` **must** point to ``geant4/geant4-dev`` and
``geant4-MAJOR-MINOR_patches`` as appropriate for the version you are submitting the patch for.


Reference
=====
Quick Start
-----

Creating a Fork of ``geant4-dev``
^^^^^
1. Go to https://gitlab.cern.ch/geant4/geant4-dev
2. Click on ``Fork``
3. Select your personal GitLab namespace to fork into, if required
4. Your fork will now be created at ``https://gitlab.cern.ch/[YOUR_CERN_ID]/geant4-dev``

Creating a Local Clone
^^^^^
.. code-block:: console

  $ git clone ssh://git@gitlab.cern.ch:7999/[YOUR_CERN_ID]/geant4-dev.git
  ...
  $ cd geant4-dev
  $ git remote add upstream ssh://git@gitlab.cern.ch:7999/geant4/geant4-dev.git
  $ git fetch upstream
  $ git rebase upstream/master
  $ ./GitUtilities/SetupForDevelopment.sh

Starting a Topic Branch
^^^^^
.. code-block:: console

  $ cd geant4-dev
  $ git checkout master
  $ git fetch upstream
  $ git rebase upstream/master
  $ git push origin master
  $ git checkout -b name-of-topic upstream/master

Developing a Topic Branch Locally
^^^^^
.. code-block:: console

  $ git checkout name-of-topic
  $ git fetch upstream
  $ git rebase upstream/master
  ...
  $ ... edit code ...
  $ ... configure, build, test ...
  $ ... fix issues ...
  ...
  $ git add <list of files changed>
  $ git commit
  ...
  ... repeat edit/test/fix cycle
  ... repeat add/commit
  ...
  ... update History file
  $ git add path/to/History
  $ git commit


Publishing a Topic Branch on Your Fork
^^^^^
.. code-block:: console

  $ git checkout name-of-topic
  ... publishing for first time
  $ git push -u origin name-of-topic
  ...
  ... after initial developments/rebases
  $ git fetch upstream
  $ git rebase upstream/master
  $ git push --force origin name-of-topic

Submitting a Merge Request to Geant4
^^^^^
1. Ensure your Topic Branch is published and up to date on your fork
2. Go to ``https://gitlab.cern.ch/[YOUR_CERN_ID]/geant4-dev/merge_requests/new``
3. Set Source branch to your fork and the Topic Branch being submitted
4. Set Target branch to ``geant4/geant4-dev`` and ``master``
5. Fill out form as described in `Submitting your Topic Branch`_

Updating a Merge Request
^^^^^
Just follow the same steps as in `Developing a Topic Branch Locally`_ and
`Publishing a Topic Branch on Your Fork`_. GitLab watches your
Topic Branch for changes and will update the Merge Request accordingly

Cleaning up a Topic Branch
^^^^^
.. code-block:: console

  $ git checkout master
  ...
  $ git fetch upstream
  ...
  $ git rebase upstream/master
  ...
  $ git fetch origin --prune
  ...
  $ git push origin master
  ...
  $ git push origin --tags
  ...
  $ git branch --delete name-of-topic
  ...
  ... if deleteing an un-merged Topic Branch
  $ git branch --delete --force name-of-topic


Additional Git Notes
-----
Git Hooks and Aliases Installed by SetupForDevelopment.sh
^^^^^

Running the ``GitUtilities/SetupForDevelopment.sh`` program configures your
Local Clone of ``geant4-dev`` with a few basic aliases and hooks to assist
with the Geant4 Git workflow. After the initial run, the program monitors itself
for updates during the commit process and will prompt you to rerun it when a
new version is available.

Aliases for ``git`` work much like Shell aliases in providing shorthands for
commonly used, or lengthy, base commands. In ``git``, aliases are used just
like standard git commands such as ``add``, ``commit``, etc. ``SetupForDevelopment.sh``
creates the following aliases for use with the documented workflow:

- ``g4.log-tree``

  - Prints commit log of current branch in oneline, decorated format, showing branch structure

- ``g4.log-stat``

  - Prints commit log of current branch in oneline, decorated format, additionally
    printing which files the commit modified together with counts of lines added/removed

- ``g4.log-merges``

  - Prints oneline format log of *merge commits* on the local ``master`` branch
  - Thus provides a log of all *Merge Requests integrated onto master*
  - For accurate reporting, ensure your Local Clone's ``master`` branch is regularly
    rebased against that of the authoritative ``geant4-dev`` repository on GitLab.

- ``g4.show-origin``

  - Reports currently available branches on your Local Clone's ``origin`` (i.e. your fork)

Note that all aliases are prefixed with ``g4.`` to avoid clashing with any
you may have setup of the same basename. Most systems with command completion
support Git aliases, so typing ``git g4<TAB>`` (Unix/Bash here) should give a list
of available ``g4`` aliases. Requests for additional aliases are welcome, subject to
review of general usefulness.

``SetupForDevelopment.sh`` also installs two client side hooks to assist in the
commit process:

- ``pre-commit``: Runs immediately on ``git commit`` and checks:

  - If ``SetupForDevelopment.sh`` is out-of-date, the hook will fail
    and ask you to rerun ``SetupForDevelopment.sh`` to update the installed
    aliases and hooks
  - If the commit contains a file more than 2MB (2e6B) in size, the hook
    will fail and report which file is oversize

- ``commit-msg``: Runs immediately after the commit message has been written and checks
  that the message:

  - Is not empty
  - Is not pure whitespace
  - Has at least two words

  If these conditions are not met, the hook will fail and report which condition has
  not been met

As these checks are enforced in the authoritative repository by GitLab and Geant4Robot,
they will ensure your Topic Branch and subsequent Merge Request will pass these checks!
As with aliases, requests for additional checks are welcome subject to review of usefulness.

Git and Case-Insensitive Filesystems
^^^^^
If you are running on a system with a case-insensitive filesystem (e.g.
macOS HFS+, APFS, Windows FAT), you may encounter issues if you need to
rename a file, changing only the case of one or more letters, e.g. ``file.cc`` to ``fIle.cc``.

Any need for case-only filename changes should be *extremely rare* if the
Geant4 coding guidelines are followed. If a case-only rename is required,
then you must either:

- Perform the rename on a *case-sensitive* filesystem, such as most Linux filesystems

- Perform the rename in two stages::

  $ git mv fname fname.tmp
  $ git mv fname.tmp Fname
  $ git commit


Synchronizing Forks with Upstream via Mirroring
-----
As noted in the main guide, your fork is an independent copy of the authoritative Geant4 git repository
and is not updated automatically with any changes made such as merges of
accepted Merge Requests to the ``master`` branch or new tags.

The git workflow naturally includes the steps needed to synchronize your fork via fetching
from the authoritative repository to your local clone, rebasing branches, and pushing branches and
tags to your fork. Whilst this is the recommended and least error-prone method
of synchronisation, GitLab provides a mechanism to set up your fork to regularly pull changes from
the authoritative repository. Please note that both GitLab and CERN **do not**
`recommend or support this functionality for fork synchronization <https://cern.service-now.com/service-portal/article.do?n=KB0004419>`_,
but if you want to set it up bearing those caveats in mind, then it is
configured as follows:

1. From the GitLab front page for your fork, go to the menu on the right hand side and
   select ``Settings -> Repository``
2. Click the ``Expand`` button next to ``Mirroring Repositories``
3. Under the ``Mirror a repository`` section:

   - Set ``Git repository URL`` to ``ssh://git@gitlab.cern.ch:7999/geant4/geant4-dev.git``
   - Set ``Mirror direction`` to ``Pull``
   - Ensure you click the ``Detect host keys`` button
   - Set ``Authentication method`` to ``SSH public key``, and take a copy of the
     public SSH key that is generated

     - Go to your GitLab Profile and go to ``SSH Keys``
     - Add the public SSH key that you copied above as a new key, naming it something like "Geant4 Pull Mirror Key"

   - Back on the mirroring page, ``Mirror user`` should be your account
   - Tick the ``Only mirror protected branches`` box
   - Click on ``Mirror repository``

GitLab should now schedule an **hourly** job that mirrors changes to the main repository's
``master`` branch, release branches, and tags, to your fork. GitLab will report the
mirroring, and its status, on the front page of your fork as "Mirrored from ssh://*****@gitlab.cern.ch:7999/geant4/geant4-dev.git".
You *may* see errors reported on, or just after, creating the mirror about "No ECDSA host key... Host key verification failed".
Should this happen, go back to ``Settings -> Repository -> Mirroring Repositories`` and
click on the "Update now" button for the entry you created in ``Mirrored repositories``.
If a failure is *still* reported, wait and recheck it again the next morning as GitLab
seems to require a few repeats to synchronize settings and keys. In the worst case,
you simply delete the mirroring operation by clicking on the read trashbin icon next
to its entry in ``Settings -> Repository -> Mirroring Repositories -> Mirrored Repositories``.

Resetting your Fork's ``master`` branch
-----
If you make commits directly on the ``master`` branch of your local clone your ``master`` 
branch is now *diverged* from that of the authoritative Geant4 repository. If you further
push the ``master`` to your Fork, this also becomes diverged, and any mirroring you have set
up can fail. More importantly, you cannot create reliable Topic Branches from the diverged 
``master`` as the divergence will increase over time and thus they *cannot* be integrated through
Merge Requests (either via conflicts, or rejection by the Shifter). Generally a divergent ``master``
can be identified by running

.. code-block:: console

  $ git fetch upstream
  ... information on fetch ...
  $ git checkout master
  ...
  $ git rebase upstream/master
  ... information on rebasing master ...
  $ git log -n1 --oneline --decorate

If the branch is *not* divergent, then the ``log`` command will show something like:

.. code-block:: console

  $ git log -n1 --oneline --decorate
  2da4779bbe (HEAD -> master, upstream/master) Merge topic 'PIXE'
  
In other words, the tip of the Local Clone's ``master`` is identical to that on the ``upstream`` remote
(which we take here to be the authoritative Geant4 repository). However, if after the rebase it displays,
for example:

.. code-block:: console

  $ git log -n1 --oneline --decorate
  430c36f40d (HEAD -> master) divergent commit

then the ``master`` is now *divergent* as the commit is *not* the same as that on the ``upstream`` remote.
This can be confirmed by running

.. code-block:: console

  $ git fetch upstream
  ...
  $ git log -n1 --oneline --decorate upstream/master
  2da4779bbe (upstream/master) Merge topic 'PIXE'

Whilst you should *never* actively commit to the ``master`` branch on your Local Clone/Fork, mistakes
can happen and so it is important to know how to reset the state of your ``master`` branch. If you have
identified that your Local Clone's ``master`` has diverged, then the first step is to fetch from the
authoritative repository and find the commit at its head

.. code-block:: console

  $ git fetch upstream
  ...
  $ git log -n1 --oneline --decorate upstream/master
  2da4779bbe (upstream/master) Merge topic 'PIXE'

As usual we have assumed that the ``upstream`` remote points to the authoritative Geant4 repository.
In this example we see that the commit ``2da4779bbe`` is at the head, and we can reset our Local Clone's
``master`` to this be doing

.. code-block:: console

  $ git checkout master
  $ git reset --hard 2da4779bbe
  HEAD is now at 2da4779bbe Merge topic 'PIXE'
  $ git log -n1 --oneline --decorate
  2da4779bbe (HEAD -> master, upstream/master) Merge topic 'PIXE'

Note that this operation **will remove any commits you have made to the** ``master``. If these commits
represent work you want to keep (e.g. you forgot to checkout a Topic Branch), then you should move them
onto a branch before the reset. For example

.. code-block:: console

  $ git log -n1 --oneline --decorate
  430c36f40d (HEAD -> master) divergent commit
  $ git checkout -b mywork
  Switched to a new branch 'mywork'
  $ git checkout master
  ... then reset as above ...

The last step is to try pushing the reset ``master`` branch to your Fork with

.. code-block:: console

  $ git push -f origin master

If this succeeds, then everything is fine. However if you see an error message like:

.. code-block:: console

  Total 0 (delta 0), reused 0 (delta 0)
  remote: GitLab: You are not allowed to force push code to a protected branch on this project.
  To ssh://gitlab.cern.ch:7999/[YOUR_CERN_ID]/geant4-dev.git
   ! [remote rejected]       master -> master (pre-receive hook declined)
  error: failed to push some refs to 'ssh://git@gitlab.cern.ch:7999/[YOUR_CERN_ID]/geant4-dev.git'

then your Fork has also diverged. By default, your Fork's ``master`` branch is "protected" to
largely prevent divergences, but it cannot protect against all cases. To resolve the above error
and resynchronize your Fork's ``master``

1. From the GitLab front page for your fork, go to the menu on the left hand side and select ``Settings -> Repository``
2. Find the ``Protected Branches`` section and click ``Expand``
3. In this section, find the ``Branches`` frame and entry for ``master``
4. In the entry for ``master``, click on ``Unprotect`` and then ``OK`` in the popup to confirm
5. Back in your Local Clone, run the ``push`` again and it should now succeed

   .. code-block:: console

     $ git push -f origin master
     Total 0 (delta 0), reused 0 (delta 0)
     To ssh://gitlab.cern.ch:7999/[YOUR_CERN_ID]/geant4-dev.git
     + <badcommit>...2da4779bbe master -> master (forced update)

6. Back on the GitLab ``Protected Branches`` section, find the ``Protect a Branch`` form
7. Select ``master`` in the ``Branch:`` entry
8. Select ``Maintainer`` in ``Allowed to push`` and ``Allowed to merge``
9. Click on ``Protect``
10. Confirm that ``master`` appears in the ``Branches`` section


.. Geant4 References
.. _Geant4 Collaboration: http://geant4.web.cern.ch/Collaboration
.. _Working Groups: http://geant4.web.cern.ch/Collaboration/working_groups

.. Git references
.. _Git: https://git-scm.com
.. _Pro Git Book: https://git-scm.com/book/en/v2
.. _Version Control With Git: http://swcarpentry.github.io/git-novice/

.. GitLab/CDash references
.. _CERN GitLab: https://gitlab.cern.ch/explore
.. _Geant4 CDash: http://cdash.cern.ch/index.php?project=Geant4
