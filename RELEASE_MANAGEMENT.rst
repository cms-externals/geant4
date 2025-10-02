Release and Patch Management
=====

Major/Minor Releases
-----
The ``master`` branch provides a rolling reference of the latest
approved (but not physics validated/verified) developments in Geant4.
New releases are made on a yearly basis with a Beta in late June, and
the final physics validated/verified release in the first week of
December.

1. Once all MRs for a release are merged to ``master``:

   - Create a Topic Branch from ``master`` and commit reference outputs, release notes and files
     for migrations (math functions, disclaimer, etc...)
   - Submit as a normal MR with checks and tests as usual
   - Merge to ``master`` once Nightly Testing passing

2. Tag merge commit on ``master`` of the above as a reference tag (e.g. ``geant4-10-05-ref-00``),
   and publish on the authoritative repository e.g.::

   $ git tag -a geant4-10-05-ref-00 -m "Reference tag geant4-10-05-ref-00"
   $ git push upstream geant4-10-05-ref-00

3. Create a new branch from this tag to be used for patches (e.g. ``geant4-10-05_patches``)
   and push it to the authoritative repository, e.g.::

   $ git checkout -b geant4-10-05_patches geant4-10-05-ref-00
   $ git push -u upstream geant4-10-05_patches

4. Create a new branch from this tag for release packaging (e.g. ``geant4-10-05_release``), and
   make a **single** commit (squashing if required) on this branch to

   - Remove unit/system tests
   - Update dates/release-ID in ``CMakeLists.txt``, ``source/GNUmakefile``, and
     ``source/global/management/include/G4Version.hh``

5. Tag this commit on the release branch as the release tag (e.g. ``geant4-10-05``) and
   push it to the authoritative repository, e.g.::

   $ git tag -a geant4-10-05 -m "Geant4 Release 10.5"
   $ git push upstream geant4-10-05

6. Create a release tarball from this tag, e.g.::

   $ git archive -o geant4-10-05.tar.gz geant4-10-05

7. Export the tarball to https://gitlab.cern.ch/geant4/geant4


Patch Releases
-----
Patches to public releases of Geant4 are made periodically, based on the severity
of the issues being fixed and the amount of fixes collected.

A patch to a release **only include bug-fixes**. It **cannot** include
new developments or features (unless strictly required to accomplish a fix),
and **cannot** make backward-incompatible changes to the public API of the release
being patched. See also the Tag and Release Policy document: https://cern.ch/geant4/collaboration/tag_and_release_policy

At each new release, a new separate branch ``geant4-xx-yy_patches``
gets created from the tag of the release and is subsequently used to collect
all fixes to be included in later patch(es). Any patch of a given release ``xx.yy``
is generated from such branch by removing the relevant tests/docs not supposed
to be published and then tagging the tip of the branch.


1. Developers directly propose MRs on the patch branch for a release (e.g. Topic Branch from/submitted against
   ``geant4-10-05_patches``)
2. These MRs are Continuous/Nightly tested and merged to the patch branch by the Release Manager
3. Once the patch branch is frozen for release, release notes are prepared and added
   to both patch branch and Master by the Release Manager
4. Reference outputs are committed to the patch branch through a MR
5. Create a new branch for the patched release packaging (e.g. ``geant4-10-05-patch-01_release``)
   and make a **single** commit (squashing if required) on this branch to

   - Remove unit/system tests
   - Update dates/release-ID in ``CMakeLists.txt``, ``source/GNUmakefile``, and
     ``source/global/management/include/G4Version.hh``

6. Tag this commit on the patch-release branch as the patch-release tag (e.g. ``geant4-10-05-patch-01``)
7. Create a release tarball from this tag, e.g.::

   $ git archive -o geant4-10-05.tar.gz geant4-10-05

8. Export the tarball to https://gitlab.cern.ch/geant4/geant4


Beta Releases
-------------
A Beta for the year's major/minor release is published in late June. Betas
are tagged and published in a similar way to major/minor releases, the difference
being that no patch branch is needed.

1. Once all MRs for the Beta are merged to ``master``:

   - Create a Topic Branch from ``master`` and commit reference outputs, release notes and files
     for migrations (math functions, disclaimer, etc...)
   - Submit as a normal MR with checks and tests as usual
   - Merge to ``master`` once Nightly Testing passing

2. Tag merge commit on ``master`` of above as a reference tag and push to the authoritative repository.
   This is invariably June's reference tag, ``geant4-10-05-ref-06``), e.g.::

   $ git tag -a geant4-10-05-ref-06 -m "Reference tag geant4-10-05-ref-06"
   $ git push geant4-dev geant4-10-05-ref-06

3. Create a new branch from this tag for beta packaging (e.g. ``geant4-10-06-beta-01_release``), and
   make a *single* commit (squashing if required) on this branch to

   - Remove unit/system tests
   - Update dates/release-ID in ``CMakeLists.txt``, ``source/GNUmakefile``, and
     ``source/global/management/include/G4Version.hh``

5. Tag this commit on the beta release branch as the beta tag (e.g. ``geant4-10-06-beta-01``)
   and push it to the authoritative repository::

   $ git tag -a geant4-10-06-beta-01 -m "Beta Release geant4-10-06-beta-01"
   $ git push upstream geant4-10-06-beta-01

6. Create a release tarball from this tag, e.g.::

   $ git archive -o geant4-10-06-beta-01.tar.gz geant4-10-06-beta-01

7. Export the tarball to https://gitlab.cern.ch/geant4/geant4

