#!/usr/bin/env python

from __future__ import print_function


# for environment
import os
import re
from shutil import rmtree
import subprocess as sb

#to get Windows version 
import platform

# from sys import stdout
import sys
from string import whitespace


def printEnv():
    for key in os.environ.keys():
        print("{0:30} {1}".format(key, os.environ[key]))


def touch(fname, times=None):
    with open(fname, "a"):
        os.utime(fname, times)


def getBuildOptionsDirPath():
    # get a list of build options
    dir_bo = ""
    try:
        # split() removes all spaces, join the result, and split at ,
        build_options = "".join(os.environ["BUILDOPTIONS"].split()).split(",")
        for opt in build_options:
            if opt:
                dir_bo = dir_bo + "/" + opt
    except KeyError:
        pass
    except TypeError:
        print("TypeError: BUILDOPTIONS  ", os.environ["BUILDOPTIONS"])
        raise
    return dir_bo


def SetDirToShortPath(skip_dirs, verbose=False):
    path_list = os.getcwd().split(os.sep)
    os.chdir(path_list[0] + os.sep + path_list[1])
    if verbose:
        print("initial dir ", os.getcwd())
    for dir in path_list[2:]:
        if verbose:
            print("checking-1 ", dir)
        dir = dir.replace("@", "_")
        # print "checking-2 %s" % dir
        if dir not in skip_dirs:
            if dir not in os.listdir("."):
                os.mkdir(dir)
            os.chdir(dir)
            if verbose:
                print("CWD:  {0}\n".format(os.getcwd()))
        else:
            if verbose:
                print("skipping : {0}\n".format(dir))

    if verbose:
        print("final CWD:  {0}\n".format(os.getcwd()))


def CreateBuildDir(dir_bo):
    path = "./" + os.environ["MODE"] + "/" + os.environ["BUILDTYPE"]
    if dir_bo:
        path += "/" + dir_bo

    # print("constructed path: {0}\n".format(path))
    try:
        os.makedirs(path)
    except OSError:
        pass
        if not os.path.isdir(path):
            raise

    os.chdir(path)


def deleteOldBuild(dir_list):
    for dir in dir_list:
        if os.path.exists(dir):
            if os.path.isdir(dir):
                print("Delete directory {0}\n".format(dir))
                # rmtree(dir,True)           # True to ignore errors occuring in .svn dir
                sb.call(["rmdir", "/S", "/Q", dir], shell=True)
            else:
                print("Error: {0} exists, but is not a directory\n".format(dir))
                exit(1)


def purge(dir, pattern):
    for f in os.listdir(dir):
        if re.search(pattern, f):
            os.remove(os.path.join(dir, f))


def doBuild():
    sys.stdout.flush()

    arg = (
        os.environ["WORKSPACE"]
        + os.sep
        + "scripts"
        + os.sep
        + "tests"
        + os.sep
        + "tools"
        + os.sep
        + "ctest"
        + os.sep
        + "g4-win" + platform.uname().release + "-"
        + os.environ["COMPILER"].lower()
        + ".bat"
    )
    dircmd = ["dir", arg]
    sb.call(dircmd, shell=True)

    cmd = "call " + arg
    print("cmd={0}\n".format(cmd))
    sys.stdout.flush()
    sb.call(cmd, shell=True)
