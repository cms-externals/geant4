#!/usr/bin/env python

# start windows build, leave WORKDIR as is, pipeline has created a short Path already.


import os
import g4jk_win_helpers as h

cwd = os.getcwd()
print("CWD:  %s", cwd)

h.touch("controlfile")

h.doBuild()
