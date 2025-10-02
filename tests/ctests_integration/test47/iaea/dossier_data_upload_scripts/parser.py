# command line parser for legacy validation @ FERMILAB

import argparse

def getargs(usage="./run_tst.py <options>"):
   parser = argparse.ArgumentParser( description = "IAEA data parser", usage=usage )
   parser.add_argument( "--path-to-datafile", action="store", dest="datafile" )
   return parser.parse_args();
