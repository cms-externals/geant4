#!/usr/bin/env python

"""
put comments here !

"""

import os
import sys
import parser

npoints = []
angle = []
ekin = []
dekin = []
xsec = []
err = []

tgt_dict = {
   "al" : "Al",
   "c"  : "C",
   "fe" : "Fe",
   "in" : "In"
}
   
def getdatablock(lines):

   del npoints[:]
   del angle[:]
   del ekin[:]
   del dekin[:]
   del xsec[:]
   del err[:]

   lcounter = 0
   ldata = 0
   npt = 0
   
   for line in lines:
      line = line.strip(); # strips the "\n" at teh EoL
      
      if 'END' in line:
         lcounter += 1
	 return lcounter

      if 'double' in line:
         lcounter += 1
	 tokens = line.split()
	 npt = int(tokens[1])
	 npoints.append(int(tokens[1]))
	 # ---> angle.append(float(tokens[2]))
	 angle.append(tokens[2])
	 # print angle
         continue
      
      ldata += 1
      tokens = line.split()
      # ekin.append(float(tokens[0]))
      # dekin.append(float(tokens[1]))
      # xsec.append(float(tokens[2]))
      # err.append(float(tokens[3]))
      ekin.append(tokens[0])
      dekin.append(tokens[1])
      xsec.append(tokens[2])
      err.append(tokens[3])
      lcounter += 1
      
      if ldata == npt:
         return lcounter

if __name__ == "__main__":
   
   args = parser.getargs();   
   print 'Path to datafile: ', args.datafile
   
   dfile = open( args.datafile, 'r' )
   lines = dfile.readlines()   
   
   oname = os.path.basename( os.path.dirname( args.datafile ) )
   
   rea, tgt, energy = oname.split("_")
      
   oname = oname + '.json'
   print 'Output: ', oname
   
   ofile = open( oname, 'w' ) 
   
   # print ' Number of lines in the datafile = ', len(lines) 
   
   print >>ofile, '['
   
   counter = 0
   while ( counter < len(lines) ):
      nrecords = getdatablock( lines[counter:] )
      # if  nrecords > 0:
      if  nrecords > 1:
         print >>ofile, '  {'
         print >>ofile, '    "modtime": null,'
         print >>ofile, '    "accesskw": "public",'
	 print >>ofile, '    "versiontagkw": "NA",'
         print >>ofile, '    "observablekw": "differential cross section dsig/dO dT", '
         print >>ofile, '    "secondarykw": "neutron",'   # make it configurable !!!
         print >>ofile, '    "trid": null,'
         print >>ofile, '    "testkw": null,'
         # print >>ofile, '    "beamkw": "KEK pi2 3.0 Gev Protons",' # make it configurable !!!
         print >>ofile, '    "beamkw": "KEK pi2', float(energy)/1000., 'GeV Protons",' 
         #
	 # NOTE: need to form a string while adding a value, then print it out to the output file
	 #       otherwise print will add extra white spaces arounf the value (tgt.upper(), that is)
	 #
	 # out_tgt = '    "targetkw": "' + tgt.upper() + '",'
	 out_tgt = '    "targetkw": "' + tgt_dict[tgt] + '",'
	 print >>ofile, out_tgt
         print >>ofile, '    "referencekw": "1501600",'            # make it configurable !!!... actually, not; leave as is
         print >>ofile, '    "reactionkw": "particle production",'
         print >>ofile, '    "modelkw": "NA",'
         print >>ofile, '    "scoreskw": "NA",'
         print >>ofile, '    "mctoolkw": "Experiment",'
         print >>ofile, '    "imageblobslnk": 0,'
         print >>ofile, '    "parnames": ["theta [deg] of secondary particle"],'
	 out_parvalues = '    "parvalues": ["' + angle[0] + '"],'
	 print >>ofile, out_parvalues
	 print >>ofile, '    "datatable": {'
	 print >>ofile, '      "dtid": null,'
	 print >>ofile, '      "datatypeskw": "1D Datapoint set with total stat+sys errors in y",'
	 out_title = '      "title": "proton + C, E=' + str(float(energy)/1000.) + 'GeV, theta=' + angle[0] + ' [deg]",'
	 print >>ofile, out_title
	 print >>ofile, '      "axisTitle": [ "Kinetic energy of secondary neutron [MeV]", "d#sigma / dE d#Theta [mb/srad/MeV]" ],'
	 print >>ofile, '      "nbins": [],'
	 out_npoints = '      "npoints": ' + str(len(ekin)) + ','
	 print >>ofile, out_npoints
	 print >>ofile, '      "binMin": [],'
	 print >>ofile, '      "binMax": [],'
	 # first, X (EKin), then Y (XSec)
	 print >>ofile, '      "val": ['
	 for i in range(len(ekin)):
	    print >>ofile, '        ', ekin[i], ','
	 for i in range(len(xsec)-1):
	    print >>ofile, '        ', xsec[i], ','
	 print >>ofile, '        ', xsec[len(xsec)-1]
	 print >>ofile, '      ],'
	 # half-errors (stat+sys) "minus"
	 print >>ofile, '      "errStatMinus": ['
	 for i in range(len(dekin)):
	    print >>ofile, '        ', float(dekin[i])/2., ','
	 for i in range(len(err)-1):
	    print >>ofile, '        ', float(err[i])/2., ','
	 print >>ofile, '        ', float(err[len(err)-1])/2.
	 print >>ofile, '      ],'
	 # half-errors (stat+sys) "plus"
	 print >>ofile, '      "errStatPlus": ['
	 for i in range(len(dekin)):
	    print >>ofile, '        ', float(dekin[i])/2., ','
	 for i in range(len(err)-1):
	    print >>ofile, '        ', float(err[i])/2., ','
	 print >>ofile, '        ', float(err[len(err)-1])/2.
	 print >>ofile, '      ],'
	 print >>ofile, '      "errSysMinus": [],'
	 print >>ofile, '      "errSysPlus": []'
	 print >>ofile, '    }'
         print >>ofile, '  },'
      counter += nrecords
   
   print >>ofile, ']'

   dfile.close()

