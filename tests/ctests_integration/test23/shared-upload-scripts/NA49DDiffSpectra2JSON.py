#!/usr/bin/env python

"""
put comments here !

"""

import urllib

url = "https://spshadrons.web.cern.ch/spshadrons/data/xsect_" # piplus.txt"

url_pH = {
   'piplus'      : 'piplus',
   'piminus'     : 'piminus',
   'proton'      : 'proton_pp',
   'antiproton'  : 'antiproton_pp',
   'kplus'       : 'kaplus_pp',
   'kminus'      : 'kaminus_pp'
}

url_pC = {
   'piplus'      : 'piplus_pC',
   'piminus'     : 'piminus_pC',
   'proton'      : 'proton_pC',
   'antiproton'  : 'antiproton_pC'
}

partname = {
   'proton' : 'proton',
   'antiproton' : 'antiproton',
   'neutron' : 'neutron',
   'piplus' : 'pi+',
   'piminus' : 'pi-',
   'kplus' : 'K+',
   'kminus' : 'K-'
}

partid = {
   'proton'     : '2212',
   'antiproton' : '-2212',
   'neutron'    : '2112',
   'piplus'     : '211',
   'piminus'    : '-211',
   'kplus'      : '321',
   'kminus'     : '-321' 
}

xf = []
pt = []
xsec = []
err = []

def getdatablock(lines):

   del xf[:]
   del pt[:]
   del xsec[:]
   del err[:]
   
   comment = 0
   
   for line in lines:
      line = line.strip(); # strips the "\n" at the EoL
      if '#' in line :
         comment += 1
	 continue
      if len(line) == 0:
	 return comment # xf, pt, xsec, err 
      columns = line.split()
      xf.append(columns[0])
      pt.append(columns[1])
      xsec.append(columns[2])
      err.append(columns[3])
   
   return comment  
    

def dumpspectrum(ofile):

   ofile.write("\n   \"val\":\n   [")   
   for ipt in range(len(pt)):
      ofile.write(pt[ipt] + ", ")   
   for ix in range(len(xsec)-1):
      ofile.write(xsec[ix] + ", ")   
   ofile.write(xsec[-1] + "],")
   
   ofile.write("\n   \"errStatPlus\":\n   [")
   for ipt in range(len(pt)):
      ofile.write("0., ")
   for ier in range(len(err)-1):
      estat = float(err[ier])
      estat *= 0.5
      ofile.write("%f, " % estat)
   estat = float(err[-1])
   estat *= 0.5
   ofile.write("%f]," % estat)
   ofile.write("\n   \"errStatMinus\":\n   [")
   for ipt in range(len(pt)):
      ofile.write("0., ")
   for ier in range(len(err)-1):
      estat = float(err[ier])
      estat *= 0.5
      ofile.write("%f, " % estat)
   estat = float(err[-1])
   estat *= 0.5
   ofile.write("%f]," % estat)
   
   ofile.write("\n   \"errSysPlus\":\n   [")
   for ipt in range(len(pt)):
      ofile.write("0., ")
   for ier in range(len(xsec)-1):
      esys = float(xsec[ier])
      esys *= 0.038
      ofile.write("%f, " % esys)
   esys = float(xsec[-1])
   estat *= 0.038
   ofile.write("%f]," % esys)
   ofile.write("\n   \"errSysMinus\":\n   [")
   for ipt in range(len(pt)):
      ofile.write("0., ")
   for ier in range(len(xsec)-1):
      esys = float(xsec[ier])
      esys *= 0.038
      ofile.write("%f, " % esys)
   esys = float(xsec[-1])
   estat *= 0.038
   ofile.write("%f]," % esys)
   
   ofile.write("\n   \"binMin\":[], \"binmax\":[]")
   
         
def dumpdatablock(ofile,key,tg,tglnk):
      
   ofile.write( "{" )
   ofile.write("\n   \"trid\":1,\"testlnk\":0,\"referencelnk\":53,\"mcdetaillnk\":1," )
   ofile.write("\n   \"beamlnk\":7,\"targetlnk\":" + tglnk + ",\"observablelnk\":9,\"secondarylnk\":" + partid[key] + ",\"reactionlnk\":1," ) 
   ofile.write("\n   \"datatable\":{" )
   ofile.write("\n   \"dtid\":1,\"datatypeslnk\":1000," )
   ofile.write("\n   \"title\":\"Production of " + partname[key] + " in proton-" + tg + " interactions at 158GeV/c, at xF=" 
              + xf[0] + "\",\"npoints\":%d," % len(xf) )
   ofile.write("\n   \"nbins\":[],\"axisTitle\":[\"E(d3sigma/dp3) [mb/(GeV**2/c**3)]\",\"pT [GeV/c]\"],")
   dumpspectrum(ofile)
   ofile.write("\n   },")
   ofile.write("\n   \"imageblobslnk\":NULL,\"scoreslnk\":2,\"accesslnk\":1,")
   ofile.write("\"parnames\":[\"xF of secondary particle\"],\"parvalues\":[\"" + xf[0] + "\"]")
   ofile.write("\n},\n")


if __name__ == "__main__":

   # data = urllib.urlopen( url ).read();
   # print data
   # print '=' * 20

   for key in url_pH.iterkeys():
      lines = urllib.urlopen( url + url_pH[key] + ".txt" ).readlines()      
      counter = 0
      ofile = open( 'NA49-proton-H-' + key + '-ddiff.json', 'w' )
      print >> ofile, "{\"ResultList\":[" # using print instead of write just to see how it works...looks like I don't need \n at the end...      
      while ( counter < len(lines) ):
         comment = getdatablock( lines[counter:] )
         dumpdatablock(ofile,key,"H","1")
         counter += ( len(xf) + comment + 1 )   
      print >> ofile, "\n]}\n"    
      ofile.close()  

   for key in url_pC.iterkeys():
      lines = urllib.urlopen( url + url_pC[key] + ".txt" ).readlines()      
      counter = 0
      ofile = open( 'NA49-proton-C-' + key + '-ddiff.json', 'w' )
      print >> ofile, "{\"ResultList\":[" # using print instead of write just to see how it works...looks like I don't need \n at the end...      
      while ( counter < len(lines) ):
         comment = getdatablock( lines[counter:] )
         dumpdatablock(ofile,key,"C","6")
         counter += ( len(xf) + comment + 1 )   
      print >> ofile, "\n]}\n"    
      ofile.close()  
   

