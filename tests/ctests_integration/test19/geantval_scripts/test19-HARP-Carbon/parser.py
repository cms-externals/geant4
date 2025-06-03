#!/usr/bin/env python

# class for sc test parsing
import os
from gts.utils import getJSON
from gts.BaseParser import BaseParser
from ROOT import gROOT, TFile, TTree
from ROOT import gDirectory
gROOT.SetBatch(True)
gROOT.ProcessLine("gErrorIgnoreLevel = 5000;")

class Test (BaseParser):
   TEST="test19-HARP-Carbon"

   def parse(self, jobs):

        mydict = { 'FW_0' : '0.05-0.1 rad',
                   'FW_1' : '0.1-0.15 rad',
                   'FW_2' : '0.15-0.2 rad',
                   'FW_3' : '0.2-0.25 rad',
                   'LA_0' : '0.35-0.55 rad',
                   'LA_1' : '0.55-0.75 rad',
                   'LA_2' : '0.75-0.95 rad',
                   'LA_3' : '0.95-1.15 rad',
                   'LA_4' : '1.15-1.35 rad',
                   'LA_5' : '1.35-1.55 rad',
                   'LA_6' : '1.55-1.75 rad',
                   'LA_7' : '1.75-1.95 rad',
                   'LA_8' : '1.95-2.15 rad'
	         }
        inspireid = -1

        # should be just one job for each job
        job = jobs[0]
        print("path = ", job["path"])
	
        material = job["MATERIAL"]
        if 'G4_' in job["MATERIAL"]:
           material = (job["MATERIAL"]).split("_")[-1] # [-1] means counting from right

        energy = job["ENERGY"]
        if 'MeV' in job["ENERGY_UNIT"]:
           en = float(energy) / 1000.
           energy = str(en)
#        if '.' in energy:
#           energy = energy.split(".")[0]

        beam = job["PARTICLE"]
        if 'pi+' in job["PARTICLE"]:
           beam = 'piplus'
        if 'pi-' in job["PARTICLE"]:
           beam = 'piminus'
	
        detector = job["DETECTOR"].lower()
        print("detector = ", detector)

        fname = detector + "-histo/" + beam + material + energy + "GeV" + job["GENERATOR"] + ".root"	
        filename = os.path.join(job["path"], fname)
        print("filename = ", filename)

        hfile = TFile(filename)

        for h in hfile.GetListOfKeys():
           if 'NSec' in h.GetName():
              print("Skip histo : ", h.GetName())
              continue
           if 'proton' in h.GetName():
              print("Skip histo : ", h.GetName())
              continue
           if 'mom_' in h.GetName():
              print("Skip histo : ", h.GetName())
              continue
           if 'proton' in beam:
              if 'FW' in h.GetName():
                 inspireid = 826544
              if 'LA' in h.GetName():
                 inspireid = 786183
           if 'piplus' in beam or 'piminus' in beam:
              if 'FW' in h.GetName():
                 inspireid = 813159
              if 'LA' in h.GetName():
                 inspireid = 825244
           if 'piplus' in h.GetName():
              secondary = "pi+"
           if 'piminus' in h.GetName():
              secondary = "pi-"
           hh = hfile.Get(h.GetName())
           firstNonZeroBin = hh.FindFirstBinAbove(0)
           lastNonZeroBin  = hh.FindLastBinAbove(0)

           nBins = []
           binEdgeLow = []
           binEdgeHigh = []
           binContent = []
           yStatErrorPlus = []
           yStatErrorMinus = []
           beamenergies = []
           beamenergies.append(en)

           for x in range(firstNonZeroBin, lastNonZeroBin+1):
               binContent.append(round(hh.GetBinContent(x),3))
               binEdgeLow.append(round(hh.GetBinLowEdge(x),2)) # it's better/safer to round it to the 2nd digit after decimal point
               binEdgeHigh.append(round((hh.GetBinLowEdge(x)+hh.GetBinWidth(x)),2))
               #
               # NOTE: don't use BinErrorLow/Up as they return BinError !!!
               #
               yStatErrorPlus.append(round(hh.GetBinError(x)/2.,3))
               yStatErrorMinus.append(round(hh.GetBinError(x)/2.,3))

           nBins.append(len(binContent))

           hname = hh.GetName()
           ptag = hname[hname.find("_")+1:]
           parameters = []
           parameters_test = { "names": "THETA", "values": mydict[ptag] }
           parameters.append(parameters_test)

           htitle = job["PARTICLE"] + " + " + material + " -> " + secondary + " + X" 

           observableName = "D2(SIG)/DP/DTHETA"
           xAxisTitle = "p (GeV/c)"
           yAxisTitle = "D2(SIG)/DP/OTHETA (mb/rad/GeV/c)"

#           print("Processing histo : ", h.GetName())
#           print("job[PARTICLE] = ", job["PARTICLE"])
#           print("secondary = ", secondary) 
#           print("parameters_test = ", parameters_test)
#           print("theta = ", mydict[ptag])

           if 'proton' in job["PARTICLE"]:
              if '0.05-0.1' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5.0 ]
                    binContentExp = [ 0.13, 0.16, 0.108, 0.071, 0.032, 0.024, 0.025, 0.008 ]
                    binHalfErrorExp = [ 0.01, 0.01, 0.0055, 0.004, 0.0035, 0.0025, 0.0025, 0.001 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5.0 ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5.0, 6.5 ]
                    binContentExp = [ 0.18, 0.2, 0.25, 0.185, 0.135, 0.06, 0.03, 0.021, 0.004 ]
                    binHalfErrorExp = [  0.01, 0.01, 0.01, 0.007, 0.0065, 0.005, 0.0035, 0.002, 0.0005 ]
              elif '0.1-0.15' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5.0 ]
                    binContentExp = [ 0.2, 0.144, 0.08, 0.05, 0.028, 0.015, 0.009, 0.003 ]
                    binHalfErrorExp = [ 0.015, 0.0075, 0.0045, 0.0035, 0.0025, 0.0015, 0.001, 0.0005 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5.0 ]
                    binContentExp = [ 0.21, 0.22, 0.2, 0.123, 0.049, 0.027, 0.012, 0.004 ]
                    binHalfErrorExp = [ 0.01, 0.01, 0.01, 0.007, 0.0035, 0.003, 0.0015, 0.0005 ]
              elif '0.15-0.2' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [7]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5 ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4. ]
                    binContentExp = [ 0.16, 0.092, 0.059, 0.036, 0.025, 0.008, 0.005  ]
                    binHalfErrorExp = [ 0.01, 0.006, 0.0045, 0.003, 0.0025, 0.001, 0.001 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5.0 ]
                    binContentExp = [ 0.22, 0.19, 0.108, 0.06, 0.032, 0.019, 0.007, 0.002 ]
                    binHalfErrorExp = [ 0.015, 0.01, 0.0055, 0.0045, 0.0025, 0.002, 0.001, 0.0005 ]
              elif '0.2-0.25' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [7]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5 ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4. ]
                    binContentExp = [ 0.13, 0.096, 0.08, 0.024, 0.009, 0.003, 0.001  ]
                    binHalfErrorExp = [ 0.01, 0.007, 0.01, 0.0035, 0.0015, 0.001, 0.0005 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5.0 ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5.0, 6.5 ]
                    binContentExp = [ 0.17, 0.12, 0.076, 0.055, 0.032, 0.017, 0.007, 0.003, 0.001 ]
                    binHalfErrorExp = [ 0.01, 0.01, 0.0065, 0.005, 0.0035, 0.0025, 0.0015, 0.001, 0.0005 ]
              elif '0.35-0.55' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [10]
                    binEdgeLowExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binEdgeHighExp = [ 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8 ]
                    binContentExp = [ 0.12, 0.16, 0.19, 0.2, 0.2, 0.19, 0.18, 0.18, 0.18, 0.16 ]
                    binHalfErrorExp = [ 0.015, 0.01, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.01 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [10]
                    binEdgeLowExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binEdgeHighExp = [ 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8 ]
                    binContentExp = [ 0.15, 0.2, 0.28, 0.3, 0.33, 0.35, 0.36, 0.36, 0.33, 0.26 ]
                    binHalfErrorExp = [ 0.015, 0.01, 0.015, 0.01, 0.015, 0.01, 0.01, 0.01, 0.015, 0.02 ]
              elif '0.55-0.75' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [11]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8 ]
                    binContentExp = [ 0.06, 0.15, 0.17, 0.18, 0.17, 0.17, 0.17, 0.15, 0.14, 0.12, 0.1 ]
                    binHalfErrorExp = [ 0.01, 0.01, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [11]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8 ]
                    binContentExp = [ 0.1, 0.19, 0.26, 0.31, 0.29, 0.3, 0.28, 0.26, 0.22, 0.17, 0.12]
                    binHalfErrorExp = [ 0.015, 0.01, 0.01, 0.015, 0.01, 0.01, 0.005, 0.005, 0.005, 0.01, 0.01]
              elif '0.75-0.95' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [10]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binContentExp = [ 0.07, 0.17, 0.17, 0.18, 0.15, 0.15, 0.13, 0.11, 0.09, 0.07 ]
                    binHalfErrorExp = [ 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [10]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binContentExp = [ 0.12, 0.21, 0.24, 0.25, 0.24, 0.22, 0.18, 0.16, 0.12, 0.07 ]
                    binHalfErrorExp = [ 0.01, 0.01, 0.01, 0.01, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005 ]
              elif '0.95-1.15' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6 ]
                    binContentExp = [ 0.1, 0.15, 0.14, 0.13, 0.12, 0.1, 0.08, 0.07, 0.05 ]
                    binHalfErrorExp = [ 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6 ]
                    binContentExp = [ 0.12, 0.23, 0.24, 0.22, 0.16, 0.13, 0.11, 0.09, 0.06 ]
                    binHalfErrorExp = [ 0.01, 0.01, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005 ]
              elif '1.15-1.35' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [ 0.1, 0.15, 0.14, 0.12, 0.09, 0.07, 0.06, 0.04 ]
                    binHalfErrorExp = [ 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [ 0.14, 0.21, 0.2, 0.15, 0.12, 0.08, 0.06, 0.04 ]
                    binHalfErrorExp = [ 0.01, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005 ]
              elif '1.35-1.55' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [ 0.1, 0.14, 0.11, 0.08, 0.06, 0.05, 0.04, 0.03 ]
                    binHalfErrorExp = [ 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [ 0.14, 0.2, 0.15, 0.1, 0.08, 0.06, 0.04, 0.02 ]
                    binHalfErrorExp = [ 0.01, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005 ]
              elif '1.55-1.75' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [ 0.09, 0.12, 0.09, 0.06, 0.04, 0.04, 0.03, 0.02 ]
                    binHalfErrorExp = [ 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [ 0.13, 0.16, 0.12, 0.07, 0.05, 0.03, 0.02, 0.01 ]
                    binHalfErrorExp = [ 0.01, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005 ]
              elif '1.75-1.95' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [ 0.09, 0.1, 0.07, 0.05, 0.03, 0.02, 0.02, 0.01 ]
                    binHalfErrorExp = [ 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [7]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binContentExp = [ 0.12, 0.14, 0.1, 0.06, 0.03, 0.02, 0.01 ]
                    binHalfErrorExp = [ 0.01, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005 ]
              elif '1.95-2.15' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [6]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4 ]
                    binContentExp = [ 0.07, 0.08, 0.06, 0.04, 0.02, 0.01 ]
                    binHalfErrorExp = [ 0.005, 0.005, 0.005, 0.005, 0.005, 0.005 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [6]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4 ]
                    binContentExp = [ 0.11, 0.11, 0.06, 0.04, 0.02, 0.01 ]
                    binHalfErrorExp = [ 0.005, 0.005, 0.005, 0.005, 0.005, 0.005 ]
           elif 'pi+' in job["PARTICLE"]:
              if '0.05-0.1' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6.5 ]
                    binContentExp = [ 0.14, 0.15, 0.15, 0.131, 0.113, 0.081, 0.104, 0.05, 0.026 ]
                    binHalfErrorExp = [ 0.01, 0.01, 0.01, 0.0075, 0.0065, 0.0045, 0.0055, 0.0025, 0.0015  ]
                 elif 'pi+' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6.5 ]
                    binContentExp = [ 0.18, 0.21, 0.19, 0.19, 0.2, 0.19, 0.14, 0.196, 0.17 ]
                    binHalfErrorExp = [ 0.015, 0.01, 0.01, 0.015, 0.015, 0.01, 0.01, 0.007, 0.0015 ]
              elif '0.1-0.15' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6.5 ]
                    binContentExp = [ 0.21, 0.15, 0.117, 0.12, 0.056, 0.06, 0.037, 0.022, 0.006 ]
                    binHalfErrorExp = [ 0.015, 0.01, 0.0075, 0.01, 0.004, 0.004, 0.003, 0.002, 0.0005 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6.5 ]
                    binContentExp = [ 0.18, 0.21, 0.25, 0.21, 0.159, 0.11, 0.097, 0.065, 0.022 ]
                    binHalfErrorExp = [ 0.015, 0.015, 0.015, 0.04, 0.088, 0.01, 0.0055, 0.0035, 0.0015 ]
              elif '0.15-0.2' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6.5 ]
                    binContentExp = [ 0.2, 0.13, 0.12, 0.095, 0.048, 0.027, 0.017, 0.005, 0.001 ]
                    binHalfErrorExp = [ 0.015, 0.01, 0.01, 0.007, 0.0045, 0.003, 0.002, 0.001, 0.0005 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6.5 ]
                    binContentExp = [ 0.2, 0.18, 0.17, 0.13, 0.08, 0.076, 0.043, 0.025, 0.005 ]
                    binHalfErrorExp = [ 0.02, 0.01, 0.01, 0.01, 0.0065, 0.006, 0.004, 0.0025, 0.001 ]
              elif '0.2-0.25' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5. ]
                    binContentExp = [ 0.15, 0.16, 0.13, 0.052, 0.02, 0.011, 0.006, 0.002 ]
                    binHalfErrorExp = [ 0.015, 0.015, 0.01, 0.0065, 0.003, 0.002, 0.0015, 0.0015 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6.5 ]
                    binContentExp = [ 0.19, 0.13, 0.13, 0.09, 0.055, 0.04, 0.032, 0.023, 0.007 ]
                    binHalfErrorExp = [ 0.015, 0.01, 0.01, 0.01, 0.02, 0.0055, 0.0045, 0.004, 0.002  ]
              elif '0.35-0.55' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [10]
                    binEdgeLowExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binEdgeHighExp = [ 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8 ]
                    binContentExp = [ 0.137, 0.154, 0.179, 0.205, 0.187, 0.175, 0.17, 0.184, 0.186, 0.17 ]
                    binHalfErrorExp = [ 0.012, 0.008, 0.008, 0.0075, 0.005, 0.005, 0.0045, 0.005, 0.0065, 0.009 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [10]
                    binEdgeLowExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binEdgeHighExp = [ 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8 ]
                    binContentExp = [ 0.181, 0.214, 0.297, 0.319, 0.336, 0.41, 0.399, 0.414, 0.412, 0.367 ]
                    binHalfErrorExp = [ 0.0135, 0.011, 0.0135, 0.0105, 0.014, 0.0165, 0.0085, 0.012, 0.0165, 0.0225  ]
              elif '0.55-0.75' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [11]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8 ]
                    binContentExp = [ 0.082, 0.136, 0.146, 0.174, 0.172, 0.16, 0.149, 0.146, 0.136, 0.114, 0.1 ]
                    binHalfErrorExp = [ 0.0135, 0.006, 0.006, 0.0075, 0.005, 0.0045, 0.004, 0.004, 0.004, 0.0045, 0.006 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [11]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8 ]
                    binContentExp = [ 0.117, 0.185, 0.249, 0.319, 0.321, 0.311, 0.325, 0.329, 0.274, 0.21, 0.146 ]
                    binHalfErrorExp = [ 0.014, 0.0075, 0.011, 0.0135, 0.01, 0.0095, 0.0085, 0.008, 0.009, 0.012, 0.014 ]
              elif '0.75-0.95' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [10]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binContentExp = [ 0.066, 0.144, 0.157, 0.141, 0.126, 0.126, 0.125, 0.103, 0.082, 0.072 ]
                    binHalfErrorExp = [ 0.0065, 0.0065, 0.006, 0.004, 0.004, 0.0045, 0.0035, 0.0035, 0.0025, 0.003 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [10]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binContentExp = [ 0.124, 0.219, 0.267, 0.288, 0.245, 0.241, 0.209, 0.181, 0.141, 0.091 ]
                    binHalfErrorExp = [ 0.0095, 0.0095, 0.0105, 0.008, 0.007, 0.006, 0.005, 0.0045, 0.0055, 0.0065 ]
              elif '0.95-1.15' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6 ]
                    binContentExp = [ 0.076, 0.131, 0.133, 0.113, 0.092, 0.08, 0.07, 0.059, 0.048 ]
                    binHalfErrorExp = [ 0.006, 0.005, 0.005, 0.0035, 0.003, 0.0025, 0.002, 0.002, 0.0015 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6 ]
                    binContentExp = [ 0.128, 0.232, 0.223, 0.207, 0.183, 0.155, 0.122, 0.095, 0.068 ]
                    binHalfErrorExp = [ 0.0085, 0.008, 0.007, 0.0065, 0.005, 0.0045, 0.0035, 0.0035, 0.0035 ]
              elif '1.15-1.35' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [ 0.069, 0.129, 0.088, 0.092, 0.077, 0.056, 0.048, 0.041 ]
                    binHalfErrorExp = [ 0.0045, 0.005, 0.003, 0.0035,0.0025, 0.002, 0.0015, 0.0015 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [ 0.134, 0.219, 0.202, 0.151, 0.111, 0.088, 0.075, 0.057 ]
                    binHalfErrorExp = [ 0.0095, 0.008, 0.006, 0.0045, 0.004, 0.0025, 0.0025, 0.003 ]
              elif '1.35-1.55' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [ 0.083, 0.114, 0.096, 0.073, 0.051, 0.041, 0.035, 0.026 ]
                    binHalfErrorExp = [ 0.005, 0.0045, 0.0035, 0.003, 0.002, 0.0015, 0.0015, 0.0015  ]
                 elif 'pi+' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [ 0.155, 0.212, 0.163, 0.11, 0.083, 0.061, 0.045, 0.031 ]
                    binHalfErrorExp = [ 0.0095, 0.0075, 0.0055, 0.004, 0.003, 0.002, 0.002, 0.002 ]
              elif '1.55-1.75' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [ 0.057, 0.087, 0.067, 0.062, 0.043, 0.029, 0.021, 0.017 ]
                    binHalfErrorExp = [ 0.004, 0.004, 0.0025, 0.0025, 0.0025, 0.0015, 0.001, 0.001 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [ 0.141, 0.152, 0.116, 0.083, 0.063, 0.041, 0.028, 0.017 ]
                    binHalfErrorExp = [ 0.0085, 0.005, 0.004, 0.003, 0.0025, 0.002, 0.002, 0.0015  ]
              elif '1.75-1.95' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [ 0.05, 0.069, 0.049, 0.037, 0.027, 0.021, 0.014, 0.011 ]
                    binHalfErrorExp = [ 0.003, 0.0035, 0.002, 0.002, 0.0015, 0.001, 0.001, 0.001 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [ 0.105, 0.129, 0.098, 0.054, 0.043, 0.027, 0.017, 0.01 ]
                    binHalfErrorExp = [ 0.0075, 0.0045, 0.0035, 0.0025, 0.002, 0.0015, 0.001, 0.001 ]
              elif '1.95-2.15' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [ 0.054, 0.054, 0.039, 0.025, 0.021, 0.014, 0.009, 0.006 ]
                    binHalfErrorExp = [ 0.0035, 0.0025, 0.002, 0.0015, 0.0015, 0.001, 0.0005, 0.0005 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [ 0.097, 0.115, 0.068, 0.036, 0.026, 0.017, 0.009, 0.004 ]
                    binHalfErrorExp = [ 0.0065, 0.004, 0.003, 0.0025, 0.0015, 0.0015, 0.001, 0.0005 ]
           elif 'pi-' in job["PARTICLE"]:
              if '0.05-0.1' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6.5 ]
                    binContentExp = [ 0.27, 0.23, 0.25, 0.24, 0.222, 0.183, 0.177, 0.172, 0.116 ]
                    binHalfErrorExp = [ 0.015, 0.01, 0.01, 0.01, 0.0065, 0.007, 0.007, 0.0055, 0.004 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6.5 ]
                    binContentExp = [ 0.14, 0.204, 0.174, 0.126, 0.094, 0.094, 0.068, 0.05, 0.023 ]
                    binHalfErrorExp = [ 0.01, 0.007, 0.0055, 0.0055, 0.005, 0.005, 0.0045, 0.0025, 0.0015 ]
              elif '0.1-0.15' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6.5 ]
                    binContentExp = [ 0.25, 0.24, 0.21, 0.2, 0.167, 0.097, 0.076, 0.046, 0.02 ]
                    binHalfErrorExp = [ 0.015, 0.01, 0.01, 0.01, 0.007, 0.005, 0.004, 0.003, 0.0015 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6.5 ]
                    binContentExp = [ 0.19, 0.167, 0.146, 0.13, 0.086, 0.053, 0.05, 0.019, 0.003 ]
                    binHalfErrorExp = [ 0.01, 0.0075, 0.0065, 0.006, 0.005, 0.0035, 0.003, 0.002, 0.0005 ]
              elif '0.15-0.2' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6.5 ]
                    binContentExp = [ 0.22, 0.23, 0.18, 0.12, 0.077, 0.042, 0.029, 0.017, 0.005 ]
                    binHalfErrorExp = [ 0.015, 0.01, 0.01, 0.0065, 0.004, 0.003, 0.002, 0.0015, 0.0005 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5. ]
                    binContentExp = [ 0.19, 0.13, 0.115, 0.091, 0.05, 0.031, 0.016, 0.009 ]
                    binHalfErrorExp = [ 0.01, 0.0065, 0.0055, 0.005, 0.0035, 0.0025, 0.0015, 0.001 ]
              elif '0.2-0.25' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6.5 ]
                    binContentExp = [ 0.15, 0.099, 0.076, 0.07, 0.053, 0.029, 0.02, 0.014, 0.003 ]
                    binHalfErrorExp = [ 0.01, 0.0065, 0.007, 0.006, 0.0045, 0.0025, 0.002, 0.0015, 0.0005 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5. ]
                    binContentExp = [ 0.14, 0.19, 0.12, 0.069, 0.036, 0.011, 0.002, 0.001 ]
                    binHalfErrorExp = [ 0.01, 0.01, 0.01, 0.006, 0.004, 0.0025, 0.001, 0.0005 ]
              elif '0.35-0.55' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [10]
                    binEdgeLowExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binEdgeHighExp = [ 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8 ]
                    binContentExp = [ 0.173, 0.219, 0.267, 0.272, 0.275, 0.301, 0.309, 0.317, 0.315, 0.285 ]
                    binHalfErrorExp = [ 0.013, 0.009, 0.011, 0.007, 0.009, 0.01, 0.0075, 0.009, 0.012, 0.0165 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [10]
                    binEdgeLowExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binEdgeHighExp = [ 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8 ]
                    binContentExp = [ 0.142, 0.168, 0.206, 0.224, 0.24, 0.255, 0.258, 0.272, 0.273, 0.227 ]
                    binHalfErrorExp = [ 0.011, 0.007, 0.0085, 0.0075, 0.0085, 0.0065, 0.006, 0.0075, 0.0125, 0.0165 ]
              elif '0.55-0.75' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [11]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8 ]
                    binContentExp = [ 0.113, 0.193, 0.241, 0.261, 0.268, 0.247, 0.24, 0.234, 0.226, 0.195, 0.16 ]
                    binHalfErrorExp = [ 0.0125, 0.008, 0.0075, 0.009, 0.0075, 0.006, 0.0055, 0.005, 0.006, 0.0085, 0.011 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [11]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8 ]
                    binContentExp = [ 0.084, 0.147, 0.183, 0.211, 0.222, 0.207, 0.223, 0.218, 0.196, 0.149, 0.094 ]
                    binHalfErrorExp = [ 0.0115, 0.006, 0.0065, 0.0085, 0.0065, 0.0045, 0.007, 0.005, 0.0065, 0.01, 0.0095 ]
              elif '0.75-0.95' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [10]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binContentExp = [ 0.099, 0.214, 0.228, 0.217, 0.202, 0.179, 0.164, 0.158, 0.136, 0.101 ]
                    binHalfErrorExp = [ 0.0075, 0.008, 0.0065, 0.006, 0.0045, 0.004, 0.0035, 0.0035, 0.004, 0.0055 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [10]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binContentExp = [ 0.101, 0.154, 0.176, 0.185, 0.174, 0.162, 0.146, 0.138, 0.11, 0.071 ]
                    binHalfErrorExp = [ 0.0075, 0.005, 0.0065, 0.005, 0.005, 0.0035, 0.0035, 0.0035, 0.0045, 0.0055 ]
              elif '0.95-1.15' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6 ]
                    binContentExp = [ 0.129, 0.202, 0.195, 0.178, 0.145, 0.121, 0.106, 0.095, 0.076 ]
                    binHalfErrorExp = [ 0.0075, 0.006, 0.005, 0.004, 0.0035, 0.0025, 0.002, 0.002, 0.0025 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6 ]
                    binContentExp = [ 0.093, 0.146, 0.152, 0.151, 0.128, 0.111, 0.089, 0.069, 0.048 ]
                    binHalfErrorExp = [ 0.006, 0.0055, 0.0045, 0.004, 0.003, 0.0025, 0.0025, 0.0025, 0.0025 ]
              elif '1.15-1.35' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [ 0.14, 0.179, 0.168, 0.129, 0.108, 0.088, 0.067, 0.056 ]
                    binHalfErrorExp = [ 0.007, 0.005, 0.004, 0.0035, 0.0025, 0.0025, 0.002, 0.002 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [ 0.086, 0.145, 0.14, 0.109, 0.088, 0.071, 0.058, 0.042 ]
                    binHalfErrorExp = [ 0.006, 0.006, 0.0035, 0.003, 0.0025, 0.002, 0.002, 0.0025 ]
              elif '1.35-1.55' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [ 0.125, 0.165, 0.139, 0.102, 0.08, 0.061, 0.047, 0.035 ]
                    binHalfErrorExp = [ 0.006, 0.005, 0.0035, 0.003, 0.0025, 0.002, 0.0015, 0.001 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [ 0.09, 0.121, 0.109, 0.081, 0.062, 0.046, 0.035, 0.025 ]
                    binHalfErrorExp = [ 0.0055, 0.0035, 0.003, 0.0025, 0.002, 0.0015, 0.0015, 0.0015 ]
              elif '1.55-1.75' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [ 0.103, 0.14, 0.109, 0.09, 0.065, 0.04, 0.028, 0.022 ]
                    binHalfErrorExp = [ 0.0055, 0.004, 0.003, 0.0025, 0.0025, 0.002, 0.001, 0.001 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [ 0.082, 0.101, 0.085, 0.066, 0.043, 0.03, 0.019, 0.011 ]
                    binHalfErrorExp = [ 0.005, 0.0035, 0.0025, 0.002, 0.0015, 0.0015, 0.001, 0.001 ]
              elif '1.75-1.95' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [ 0.097, 0.123, 0.09, 0.052, 0.042, 0.029, 0.021, 0.016 ]
                    binHalfErrorExp = [ 0.0055, 0.0035, 0.003, 0.002, 0.0015, 0.0015, 0.001, 0.001 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [ 0.074, 0.094, 0.067, 0.044, 0.027, 0.019, 0.013, 0.008 ]
                    binHalfErrorExp = [ 0.0045, 0.003, 0.0025, 0.0015, 0.001, 0.001, 0.001, 0.001 ]
              elif '1.95-2.15' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [ 0.084, 0.104, 0.063, 0.044, 0.027, 0.02, 0.015, 0.01 ]
                    binHalfErrorExp = [ 0.004, 0.003, 0.002, 0.0015, 0.0015, 0.0005, 0.001, 0.0005 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [ 0.056, 0.076, 0.047, 0.026, 0.022, 0.013, 0.008, 0.004 ]
                    binHalfErrorExp = [ 0.003, 0.003, 0.0025, 0.001, 0.001, 0.001, 0.0005, 0.0005 ]

# NOW NEED TO MILTIPLY EACH BIN CONTENT?ERROR by 1000. (something with units)
           binContentExp = [ 1000.*i for i in binContentExp ]
           binHalfErrorExp = [ 1000.*i for i in binHalfErrorExp ]

#           print("inspireid = ",inspireid)

           rjsonexp = getJSON( job, "histogram",
                            inspireId=inspireid,
                            mctool_name="experiment",
                            mctool_version="experiment",
                            mctool_model="experiment",
                            testName="experiment",
                            plotType="TH1",
                            observableName=observableName,
                            reaction="particle production", # otherwise it'll write the default which is "reaction name"
                            targetName=material,
                            beamParticle=job["PARTICLE"],
                            beamEnergies=beamenergies,
                            secondaryParticle=secondary,
                            parameters=parameters,
                            nBins=nBinsExp,
                            binContent=binContentExp,
                            binEdgeLow=binEdgeLowExp,
                            binEdgeHigh=binEdgeHighExp,
                            yStatErrorsPlus=binHalfErrorExp,
                            yStatErrorsMinus=binHalfErrorExp,
                            xAxisName=xAxisTitle,
                            yAxisName=yAxisTitle,
                            title=htitle
                          )

           rjson = getJSON( job, "histogram",
                            mctool_name="GEANT4",
                            mctool_version=job["VERSION"],
                            mctool_model=job["GENERATOR"],
                            testName="test19-HARP-Carbon",
                            plotType="TH1",
                            observableName=observableName,
                            reaction="particle production", # otherwise it'll write the default which is "reaction name"
                            targetName=material,
                            beamParticle=job["PARTICLE"],
                            beamEnergies=beamenergies,
                            secondaryParticle=secondary,
                            parameters=parameters,
                            nBins=nBins,
                            binContent=binContent,
                            binEdgeLow=binEdgeLow,
                            binEdgeHigh=binEdgeHigh,
                            yStatErrorsPlus=yStatErrorPlus,
                            yStatErrorsMinus=yStatErrorMinus,
                            xAxisName=xAxisTitle,
                            yAxisName=yAxisTitle,
                            title=htitle
                          )

# --> already generated -->           yield rjsonexp
           yield rjson
