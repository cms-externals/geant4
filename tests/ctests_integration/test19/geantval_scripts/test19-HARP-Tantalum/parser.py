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
   TEST="test19-HARP-Tantalum"

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

           if 'proton' in job["PARTICLE"]:
              if '0.05-0.1' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5. ]
                    binContentExp = [       680.0,
      600.0,
      420.0,
      270.0,
      120.0,
      46.0,
      90.0,
      17.0 ]
                    binHalfErrorExp = [       50.0,
      40.0,
      30.0,
      20.0,
      15.0,
      7.5,
      10.0,
      5.5 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5.  ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6.5 ]
                    binContentExp = [       870.0,
      950.0,
      800.0,
      590.0,
      340.0,
      150.0,
      50.0,
      50.0,
      7.0 ]
                    binHalfErrorExp = [       60.0,
      450.0,
      400.0,
      350.0,
      300.0,
      250.0,
      150.0,
      100.0,
      15.0 ]
              elif '0.1-0.15' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6.5 ]
                    binContentExp = [       1180.0,
      680.0,
      410.0,
      170.0,
      100.0,
      70.0,
      31.0,
      14.0,
      1.0 ]
                    binHalfErrorExp = [       75.0,
      40.0,
      25.0,
      15.0,
      10.0,
      10.0,
      5.0,
      3.0,
      0.5 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6.5 ]
                    binContentExp = [       930.0,
      820.0,
      820.0,
      400.0,
      150.0,
      90.0,
      20.0,
      9.0,
      2.0 ]
                    binHalfErrorExp = [       55.0,
      50.0,
      450.0,
      300.0,
      200.0,
      150.0,
      100.0,
      20.0,
      5.0 ]
              elif '0.15-0.2' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [6]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5 ]
                    binContentExp = [       860.0,
      550.0,
      200.0,
      160.0,
      60.0,
      31.0 ]
                    binHalfErrorExp = [       60.0,
      35.0,
      20.0,
      15.0,
      10.0,
      7.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6.5 ]
                    binContentExp = [       930.0,
      770.0,
      380.0,
      200.0,
      120.0,
      60.0,
      16.0,
      7.0,
      2.0 ]
                    binHalfErrorExp = [       60.0,
      450.0,
      300.0,
      200.0,
      150.0,
      100.0,
      45.0,
      20.0,
      10.0 ]
              elif '0.2-0.25' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [5]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5 ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3. ]
                    binContentExp = [       760.0,
      320.0,
      140.0,
      60.0,
      21.0 ]
                    binHalfErrorExp = [       55.0,
      30.0,
      20.0,
      10.0,
      5.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6.5 ]
                    binContentExp = [       630.0,
      460.0,
      200.0,
      240.0,
      90.0,
      27.0,
      16.0,
      11.0,
      2.0 ]
                    binHalfErrorExp = [       50.0,
      400.0,
      250.0,
      250.0,
      150.0,
      75.0,
      65.0,
      50.0,
      25.0 ]
              elif '0.35-0.55' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [10]
                    binEdgeLowExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binEdgeHighExp = [ 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8 ]
                    binContentExp = [       1310.0,
      1590.0,
      1770.0,
      1880.0,
      1670.0,
      1490.0,
      1420.0,
      1280.0,
      1150.0,
      920.0 ]
                    binHalfErrorExp = [       155.0,
      80.0,
      70.0,
      55.0,
      40.0,
      35.0,
      30.0,
      35.0,
      40.0,
      55.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [10]
                    binEdgeLowExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binEdgeHighExp = [ 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8 ]
                    binContentExp = [       1160.0,
      1540.0,
      2000.0,
      2130.0,
      2120.0,
      2030.0,
      2010.0,
      2020.0,
      1780.0,
      1220.0 ]
                    binHalfErrorExp = [       150.0,
      95.0,
      70.0,
      80.0,
      50.0,
      40.0,
      65.0,
      55.0,
      90.0,
      100.0 ]
              elif '0.55-0.75' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [11]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8 ]
                    binContentExp = [       1180.0,
      1740.0,
      1990.0,
      1910.0,
      1570.0,
      1430.0,
      1330.0,
      1180.0,
      1050.0,
      840.0,
      640.0 ]
                    binHalfErrorExp = [       195.0,
      90.0,
      70.0,
      50.0,
      40.0,
      35.0,
      30.0,
      25.0,
      25.0,
      40.0,
      40.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [11]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8 ]
                    binContentExp = [       980.0,
      1790.0,
      2100.0,
      2130.0,
      2170.0,
      1910.0,
      1870.0,
      1650.0,
      1410.0,
      950.0,
      580.0 ]
                    binHalfErrorExp = [       185.0,
      95.0,
      75.0,
      65.0,
      45.0,
      55.0,
      45.0,
      45.0,
      50.0,
      65.0,
      55.0 ]
              elif '0.75-0.95' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [10]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binContentExp = [       1420.0,
      2040.0,
      1870.0,
      1680.0,
      1330.0,
      1220.0,
      1010.0,
      880.0,
      710.0,
      560.0 ]
                    binHalfErrorExp = [       160.0,
      65.0,
      50.0,
      40.0,
      30.0,
      30.0,
      25.0,
      20.0,
      20.0,
      25.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [10]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binContentExp = [       1130.0,
      2060.0,
      2150.0,
      1920.0,
      1780.0,
      1530.0,
      1330.0,
      1140.0,
      820.0,
      480.0 ]
                    binHalfErrorExp = [       150.0,
      60.0,
      70.0,
      45.0,
      40.0,
      40.0,
      30.0,
      30.0,
      35.0,
      35.0 ]
              elif '0.95-1.15' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6 ]
                    binContentExp = [       1860.0,
      2130.0,
      1690.0,
      1420.0,
      1100.0,
      890.0,
      750.0,
      600.0,
      460.0 ]
                    binHalfErrorExp = [       145.0,
      65.0,
      45.0,
      40.0,
      30.0,
      20.0,
      15.0,
      15.0,
      15.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6 ]
                    binContentExp = [       1400.0,
      2060.0,
      1980.0,
      1670.0,
      1300.0,
      1020.0,
      820.0,
      650.0,
      380.0 ]
                    binHalfErrorExp = [       135.0,
      65.0,
      55.0,
      40.0,
      40.0,
      35.0,
      25.0,
      20.0,
      20.0 ]
              elif '1.15-1.35' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [       2110.0,
      2070.0,
      1580.0,
      1200.0,
      850.0,
      670.0,
      510.0,
      410.0 ]
                    binHalfErrorExp = [       180.0,
      70.0,
      45.0,
      40.0,
      25.0,
      20.0,
      15.0,
      15.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [       1640.0,
      2030.0,
      1840.0,
      1290.0,
      900.0,
      660.0,
      490.0,
      360.0 ]
                    binHalfErrorExp = [       85.0,
      50.0,
      40.0,
      30.0,
      20.0,
      15.0,
      15.0 ]
              elif '1.35-1.55' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [       2370.0,
      1870.0,
      1320.0,
      900.0,
      600.0,
      450.0,
      340.0,
      270.0 ]
                    binHalfErrorExp = [       235.0,
      90.0,
      50.0,
      35.0,
      20.0,
      15.0,
      10.0,
      10.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [       1610.0,
      1910.0,
      1640.0,
      1060.0,
      750.0,
      500.0,
      320.0,
      220.0 ]
                    binHalfErrorExp = [       175.0,
      105.0,
      55.0,
      35.0,
      25.0,
      20.0,
      15.0,
      10.0 ]
              elif '1.55-1.75' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [       1840.0,
      1620.0,
      1110.0,
      670.0,
      430.0,
      300.0,
      210.0,
      150.0 ]
                    binHalfErrorExp = [       190.0,
      80.0,
      40.0,
      25.0,
      20.0,
      10.0,
      5.0,
      5.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [7]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binContentExp = [       1570.0,
      1690.0,
      1360.0,
      740.0,
      480.0,
      320.0,
      210.0 ]
                    binHalfErrorExp = [       185.0,
      85.0,
      45.0,
      30.0,
      15.0,
      15.0,
      10.0 ]
              elif '1.75-1.95' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [       1540.0,
      1270.0,
      810.0,
      510.0,
      280.0,
      210.0,
      180.0,
      120.0 ]
                    binHalfErrorExp = [       115.0,
      40.0,
      25.0,
      20.0,
      15.0,
      5.0,
      5.0,
      5.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [7]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binContentExp = [       1280.0,
      1390.0,
      970.0,
      540.0,
      320.0,
      190.0,
      120.0 ]
                    binHalfErrorExp = [       100.0,
      40.0,
      30.0,
      25.0,
      15.0,
      10.0,
      5.0 ]
              elif '1.95-2.15' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [       1430.0,
      990.0,
      570.0,
      360.0,
      190.0,
      150.0,
      110.0,
      70.0 ]
                    binHalfErrorExp = [       85.0,
      25.0,
      15.0,
      15.0,
      10.0,
      5.0,
      5.0,
      5.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [       1040.0,
      1090.0,
      640.0,
      340.0,
      190.0,
      110.0,
      50.0,
      30.0 ]
                    binHalfErrorExp = [       80.0,
      25.0,
      25.0,
      20.0,
      10.0,
      10.0,
      5.0,
      5.0 ]
           if 'pi+' in job["PARTICLE"]:
              if '0.05-0.1' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0 ]
                    binEdgeHighExp = [ 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.5 ]
                    binContentExp = [       790.0,
      900.0,
      770.0,
      530.0,
      560.0,
      410.0,
      300.0,
      240.0,
      110.0 ]
                    binHalfErrorExp = [       65.0,
      60.0,
      45.0,
      35.0,
      35.0,
      25.0,
      20.0,
      15.0,
      7.5 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0 ]
                    binEdgeHighExp = [ 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.5 ]
                    binContentExp = [       1060.0,
      1030.0,
      630.0,
      710.0,
      730.0,
      610.0,
      590.0,
      810.0,
      680.0 ]
                    binHalfErrorExp = [       80.0,
      60.0,
      45.0,
      60.0,
      45.0,
      40.0,
      40.0,
      30.0,
      50.0 ]
              elif '0.1-0.15' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0 ]
                    binEdgeHighExp = [ 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.5 ]
                    binContentExp = [       1190.0,
      780.0,
      750.0,
      480.0,
      430.0,
      270.0,
      170.0,
      110.0,
      28.0 ]
                    binHalfErrorExp = [       90.0,
      50.0,
      45.0,
      35.0,
      30.0,
      20.0,
      15.0,
      10.0,
      4.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0 ]
                    binEdgeHighExp = [ 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.5 ]
                    binContentExp = [       1380.0,
      880.0,
      940.0,
      770.0,
      580.0,
      390.0,
      450.0,
      250.0,
      78.0 ]
                    binHalfErrorExp = [       90.0,
      60.0,
      70.0,
      75.0,
      75.0,
      30.0,
      30.0,
      20.0,
      10.0 ]
              elif '0.15-0.2' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0 ]
                    binEdgeHighExp = [ 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.5 ]
                    binContentExp = [       1230.0,
      750.0,
      600.0,
      490.0,
      190.0,
      140.0,
      50.0,
      24.0,
      3.0 ]
                    binHalfErrorExp = [       90.0,
      55.0,
      45.0,
      40.0,
      20.0,
      15.0,
      10.0,
      5.0,
      1.5 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0 ]
                    binEdgeHighExp = [ 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.5 ]
                    binContentExp = [       1380.0,
      700.0,
      560.0,
      510.0,
      280.0,
      250.0,
      210.0,
      110.0,
      29.0 ]
                    binHalfErrorExp = [       95.0,
      50.0,
      40.0,
      40.0,
      25.0,
      20.0,
      20.0,
      15.0,
      6.0 ]
              elif '0.2-0.25' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [6]
                    binEdgeLowExp = [ 0.5, 1.0, 1.5, 2.0, 2.5, 3.0 ]
                    binEdgeHighExp = [ 1.0, 1.5, 2.0, 2.5, 3.0, 3.5 ]
                    binContentExp = [      1070.0,
      720.0,
      350.0,
      120.0,
      26.0,
      16.0 ]
                    binHalfErrorExp = [       90.0,
      65.0,
      40.0,
      20.0,
      6.0,
      6.5 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0 ]
                    binEdgeHighExp = [ 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.5 ]
                    binContentExp = [       690.0,
      500.0,
      370.0,
      390.0,
      230.0,
      180.0,
      100.0,
      50.0,
      20.0 ]
                    binHalfErrorExp = [       65.0,
      55.0,
      45.0,
      45.0,
      35.0,
      30.0,
      15.0,
      10.0,
      10.0 ]
              elif '0.35-0.55' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [10]
                    binEdgeLowExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binEdgeHighExp = [ 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8 ]
                    binContentExp = [       1310.0,
      1650.0,
      2080.0,
      2130.0,
      1700.0,
      1570.0,
      1570.0,
      1440.0,
      1360.0,
      1200.0 ]
                    binHalfErrorExp = [       270.0,
      185.0,
      95.0,
      65.0,
      50.0,
      45.0,
      40.0,
      40.0,
      50.0,
      75.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [10]
                    binEdgeLowExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binEdgeHighExp = [ 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8 ]
                    binContentExp = [       1280.0,
      1730.0,
      2380.0,
      2600.0,
      2410.0,
      2410.0,
      2530.0,
      2740.0,
      2330.0,
      1780.0 ]
                    binHalfErrorExp = [       280.0,
      185.0,
      115.0,
      100.0,
      55.0,
      60.0,
      90.0,
      75.0,
      125.0,
      120.0 ]
              elif '0.55-0.75' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [11]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8 ]
                    binContentExp = [       950.0,
      1690.0,
      2200.0,
      1850.0,
      1610.0,
      1560.0,
      1510.0,
      1320.0,
      1190.0,
      980.0,
      760.0 ]
                    binHalfErrorExp = [       275.0,
      190.0,
      85.0,
      55.0,
      45.0,
      45.0,
      40.0,
      35.0,
      30.0,
      40.0,
      50.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [11]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8 ]
                    binContentExp = [       600.0,
      1630.0,
      2360.0,
      2620.0,
      2670.0,
      2310.0,
      2180.0,
      2190.0,
      1970.0,
      1450.0,
      930.0 ]
                    binHalfErrorExp = [       220.0,
      195.0,
      130.0,
      75.0,
      70.0,
      55.0,
      70.0,
      55.0,
      65.0,
      85.0,
      80.0 ]
              elif '0.75-0.95' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [10]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binContentExp = [       1460.0,
      1960.0,
      2020.0,
      1590.0,
      1500.0,
      1280.0,
      1160.0,
      1030.0,
      810.0,
      640.0 ]
                    binHalfErrorExp = [       220.0,
      100.0,
      60.0,
      45.0,
      40.0,
      30.0,
      30.0,
      25.0,
      25.0,
      25.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [10]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binContentExp = [       1020.0,
      2300.0,
      2490.0,
      2310.0,
      2180.0,
      1900.0,
      1680.0,
      1450.0,
      1090.0,
      670.0 ]
                    binHalfErrorExp = [       215.0,
      105.0,
      90.0,
      60.0,
      55.0,
      45.0,
      40.0,
      40.0,
      40.0,
      45.0 ]
              elif '0.95-1.15' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6 ]
                    binContentExp = [       2020.0,
      2020.0,
      1810.0,
      1520.0,
      1190.0,
      930.0,
      820.0,
      690.0,
      540.0 ]
                    binHalfErrorExp = [       170.0,
      75.0,
      60.0,
      45.0,
      35.0,
      25.0,
      20.0,
      20.0,
      20.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6 ]
                    binContentExp = [       1560.0,
      2260.0,
      2340.0,
      1830.0,
      1540.0,
      1360.0,
      1120.0,
      870.0,
      570.0 ]
                    binHalfErrorExp = [       165.0,
      80.0,
      70.0,
      50.0,
      45.0,
      40.0,
      35.0,
      30.0,
      25.0 ]
              elif '1.15-1.35' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [       2290.0,
      2020.0,
      1520.0,
      1110.0,
      960.0,
      760.0,
      580.0,
      480.0 ]
                    binHalfErrorExp = [       195.0,
      75.0,
      50.0,
      35.0,
      30.0,
      25.0,
      20.0,
      15.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [       1830.0,
      2220.0,
      2000.0,
      1490.0,
      1090.0,
      910.0,
      710.0,
      530.0 ]
                    binHalfErrorExp = [       165.0,
      100.0,
      65.0,
      50.0,
      35.0,
      25.0,
      20.0,
      20.0 ]
              elif '1.35-1.55' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [       2630.0,
      1980.0,
      1300.0,
      870.0,
      690.0,
      560.0,
      420.0,
      320.0 ]
                    binHalfErrorExp = [       250.0,
      90.0,
      50.0,
      35.0,
      25.0,
      20.0,
      15.0,
      15.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [       1870.0,
      2160.0,
      1830.0,
      1140.0,
      840.0,
      630.0,
      470.0,
      350.0 ]
                    binHalfErrorExp = [       215.0,
      120.0,
      60.0,
      40.0,
      30.0,
      20.0,
      15.0,
      15.0 ]
              elif '1.55-1.75' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [       1990.0,
      1700.0,
      1070.0,
      730.0,
      570.0,
      410.0,
      290.0,
      200.0 ]
                    binHalfErrorExp = [       205.0,
      90.0,
      45.0,
      30.0,
      20.0,
      15.0,
      10.0,
      10.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [       1760.0,
      1920.0,
      1510.0,
      870.0,
      630.0,
      450.0,
      280.0,
      180.0 ]
                    binHalfErrorExp = [       215.0,
      100.0,
      55.0,
      30.0,
      20.0,
      20.0,
      15.0,
      10.0 ]
              elif '1.75-1.95' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [       1550.0,
      1330.0,
      850.0,
      560.0,
      330.0,
      260.0,
      190.0,
      140.0 ]
                    binHalfErrorExp = [       120.0,
      45.0,
      25.0,
      25.0,
      15.0,
      10.0,
      10.0,
      5.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [       1440.0,
      1550.0,
      1070.0,
      600.0,
      390.0,
      260.0,
      170.0,
      100.0 ]
                    binHalfErrorExp = [       110.0,
      55.0,
      35.0,
      25.0,
      15.0,
      10.0,
      10.0,
      10.0 ]
              elif '1.95-2.15' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [       1360.0,
      1040.0,
      600.0,
      420.0,
      250.0,
      190.0,
      150.0,
      120.0 ]
                    binHalfErrorExp = [       90.0,
      30.0,
      20.0,
      20.0,
      10.0,
      10.0,
      10.0,
      10.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [       1090.0,
      1240.0,
      720.0,
      450.0,
      300.0,
      190.0,
      100.0,
      50.0 ]
                    binHalfErrorExp = [       80.0,
      30.0,
      25.0,
      20.0,
      10.0,
      15.0,
      10.0,
      5.0 ]
           if 'pi-' in job["PARTICLE"]:
              if '0.05-0.1' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6.5 ]
                    binContentExp = [       1480.0,
      1100.0,
      1220.0,
      1230.0,
      1000.0,
      810.0,
      690.0,
      670.0,
      490.0 ]
                    binHalfErrorExp = [       65.0,
      45.0,
      40.0,
      50.0,
      35.0,
      40.0,
      35.0,
      30.0,
      25.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6.5 ]
                    binContentExp = [       1060.0,
      1270.0,
      1040.0,
      700.0,
      360.0,
      370.0,
      290.0,
      180.0,
      110.0 ]
                    binHalfErrorExp = [       55.0,
      45.0,
      35.0,
      35.0,
      30.0,
      30.0,
      25.0,
      15.0,
      10.0 ]
              if '0.1-0.15' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6.5 ]
                    binContentExp = [       1800.0,
      1140.0,
      990.0,
      910.0,
      670.0,
      380.0,
      310.0,
      210.0,
      95.0 ]
                    binHalfErrorExp = [       100.0,
      50.0,
      45.0,
      40.0,
      30.0,
      25.0,
      20.0,
      15.0,
      7.5 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6.5 ]
                    binContentExp = [       1370.0,
      920.0,
      680.0,
      580.0,
      420.0,
      220.0,
      220.0,
      80.0,
      19.0 ]
                    binHalfErrorExp = [       70.0,
      45.0,
      35.0,
      30.0,
      25.0,
      20.0,
      15.0,
      10.0,
      3.5 ]
              if '0.15-0.2' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6.5 ]
                    binContentExp = [       1430.0,
      1310.0,
      880.0,
      530.0,
      310.0,
      170.0,
      160.0,
      79.0,
      34.0 ]
                    binHalfErrorExp = [       85.0,
      60.0,
      45.0,
      30.0,
      20.0,
      15.0,
      10.0,
      7.0,
      3.5 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6.5 ]
                    binContentExp = [       1220.0,
      640.0,
      510.0,
      410.0,
      270.0,
      160.0,
      100.0,
      24.0,
      4.0 ]
                    binHalfErrorExp = [       70.0,
      35.0,
      30.0,
      25.0,
      20.0,
      15.0,
      10.0,
      3.5,
      1.0 ]
              if '0.2-0.25' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6.5 ]
                    binContentExp = [       1100.0,
      400.0,
      390.0,
      280.0,
      210.0,
      150.0,
      80.0,
      60.0,
      25.0 ]
                    binHalfErrorExp = [       70.0,
      30.0,
      40.0,
      25.0,
      20.0,
      15.0,
      10.0,
      10.0,
      5.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 5. ]
                    binEdgeHighExp = [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 5., 6.5 ]
                    binContentExp = [       1010.0,
      900.0,
      790.0,
      420.0,
      220.0,
      50.0,
      16.0,
      7.0,
      2.0 ]
                    binHalfErrorExp = [       70.0,
      60.0,
      60.0,
      40.0,
      25.0,
      10.0,
      6.0,
      3.5,
      1.0 ]
              if '0.35-0.55' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [10]
                    binEdgeLowExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binEdgeHighExp = [ 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8 ]
                    binContentExp = [       1710.0,
      2020.0,
      2400.0,
      2430.0,
      2290.0,
      2160.0,
      2140.0,
      2040.0,
      1900.0,
      1620.0 ]
                    binHalfErrorExp = [       285.0,
      175.0,
      110.0,
      75.0,
      55.0,
      50.0,
      45.0,
      50.0,
      65.0,
      95.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [10]
                    binEdgeLowExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binEdgeHighExp = [ 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8 ]
                    binContentExp = [       1230.0,
      1280.0,
      1910.0,
      1940.0,
      1930.0,
      1770.0,
      1880.0,
      1780.0,
      1590.0,
      1240.0 ]
                    binHalfErrorExp = [       255.0,
      180.0,
      100.0,
      60.0,
      55.0,
      40.0,
      50.0,
      45.0,
      70.0,
      95.0 ]
              if '0.55-0.75' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [11]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8 ]
                    binContentExp = [       1250.0,
      2110.0,
      2440.0,
      2480.0,
      2330.0,
      2250.0,
      1980.0,
      1780.0,
      1580.0,
      1310.0,
      990.0 ]
                    binHalfErrorExp = [       270.0,
      165.0,
      95.0,
      70.0,
      70.0,
      50.0,
      40.0,
      35.0,
      40.0,
      55.0,
      65.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [11]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8 ]
                    binContentExp = [       830.0,
      1330.0,
      1940.0,
      1950.0,
      1950.0,
      1790.0,
      1580.0,
      1550.0,
      1330.0,
      1000.0,
      650.0 ]
                    binHalfErrorExp = [       245.0,
      160.0,
      100.0,
      60.0,
      65.0,
      40.0,
      35.0,
      35.0,
      45.0,
      60.0,
      60.0 ]
              if '0.75-0.95' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [10]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binContentExp = [       1850.0,
      2630.0,
      2520.0,
      2200.0,
      1920.0,
      1710.0,
      1450.0,
      1210.0,
      950.0,
      760.0 ]
                    binHalfErrorExp = [       215.0,
      105.0,
      75.0,
      55.0,
      45.0,
      40.0,
      35.0,
      30.0,
      25.0,
      30.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [10]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7 ]
                    binContentExp = [       910.0,
      1700.0,
      2030.0,
      1900.0,
      1580.0,
      1360.0,
      1270.0,
      1120.0,
      860.0,
      550.0 ]
                    binHalfErrorExp = [       175.0,
      105.0,
      60.0,
      50.0,
      35.0,
      35.0,
      30.0,
      30.0,
      30.0,
      40.0 ]
              if '0.95-1.15' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6 ]
                    binContentExp = [       2390.0,
      2850.0,
      2400.0,
      1910.0,
      1480.0,
      1190.0,
      1030.0,
      870.0,
      670.0 ]
                    binHalfErrorExp = [       195.0,
      85.0,
      60.0,
      50.0,
      40.0,
      30.0,
      25.0,
      20.0,
      20.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [9]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6 ]
                    binContentExp = [       1190.0,
      1810.0,
      1820.0,
      1580.0,
      1330.0,
      1070.0,
      840.0,
      640.0,
      410.0 ]
                    binHalfErrorExp = [       125.0,
      75.0,
      50.0,
      45.0,
      35.0,
      25.0,
      25.0,
      25.0,
      20.0 ]
              if '1.15-1.35' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [       2580.0,
      2630.0,
      2000.0,
      1540.0,
      1120.0,
      880.0,
      700.0,
      580.0 ]
                    binHalfErrorExp = [       235.0,
      95.0,
      60.0,
      45.0,
      35.0,
      25.0,
      20.0,
      15.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [       1310.0,
      1720.0,
      1550.0,
      1280.0,
      930.0,
      710.0,
      560.0,
      420.0 ]
                    binHalfErrorExp = [       125.0,
      65.0,
      45.0,
      35.0,
      30.0,
      20.0,
      15.0,
      20.0 ]
              if '1.35-1.55' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [       2570.0,
      2290.0,
      1690.0,
      1160.0,
      800.0,
      630.0,
      500.0,
      370.0 ]
                    binHalfErrorExp = [       245.0,
      110.0,
      60.0,
      45.0,
      30.0,
      20.0,
      15.0,
      15.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [       1100.0,
      1560.0,
      1380.0,
      940.0,
      680.0,
      490.0,
      350.0,
      230.0 ]
                    binHalfErrorExp = [       130.0,
      90.0,
      50.0,
      35.0,
      25.0,
      15.0,
      15.0,
      10.0 ]
              if '1.55-1.75' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [       2580.0,
      1960.0,
      1390.0,
      890.0,
      620.0,
      430.0,
      320.0,
      230.0 ]
                    binHalfErrorExp = [       250.0,
      85.0,
      50.0,
      35.0,
      25.0,
      15.0,
      10.0,
      10.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [       1190.0,
      1390.0,
      1220.0,
      730.0,
      510.0,
      350.0,
      230.0,
      140.0 ]
                    binHalfErrorExp = [       125.0,
      65.0,
      40.0,
      25.0,
      20.0,
      15.0,
      10.0,
      10.0 ]
              if '1.75-1.95' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [       2020.0,
      1560.0,
      1010.0,
      570.0,
      420.0,
      320.0,
      230.0,
      170.0 ]
                    binHalfErrorExp = [       145.0,
      50.0,
      30.0,
      25.0,
      15.0,
      10.0,
      10.0,
      5.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [       1200.0,
      1170.0,
      830.0,
      470.0,
      320.0,
      210.0,
      140.0,
      80.0 ]
                    binHalfErrorExp = [       95.0,
      35.0,
      25.0,
      20.0,
      10.0,
      10.0,
      10.0,
      5.0 ]
              if '1.95-2.15' in mydict[ptag]: # parameters_test:
                 if 'pi-' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [       1700.0,
      1260.0,
      650.0,
      380.0,
      280.0,
      230.0,
      160.0,
      120.0 ]
                    binHalfErrorExp = [       120.0,
      35.0,
      25.0,
      15.0,
      10.0,
      10.0,
      5.0,
      5.0 ]
                 elif 'pi+' in secondary:
                    nBinsExp = [8]
                    binEdgeLowExp = [ 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45 ]
                    binEdgeHighExp = [ 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5 ]
                    binContentExp = [       930.0,
      860.0,
      590.0,
      310.0,
      180.0,
      120.0,
      100.0,
      50.0 ]
                    binHalfErrorExp = [       70.0,
      20.0,
      20.0,
      15.0,
      10.0,
      5.0,
      5.0,
      5.0 ]


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
                            testName="test19-HARP-Tantalum",
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
