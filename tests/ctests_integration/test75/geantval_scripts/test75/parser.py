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
    TEST = "test75"

    def parse(self, jobs):
        # should be just one job for each job
        job = jobs[0]
        print("path = ", job["path"])
	
        material = job["MATERIAL"]
        if 'G4_' in job["MATERIAL"]:
           material = (job["MATERIAL"]).split("_")[-1] # [-1] means counting from right

        energy = job["ENERGY"]
        if not 'MeV' in job["ENERGY_UNIT"]:
           en = float(energy) * 1000.
           energy = str(en)
        if '.' in energy:
           energy = energy.split(".")[0]

        inspireid = -1
        if '300' in energy:
           inspireid = 180991
        elif '668' in energy:
           inspireid = 140013

        fname = job["PARTICLE"] + energy + "MeV" + job["MATERIAL"] + job["GENERATOR"] + ".root"
        filename = os.path.join(job["path"], fname)
# -->         print("filename = ", filename)

        hfile = _open_file(filename)

        # FIXME !!! Also need to skip those theta that we don't have data for
        for h in hfile.GetListOfKeys():
           if '300' in energy:
              if 'pip' in h.GetName() or 'pim' in h.GetName():
                 print("Skip histo ", h.GetName())
                 continue
              if '60' in h.GetName() or '72' in h.GetName() or '84' in h.GetName():
                 print("Skip histo ", h.GetName()) 
                 continue
              if 'proKE' in h.GetName() or 'Cos' in h.GetName():
                 print("Skip histo ", h.GetName())
                 continue
              secondary = 'proton'
           elif '668' in energy:
              if 'Pb' in material and '28' in h.GetName(): 
                 print("Skip histo ", h.GetName())
                 continue
              if 'pip' in h.GetName():
                 secondary = "pi+"
              elif 'pim' in h.GetName():
                 secondary = "pi-"
              else:
                 print("Skip histo ", h.GetName()) 
                 continue
           else:
              print(" Unknown study case; skip")
              continue

           parameter_name = []
           parameter_name.append("THETA")
           parameter_value = []
#           parvalue = filter(str.isdigit, h.GetName())
           parvalue = ""
           for c in h.GetName():
              if c.isdigit():
                 parvalue += c
           if '28' in parvalue:
              parvalue += ".4"
           if '44' in parvalue:
              parvalue += ".2" 
           parameter_value.append( parvalue )
   
           hh = hfile.Get( h.GetName() );
   	   
           firstNonZeroBin = hh.FindFirstBinAbove(0)
           lastNonZeroBin  = hh.FindLastBinAbove(0)

           nBins = []
           binEdgeLow = []
           binEdgeHigh = []
           binContent = []
           yStatErrorPlus = []
           yStatErrorMinus = []
           beamenergies = []
	   en_to_json = float(energy)
           if 'MeV' in job["ENERGY_UNIT"]:
              en_to_json = en_to_json / 1000.
           # beamenergies.append(float(job['ENERGY']))
           # beamenergies.append(float(energy))
           beamenergies.append(en_to_json)

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
    
           parameters = []
           parameters_test = {"names": parameter_name[0], "values": parameter_value[0] + " deg"}
           parameters.append(parameters_test)

           htitle = job["PARTICLE"] + " + " + material + " -> " + secondary + " + X" 
   
           observableName = "D2(SIG)/DP/DOMEGA"
           xAxisTitle = "p (MeV/c)"
           yAxisTitle = "D2(SIG)/DP/OMEGA (mub/sr/MeV)"
           if 'proton' in secondary:
              observableName = "D2(SIG)/DE/DOMEGA"
              xAxisTitle = "EKin (MeV)"
              yAxisTitle = "D2(SIG)/DE/DOMEGA (mub/sr/MeV)"
    
           if '300' in energy:
              if '45' in parvalue:
                 nBinsExp = [17]
                 binEdgeLowExp = [44.65,  49.35,  54.3,  59.7,  65.6,  72.1,  79.2,  86.95,  
		                  95.45,  104.75,  114.9,  126,  138.1,  151.25,  165.55,  181.15,  198.1]
                 binEdgeHighExp = [49.35,  54.3,  59.7,  65.6,  72.1,  79.2,  86.95,  95.45,  
		                   104.75,  114.9,  126,  138.1,  151.25,  165.55,  181.15,  198.1,  215.7]
                 binContentExp = [29,  30,  25,  33,  23,  15.5,  13.7,  10.7,
		                  11.2,  11,  8.9,  7.4,   7.1,  5.5,  4,  3.4,  2.5]
                 binHalfErrorExp = [3.16705,  2.25887,  2.18303,  2.30955,  2.15593,  0.966595,  
		                    1.01976,  0.665395,  1.0277,  0.488083,  0.999766,  0.519212,  
				    0.389554,  0.244041,  0.331059,  0.191471,  0.132877]
              elif '90' in parvalue:
                 binEdgeLowExp = [44.65,  49.35,  54.3,  59.7,  65.6,  72.1,  79.2,  86.95,  
		                  95.45,  104.75,  114.9,  126,  138.1,  151.25,  165.55,  181.15,  198.1]
                 binEdgeHighExp = [49.35,  54.3,  59.7,  65.6,  72.1,  79.2,  86.95,  95.45,  
		                   104.75,  114.9,  126,  138.1,  151.25,  165.55,  181.15,  198.1,  215.7]
                 binContentExp = [19.8,  14,  14,  11.6,  10.8,  10.4,  7.7,  6.8,  
		                  8.6,  4.9,  3.1,  3.3,  2.4,  1.9,  1.19,  0.65,  0.27]
                 binHalfErrorExp = [1.05842,  0.774661,  1.43618,  0.987338,  1.16314,  0.578788,  
		                    0.568006,  0.382941,  0.949,  0.576118,  0.414454,  0.152775,  
				    0.130599,  0.10403,  0.0544033,  0.0302913,  0.0137587]
              elif '135' in parvalue:
                 nBinsExp = [16]
                 binEdgeLowExp = [44.65,  49.35,  54.3,  59.7,  65.6,  72.1,  79.2,  86.95,  
		                  95.45,  104.75,  114.9,  126,  138.1,  151.25,  165.55,  189.95]
                 binEdgeHighExp = [49.35,  54.3,  59.7,  65.6,  72.1,  79.2,  86.95,  95.45,  
		                   104.75,  114.9,  126,  138.1,  151.25,  165.55,  189.95,  223.85]
                 binContentExp = [11.6,  9.8,  10.5,  8.4,  9.9,  8,  8.5,  5.06,  
		                  3.31,  3.04,  2.28,  1.42,  1.04,  0.51,  0.288,  0.039]
                 binHalfErrorExp = []
              elif '150' in parvalue:
                 nBinsExp = [5]
                 binEdgeLowExp = [42.05,  51.95,  62.8,  75.8,  91.4]
                 binEdgeHighExp = [51.95,  62.8,  75.8,  91.4,  108.4]
                 binContentExp = [12.4,  7.2,  6.2,  4.9,  2.5]
                 binHalfErrorExp = [0.781573,  0.35497,  0.295108,  0.227843,  0.125524]
           if '668' in energy:
              if 'Cu' in material:
                 if 'pi-' in secondary:
                    if '28' in parvalue:
                       nBinsExp = [17]
                       binEdgeLowExp = [422,  436,  450,  464.5,  478.5,  488.5,  500,  516.5,  
		                        531.5,  543.5,  558,  576,  595,  612,  619.5,  629.5,  649]
                       binEdgeHighExp = [436,  450,  464.5,  478.5,  488.5,  500,  516.5,  531.5,  
		                         543.5,  558,  576,  595,  612,  619.5,  629.5,  649,  669]
                       binContentExp = [0.84,  0.17,  0.8,  0.86,  1.26,  1,  1.54,  2.05,  
		                        1.42,  1.65,  2.44,  1.98,  1.37,  0.82,  0.76,  0.36,  0.24]
                       binHalfErrorExp = [0.31726,  0.295099,  0.267434,  0.248038,  0.266111,  0.318198,  
		                          0.254612,  0.247811,  0.214728,  0.275205,  0.210846,  0.178784,  
					  0.139376,  0.0972708,  0.096279,  0.0963714,  0.0511531]
                    elif '44' in parvalue:
                       nBinsExp = [15]
                       binEdgeLowExp = [423.5,  436.5,  450.5,  465.5,  474.5,  480.5,  488,  499,  
		                        515.5,  531,  543.5,  557.5,  575,  594,  611.5]
                       binEdgeHighExp = [436.5,  450.5,  465.5,  474.5,  480.5,  488,  499,  515.5,  
		                         531,  543.5,  557.5,  575 , 594,  611.5,  626.5]
                       binContentExp = [0.51,  1.19,  1.15,  0.79,  0.81,  1.03,  1.61,  0.93,  
		                        1.07,  1.03,  0.96,  0.69,  0.48,  0.31,  0.13]
                       binHalfErrorExp = [ 0.146805,  0.145233,  0.135289,  0.115602,  0.125414,  0.152228,  
		                           0.131716,  0.10381,  0.0891259,  0.105704,  0.0865519,  0.0675581,  
					   0.0454594,  0.0286287,  0.0305651]
                 elif 'pi+' in secondary:
                    if '28' in parvalue:   
                       nBinsExp = [17]
                       binEdgeLowExp = [422,  436,  450,  465,  479,  488,  499,  515.5,  
		                        531,  543.5,  557.5,  575,  594,  611,  619,  629,  648]
                       binEdgeHighExp = [436,  450,  465,  479,  488,  499,  515.5,  531,  
		                         543.5,  557.5,  575,  594,  611,  619,  629,  648,  668]
                       binContentExp = [1.16,  1.05,  0.78,  0.97,  0.95,  1.37,  1.01,  1.98,  
		                        2.02,  1.61,  1.52,  1.47,  1.31,  1.32,  0.75,  0.25,  0.14]
                       binHalfErrorExp = [0.334103,  0.31358,  0.287153,  0.263639,  0.362529,  0.276949,  
		                          0.24918,  0.242,  0.318257 , 0.171096,  0.138125,  0.115762,  
					  0.0915156,  0.103698,  0.0732398,  0.0367636,  0.0257816]
                    elif '44' in parvalue:   
                       nBinsExp = [15]
                       binEdgeLowExp = [423,  437,  451,  465.5,  475,  481.5,  489,  500,  
		                        516.5,  532,  544,  557.5,  575,  594,  612]
                       binEdgeHighExp = [437,  451,  465.5,  475,  481.5,  489,  500,  516.5,  
		                         532,  544,  557.5,  575,  594,  612,  628]
                       binContentExp = [0.92,  1,  0.87,  1.13,  1.15,  1.26,  1.11,  0.92,  
		                        1.33,  1.19,  1.09,  1.08,  0.79,  0.42,  0.17]
                       binHalfErrorExp = [0.155608,  0.142302,  0.130987,  0.125741,  0.130683,  0.160359,  
		                          0.12081,  0.0990654,  0.0999101,  0.113435,  0.106915,  0.0893698,  
					  0.0654889,  0.039777,  0.03096]
              elif 'Pb' in material:
                 if 'pi-' in secondary:
                    nBinsExp = [15]
                    binEdgeLowExp = [423,  437,  451,  465.5,  475,  481,  488.5,  500,  
		                     516.5,  531.5,  544,  559,  577,  596,  613.5]
                    binEdgeHighExp = [437,  451,  465.5,  475,  481,  488.5,  500,  516.5,  
		                      531.5,  544,  559,  577,  596,  613.5,  628.5]
                    binContentExp = [2.27,  1.51,  2.4,  2.38,  2.7,  2.02,  2.43,  3.28,  
		                     3.07,  2.64,  2.6,  1.58,  1.43,  0.72,  0.2]
                    binHalfErrorExp = [0.446833,  0.420526,  0.399861,  0.351704,  0.300645,  0.444396,  
		                       0.263737,  0.264926,  0.226904,  0.245639,  0.236038,  0.175086,  
				       0.13178,  0.0771347,  0.0904489]  
                 elif 'pi+' in secondary:
                    nBinsExp = [15]
                    binEdgeLowExp = [422,  436,  450,  464.5,  474.5,  481,  488.5,  500,  
		                     516.5,  532,  544.5,  558.5,  576.5,  595.5,  613]
                    binEdgeHighExp = [436,  450,  464.5,  474.5,  481,  488.5,  500,  516.5,  
		                      532,  544.5,  558.5,  576.5,  595.5,  613,  629]
                    binContentExp = [2.32,  2.27,  1.55,  2.41,  1.39,  1.63,  3.05,  2.82,  
		                     2.11,  2.27,  2.01,  1.52,  1.61,  1.26,  0.35]
                    binHalfErrorExp = [0.341357,  0.312185,  0.293411,  0.277104,  0.282024,  0.333174,  
		                       0.28959,  0.253976,  0.207944,  0.260834,  0.224067,  0.192558,  
				       0.144478,  0.0980555,  0.0864469]  

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
                            testName="test75",
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
     
# --> exp data already in geant-val -->           yield rjsonexp
           yield rjson

# =========================================
#
# root functions
#


def _open_file(file_path):
    # Open as ROOT file
    try:
        f = TFile.Open(file_path)
        if f and f.IsOpen():
            if f.GetListOfKeys().GetEntries() == 0:
                return None
            return f
        else:
            return None
    except Exception as e:
        print("Exception: ", str(e))
        # ROOT related exception
        return None


def _close_file(root_file):
    try:
        root_file.Close()
        del root_file
    except:
        # ROOT related exception
        return None
