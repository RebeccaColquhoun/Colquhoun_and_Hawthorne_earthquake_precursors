import numpy as np
import pickle
import os

files = os.listdir(PATH)

def loadP(eqFolder):
    with open(rootPath+eqFolder+"/P.pkl", "rb") as input: 
        P = pickle.load(input)
    return P

for i in range(0, len(files)):
    print(i)
    eqName = files[i][0:17]
    P = loadP(eqName)

    sdComp = np.std(P.Cp['Cpcomp'])
    meanComp = np.mean(P.Cp['Cpcomp'])
    sdStat = np.std(P.Cp['Cpstat'])
    meanStat = np.mean(P.Cp['Cpstat'])

    peakLocations2Comp = []
    peakPC2Comp = []
    peakLocations3Comp = []
    peakPC3Comp = []
    peakLocations4Comp = []
    peakPC4Comp = []
    peakLocations2 = []
    peakPC2 = []
    peakLocations3 = []
    peakPC3 = []
    peakLocations4 = []
    peakPC4 = []
    for j in range(0, len(P.Cp['Cpstat'])):
        if P.Cp['Cpstat'][j]>(2*sdStat+meanStat):
            peakLocations2.append(j)
            peakPC2.append(P.Cp['Cpstat'][j])
        if P.Cp['Cpstat'][j]>(3*sdStat+meanStat):
            peakLocations3.append(j)
            peakPC3.append(P.Cp['Cpstat'][j])
        if P.Cp['Cpstat'][j]>(4*sdStat+meanStat):
            peakLocations4.append(j)
            peakPC4.append(P.Cp['Cpstat'][j])
    for j in range(0, len(P.Cp['Cpcomp'])):
        if P.Cp['Cpcomp'][j]>(2*sdComp+meanComp):
            peakLocations2Comp.append(j)
            peakPC2Comp.append(P.Cp['Cpcomp'][j])
        if P.Cp['Cpcomp'][j]>(3*sdComp+meanComp):
            peakLocations3Comp.append(j)
            peakPC3Comp.append(P.Cp['Cpcomp'][j])
        if P.Cp['Cpcomp'][j]>(4*sdComp+meanComp):
            peakLocations4Comp.append(j)
            peakPC4Comp.append(P.Cp['Cpcomp'][j])
    np.save('/home/earthquakes1/homes/Rebecca/NZPaperData/detectionsAll/'+eqName+'detections.npy',[peakLocations2, peakPC2, peakLocations3, peakPC3, peakLocations4, peakPC4, peakLocations2Comp, peakPC2Comp, peakLocations3Comp, peakPC3Comp, peakLocations4Comp, peakPC4Comp])
