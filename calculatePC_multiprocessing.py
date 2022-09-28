#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 13:43:10 2020

@author: Rebecca
"""
import NZfunctions
import simpleexample3
import numpy as np
import matplotlib
import pickle
from obspy.clients.fdsn import Client
from multiprocessing import Pool

matplotlib.use('Agg')


'''
runs makeTemplate2.NZLoadData() which finds data ('HH' channel only), catalog of events, the number of earthquakes in the time window, and list of folders containing seismograms (not empty folders)

for each earthquake, find the templates (Trace)

runs calcPC on data (seismograms) and templates (Trace). This saves results.pdf in path, and returns P, t and location of peaks. 
saves P as a .npy

'''

client = Client("IRIS")
rootPath ='/home/earthquakes1/homes/Rebecca/NZPaperData/' 

cat,numEQ = NZfunctions.NZmakeCat() 
eqFolders = NZfunctions.NZopenEQList()

failedPC = []
sequence2 = []
sequence3 = []
sequence4 = []
manyStaList = []
#eventTimes, eqLong, eqLat, eqDepth, sta = NZfunctions.NZEQDataLists(cat) 
CpcompList = []
CpstatList = []

def runPC(j):
# =============================================================================
#     stationType = 0
#     numSta = 0
#     for i in range(0,len(sta[j])):
#         numSta = numSta + len(sta[j][i])
# =============================================================================
    print(eqFolders[j])
    print(j)
    try: #completed commands 
        print('try')
        data = NZfunctions.NZLoadData(eqFolders[j],'sac')
        templates = NZfunctions.NZLoadTempl(eqFolders[j])
    except:
        print('except')
        failedPC.append(j)
    else:
        print('path')
        path = '/home/earthquakes1/homes/Rebecca/NZPaperData/PCResults/'+str(eqFolders[j])+'.pdf'
        print('simple exmaple')
        P, t, peakLocations2, peakLocations3, peakLocations4,  seqLocations2, seqLocations3, seqLocations4 = simpleexample3.calcPC(data,templates,path, eqFolders[j], False)
        # print('del')
        # delattr(P, 'd')
        with open('/home/earthquakes1/homes/Rebecca/NZPaperData/'+eqFolders[j]+'/P.pkl', 'wb') as output: 
            pickle.dump(P, output, pickle.HIGHEST_PROTOCOL) 
        print('comp')
        CpcompList.extend(P.Cp['Cpcomp'])
        print('stat')
        CpstatList.extend(P.Cp['Cpstat'])
    print(j)
    
    
    
"""~~~main code~~~"""
j = np.arange(0,len(eqFolders))
if __name__=='__main__':
    with Pool(5) as p:
        p.map(runPC,j)


        
