#!/usr/bin/env python3.6

# Import externales
import numpy as np
import datetime
import obspy
import os
import copy
import sys
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import databseis
import seisproc
import scipy as sp
import csv


# Import internals
from PhaseCoherence import PhaseCoherence as PC


##########################################################################################################
def computePC(template,data,wintemp,buftemp,tlook,wlenlook,blim,reftemp,shifts,reflook):
    '''
    Compute one PC per hour
    '''

 

    # Initialise Phsae coherence object
    P = PC.PhaseCoherence('test',template,data=data) #template can be anything as long as includes wanted data and t pick
    #for i in range(0, len(P.template)):
        #print('t3', P.template[i].stats.t3)
        #P.template[i].stats.t3=3.
        #print(P.template[i].stats.t3)
#    import code
#    code.interact(local=locals())
    # Prepare data for computation
    P.PrepareData()
#    import code
#    code.interact(local=locals())
    # Define few extra parameters
    shtemp='t3'; shgrid=None 
    shtry=shifts;
    shlook='t3'
    #temp = template
    #look = datat to look at
    # Set parameters in objects
    P.setParams(reftemp, shtemp, wintemp,buftemp, reflook,shlook, wlenlook, tlook, blim, shtry, shgrid)
#    import code
#    code.interact(local=locals())
    #P.setParams(reftemp, shtemp, wintemp,buftemp, reflook,shlook, wlenlook, tlook, blim band limit, shtry, shgrid)
    # Make cross-correlation
    P.crosscorrDT()
#    import code
#    code.interact(local=locals())
    # Make tapering
    Mw  = int(np.round(P.params['wlenlook']/P.dtim)) # length of window
    tap = np.hamming(Mw).reshape(Mw,1)
#    import code
#    code.interact(local=locals())    
    # Taper firtst cross-correlation 
    P.taperCrosscorr(taper=tap)
#    import code
#    code.interact(local=locals())
    # Compute phase coherence
    P.computeCp()
#    import code
#    code.interact(local=locals())
    time = P.params['tlook']
#    import code
#    code.interact(local=locals())

    return P,time

def findPeaks(y,threshold): # uses sliding median window
    peakLoc = []
    med_filt_tr = sp.signal.medfilt(abs(y), kernel_size=5)
    spikes = abs((y - med_filt_tr)) > threshold
    for i in range(0,len(spikes)):
        if spikes[i]==True:
            peakLoc.append(i)
    return peakLoc
            
##########################################################################################################
def calcPC(D,T,path,eqName, plot):
    '''
    runs the phase coherence code for one earthquake. 
    paramaters:
        D: data stream, all seismograms of that eq
        T: template, template seismograms
        path: path to save results file to
    
    Results are stored in a dictionnary
    P.Cp['Cpstat'] for inter-station phase coherence
    P.Cp['Cpcomp'] for inter-component phase coherence
    '''
    # Define dates of interests
    #t1 = obspy.UTCDateTime(2011,7,26,17,12) #EQ time 2011-07-26 17:42:13
    #t2 = obspy.UTCDateTime(2011,7,26,18,12)

    # frequency range
    blim=[1,10]

    # window for template, 
    wintemp = [-1,2] #make -1, 2

    # we'll buffer by buftemp on either side of the template
    # the template tapers to zero with a cosine taper within the buffer
    buftemp = 0.3

    # times of window centers: every 6 seconds
    # relative to reflook + shtemp + shtry
    wlenlook=1. # Windows where PC is computed are wlenlook seconds long, could make 2
    dtl = 0.2 # Windows are separated by dtl seconds
    trange=[-3500,3600] # We look at 3600. seconds of data, relative to t3
    shift = []
    tlook = np.arange(-1800,1800,1.) 
    #tlook = [0]
    #tlook = np.arange(trange[0],trange[1],dtl)
    #print(wintemp, buftemp, wlenlook, dtl, tlook, trange)
    # need to highpass to avoid aliasing
    hpfilt=np.minimum(np.diff(wintemp)[0],wlenlook)
    hpfilt=3/hpfilt
    
    # Get data waveform
    #Data = obspy.read('/home/earthquakes1/homes/Rebecca/TaperedSynthetics/*.SAC') 
    #D = Data#[9:25]
    toPlot = D[0].copy()
    D.merge()   
    #D.trim(t1,t2) # Select time of interest
    starttime = D[0].stats.starttime.datetime
    Dtimes = toPlot.times()
    databseis.copyfromsacheader(toPlot)
    Ddates = Dtimes - toPlot.stats.t3
    # Get template 
    #Trace = obspy.read('/home/earthquakes1/homes/Rebecca/Template/uniformTemplates/*.SAC')
    #T = Trace#[9:25]
    
    #tr.stats.t3=10.
    
    # Get time shifts based on arrival time differences
    shifts = {}
    lini = 99999999.
    databseis.copyfromsacheader(T)
    databseis.copyfromsacheader(D)
    for tr in T:
        nid = tr.stats.network+'.'+tr.stats.station+'.'+tr.stats.channel[-1]
        #tr.stats.t3 = 3.
        shifts[nid] = tr.stats.t3 
        if float(tr.stats.t3 )<float(lini):
            lini=tr.stats.t3 
        #tr.stats.t3=tr.stats.sac['kt3']
    # Make sure smallest time-shift is 0.
    for k in shifts.keys():
        #shifts[k] = float(shifts[k])-float(lini)
        shifts[k] = 0.

    
# =============================================================================
#     toRemoveT = []
#     for t in range(0, len(T)): 
#         if hasattr(T[t].data,'mask')==True: 
#             toRemoveT.append(T[t]) 
#     for r in toRemoveT:
#         T.remove(r) 
#         
#     
#     toRemoveD = []
#     for d in range(0, len(D)): 
#         if hasattr(D[d].data,'mask')==True: 
#             toRemoveD.append(D[d])
#     for r in toRemoveD:
#         D.remove(r)
#     
#     if toRemoveD!=[] or toRemoveT!=[]:
#         print('masked removed')
#         print(toRemoveD)
#         print(toRemoveT)
# =============================================================================
    mskD=seisproc.prepfiltmask(D)
    mskT=seisproc.prepfiltmask(T)
    T.filter('bandpass',freqmin=hpfilt,freqmax=19, zerophase = True)
    D.filter('bandpass',freqmin=hpfilt,freqmax=19, zerophase = True) 
    seisproc.addfiltmask(D,mskD)
    seisproc.addfiltmask(T,mskT)
    starttime = D[0].stats.starttime.datetime
       
# =============================================================================
#     mskD = obspy.Stream()
#     mskT = obspy.Stream()
# #    import code
# #    code.interact(local=locals())
#     mskD=seisproc.prepfiltmask(D) 
#     mskT=seisproc.prepfiltmask(T)
# #    import code
# #    code.interact(local=locals())
#     # Filter that 
#     T.filter('bandpass',freqmin=hpfilt,freqmax=19)
#     D.filter('bandpass',freqmin=hpfilt,freqmax=19)
# #    import code
# #    code.interact(local=locals())
#     seisproc.addfiltmask(D,mskD)
#     seisproc.addfiltmask(T,mskT)
#     print('I finished filtering')
# =============================================================================
    
    # Resample data and template
    #T.interpolate(sampling_rate=40,method="lanczos") 
    #T.resample(sampling_rate=40,no_filter=True)
    #D.interpolate(sampling_rate=40,method="lanczos") 
    #D.resample(sampling_rate=40,no_filter=True)

    
    # Compute PC
    P,t   = computePC(T,D,wintemp,buftemp,tlook,wlenlook,blim,None,shifts,None)
    threshold =  max(P.Cp['Cpstat'][0:200])*2
    #print(threshold, len(D))
    #peakLocations = findPeaks(P.Cp['Cpstat'],threshold)
    sd = np.std(P.Cp['Cpstat'])
    mean = np.mean(P.Cp['Cpstat'])
    peakLocations2 = []
    peakPC2 = []  
    peakLocations3 = []
    peakPC3 = []  
    peakLocations4 = []
    peakPC4 = []  
    for i in range(0, len(P.Cp['Cpstat'])):
        if P.Cp['Cpstat'][i]>(2*sd + mean):
            peakLocations2.append(i)
            peakPC2.append(P.Cp['Cpstat'][i])
        if P.Cp['Cpstat'][i]>(3*sd + mean):
            peakLocations3.append(i)
            peakPC3.append(P.Cp['Cpstat'][i])
        if P.Cp['Cpstat'][i]>(4*sd+mean):
            peakLocations4.append(i)
            peakPC4.append(P.Cp['Cpstat'][i])
    seqLocations2 = []
    peakPCseq2 = []     
    seqLocations3 = []
    peakPCseq3 = []  
    seqLocations4 = []
    peakPCseq4 = []  
    for k in range(5, len(P.Cp['Cpstat'])-6): 
        seq2 = False
        seq3 = False
        seq4 = False

        for j in range(k-5, k):
            #print('test sequence')
            if P.Cp['Cpstat'][j]>2*sd and P.Cp['Cpstat'][i]>2*sd and j!=i: 
                seq2 = True
                #print(i,j)
                
            if P.Cp['Cpstat'][j]>3*sd and P.Cp['Cpstat'][i]>3*sd and j!=i: 
                seq3 = True
                #print(i,j)
                
            if P.Cp['Cpstat'][j]>4*sd and P.Cp['Cpstat'][i]>4*sd and j!=i: 
                seq4 = True
                #print(i,j)
                
        if seq2 == True:
            seqLocations2.append(k)
            peakPCseq2.append(P.Cp['Cpstat'][k])
        if seq3 == True:
            seqLocations3.append(k)
            peakPCseq3.append(P.Cp['Cpstat'][k])
        if seq4 == True:
            seqLocations4.append(k)
            peakPCseq4.append(P.Cp['Cpstat'][k])
            
    sdComp = np.std(P.Cp['Cpcomp'])
    meanComp = np.mean(P.Cp['Cpcomp'])
    peakLocations2Comp = []
    peakPC2Comp = []  
    peakLocations3Comp = []
    peakPC3Comp = []  
    peakLocations4Comp = []
    peakPC4Comp = []  
    for i in range(0, len(P.Cp['Cpcomp'])):
        if P.Cp['Cpcomp'][i]>(2*sd+meanComp):
            peakLocations2Comp.append(i)
            peakPC2Comp.append(P.Cp['Cpcomp'][i])
        if P.Cp['Cpcomp'][i]>(3*sd+meanComp):
            peakLocations3Comp.append(i)
            peakPC3Comp.append(P.Cp['Cpcomp'][i])
        if P.Cp['Cpcomp'][i]>(4*sd+meanComp):
            peakLocations4Comp.append(i)
            peakPC4Comp.append(P.Cp['Cpcomp'][i])
            
    seqLocations2Comp = []
    peakPCseq2Comp = []     
    seqLocations3Comp = []
    peakPCseq3Comp = []  
    seqLocations4Comp = []
    peakPCseq4Comp = []  
    for k in range(5, len(P.Cp['Cpcomp'])-6): 
        seq2Comp = False
        seq3Comp = False
        seq4Comp = False

        for j in range(k-5, k):
            #print('test sequence')
            if P.Cp['Cpcomp'][j]>2*sd and P.Cp['Cpcomp'][i]>2*sd and j!=i: 
                seq2Comp = True
                #print(i,j)
                
            if P.Cp['Cpcomp'][j]>3*sd and P.Cp['Cpcomp'][i]>3*sd and j!=i: 
                seq3Comp = True
                #print(i,j)
                
            if P.Cp['Cpcomp'][j]>4*sd and P.Cp['Cpcomp'][i]>4*sd and j!=i: 
                seq4Comp = True
                #print(i,j)
                
        if seq2Comp == True:
            seqLocations2Comp.append(k)
            peakPCseq2Comp.append(P.Cp['Cpcomp'][k])
        if seq3 == True:
            seqLocations3Comp.append(k)
            peakPCseq3Comp.append(P.Cp['Cpcomp'][k])
        if seq4 == True:
            seqLocations4Comp.append(k)
            peakPCseq4Comp.append(P.Cp['Cpcomp'][k])
    #import code
    ## Plot results
    #Tdates = [datetime.timedelta(seconds=tt) + D[0].stats.starttime.datetime for tt in t] 
    #Ddates = [datetime.timedelta(seconds=tt) + D[0].stats.starttime.datetime for tt in Dtimes] 
    
    
    #problem with D[0].times() still being masked
# =============================================================================
#     f,ax1 = plt.subplots()
#     ax2 = ax1.twinx()
#     ax1.plot(Ddates,D[0].data,'k')
#     ax1.set_ylabel('Data (counts)')
#     ax2.plot(P.params['tlook'],P.Cp['Cpstat'],'r')
#     ax2.spines['right'].set_edgecolor('r')
# 
#     ax1.set_xlabel('Time')
#     ax2.set_ylabel('Phase coherence')
#     ax2.yaxis.label.set_color('r')
#     ax2.tick_params(axis='y',colors='r')
#     f.savefig(path,bbox_inches='tight')
# =============================================================================
    plt.ioff()
    figure, ax = plt.subplots(3, sharex=True)
    figure.suptitle('Seismograph and Phase Coherence Plot \n of an Earthquake')
    #import code
    #code.interact(local=locals())
    #DataPlot = obspy.read('/home/earthquakes1/homes/Rebecca/TaperedSynthetics/*.SAC') 
    #databseis.copyfromsacheader(DataPlot)
    #DPlotdates = DataPlot[0].times() - DataPlot[0].stats.t3
    #print(DPlotdates)
    #print(DataPlot[0].data)
    gain = 2.51658E+03
    #toPlot = DataPlot[0].data/gain
    ax[0].plot(Ddates,toPlot.data/gain,'k')
    #ax[0].ticklabel_format(style='sci')
    #ax[0].set_xlabel('Time')
    
    ax[1].plot([-1800,1800],[1*sd,1*sd])
    ax[1].plot([-1800,1800],[2*sd,2*sd])
    ax[1].plot([-1800,1800],[3*sd,3*sd])
    ax[1].plot([-1800,1800],[4*sd,4*sd])
    ax[1].plot(P.params['tlook'],P.Cp['Cpstat'],'r')
    ax[1].set_xlabel('Time relative to earthquake P wave arrival (s)')
    ax[1].set_ylabel('Interstation \n Phase coherence')
    ax[2].plot([-1800,1800],[1*sdComp,1*sdComp])
    ax[2].plot([-1800,1800],[2*sdComp,2*sdComp])
    ax[2].plot([-1800,1800],[3*sdComp,3*sdComp])
    ax[2].plot([-1800,1800],[4*sdComp,4*sdComp])
    ax[2].plot(P.params['tlook'],P.Cp['Cpcomp'],'b')
    ax[2].set_xlabel('Time relative to earthquake P wave arrival (s)')
    ax[2].set_ylabel('Intercomponent \n Phase coherence')
    #ax[0].set_ylim(-10000,10000)
    ax[0].set_ylabel('Displacement Rate \n ($×10^{-6}$ms$^{-1}$)')
# =============================================================================
#     ax[0].set_ylim(-10000,10000)
#     ax[0].set_ylabel('Displacement Rate \n ($× 10^{-5} $ms$^{-1}$)')
# =============================================================================
    #import code
    #code.interact(local=locals())
    if plot ==True:
        plt.show()
    figure.savefig(path,bbox_inches='tight')
    #plt.show()
    plt.close()
    toWrite = [eqName, toPlot.stats.starttime, toPlot.stats.t3, tlook, sd, peakLocations2, peakPC2, peakLocations3, peakPC3,peakLocations4, peakPC4, seqLocations2, peakPCseq2,seqLocations3, peakPCseq3,seqLocations4, peakPCseq4, peakLocations2Comp, peakPC2Comp, peakLocations3Comp, peakPC3Comp,peakLocations4Comp, peakPC4Comp, seqLocations2Comp, peakPCseq2Comp,seqLocations3Comp, peakPCseq3Comp,seqLocations4Comp, peakPCseq4Comp]
    

    np.save('/home/earthquakes1/homes/Rebecca/NZPaperData/attributes/'+eqName+'attributes.npy',toWrite)
    print('saved to /home/earthquakes1/homes/Rebecca/NZPaperData/attributes/'+eqName+'attributes.npy')


    

    
# =============================================================================
# import databseis
# Data = obspy.read('/home/earthquakes1/homes/Rebecca/TaperedSynthetics/*.SAC')                                                                                               
# D = Data                                                                                                                                                                    
# Dtimes = D[0].times()  
# databseis.copyfromsacheader(D)                                                                                                                                                    
# Ddates = Dtimes - D[0].stats.t3  
# 
# plt.plot(P.params['tlook'],P.Cp['Cpstat'],'r')  
# plt.xlabel('Time')
# plt.ylabel('Phase coherence')
# plt.savefig('SynthPCResults.pdf',bbox_inches='tight')
# plt.show()
# plt.plot(Ddates,D[0].data,'k')   
# plt.xlabel('Time')
# plt.ylabel('Data Counts') 
# plt.savefig('SynthSeismogramResults.pdf',bbox_inches='tight')
# plt.show()                                                                                                                            
# =============================================================================
    return P, t, peakLocations2, peakLocations3, peakLocations4,  seqLocations2, seqLocations3, seqLocations4


def makeSynth(data553, data1005):
    synth = data553.copy()       
    synth[0].data[260000:290000]=data1005[0].data[360000:390000].copy()
    synth[1].data[260000:290000]=data1005[1].data[360000:390000].copy()
    synth[2].data[260000:290000]=data1005[2].data[360000:390000].copy()
    synth[3].data[260000:290000]=data1005[3].data[360000:390000].copy()
    synth[4].data[260000:290000]=data1005[4].data[360000:390000].copy()
    synth[5].data[260000:290000]=data1005[5].data[360000:390000].copy()
    synth[6].data[260000:290000]=data1005[6].data[360000:390000].copy()
    synth[7].data[260000:290000]=data1005[7].data[360000:390000].copy()
    synth[8].data[260000:290000]=data1005[8].data[360000:390000].copy()
    synth[9].data[260000:290000]=data1005[12].data[360000:390000].copy()
    synth[10].data[260000:290000]=data1005[13].data[360000:390000].copy()
    synth[11].data[260000:290000]=data1005[14].data[360000:390000].copy()
    return synth
def makeSynth2(data553, data1005,a):
    synth = data553.copy()       
    synth[0].data[260000:290000]=a*data553[0].data[360000:390000].copy()
    synth[1].data[260000:290000]=a*data553[1].data[360000:390000].copy()
    synth[2].data[260000:290000]=a*data553[2].data[360000:390000].copy()
    synth[3].data[260000:290000]=a*data553[3].data[360000:390000].copy()
    synth[4].data[260000:290000]=a*data553[4].data[360000:390000].copy()
    synth[5].data[260000:290000]=a*data553[5].data[360000:390000].copy()
    synth[6].data[260000:290000]=a*data553[6].data[360000:390000].copy()
    synth[7].data[260000:290000]=a*data553[7].data[360000:390000].copy()
    synth[8].data[260000:290000]=a*data553[8].data[360000:390000].copy()
    synth[9].data[260000:290000]=a*data553[9].data[360000:390000].copy()
    synth[10].data[260000:290000]=a*data553[10].data[360000:390000].copy()
    synth[11].data[260000:290000]=a*data553[11].data[360000:390000].copy()
    return synth
    
