import obspy
import databseis
import plotseis
import seisproc
import os
import numpy as np

from obspy.geodetics import locations2degrees as loc2deg
import obspy.signal.filter
from obspy.signal.trigger import classic_sta_lta 
#cd 
from obspy.core.event import read_events 
from obspy.taup import TauPyModel 
model = TauPyModel(model="ak135")
from obspy.clients.fdsn import Client 
client = Client("IRIS")                                                                                                                                                                                                                               


rootPath ='/home/earthquakes1/homes/Rebecca/NZData2/'
    
t1 = obspy.UTCDateTime("2000-01-01")
t2 = obspy.UTCDateTime("2020-01-01")
  
def cutTemplate(st, traceNo, templStart,templEnd, eventTime):
    '''
    makes a copy and cuts it to only window of interest. then plots it 
    args: 
        st: seismograms (stream) MUST BE A STREAM NOT A TRACE
        traceNo: which trace of st to cut
        templStart: time to cut before EQ
        templEnd: time to cut after EQ
        eventTime:time of earthquake
    returns:
        templ: cut st[traceNo]
    '''
    trace = st[traceNo].copy()
    templ = seisproc.trimandshift(trace,(eventTime+templStart), (eventTime+templEnd))

    return templ 

def plotTempl(template):    
    template.plot()
    

def arrivalTimeCalculator(sta, eqLat,eqLong,eqDepth):
    
    stationLat = sta.latitude
    stationLong = sta.longitude
    
    angDist = loc2deg(eqLat,eqLong,stationLat,stationLong)
    arrivals =  model.get_travel_times(source_depth_in_km = eqDepth, distance_in_degree = angDist, phase_list = ["p","P"]) #if change model need to check if need to change phase list particularly p/P
    arrivalTime =  arrivals[0].time #if only interested in 1 phase
    return arrivalTime

def NZLoadData(pathToEQ, subfolder):
    ''' 
    Loads seismograms for 1 EQ from memory discarding any which are not channel HH*. 
    args: 
        pathToEQ - takes one item in EQ folders
        subfolder = data source to load from e.g. sac, raw
    Returns:
        data: obspy Stream with all seismograms with HH channel of that eq
        cat: EQ catalog between t1 and t2
    '''
 
    path = rootPath + pathToEQ + '/'+subfolder+'/*'
    st = obspy.read(path)
    data = obspy.Stream() 
    for i in range(0,len(st)): 
        if st[i].stats['channel'][0:2]=='HH':# and st[i].stats['network']=='NZ':
            tr = st[i] 
            if hasattr(tr.data,'mask')==False:
                data = data + tr
                        
    return data

 
def NZLoadTempl(pathToEQ):
    ''' 
    Loads seismograms for 1 EQ from memory discarding any which are not channel HH*. 
    args: 
        pathToEQ - takes one item in EQ folders
    Returns:
        data: obspy Stream with all seismograms with HH channel of that eq
        cat: EQ catalog between t1 and t2
    '''
 
    path = rootPath + pathToEQ + '/Templates/*'
    t = obspy.read(path)
    templ = obspy.Stream()
    for i in range(0,len(t)): 
        if t[i].stats['network']=='NZ':
            templ = templ + t[i]
    
                        
    return templ

def NZmakeCat():
    '''      
    Downloads EQ catalog between t1 and t2
    Returns:
    data: obspy Stream with all seismograms with HH channel
    cat: EQ catalog between t1 and t2
    '''
    cat = read_events('/home/earthquakes1/homes/Rebecca/NZData2/EVENTS-INFO/catalog.ml')
    #obspyDMT --datapath NZDataJuneJuly2011-2 --event_catalog IRIS --data_source IRIS --event_rect 175/180/-41/-37 --min_mag 4.0 --min_date 2011-06-01 --max_date 2011-08-01  --max_epi 4 --net "NZ", "II", 'IU' 
    noEQ=len(cat)                            
    return cat, noEQ

def NZmakeFolderList():
    #makes a complete list. doesn't exclude those folders where there is no STALTA pick. these are excluded when running the makeTemplate command
    subfolders = os.listdir(rootPath)
    eqFolders=[]
    empty = []
    for k in range(0,len(subfolders)):
        if subfolders[k][0]=='2':
            dirContents = os.listdir(rootPath+subfolders[k])
            eqFolders.append(subfolders[k])
        else:
            print('failed', subfolders[k])
    with open('eqFoldersList2.txt', 'w') as f:
        for item in eqFolders:
            f.write("%s\n" % item)
    eqFolders.sort(reverse=True)
    return eqFolders, subfolders, empty      
def NZopenEQList():
    eqFolders = [line.rstrip('\n') for line in open('eqFoldersList2.txt')]
    eqFolders.sort(reverse=True)
    return eqFolders

def NZSortEQData(eventList):
    '''
    Creates lists where each element is info on a specific earthquake
    Args:
        eventList: catalog of all events
        
    Returns:
        eventTimes: list of times of each earthquake
        eqLong: list of longitudes of each earthquake
        eqLat: list of latitudes of each earthquake
        eqDepth: list of depths of each earthquake
        stations: list of stations within 4 degrees of each earthquake
    '''    

    eventTimes=[]
    earthquakeLong = []
    earthquakeLat =[]
    earthquakeDepth=[]
    stations = []
    
    for i in range(0, len(eventList)):
        print(i) #for each earthquake
        eventTimes.append(eventList[i].__dict__['origins'][0].time)
        earthquakeLong.append(eventList[i].__dict__['origins'][0].longitude)
        earthquakeLat.append(eventList[i].__dict__['origins'][0].latitude)
        earthquakeDepth.append((eventList[i].__dict__['origins'][0].depth)/1000)
        stations.append(client.get_stations(network='NZ', channel = 'HH*', starttime=t1, endtime=t2, latitude = earthquakeLat[i], longitude = earthquakeLong[i], maxradius = 4))

        
        
    
    return eventTimes, earthquakeLong, earthquakeLat, earthquakeDepth, stations
def CASortEQData(eventList):
    '''
    Creates lists where each element is info on a specific earthquake
    Args:
        eventList: catalog of all events
        
    Returns:
        eventTimes: list of times of each earthquake
        eqLong: list of longitudes of each earthquake
        eqLat: list of latitudes of each earthquake
        eqDepth: list of depths of each earthquake
        stations: list of stations within 4 degrees of each earthquake
    '''    

    eventTimes=[]
    earthquakeLong = []
    earthquakeLat =[]
    earthquakeDepth=[]
    stations = []
    t1 = obspy.UTCDateTime("2016-03-11")
    t2 = obspy.UTCDateTime("2016-03-13")
    
    for i in range(0, len(eventList)): #for each earthquake
        eventTimes.append(eventList[i].__dict__['origins'][0].time)
        earthquakeLong.append(eventList[i].__dict__['origins'][0].longitude)
        earthquakeLat.append(eventList[i].__dict__['origins'][0].latitude)
        earthquakeDepth.append((eventList[i].__dict__['origins'][0].depth)/1000)
        stations.append(client.get_stations( channel = "HH*", starttime=t1, endtime=t2, latitude = earthquakeLat[i], longitude = earthquakeLong[i], maxradius = 2))
    return eventTimes, earthquakeLong, earthquakeLat, earthquakeDepth, stations 


def STALTA(numEQ, stream):
    trigger = [] 
    exceeds = []
    #earthquakes
    data = stream 

    for i in range(0, len(data),3):#stations 
        exceeded = False
        trigger.append([]) #at each station
        for k in range(0,3): #components
            trace = data[i+k] 
            databseis.copyfromsacheader(trace) 
            short = cutTemplate(stream, k, -150, 150, trace.stats.starttime+trace.stats.t3) 
            trigger[int(i/3)].append([])
            df = trace.stats.sampling_rate
            cft = classic_sta_lta(short.data, int(1 * df), int(10 * df)) 
            for l in range(0,len(cft)):  
                if cft[l]>6 and trigger[int(i/3)][k] == []:  
                    trigger[int(i/3)][k].append(round((l/100)-150,2))
                    if exceeded == False:
                        exceeds.append([i, k])
                        exceeded = True
                             
             
                             
                             
                    
        
def NZTemplate(traces,eventTimes, eqLong, eqLat, eqDepth, sta, eqFolder):

    '''
    Args: 
        traces: all seismograms
        eventList: catalog of events
    '''

    traceNames=[]  
    arrivalTime=[] 
    exceeds = [] 
    blankSTALTA = []
    partlyBlankSTALTA = []


    tracesOfOneEQ = traces

    for i in range(0, len(tracesOfOneEQ), 3): #for each station
        traceNames.append(tracesOfOneEQ[i].stats['network']+'.'+tracesOfOneEQ[i].stats['station']) #gradually compiles list of station names.
        if (len(tracesOfOneEQ[i]) == 720000) and (len(tracesOfOneEQ[i+1]) == 720000)  and (len(tracesOfOneEQ[i+2]) == 720000):
            '''station list'''
             
            stationType = 0
            stationNumber = 0
            if len(sta) ==1:
                stationNumber = int(i/3)
            elif len(sta) ==2:
                if int(i/3)>len(sta[0]):
                    stationType = 1
                    stationNumber = int(i/3)-len(sta[0])
            elif len(sta) ==3:
                if int(i/3)>len(sta[0]) and int(i/3)<len(sta[1]):
                    stationType = 1
                    stationNumber = int(i/3)-len(sta[0])
                elif int(i/3)>len(sta[1]):
                    stationType = 2
                    stationNumber = int(i/3)-len(sta[0])-len(sta[1])
                    
            TravelTime = arrivalTimeCalculator(sta[stationType][stationNumber],eqLat,eqLong,eqDepth) #work out travel time for this eq to this station

    
            arrivalTime = eventTimes+TravelTime #arrival time of this eq at this station

            '''load data and resave as SAC'''
            E = tracesOfOneEQ[i]
            N = tracesOfOneEQ[i+1]
            Z = tracesOfOneEQ[i+2]
            SACstream=obspy.core.Stream(traces=[N,E,Z])
            
            for tr in SACstream:
                if isinstance(tr.data, np.ma.masked_array):
                    tr.data = tr.data.filled()
            print(traceNames[int(i/3)]+ 'saved to hold data')
            SACstream.write("/home/earthquakes1/homes/Rebecca/NZData2/holdData/"+traceNames[int(i/3)]+".test.SAC", format = "SAC")
            
            
            '''load SAC data'''
            sacData = obspy.read("/home/earthquakes1/homes/Rebecca/NZData2/holdData/*")
            
            filename = ''
            #empties hold data
            for filename in os.listdir("/home/earthquakes1/homes/Rebecca/NZData2/holdData/"):
                os.unlink("/home/earthquakes1/homes/Rebecca/NZData2/holdData/"+filename)
                
            startTime = E.stats['starttime']
            
            '''filter and find and set kt3'''
            N = sacData[1]
            N.data = obspy.signal.filter.highpass(N.data, 1.5, corners = 1, df = 100)  
            E = sacData[0]
            E.data = obspy.signal.filter.highpass(E.data, 1.5, corners = 1, df = 100)  
            Z = sacData[2]
            Z.data = obspy.signal.filter.highpass(Z.data, 1.5, corners = 1, df = 100)
    
                        
            forSTA=obspy.core.Stream(traces=[N,E,Z])
            
            
            trigger = []
            for k in range(0,3): #components
                 trace = forSTA[k]
                 df = trace.stats.sampling_rate
                 short = cutTemplate(forSTA, k, -150, 150, arrivalTime)
                 cft = classic_sta_lta(short.data, int(1 * df), int(10 * df))
                 trigger.append([])
                 for l in range(0,len(cft)):  
                     if cft[l]>6 and trigger[k] == []:  
                         trigger[k].append(round((l/100)-150,2))
                         exceeds.append(traceNames[int(i/3)])
        
            if trigger[0]==[] and trigger[1]==[] and trigger[2]==[]:
                blankSTALTA.append(traceNames[int(i/3)])
                continue
            elif trigger[0]==[] and trigger[1]!=[] and trigger[2]!=[] : 
                trigger[0].append(min(trigger[1][0], trigger[2][0]))
                partlyBlankSTALTA.append(traceNames[int(i/3)])
            elif trigger[1]==[] and trigger[0]!=[] and trigger[2]!=[]: 
                trigger[1].append(min(trigger[0][0], trigger[2][0]))
                partlyBlankSTALTA.append(traceNames[int(i/3)])
            elif trigger[2]==[] and trigger[1]!=[] and trigger[0]!=[]:
                trigger[2].append(min(trigger[0][0], trigger[1][0]))
                partlyBlankSTALTA.append(traceNames[int(i/3)])
            elif trigger[0]==[] and trigger[1]==[] and trigger[2]!=[] : 
                trigger[0].append(trigger[2][0])
                trigger[1].append(trigger[2][0])
                partlyBlankSTALTA.append(traceNames[int(i/3)])
            elif trigger[0]==[] and trigger[1]!=[] and trigger[2]==[]: 
                trigger[0].append(trigger[1][0])
                trigger[2].append(trigger[1][0])
                partlyBlankSTALTA.append(traceNames[int(i/3)])
            elif trigger[0]!=[] and trigger[1]==[] and trigger[2]==[]:
                trigger[1].append(trigger[0][0])
                trigger[2].append(trigger[0][0])
                partlyBlankSTALTA.append(traceNames[int(i/3)])           
            
            Nt3 = arrivalTime-N.stats.starttime
            
            Et3 = arrivalTime-E.stats.starttime
            Zt3 = arrivalTime-Z.stats.starttime
            
            N.stats.sac.__dict__.update({'kt3' : trigger[0][0]+Nt3})
            E.stats.sac.__dict__.update({'kt3' : trigger[1][0]+Et3})  
            Z.stats.sac.__dict__.update({'kt3' : trigger[2][0]+Zt3})
    
            
            streamT3=obspy.core.Stream(traces=[N,E,Z])
            if not os.path.exists("/home/earthquakes1/homes/Rebecca/NZData2/"+eqFolder+"/sac/"):
                os.makedirs("/home/earthquakes1/homes/Rebecca/NZData2/"+eqFolder+"/sac/")
            path = "/home/earthquakes1/homes/Rebecca/NZData2/"+eqFolder+"/sac/"+traceNames[int(i/3)]+".SAC"
            streamT3.write(path, format = "SAC")
            

    # =============================================================================
    #             MAKES AND WRITES TEMPLATES
    # =============================================================================
            templateE = cutTemplate(sacData, 0, -20, 20, (E.stats.starttime+E.stats.sac['kt3']))
            templateE.stats.sac.__dict__.update({'kt3' : 23.})
            templateN = cutTemplate(sacData, 1, -20, 20, (N.stats.starttime+N.stats.sac['kt3']))
            templateN.stats.sac.__dict__.update({'kt3' : 23.})
            templateZ = cutTemplate(sacData, 2, -20, 20, (Z.stats.starttime+Z.stats.sac['kt3']))
            templateZ.stats.sac.__dict__.update({'kt3' : 23.})
            
            stream=obspy.core.Stream(traces=[templateN,templateE,templateZ])

            if not os.path.exists("/home/earthquakes1/homes/Rebecca/NZData2/"+eqFolder+"/Templates/"):
                os.makedirs("/home/earthquakes1/homes/Rebecca/NZData2/"+eqFolder+"/Templates/")

            stream.write("/home/earthquakes1/homes/Rebecca/NZData2/"+eqFolder+"/Templates/"+traceNames[int(i/3)]+".Template.SAC", format = "SAC") #saves Templates

        else:
            print('too short')
    if blankSTALTA != [] and partlyBlankSTALTA !=[]:
        return blankSTALTA, partlyBlankSTALTA
    elif blankSTALTA == []:
        return partlyBlankSTALTA
    elif partlyBlankSTALTA ==[]:
        return blankSTALTA


def openAttributesFile(eqName):
    att = np.load('/home/earthquakes1/homes/Rebecca/NZData2/'+eqName+'/attributes.npy', allow_pickle = True)
    return att
