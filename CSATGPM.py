from pyorbital.orbital import Orbital
from pyorbital.tlefile import Tle
from datetime import datetime
from datetime import timedelta
import calendar
import shapely as shp
import numpy.ma as ma
import shapely as shp
import netCDF4
from haversine import haversine as hs
from sklearn.neighbors import NearestNeighbors as sknn
import h5py
import pandas as pd
import geopandas as gpd
from datetime import timedelta
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from pyhdf.SD import SD, SDC
from pyhdf import HDF, VS, V
from datetime import timedelta
import re
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib import colormaps
import numpy.ma as ma

def CSATGPM():

    ## USER INPUT

    # Provide desired time range for coincidences
    # to get coincidences for all months of provided year, leave months as empty list
    # to get coincidences for all days of provided months, leave days as None
    YEARS = []
    MONTHS = []
    DAYS = None

    # If coincidence list for desired time range is already stored in file, set to False and provide filepath
    NEED_COINCIDENCES = True
    COINCIDENCE_FILEPATH = None

    # Provide paths to downloaded files
    GMI_PATH = 'GMI_files'
    DPR_PATH = 'DPR_files'
    CS_PATH = 'CS_files'

    # Provide filepath to save CSATGPM files
    CSATGPM_PATH = ''


    CS_TLEfile = r"TLEs\CloudSat_All_TLEs.txt"
    GPM_TLEfile = r"TLEs\GPM_TLEs.txt"
    PYORBITAL_CONFIG_PATH = r'\path\platforms'

    # 


    ## 


    if NEED_COINCIDENCES:
        coincidence_list = getCoincidences(YEARS, MONTHS, DAYS, CS_TLEfile, GPM_TLEfile, PYORBITAL_CONFIG_PATH)
    
    else:
        print()
        #TODO: write code to read coincidences from file
    '''
    for coincidence in coincidence_list:
        CS_time = coincidence[3]
        GPM_time = coincidence[7]
        GMIfiles, DPRfiles, CSfiles = findFiles(CS_time, GPM_time, GMI_PATH, DPR_PATH, CS_PATH)
        map(GPM_time, CS_time, GMIfiles, CSfiles, DPRfiles, False, 1)
    ''' 





def getCoincidences(years, months, singleDay, GMI_TLEs, CS_TLEs, PlatformPath):

    # reference files
    CS_TLEfile = CS_TLEs
    GPM_TLEfile = GMI_TLEs
    PYORBITAL_CONFIG_PATH = PlatformPath
    


    RANGE = 180*5 #interval within which crossings will be considered overlaps - number of seconds in 15 minutes
    plus_1s = timedelta(seconds=1)

    #initialize list to hold CSAT,GPM pairs that have coincidences
    coincidences = []
    GPMinterval = [] #initialize list that will store GPM locations within interval (+- 15 minutes from Cloudsat)
    lenGPMinterval = 0 #number of orbits stored in GPM interval list
    
    # find numnber of minutes in the desired time range, initialize current time datetime object at start time
    year = years[0]
    month = months[0]
    if singleDay != None:
        days = 1
        ct = datetime(year,month,singleDay)
    else:
        days = 0
        for year in years:
            for month in months:
                days += (calendar.monthrange(year,month)[1])
        ct = datetime(year,month, 1)
    seconds = days * 24 * 60 * 60 # num days x 24 hours x 60 minutes x 60 seconds


    # read relevant TLEs from TLE files - returns dictionary of dates and corresponding TLEs 
    #print('Getting TLEs.')
    startTime = ct
    endTime = datetime(years[-1],months[-1],(calendar.monthrange(years[-1],months[-1])[1]))
    GPM_TLEs = readTLEs('GPM-Core',GPM_TLEfile,startTime,endTime)
    CS_TLEs = readTLEs('CloudSat',CS_TLEfile,startTime,endTime)
    


    ## populate list of GPM orbits within interval for the 15 minutes behind and ahead of starting time

    #print("Populating interval.")
    
    for i in range(RANGE): # go back 15 minutes in time
        ct -= plus_1s
    # access GPM TLE in order to create += 15 minute interval
    [currentGPMTLE, GPM_TLEindex, nextTimeGPM] = getTLE(ct, GPM_TLEs, None)
    

    GPM = Orbital('GPM-Core', None, currentGPMTLE.line1, currentGPMTLE.line2)
    for i in range(RANGE*2): # iterate through the half hour interval
        if timeTLE(ct) >= nextTimeGPM:
            [currentGPMTLE, GPM_TLEindex, nextTimeGPM] = getTLE(ct, GPM_TLEs, GPM_TLEindex)
            GPM = Orbital('GPM-Core', None, currentGPMTLE.line1, currentGPMTLE.line2)
        GPMorbit = GPM.get_lonlatalt(ct)
        GPMinterval.append([GPMorbit[0],GPMorbit[1],ct,GPM.get_orbit_number(ct)])
        ct += plus_1s


    if (len(GPMinterval) != 60*15*2): # check to make sure that interval of GPM locations populated correctly - 1 per second for +- 15 minutes
        print("wrong length of interval")
        return


    # reset time to starting time
    if singleDay != None:
        ct = datetime(year,month,day=singleDay)
    else:
        ct = datetime(year,month, 1)

    # initialize datetime object to hold the time 15 minutes ahead of current time (for maintaining GPM interval)
    ct_p15 = ct + timedelta(0,15*60)

    # access TLEs and initialize TLE files
    [currentGPMTLE, GPM_TLEindex, nextTimeGPM] = getTLE(ct, GPM_TLEs, None)
    #print(currentGPMTLE)
    GPM = Orbital('GPM-Core', None, currentGPMTLE.line1, currentGPMTLE.line2)
    [currentCSTLE, CS_TLEindex, nextTimeCS] = getTLE(ct, CS_TLEs, None)
    CS = Orbital('CloudSat', None, currentCSTLE.line1, currentCSTLE.line2)

    lenGPMinterval = len(GPMinterval)

    '''
    breaktime = 1296000 # arbitrary time to end early for testing
    if day != None:
        breaktime = 24 * 60 * 60
    '''


    recentCoincidence = False
    i_sinceCoincidence = 0


    #points = []
    #for GPM_O in GPMinterval:
        #points.append((GPM_O[0],GPM_O[1]))


    for i in range(seconds):

        #if i%(24*1800*2) == 0:
            #print("Current time: " + str(ct)) #print out time every day for testing purposes
            #print(len(coincidences))


        ## Maintain GPM +- 15 minute interval
        GPMorbit = GPM.get_lonlatalt(ct_p15)

        GPMinterval.append([GPMorbit[0],GPMorbit[1],ct_p15, GPM.get_orbit_number(ct_p15)])
        lenGPMinterval+=1
        #points.pop(0)
        #points.append((GPMorbit[0],GPMorbit[1]))

        if lenGPMinterval > (2*RANGE): #maintain +- 15 minute interval
            GPMinterval.pop(0)
            lenGPMinterval -= 1
        

        #add if statement for everything below - only check if time falls within cloudsat daytime orbit

        CSorbit = CS.get_lonlatalt(ct)

        

        



        
        if recentCoincidence == True:
            i_sinceCoincidence += 1
            if i_sinceCoincidence > 60*10: # assume there will never be two separate coincidences within 1 minute
                recentCoincidence = False
                i_sinceCoincidence = 0
        
        if recentCoincidence == False:

            
           
            #if (abs(GPMinterval[lenGPMinterval//2][0]-CSorbit[0]) < 60) and (abs(GPMinterval[lenGPMinterval//2][1]-CSorbit[1]) < 60):

                for k in range(0,lenGPMinterval,30): # loop through all GPM data points within 15 minutes in either direction of the CS scan
                    
                    GPM_O = GPMinterval[k]
                    # check to see if the lat/lon coordinates are close to one another
                    if abs(GPM_O[0]-CSorbit[0]) < 12:
                        if abs(GPM_O[1]-CSorbit[1]) < 12:

                            #print('checking all 30s')
                            #print(ct)
                            
                            for h in range(max(0,k-15),min(k+15,lenGPMinterval), 5):

                                GPM_O = GPMinterval[h]

                                if abs(GPM_O[0]-CSorbit[0]) < 6:
                                    if abs(GPM_O[1]-CSorbit[1]) < 6:
                                        
                                        for m in range(max(0,h-3),min(k+3,lenGPMinterval)):
                                        
                                            GPM_O = GPMinterval[m]

                                            if abs(GPM_O[0]-CSorbit[0]) < .05:
                                                if abs(GPM_O[1]-CSorbit[1]) < .05:
                                                    print("coincidence")
                                                    recentCoincidence = True
                                                    print(GPM_O)
                                                    # CS lat, CS lon, CS orbit number, CS time, GPM lat, GPM lon, GPM orbit number, GPM time, time offset
                                                    coincidences.append([CSorbit[0],CSorbit[1], CS.get_orbit_number(ct), ct, GPM_O[0], GPM_O[1], GPM_O[3], GPM_O[2], ct-GPM_O[2]])
                                
                                            if recentCoincidence:
                                                #print('breaking')
                                                break
                                        
                                if recentCoincidence:
                                    #print('breaking')
                                    break
                    
                    if recentCoincidence:
                        break
                                        
        


        # increment current time and time 15 minutes ahead by one second
        ct += plus_1s
        ct_p15 += plus_1s


        # get new TLE each time the current time passes the time of the next TLE
        if (nextTimeGPM != None) and (timeTLE(ct_p15) >= nextTimeGPM):
            #print("Getting GPM TLEs")
            [currentGPMTLE, GPM_TLEindex, nextTimeGPM] = getTLE(ct_p15, GPM_TLEs, GPM_TLEindex)
            #print(currentGPMTLE)
            try:
                GPM = Orbital('GPM-Core', None, currentGPMTLE.line1, currentGPMTLE.line2)
            except:
                print('failed to get new TLE')
        
        if (nextTimeGPM != None) and (timeTLE(ct) >= nextTimeCS):
            #print("Getting CS TLEs")
            [currentCSTLE, CS_TLEindex, nextTimeCS] = getTLE(ct, CS_TLEs, CS_TLEindex)
            try:
                CS = Orbital('CloudSat', None, currentCSTLE.line1, currentCSTLE.line2)
            except:
                print('failed to get new TLE')
        
        #if i > breaktime:
            #return coincidences

    return coincidences




        
def readTLEs(Sat, TLEfile, startDate, endDate):
    '''
    Takes arguments of satellite, file containing TLE entries, and start and end date
    Returns dictionary containing each date within the start and end dates and its corresponding TLE
    '''

    #need to get start and end dates into correct format to match TLEs
    startDateTLE = timeTLE(startDate-timedelta(2))
    endDateTLE = timeTLE(endDate+timedelta(days=1))

    fs = open(TLEfile, 'r')

    # initialize list to hold dates and TLEs
    TLEs = []

    TLEline = fs.readline()
    #print(TLEline)

    while not TLEline[0] == '1':#skip lines until the start of a TLE
        TLEline = fs.readline()
        #print(TLEline)

    #since TLEs are in chronological order and will be searching in chronological order, can skip those before start date
    while float(TLEline[18:23]) < (float(startDateTLE)): #find first TLE

        #skip 3 lines each time since each TLE is three lines
        TLEline = fs.readline()
        while TLEline[0] != '1':
            TLEline = fs.readline()

    sdtle = 0
    
    while float(TLEline[18:23]) < (float(endDateTLE)):

        TLE1 = TLEline # first line

        TLE2 = fs.readline() # second line

        newTLE = Tle(Sat, None, TLE1, TLE2) # create TLE object

        dt = TLE1[18:26] # two digits of year, three digits for day of year, period, two digits for fractional hour\
        #print(dt)

        TLEs.append([dt,newTLE]) # add dictionary entry of TLE for date

        TLEline = fs.readline()

        while TLEline[0] != '1':
            TLEline = fs.readline()

            #print('Num TLEs:' + str(len(TLEs)))
    return TLEs




def getTLE(currentTime, TLEList, lastIndex):
    '''
    gets the closest TLE to the date being accessed from the TLE dictionary
    '''

    timeInFormat = timeTLE(currentTime)

    i = 1 #60 sec *

    currentTLE = None

    if lastIndex != None:
        Index = lastIndex
    
    elif lastIndex == None: # if no index provided for previous TLE, loop through all TLEs
        i = 0
        
        thisTLE = (float(TLEList[i][0]))
        
        while not (thisTLE <= float(timeInFormat) < float(TLEList[i+1][0])):
            i += 1
            #print(i)
            thisTLE = (float(TLEList[i][0]))
        
        currentTLE = TLEList[i][1]
        Index = i
    

    while currentTLE == None and Index + 1 < len(TLEList):
        # if still between the dates of the most recent TLE and the next TLE, return the most recent TLE
        if float(TLEList[Index][0]) <= float(timeInFormat) < float(TLEList[Index+1][0]):
            currentTLE = TLEList[Index][1]
        Index += 1
            
    if currentTLE == None:
        currentTLE = TLEList[lastIndex]
        
    #if i > 1:
        #print("no TLE found on date. days off of date: " + str(i))
    
    try:
        nextTime = TLEList[Index+1][0]
    except:
        nextTime = None

    return [currentTLE, Index, nextTime]




def timeTLE(dt):
    '''
    takes a datetime as an argument and returns a date in the format used in TLEs
    '''
    #print(dt)

    fractionalHour = dt.hour / 24
    #print(fractionalHour)
    #print(dt.strftime('%y%j'))
    return dt.strftime('%y%j') + str(fractionalHour)[1:3]

def createFile(coincidences, name):
    

    # CS lat, CS lon, CS orbit number, CS time, GPM lat, GPM lon, GPM orbit number, GPM time, time offset
    
    f = open(name, "w")

    for coincidence in coincidences:
        GPMdate = str(coincidence[7].year) + str(coincidence[7].month).zfill(2) + str(coincidence[7].day).zfill(2)
        GPMtime = str(coincidence[7].hour).zfill(2) + str(coincidence[7].minute).zfill(2) + str(coincidence[7].second).zfill(2)
        CSdate = str(coincidence[3].year) + str(coincidence[3].month).zfill(2) + str(coincidence[3].day).zfill(2)
        CStime = str(coincidence[3].hour).zfill(2) + str(coincidence[3].minute).zfill(2) + str(coincidence[3].second).zfill(2)

        f.write('' + GPMdate + '\t' + GPMtime + '\t' + str(coincidence[4]).zfill(6)[:6] + '\t' + str(coincidence[5]).zfill(6)[:6] + '\t' + str(coincidence[6]) + '\t' + CSdate + '\t' + CStime + '\t' + str(coincidence[0]).zfill(6)[:6] + '\t' + str(coincidence[1]).zfill(6)[:6] + '\t' + str(coincidence[2])+ '\t' + str(coincidence[8]) + '\n')

    f.close()

def map(GPMtime, CStime, GMIfiles, CSATfiles, DPRfiles, plot, channel):


    # get CS and GPM data for full orbits, as well as very rough approximation of overlap region
    print('read gmi')
    coinLatArray, coinLonArray, coinTbArray, TimeArray, LatArray, LonArray, TbArray, coinLatArrayS2, coinLonArrayS2, coinTbArrayS2, startIndex, endIndex = readGMI(GMIfiles, GPMtime, False)
    
    print("read DPR")
    FS_ZFactorMeasured, FS_Lat, FS_Lon, HS_ZFactorMeasured, HS_Lat, HS_Lon, CSendIndex = readDPR(DPRfiles, GPMtime, False)

    print("read CSAT")
    [coinLat, coinLon, TimeArray, coinRadarReflectivity, coinHeight, cslatitude, cslongitude] = readCS(CSATfiles, CStime, False)
    
    
    S1startCrossing, S1endCrossing, swathS1 = findCrossing(coinLat,coinLon,coinLatArray,coinLonArray)

    
    
    # create arrays that contain only CloudSat coordinates from crossing over S1 swath - reduces computations for later crossings since S1 has widest swath
    crossingLat = coinLat[S1startCrossing:S1endCrossing+1]
    crossingLon = coinLon[S1startCrossing:S1endCrossing+1]

    if plot:
        # select channel to be plotted
        channel -= 1
        Tb = TbArray[:,:,channel]
        coinTb = coinTbArray[:,:,channel]
        S1channels = ['10.65 GHz vertically-polarized Tb', '10.65 GHz horizontally-polarized Tb', '18.7 GHz vertically-polarized Tb', '18.7 GHz horizontally-polarized Tb', '23.8 GHz vertically-polarized Tb', '36.64 GHz vertically-polarized Tb', '36.64 GHz horizontally-polarized Tb', '89.0 GHz vertically-polarized Tb', '89.0 GHz horizontally-polarized Tb']
        S1channel_name = S1channels[channel]
        coinTbS2 = coinTbArrayS2[:,:,1]
        visualize(swathS1, LonArray, LatArray, Tb, cslongitude, cslatitude, crossingLon, crossingLat, S1channel_name, coinLonArray, coinLatArray, coinTb, coinLonArrayS2, coinLatArrayS2, coinTbS2, coinLon, coinLat, CStime)
    
    
    S2startCrossing, S2endCrossing, swath = findCrossing(crossingLat,crossingLon,coinLatArrayS2,coinLonArrayS2)
    
    FSstartCrossing, FSendCrossing, swath = findCrossing(crossingLat,crossingLon,FS_Lat,FS_Lon)

    HSstartCrossing, HSendCrossing, swath = findCrossing(crossingLat,crossingLon,HS_Lat,HS_Lon)
    

    # populate list of all CloudSat points that intersect with GMI swath
    CloudSatPoints = []
    for i in range(len(crossingLat)):
        CloudSatPoints.append((crossingLat[i][0], crossingLon[i][0]))
    

    # GMI points to interpolate from as well as array containing their indexes in the original array in order to later get data
    relevantS1points_coords = []
    relevantS1points_indeces = []

    relevantS2points_coords = []
    relevantS2points_indeces = []

    relevantHSpoints_coords = []
    relevantHSpoints_indeces = []

    relevantFSpoints_coords = []
    relevantFSpoints_indeces = []

    for i in range(np.shape(coinLatArray)[0]):
        for j in range(np.shape(coinLatArray)[1]):
            relevantS1points_coords.append((coinLatArray[i][j], coinLonArray[i][j]))
            relevantS1points_indeces.append((i,j,len(relevantS1points_coords)-1))

    for i in range(np.shape(coinLatArrayS2)[0]):
        for j in range(np.shape(coinLatArrayS2)[1]):
            relevantS2points_coords.append((coinLatArrayS2[i][j], coinLonArrayS2[i][j]))
            relevantS2points_indeces.append((i,j))

    for i in range(np.shape(HS_Lat)[0]):
        for j in range(np.shape(HS_Lat)[1]):
            relevantHSpoints_coords.append((HS_Lat[i][j], HS_Lon[i][j]))
            relevantHSpoints_indeces.append((i,j))
    
    for i in range(np.shape(FS_Lat)[0]):
        for j in range(np.shape(FS_Lat)[1]):
            relevantFSpoints_coords.append((FS_Lat[i][j], FS_Lon[i][j]))
            relevantFSpoints_indeces.append((i,j))

    nnS1 = sknn(n_neighbors=9)
    nnS1.fit(relevantS1points_coords)
    S1_9nn = nnS1.kneighbors(CloudSatPoints, return_distance=True)

    nnS2 = sknn(n_neighbors=9)
    nnS2.fit(relevantS2points_coords)
    S2_9nn = nnS2.kneighbors(CloudSatPoints, return_distance=True)

    nnHS = sknn(n_neighbors=9)
    nnHS.fit(relevantHSpoints_coords)
    HS_9nn = nnHS.kneighbors(CloudSatPoints, return_distance=True)

    nnFS = sknn(n_neighbors=9)
    nnFS.fit(relevantFSpoints_coords)
    FS_9nn = nnFS.kneighbors(CloudSatPoints, return_distance=True)

    numPoints = len(CloudSatPoints)    
    
    InterpS1 = np.empty((numPoints, 9)) * np.NaN
    InterpS2 = np.empty((numPoints, 4))
    InterpS2.fill(-99.0)
    InterpFS_ZFactorMeasured = np.empty((numPoints, 176, 2)) * np.NaN
    InterpHS_ZFactorMeasured = np.empty((numPoints, 88)) * np.NaN

    
    for i in range(len(CloudSatPoints)):
        
        mindist = None
        nn = 0

        for j in range(9):
            
            dist = hs(relevantS1points_coords[S1_9nn[1][i][j]],CloudSatPoints[i])
            if mindist == None or dist<mindist:
                mindist = dist
                nn = j
        
        InterpS1[i,:] = coinTbArray[relevantS1points_indeces[S1_9nn[1][i][nn]][0],relevantS1points_indeces[S1_9nn[1][i][nn]][1]]

        if S2startCrossing <= i <= S2endCrossing:

            mindist = None
            nn = 0

            for j in range(9):
                
                dist = hs(relevantS2points_coords[S2_9nn[1][i][j]],CloudSatPoints[i])
                if mindist == None or dist<mindist:
                    mindist = dist
                    nn = j
            
            InterpS2[i,:] = coinTbArrayS2[relevantS2points_indeces[S2_9nn[1][i][nn]][0],relevantS2points_indeces[S2_9nn[1][i][nn]][1]]

        if FSstartCrossing <= i <= FSendCrossing:

            mindist = None
            nn = 0

            for j in range(9):
                
                dist = hs(relevantFSpoints_coords[FS_9nn[1][i][j]],CloudSatPoints[i])
                if mindist == None or dist<mindist:
                    mindist = dist
                    nn = j
            
            InterpFS_ZFactorMeasured[i,:,:] = FS_ZFactorMeasured[relevantFSpoints_indeces[FS_9nn[1][i][nn]][0],relevantFSpoints_indeces[FS_9nn[1][i][nn]][1]]

        if HSstartCrossing <= i <= HSendCrossing:

            mindist = None
            nn = 0

            for j in range(9):
                
                dist = hs(relevantHSpoints_coords[HS_9nn[1][i][j]],CloudSatPoints[i])
                if mindist == None or dist<mindist:
                    mindist = dist
                    nn = j
            
            InterpHS_ZFactorMeasured[i,:] = HS_ZFactorMeasured[relevantHSpoints_indeces[HS_9nn[1][i][nn]][0],relevantHSpoints_indeces[HS_9nn[1][i][nn]][1]]


    # for plotting:
    
    InterpS1 = np.array(InterpS1)
    InterpS2 = np.array(InterpS2)
    InterpSwath = np.empty((numPoints,13))
    InterpSwath[:,0:4] = np.flip(InterpS2,1)
    InterpSwath[:,4:13] = np.flip(InterpS1,1)
    InterpSwath = np.flip(InterpSwath, 1)   

    InterpFS_ZFactorMeasured = InterpFS_ZFactorMeasured/100
    InterpHS_ZFactorMeasured = InterpHS_ZFactorMeasured/100
    InterpFS_ZFactorMeasured = np.ma.masked_equal(InterpFS_ZFactorMeasured, -9999.9)
    InterpHS_ZFactorMeasured = np.ma.masked_equal(InterpHS_ZFactorMeasured, -9999.9)



    #InterpS2[InterpS2 == np.NaN] = -99

    
    
     # Draw the ALONG-TRACK VERTICAL PROFILE of GMI Swath
    fig = plt.figure(figsize=(4,1))
    ax = plt.axes()
    ax.set_xlabel('along-track axis index')
    ax.set_ylabel('Channel')
    ax.set_title('Vertical Cross-Section of GMI Scan (Along-Track)')
    
    locs, labels = plt.yticks()  # Get the current locations and labels.
    plt.yticks(np.arange(0, 13, step=1))  # Set label locations.
    plt.yticks(np.arange(13), ['10V', '10H', '18V', '18H', '23V', '36V', '36H', '89V', '89H', '166V','166H', '183+-3','183+-7'])  # Set text labels.
    plt.xticks(np.arange(0,numPoints, step=numPoints/5))
    ax = ax.invert_yaxis()
    plt.autoscale(False)
    pp = plt.imshow(np.transpose(InterpSwath), cmap='jet', aspect = 10.0, interpolation='None')
    # Add a colorbar to the bottom of the plot.
    fig.subplots_adjust(bottom=0.15, left=0.06, right=0.94)
    cbar_ax = fig.add_axes([0.12, 0.11, 0.76, 0.025])  
    cbar = plt.colorbar(pp, cax=cbar_ax, orientation='horizontal')
    
    #cbar.set_label(label=pr.attrs.get('units').decode('utf-8'),size=10)
    plt.show()

    
    CloudSat_Scantimes = TimeArray[S1startCrossing:S1endCrossing]
    CloudSat_Lat = crossingLat
    CloudSat_Lon = crossingLon
    CloudSat_RadarRefl = coinRadarReflectivity[S1startCrossing:S1endCrossing]
    CloudSat_Height = coinHeight[S1startCrossing:S1endCrossing]
    CloudSat_Indexes = []
    for i in range(numPoints):
        CloudSat_Indexes.append(i)

    CloudSat_Height = CloudSat_Height.tolist()
   
    return InterpSwath, InterpS2, S2startCrossing, S2endCrossing, InterpHS_ZFactorMeasured, InterpFS_ZFactorMeasured

def readGMI(GMI_files, CoinTime, secondFile):
    '''
    Reads relevant GMI data from specified file and returns Lat, Lon, Tb, and ScanTime arrays from +- 30 minutes from coincidence time
    '''
        
    f = h5py.File(GMI_files[0],'r')

    S2 = list(f.keys())[1]

    Year = CoinTime.year
    Month = CoinTime.month
    Day = (GMI_files[0])[28:30]
    #print(Day)
    Day = int(Day)
    OriginalDay = Day
    NextDay = Day+1

    HourArray = np.array(f[S2].get('ScanTime').get('Hour'))
    MinuteArray = np.array(f[S2].get('ScanTime').get('Minute'))
    SecondArray = np.array(f[S2].get('ScanTime').get('Second'))

    TimeArray = np.empty(0)

    minus10 = CoinTime - timedelta(minutes=7.5)
    plus10 = CoinTime + timedelta(minutes=7.5)
    startIndex = None
    endIndex = None

    #print(minus10)

    if secondFile == True:
        startIndex = 0

    # find section of DPR data that falls within +- 7.5 minutes of coincidence
    for i in range(HourArray.size):
        if HourArray[i] == 0 and MinuteArray[i] == 0 and Day == OriginalDay and MinuteArray[0] != 0:
            Day = NextDay
        scantime = datetime(Year, Month, Day, HourArray[i], MinuteArray[i], SecondArray[i])
        
        if minus10 < scantime < plus10:
            TimeArray = np.append(TimeArray, scantime)
            if startIndex == None:
                startIndex = i
        elif scantime > plus10:
            endIndex = i
            #print('found end')
            break


    if startIndex == None:
        print("Error in reading GMI file: start time of overlap not within data.")
        print(scantime)
        return


    #print(startIndex)
    #print(endIndex)

    S2TbArray = np.array(f[S2].get('Tc'))

    S2coinTbArray = S2TbArray[startIndex:endIndex,:,:]

    S2LatArray = np.array(f[S2].get('Latitude'))

    


    S2coinLatArray = S2LatArray[startIndex:endIndex,:]
    S2LonArray = np.array(f[S2].get('Longitude'))

    S2coinLonArray = S2LonArray[startIndex:endIndex,:]



    S1 = list(f.keys())[0]

    S1TbArray = np.array(f[S1].get('Tc')) #changed to 1CGMI instead of 1B, so Tc instead of Tb

    S1coinTbArray = S1TbArray[startIndex:endIndex,:,:]

    S1LatArray = np.array(f[S1].get('Latitude'))

    #print(np.shape(S1LatArray))


    S1coinLatArray = S1LatArray[startIndex:endIndex,:]

    S1LonArray = np.array(f[S1].get('Longitude'))

    S1coinLonArray = S1LonArray[startIndex:endIndex,:]

    S1TbArray[S1TbArray == -9999.9] = np.NAN

    S1coinTbArray[S1coinTbArray == -9999.9] = np.NAN

    S2coinTbArray[S2coinTbArray == -9999.9] = np.NAN

    #print('GMI location at coinTime:')
    #print('' + str(S1LatArray[coin_i][100]) + ', ' + str(S1LonArray[coin_i][100]))
    #print(np.shape(S1coinLonArray))

    if len(GMI_files) > 1:
        nextGMIFile = [GMI_files[1]]
        [next_S1coinLatArray, next_S1coinLonArray, next_S1coinTbArray, next_TimeArray, next_S1LatArray, next_S1LonArray, next_S1TbArray, next_S2coinLatArray, next_S2coinLonArray, next_S2coinTbArray, nextStartIndex, nextEndIndex] = readGMI(nextGMIFile, CoinTime, True)
        
        S1coinLatArray = np.append(S1coinLatArray, next_S1coinLatArray, axis=0)
        S1coinLonArray = np.append(S1coinLonArray, next_S1coinLonArray, axis=0)
        S1coinTbArray = np.append(S1coinTbArray, next_S1coinTbArray, axis=0)
        TimeArray = np.append(TimeArray, next_TimeArray, axis=0)
        S1LatArray = np.append(S1LatArray, next_S1LatArray, axis=0)
        S1LonArray = np.append(S1LonArray, next_S1LonArray, axis=0)
        S1TbArray = np.append(S1TbArray, next_S1TbArray, axis=0)
        S2coinLatArray = np.append(S2coinLatArray, next_S2coinLatArray, axis=0)
        S2coinLonArray = np.append(S2coinLonArray, next_S2coinLonArray, axis=0)
        S2coinTbArray = np.append(S2coinTbArray, next_S2coinTbArray, axis=0)
        endIndex = nextEndIndex

    #print('shape after append: ')
    #print(np.shape(S1coinLonArray))

    if endIndex == None:
        print('Error: second GMI file needed but not provided')
        print(scantime)
        return
    
    #print(startIndex)
    #print(endIndex)
    return [S1coinLatArray, S1coinLonArray, S1coinTbArray, TimeArray, S1LatArray, S1LonArray, S1TbArray, S2coinLatArray, S2coinLonArray, S2coinTbArray, startIndex, endIndex]

def readDPR(DPRfiles, CoinTime, secondFile):
    '''
    Reads relevant DPR data from specified file and returns Lat, Lon, Tb, and ScanTime arrays from +- 30 minutes from coincidence time
    '''
    DPRfile = DPRfiles[0]
    f = h5py.File(DPRfile,'r')


    #print(len(list(f.keys())))
    FS = list(f.keys())[1]


    Year = CoinTime.year
    Month = CoinTime.month
    Day = (DPRfiles[0])[29:31]
    #print(Day)
    Day = int(Day)
    OriginalDay = Day
    NextDay = Day+1

    HourArray = np.array(f[FS].get('ScanTime').get('Hour'))
    MinuteArray = np.array(f[FS].get('ScanTime').get('Minute'))
    SecondArray = np.array(f[FS].get('ScanTime').get('Second'))

    TimeArray = np.empty(0)

    minus10 = CoinTime - timedelta(minutes=7.5)
    plus10 = CoinTime + timedelta(minutes=7.5)
    startIndex = None
    endIndex = None


    #print(minus10)

    if secondFile == True:
        startIndex = 0

    # find section of DPR data that falls within +- 7.5 minutes of coincidence
    for i in range(HourArray.size):

        if HourArray[i] == 0 and MinuteArray[i] == 0 and Day == OriginalDay and MinuteArray[0] != 0:
            Day = NextDay
        scantime = datetime(Year, Month, Day, HourArray[i], MinuteArray[i], SecondArray[i])
        
        if minus10 < scantime < plus10:
            TimeArray = np.append(TimeArray, scantime)
            if startIndex == None:
                startIndex = i
        elif scantime > plus10:
            endIndex = i
            break

    #print(startIndex)
    #print(endIndex)

    if startIndex == None:
        print("Error in reading DPR file: start time of overlap not within data.")
        return

    FS_ZFactorMeasured = np.array(f[FS].get('PRE').get('zFactorMeasured'))
    FS_Lat = np.array(f[FS].get('Latitude'))
    FS_Lon = np.array(f[FS].get('Longitude'))

    FS_ZFactorMeasured = FS_ZFactorMeasured[startIndex:endIndex,:,:]
    FS_Lat = FS_Lat[startIndex:endIndex,:]
    FS_Lon = FS_Lon[startIndex:endIndex,:]


    # other variables retrieved by old CSATGPM algorithm
    FS_localZenithAngle = (np.array(f[FS].get('PRE').get('localZenithAngle')))[startIndex:endIndex]
    FS_elevation = (np.array(f[FS].get('PRE').get('elevation')))[startIndex:endIndex]
    FS_heightZeroDeg = (np.array(f[FS].get('VER').get('heightZeroDeg')))[startIndex:endIndex]

    # extra useful variables
    FS_precipRateESurface = (np.array(f[FS].get('SLV').get('precipRateESurface')))[startIndex:endIndex]
    FS_phaseNearSurface = (np.array(f[FS].get('SLV').get('phaseNearSurface')))[startIndex:endIndex]



    #FSgroups = list(f[FS])
    #print(FSgroups)

    HS = list(f.keys())[2]
    HSgroups = list(f[HS])
    #print(HSgroups)

    HS_ZFactorMeasured = np.array(f[HS].get('PRE').get('zFactorMeasured'))
    HS_Lat = np.array(f[HS].get('Latitude'))
    HS_Lon = np.array(f[HS].get('Longitude'))

    HS_ZFactorMeasured = HS_ZFactorMeasured[startIndex:endIndex,:,:]
    HS_Lat = HS_Lat[startIndex:endIndex,:]
    HS_Lon = HS_Lon[startIndex:endIndex,:]

    if len(DPRfiles) > 1:
        nextDPRFile = [DPRfiles[1]]
        [nextFS_ZFactorMeasured, nextFS_Lat, nextFS_Lon, nextHS_ZFactorMeasured, nextHS_Lat, nextHS_Lon, nextEndIndex] = readDPR(nextDPRFile, CoinTime, True)
        FS_Lat = np.append(FS_Lat, nextFS_Lat, axis=0)
        FS_Lon = np.append(FS_Lon, nextFS_Lon, axis=0)
        FS_ZFactorMeasured = np.append(FS_ZFactorMeasured, nextFS_ZFactorMeasured, axis=0)
        HS_ZFactorMeasured = np.append(HS_ZFactorMeasured, nextHS_ZFactorMeasured, axis=0)
        HS_Lat = np.append(HS_Lat, nextHS_Lat, axis=0)
        HS_Lon = np.append(HS_Lon, nextHS_Lon, axis=0)
        
        endIndex = nextEndIndex

            
    if endIndex == None:
        print('Error: second GMI file needed but not provided')
        return


    return [FS_ZFactorMeasured, FS_Lat, FS_Lon, HS_ZFactorMeasured, HS_Lat, HS_Lon, endIndex]

def readCS(files, coinTime, secondFile):
# Open HDF4 file.
    FILE_NAME = files[0]
    
    hdf = SD(FILE_NAME, SDC.READ)

    #print(hdf.datasets())

    # Read datasets.
    DATAFIELD_NAME = 'Radar_Reflectivity'
    dset = hdf.select(DATAFIELD_NAME)
    RadarReflectivity = dset[:,:]

    #print(np.shape(RadarReflectivity))

    ht = hdf.select('Height')
    height = ht[:,:]

    #print(np.shape(height))

    h = HDF.HDF(FILE_NAME)

    vs = h.vstart()


    xid = vs.find('Latitude')

    latid = vs.attach(xid)
    latid.setfields('Latitude')
    nrecs, _, _, _, _ = latid.inquire()
    latitude = latid.read(nRec=nrecs)
    latid.detach()

    lonid = vs.attach(vs.find('Longitude'))
    lonid.setfields('Longitude')
    nrecs, _, _, _, _ = lonid.inquire()
    longitude = lonid.read(nRec=nrecs)
    lonid.detach()


    timeid = vs.attach(vs.find('Profile_time'))
    timeid.setfields('Profile_time')
    nrecs, _, _, _, _ = timeid.inquire()
    time = timeid.read(nRec=nrecs)
    #units_t =  timeid.attr('units').get()
    #longname_t = timeid.attr('long_name').get()
    timeid.detach()

    filename = (re.search(r'' + '\d+',FILE_NAME)).group()
    #print(filename)

    ##these variables are dependent on what folder file is in - need to change to be more general
    year = int(filename[0:4])
    day = int(filename[4:7])
    hour = int(filename[7:9])
    minute = int(filename[9:11])
    second = int(filename[11:13])

    TimeArray = np.empty(0)
    startTime = datetime(year,1,1) + timedelta(days = day-1, seconds=second, hours=hour,minutes=minute)
    minus5 = coinTime - timedelta(minutes=5)
    plus5 = coinTime + timedelta(minutes=5)
    startIndex = None
    endIndex = None
    
    #print(startTime)

    #print(startTime)
    #print(time[0][0])
    #print(print(time[-1][0]))
    
    # find section of CloudSat data that falls within +-5 minutes of coincidence
    for i in range(np.shape(time)[0]):
        scantime = startTime + timedelta(seconds=int((time[i][0])))
        
        if minus5 < scantime < plus5:
            TimeArray = np.append(TimeArray, scantime)

            if startIndex == None:
                startIndex = i
                
        elif scantime > plus5:
            
            endIndex = i - 1
            #print("passed plus10: ")
            #print("scantime: " + str(scantime))
            #print("plus10: " + str(plus10))
            break

    #print("startindex: " + str(startIndex))
    #print("endindex: " + str(endIndex))
    if endIndex == None:
        endIndex = -1
        if len(files) == 1:
            print('error: second CloudSat file needed but not provided')

    coinLat = latitude[startIndex:endIndex]
    coinLon = longitude[startIndex:endIndex]
    coinTime = TimeArray
    coinRadarReflectivity = RadarReflectivity[startIndex:endIndex]
    coinHeight = height[startIndex:endIndex]


    data = [coinLat, coinLon, coinTime, coinRadarReflectivity, coinHeight, latitude, longitude]

    if len(files) > 1 and endIndex == -1:

        nextFileData = readCS(files[1], coinTime, True)

        for i in range(len(nextFileData)):
            data[i] = np.append(data[i], nextFileData[i], axis=0)

    return data


def findCrossing(CloudSat_lat, CloudSat_lon, swathLat, swathLon):
    
    #create shapely polygon containing section of GMI swath near overlap area using centerpoints along edge of swath
    swathBoundary = []
    for i in range(np.shape(swathLat)[0]): #1
        swathBoundary.append([swathLon[i][0], swathLat[i][0]])
    for i in range(np.shape(swathLat)[1]): #2
        swathBoundary.append([swathLon[0][i], swathLat[0][i]])
    for i in range(np.shape(swathLat)[0]-1, -1, -1): #3
        swathBoundary.append([swathLon[i][-1], swathLat[i][-1]])
    for i in range(np.shape(swathLat)[1]-1, -1, -1): #4
        swathBoundary.append([swathLon[-1][i], swathLat[-1][i]])
    swath = shp.Polygon(swathBoundary)
    
    startCrossing = None
    # check to find where CloudSat overlap with swath begins
    for i in range(len(CloudSat_lat)):
        point = shp.Point(CloudSat_lon[i],CloudSat_lat[i])
        if point.within(swath):
            startCrossing = i
            break
    
    if startCrossing == None:
        print('error: crossing not found')
        endCrossing = None
        return startCrossing, endCrossing, swath
    
    # check to find where CloudSat overlap with swath ends
    for i in range(len(CloudSat_lat)-1, 0, -1):
        point = shp.Point(CloudSat_lon[i],CloudSat_lat[i])
        if point.within(swath):
            endCrossing = i
            break
    
    # Swath polygon used to find cloudsat points within GMI swath is created using the centerpoints of the grids along the edge of the swath.
    # This misses the cloudSat points that are within the GMI swath but on the outer half of the grids on the edge. Corrected by adding cloudSat
    # points within 0.01 degrees of the edge of the polygon.
    while True:
        point = shp.Point(CloudSat_lon[startCrossing-1],CloudSat_lat[startCrossing-1])
        if point.distance(swath) < 0.01: #fix so that it actually uses km - need to find thickness of edge grid by measuring dist to next inside pt and dividing by 2?
            startCrossing -= 1
        else:
            break
    
    while True:
        point = shp.Point(CloudSat_lon[endCrossing+1],CloudSat_lat[endCrossing+1])
        if point.distance(swath) < 0.01: # fix so that it actually uses km - need to find thickness of edge grid
            endCrossing += 1
        else:
            break
    
    return startCrossing, endCrossing, swath

def findFiles(GPMTime, CSTime, GMI_PATH, DPR_PATH, CS_PATH):
    '''
    return files for each data product at the coincidence time
    '''

    GMIfile = []
    DPRfile = []
    CSfile = []


    GPMstart = GPMTime - timedelta(minutes=15)
    GPMend = GPMTime + timedelta(minutes=15)
    GPMstartdate = GPMstart.strftime('20' + '%y%m%d')
    GPMstarthour = GPMstart.strftime('%H%M%S')
    GPMendhour = GPMend.strftime('%H%M%S')

    secondFile = False
    secondGMIFile = False

    for filename in os.scandir(GMI_PATH):

        if filename.is_file():

            if secondGMIFile == True:
                GMIfile.append(filename.name)
                secondGMIFile = False
                break
            
            startdate = (filename.name)[22:30]
            
            if startdate == GPMstartdate:

                starthour = (filename.name)[32:38]
                endhour = (filename.name)[40:46]
                
                if (int(starthour) < int(GPMstarthour)) & (int(endhour) > int(GPMendhour)):
                    GMIfile.append(filename.name)
                    break
                elif (int(starthour) < int(GPMstarthour)) & (int(endhour) < int(GPMendhour)):
                    secondGMIFile = True

    GMIorbitnum = (GMIfile[0][47:53])

    secondDPRfile = False

    for filename in os.scandir(DPR_PATH):

        if filename.is_file():

            if secondDPRfile == True:
                DPRfile.append(filename.name)
                secondDPRfile = False
                break

            orbitnum = ((filename.name)[48:54])

            if orbitnum == GMIorbitnum:
                DPRfile.append(filename.name)
                if secondFile == True:
                    secondDPRfile = True


    CSstart = CSTime - timedelta(minutes=15)
    CSend = CSTime + timedelta(minutes=15)
    CSstartdate = CSstart.strftime('20' + '%y%j')
    CSstarthour = CSstart.strftime('%H%M%S')
    CSendhour = CSend.strftime('%H%M%S')

    secondFile = False
    secondCSFile = False
    
    for filename in os.scandir(CS_PATH):
        if filename.is_file():

            if secondCSFile == True:
                CSfile.append(filename.name)
                secondCSFile = False
                break
            
            startdate = (filename.name)[0:7]
            
            
            if startdate == CSstartdate:

                starthour = (filename.name)[7:13]
                starttime = datetime.strptime((startdate+starthour)[2:-1],'%y%j%H%M%S')
                endtime = starttime + timedelta(minutes=98)
                
                if (int(starthour) < int(CSstarthour)) & (endtime > CSend):
                    CSfile.append(filename.name)
                    break
                elif (int(starthour) < int(CSstarthour)) & (endtime < CSend):
                    secondCSFile = True

    
    return [GMIfile, DPRfile, CSfile]

def visualize(swath, LonArray, LatArray, Tb, cslongitude, cslatitude, crossingLon, crossingLat, channel_name, coinLonArray, coinLatArray, coinTb, coinLonArrayS2, coinLatArrayS2, coinTbS2, coinLon, coinLat, CStime):
    print("plotting")
    
    f1 = plt.figure()
    projection = ccrs.PlateCarree(-20) # offset center of plot by 0.1 degree - for some reason, centering at 0 breaks the map
    ax = plt.axes(projection=projection)
    ax.coastlines(resolution='10m',linewidth=0.5)
    gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,
                    linewidth=0.8,color='#555555',alpha=0.5,linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    
    # plot GMI swath
    x1 = LonArray
    y1 = LatArray
    z = Tb
    zmask = ma.masked_where(np.isnan(z),z)
    mesh = ax.pcolormesh(x1,y1,zmask, transform=ccrs.PlateCarree(), cmap='jet')
    f1.colorbar(mesh,orientation='horizontal', label='Tb [K]', shrink=.8)

    # plot entire CS orbit, as well as overlap region in red
    x3 = cslongitude
    y3 = cslatitude
    plt.scatter(x3,y3, transform=ccrs.PlateCarree(), marker='.')
    plt.scatter(crossingLon,crossingLat, transform=ccrs.PlateCarree(), marker='.', color='red')

    # plot shapely swath polygon - only for demonstration purposes
    #swathx,swathy = swath.exterior.xy
    #plt.plot(swathx,swathy, transform=ccrs.PlateCarree())
   
    plt.title('GMI ' + channel_name + ' and CloudSat path, coincidence at ' + str(CStime))
    #plt.show()


    filename = '04262015_' + channel_name.replace(' ','') + '.png'
    plt.show

   

    # zoomed in plot of overlap region

    f2 = plt.figure()
    ax2 = plt.axes(projection=ccrs.PlateCarree())
    ax2.coastlines(resolution='10m',linewidth=0.5)
    gl2 = ax2.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,
                    linewidth=0.8,color='#555555',alpha=0.5,linestyle='--')
    
    # plot rough coincidence region of GMI swath
    x21 = coinLonArray
    y21 = coinLatArray
    z2 = coinTb
    z2mask = ma.masked_where(np.isnan(z2),z2)
    mesh = ax2.pcolormesh(x21,y21,z2mask,transform=ccrs.PlateCarree(), cmap='jet')
    f2.colorbar(mesh, orientation='horizontal', label='Tb [K]', shrink=.8)
   
    # plot rough coincidence region of CloudSat orbit
    x22 = coinLon
    y22 = coinLat
    ax.autoscale(False) # To avoid scatter changing limits
    ax = plt.scatter(x22,y22,transform=ccrs.PlateCarree(), marker=".")

    # plot exact overlap points of CloudSat in red
    plt.scatter(crossingLon,crossingLat, transform=ccrs.PlateCarree(), marker='.', color='red')

    # plot shapely GMI swath polygon - demonstration purposes only
    #plt.plot(swathx,swathy, transform=ccrs.PlateCarree())
    
    plt.title('GMI S1' + channel_name + ' and  CloudSat path')
    #plt.show()

    filename = '04262015_' + channel_name.replace(' ','') + 'zoomed.png'
    plt.savefig(filename)


    