import os
import sys
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use('ggplot')
from datetime import datetime

AU = 149597870700 # meters

pd.set_option('display.max_rows', 200)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 120)

rootDir          = '/home/rharvey/Documents/GS-Analysis/'
dataDir          = '/home/rharvey/data/'
inputDir         = 'SearchResult/no_overlap/'
outputDir        = 'SearchResult/selected_events/'
outputPlotDir    = 'GS_statistical_plot/'
inputListDir     = 'SearchResult/detailed_info/'
# inputShockDir    = 'GS_DataCsv/'
# inputShockFile   = 'shocks_20170303_035000.csv'
# inputHcsDir      = 'GS_DataCsv/'
# inputHcsFile     = 'Leif_Svalgaard_IMF_Sector_Boundaries_1926_2017.csv'
# inputExhaustDir  = 'GS_DataCsv/'
# inputExhaustFile = 'ACE_ExhaustList.csv'
###############################################################################

# Read Shock Data.
IPShock_DF = pd.read_csv(rootDir+inputShockDir+inputShockFile, header=0)
shockTime = IPShock_DF.apply(lambda row: datetime(row['Year'], row['Month'], row['Day'], row['Hour'], row['Minute'], row['Second']), axis=1)
IPShock_DF.insert(0, 'shockTime', shockTime) # This command is able to specify column position.
IPShock_DF.drop(['Year', 'Month', 'Day', 'Hour', 'Minute', 'Second'],inplace=True,axis=1)

# Pick only Wind events, and restrict the time range from 1996 to 2016.
IPShock_WIND_DF = IPShock_DF[IPShock_DF['Spacecraft'].str.contains('Wind')]
IPShock_WIND_1996_2016_DF = IPShock_WIND_DF[(IPShock_WIND_DF['shockTime']>=datetime(1996,1,1,0,0,0))&(IPShock_WIND_DF['shockTime']<=datetime(2016,12,31,23,59,59))]
IPShock_WIND_1996_2016_DF = IPShock_WIND_1996_2016_DF.sort_values(by='shockTime').copy()
IPShock_WIND_1996_2016_DF.reset_index(drop=True, inplace=True)
# print(IPShock_WIND_1996_2016_DF[['shockTime', 'Spacecraft']])

# Read HCS Data.
HCS_DF = pd.read_csv(rootDir+inputHcsDir+inputHcsFile, header=0)
hcsTime = HCS_DF.apply(lambda row: datetime(row['Year'], row['Month'], row['Day'], 12, 0, 0), axis=1)
HCS_DF.insert(0, 'hcsTime', hcsTime) # This command is able to specify column position.
HCS_DF.drop(['Year', 'Month', 'Day'],inplace=True,axis=1)

# Restrict the time range from 1996 to 2016.
HCS_1998_2017_DF = HCS_DF[(HCS_DF['hcsTime']>=datetime(1998,1,1,0,0,0))&(HCS_DF['hcsTime']<=datetime(2017,12,31,23,59,59))]
HCS_1998_2017_DF = HCS_1998_2017_DF.sort_values(by='hcsTime').copy()
HCS_1998_2017_DF.reset_index(drop=True, inplace=True)
# print(HCS_1996_2016_DF)


# Read HCS Data.
month2Num = {'Jan':1,'Feb':2,'Mar':3,'Apr':4,'May':5,'Jun':6,'Jul':7,'Aug':8,'Sep':9,'Oct':10,'Nov':11,'Dec':12}
Exhaust_DF = pd.read_csv(rootDir+inputExhaustDir+inputExhaustFile, header=0)
# print(month2Num[Exhaust_DF.loc[1,'Month']])
exhaustTime = Exhaust_DF.apply(lambda row: datetime(int(row['Year']), int(month2Num[row['Month']]), int(row['Day']), int(row['Hour']), int(row['Minute']), 0), axis=1)
Exhaust_DF.insert(0, 'exhaustTime', exhaustTime)
Exhaust_DF.drop(['Year', 'Month', 'Day', 'Hour', 'Minute'],inplace=True,axis=1)

# Restrict the time range from 1996 to 2016. Actually, the range is only from 1997 to 2006.
Exhaust_1996_2016_DF = Exhaust_DF[(Exhaust_DF['exhaustTime']>=datetime(1996,1,1,0,0,0))&(Exhaust_DF['exhaustTime']<=datetime(2016,12,31,23,59,59))]
Exhaust_1996_2016_DF = Exhaust_1996_2016_DF.sort_values(by='exhaustTime').copy()
Exhaust_1996_2016_DF.reset_index(drop=True, inplace=True)
# print(Exhaust_1996_2016_DF)

isVerbose = False
year = sys.argv[1]

B_threshold = ('(|B|>0)', '(|B|<5)', '(|B|>=5)')
isPlotDaysToHCS = 0
isPlotDaysToExhaust = 0
isPlot = isPlotDaysToHCS or isPlotDaysToExhaust

# Plot layout: 7 rows, 3 columns
if isPlot:
    max_row = 7
    max_column = 3

for B_lim in B_threshold:
    n_total = 0
    
    if isPlot:
        fig1,ax1 = plt.subplots(max_row, max_column, figsize=(8, 12))
        row = 0
        column = 0

    # Generate dataset with all |B| value.
    inputFileName = namestr[:-1] + probe_str +'_detailed_info.p'
    FR_selected_events_DF = pd.read_pickle(open(rootDir + inputListDir + inputFileName,'rb'))

    # Apply criteria
    if B_lim == '(|B|>0)':
        FR_selected_events_DF = FR_selected_events_DF[FR_selected_events_DF['B_abs_mean']>=0]
    elif B_lim == '(|B|<5)':
        FR_selected_events_DF = FR_selected_events_DF[FR_selected_events_DF['B_abs_mean']<5]
    elif B_lim == '(|B|>=5)':
        FR_selected_events_DF = FR_selected_events_DF[FR_selected_events_DF['B_abs_mean']>=5]
    else:
        print('Error in applying criteria!')
        exit()
    FR_selected_events_DF.reset_index(drop=True, inplace=True)

    # Get Unit vectors for flux rope frame.
    X_unitVector = np.array(FR_selected_events_DF[['X_unitVector[0]', 'X_unitVector[1]', 'X_unitVector[2]']])
    Y_unitVector = np.array(FR_selected_events_DF[['Y_unitVector[0]', 'Y_unitVector[1]', 'Y_unitVector[2]']])
    Z_unitVector = np.array(FR_selected_events_DF[['Z_unitVector[0]', 'Z_unitVector[1]', 'Z_unitVector[2]']])
    
    # Construct transformation matrices.
    Matrix_transToFRframe = np.zeros((len(FR_selected_events_DF),3,3))
    Matrix_transToFRframe[:,:,0] = X_unitVector
    Matrix_transToFRframe[:,:,1] = Y_unitVector
    Matrix_transToFRframe[:,:,2] = Z_unitVector
    
    # Get VHT and project into FR Frame.
    VHT = np.array(FR_selected_events_DF[['VHT_inGSE[0]', 'VHT_inGSE[1]', 'VHT_inGSE[2]']])
    VHT_inFR = np.zeros((len(FR_selected_events_DF),3))
    for i in range(len(VHT)):
        VHT_inFR[i,:] = VHT[i,:].dot(Matrix_transToFRframe[i,:,:])
    
    # Get time duration
    duration = np.array(FR_selected_events_DF['duration']) # Convert minutes to seconds.
    size = - VHT_inFR[:,0] * 1000.0 * duration             # Space increment along X axis. Convert km/s to m/s.
    size_inAU_array = size/AU                              # Divided by AU.
    
    # Add daysToExhaust.
    FR_selected_events_DF.insert(7, 'daysToExhaust', [np.nan]*len(FR_selected_events_DF))
    
    # Extend one month time boundary.
    if (year == year_start)&(year != year_end):
        Exhaust_currentYear_temp_DF = Exhaust_1996_2016_DF[(Exhaust_1996_2016_DF['exhaustTime']>=datetime(year,1,1,0,0,0))&(Exhaust_1996_2016_DF['exhaustTime']<=datetime(year+1,1,31,23,59,59))]
    elif (year != year_start)&(year == year_end):
        Exhaust_currentYear_temp_DF = Exhaust_1996_2016_DF[(Exhaust_1996_2016_DF['exhaustTime']>=datetime(year-1,12,1,0,0,0))&(Exhaust_1996_2016_DF['exhaustTime']<=datetime(year,12,31,23,59,59))]
    elif (year == year_start)&(year == year_end):
        Exhaust_currentYear_temp_DF = Exhaust_1996_2016_DF[(Exhaust_1996_2016_DF['exhaustTime']>=datetime(year,1,1,0,0,0))&(Exhaust_1996_2016_DF['exhaustTime']<=datetime(year,12,31,23,59,59))]
    elif (year != year_start)&(year != year_end):
        Exhaust_currentYear_temp_DF = Exhaust_1996_2016_DF[(Exhaust_1996_2016_DF['exhaustTime']>=datetime(year-1,12,1,0,0,0))&(Exhaust_1996_2016_DF['exhaustTime']<=datetime(year+1,1,31,23,59,59))]
    else:
        print('Something wrong, please check.')
    
    # If Exhaust_currentYear_temp_DF is not empty, Loop flux rope events list, calculate tempral distance of each event to Exhaust.
    if not Exhaust_currentYear_temp_DF.empty:
        for event_index in FR_selected_events_DF.index:
            if isVerbose:
                print('event_index = {}'.format(event_index))
            oneEventRecord = FR_selected_events_DF.loc[event_index] # Note the difference of iloc and loc.
            FluxRopeTurnTime = oneEventRecord['turnTime']
            Exhaust_time_list = Exhaust_currentYear_temp_DF['exhaustTime']
            
            # Time difference in days(rounded). Negative value indicate the flux rope preceded the exhaust.
            timeDiff_list = [round(x.total_seconds()/(24*3600), 0) for x in (FluxRopeTurnTime - Exhaust_time_list)]
            daysToExhaust = min(timeDiff_list, key=lambda x:abs(x))
            FR_selected_events_DF.set_value(event_index, 'daysToExhaust', daysToExhaust)
            # print(FR_selected_events_DF.loc[event_index])
    
    if isPlotDaysToExhaust:
        ax = FR_selected_events_DF['daysToExhaust'].plot.hist(bins=np.arange(-30.5,30.5), ylim=[0,400], color='#0077C8')
        ax.set_xlim([-30,30])
        ax.set_ylim([0,350]) # For B>=5
        fig = ax.get_figure()
        plt.title('Time to Reconnection Exhaust (Year '+str(year)+')')
        plt.xlabel('Day to Reconnection Exhaust')
        plt.ylabel('Flux Rope Occurrence Counts')
        fig.savefig(rootDir + '/plot_temp/TimeToExhaust_'+str(year)+'_'+B_lim+'.eps', format='eps')
        plt.close('all')
        
        # Add shock information.
        FR_selected_events_DF.insert(8, 'afterShock', [False]*len(FR_selected_events_DF))
        FR_selected_events_DF.insert(9, 'shockTheta_deg', [np.nan]*len(FR_selected_events_DF))

        # Loop shock list. Add info to FR_selected_events_DF.
        IPShock_WIND_currentYear_temp_DF = IPShock_WIND_1996_2016_DF[(IPShock_WIND_1996_2016_DF['shockTime']>=datetime(year,1,1,0,0,0))&(IPShock_WIND_1996_2016_DF['shockTime']<=datetime(year,12,31,23,59,59))]
        # IPShock_WIND_currentYear_temp_DF.reset_index(inplace=True, drop=True)
        
        # FR_selected_events_DF = FR_selected_events_DF[(FR_selected_events_DF['walenTest_slope']>=0.300) | (FR_selected_events_DF['walenTest_slope']<=-0.300)]
        # FR_selected_events_DF = FR_selected_events_DF[(FR_selected_events_DF['walenTest_r_value']>=0.800) | (FR_selected_events_DF['walenTest_r_value']<=-0.800)]
        
        # If outputDir does not exist, create it.
        if not os.path.exists(rootDir + outputDir):
            os.makedirs(rootDir + outputDir)
        
        # Save all events to pickle file.
        FR_selected_events_DF.to_pickle(rootDir + outputDir + namestr[:-1] + probe_str + B_lim + '_selected_events.p')
       
        # Pick and save the events after shock.
        FR_selected_events_DF = FR_selected_events_DF[FR_selected_events_DF['afterShock']==True]
        FR_selected_events_DF.reset_index(inplace=True, drop=True)
        FR_selected_events_DF.to_pickle(rootDir + outputDir + namestr[:-1] + probe_str + B_lim + '_selected_events_afterShock.p')
    
        # Save Plot.
        if isPlot:
            if not os.path.exists(rootDir + outputPlotDir):
                os.makedirs(rootDir + outputPlotDir)
            fig1.savefig(rootDir + outputPlotDir + 'TimeToHCS_' + B_lim + '.eps', format='eps')
            plt.close('all')

