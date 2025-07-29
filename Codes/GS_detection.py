#!/usr/local/bin/python

'''
2016-11-16
Use brute force method to detect flux rope. Window center sliding method.
1) v1.4 is able to handle multiple structure.
2) comment out print statement to accelerate.
3) improve the method to find location of turn point.
4) put true and false result into two subfolder, and change cycle counter to 10000.
2016-11-30
version 1.7
1) change turnPointOnTop = (Pt_turnPoint>(Pt_max-(Pt_max-Pt_min)*0.2)) to turnPointOnTop = (Pt_turnPoint>(Pt_max-(Pt_max-Pt_min)*0.15))

2017-01-03
version 1.8
1) use flexible smooth window when get A curve before testing turn point counts. savgol_filter_window = half_window, no less than 9.

2017-02-01
version 1.9
1) Implement VHT.

2017-02-03
version 2.0
1) Change the interval specifying method. half_window_list = [5, 10, 15, 20, 25, 30]

2017-02-07
version 2.1
1) Combine three version into one.(For macbook, blueshark, and bladerunner)
2) Remove the command to get Np_DataFrame. Np is not used in this code.
3) Create and destroy pool in inner for loop. Since we want to save result in file when we finish one inner loop(complete one window size), we need to block the processing until all parallel subproceure finish with in the inner loop, and then save file. After saving file, create a new pool and start a new inner loop.
4) Add searching time range printing in searchFluxRopeInWindow() function to indicate procedure.

2017-02-07
version 2.2
1) Add one more command line argument to specify the time range list. Command format is:
    python GS_detectFluxRope_multiprocessing_v2.2.py 1996 '((20,30),(30,40),(40,50))'

2017-02-09
version 2.3
1) In version 2.2, we used a outer for loop to iterate duration range, and a inner for loop to iterate sliding window. We create a multiprocessing pool before the beginning of inner for loop, and put all worker processes generated in inner for loop into the processes pool. When the inner for loop finish, we use pool.close() and pool.join() to wait for the worker porcesses to finish. Then save the results generated in innner for loop. After that, we delete pool object and collect memory. For the next iteration in outer for loop, we create a new pool object and do the same thing as above. This code works well on the bladerunner server in UAH, however, does not work properly on blueshark server in FIT. When running on blueshark, once the outer loop steps into the second iteration, the server becomes very slow. In first outer iteration, it takes 20 minutes to go through one month, but in the second outer iteration, it takes 2 hours or more. I think it may due to the version of multiprocessing package or python, which cannot handle mulitpy pools properly. In this version, we create only one multiprocessing pool outside the outer loop, and close it when outer for loop finish. On each end of inner loop, we use trick to block the main processe to wait for the worker processes to finish and to save the results from inner loop.

2017-02-10
version 2.4
1) Change the residue calculating formula. Divided by N within the square root.

2017-02-11
version 3.0
1) A new major version. Retrun more informations.

2017-06-10
version 3.1
1) In this version, we change the A value smoothing method. Firstly, downsample A to 20 points, then apply savgol_filter, then upsample to original data points number. In this way, we do not need to specify savgol_filter smoothing window size for different size detection windows.
2) Improve sliding window generating method.
3) Use a new duration tuple specify method. User provide the min and max duration, and window size range step, program calcuate the duration tuples.

2019-10-01
version 3.1.*
1) Apply the full expression for Pt = NKT + Bz^2/2mu0.

2020-04
version 3.1.3
1) Apply the extended GS-based equation (Teh 2018, EP&S), which has new Pt'(A'):
    A'(x,0) ~ -(1-a)By, a = MA^2
    Pt' = (1-a)^2 * Bz^2/2mu0 + (1-a)p + a(1-a)B^2/2mu0

'''

############################################ Command line format ############################################

# GS_detectFluxRope_multiprocessing_v3.1.py 1996 10 50 5
# First parameter is year, second is minimum size, third is maxmum size, fourth is size step.
# The command above generate the following size tuple:
# ((10, 15), (15, 20), (20, 25), (25, 30), (30, 35), (35, 40), (40, 45), (45, 50))
# More info can be found on Hu et al. 2018 (doi: 10.3847/1538-4365/aae57d)

#############################################################################################################


from __future__ import division # Treat integer as float.
import os
import gc
import sys
import math
import time
import pickle
import multiprocessing
import pandas as pd
import matplotlib.dates as mdates
import numpy as np
from numpy import linalg as la
import scipy as sp
from scipy import integrate, stats
from scipy.signal import savgol_filter
from datetime import datetime, timedelta
from ast import literal_eval
# import spacepy.toolbox as tb
# tb.update(leapsecs=True)
import fluxrope as FR
import warnings
warnings.filterwarnings(action="ignore", module="scipy", message="^internal gelsd")
#############################################################################################################

# Physics constants.
global mu0              # (N/A^2) magnetic constant permeability of free space vacuum permeability
global m_proton         # Proton mass [kg]
global factor_deg2rad   # Convert degree to radians
global k_Boltzmann      # Boltzmann constant, in J/K.

mu0 = 4.0 * np.pi * 1e-7     # (N/A^2)
m_proton = 1.6726219e-27     # kg
factor_deg2rad = np.pi/180.0 # radians
k_Boltzmann = 1.3806488e-23  # J/K

#############################################################################################################

# Choose root directory according to environment.
rootDir = '/home/rharvey/Documents/GS-Analysis/'
dataDir = '/home/rharvey/data/'
verbose = True
flag_smoothA = True

###############################################################################################################
# Command line argument.
time_head      = sys.argv[1]+' '+sys.argv[2]
time_tail      = sys.argv[3]+' '+sys.argv[4]
size_start_str = sys.argv[5]
size_end_str   = sys.argv[6]
size_step_str  = sys.argv[7]
overlap_str    = sys.argv[8]
namestr        = sys.argv[9]
probe_str      = sys.argv[10]

size_start = int(size_start_str)
size_end   = int(size_end_str)
size_step  = int(size_step_str)
overlap    = int(overlap_str)

if verbose:
    print(f'size_start_str={size_start_str}')
    print(f'size_end_str={size_end_str}')
    print(f'size_step_str={size_step_str}')
    print(f'overlap_str={overlap_str}')

duration_range_list = []
size_start_one_iter_temp = size_start
size_end_one_iter_temp = size_start
while size_end_one_iter_temp < size_end:
    size_end_one_iter_temp = size_start_one_iter_temp + size_step
    duration_range_list.append([size_start_one_iter_temp, size_end_one_iter_temp])
    size_start_one_iter_temp = size_end_one_iter_temp
if (duration_range_list[-1][-1] > size_end):
    duration_range_list[-1][-1] = size_end

'''
for i in range(len(duration_range_list)):
    if i==0:
        duration_range_list[i][1] += overlap
    elif i==len(duration_range_list)-1:
        duration_range_list[i][0] -= overlap
    else:
        duration_range_list[i][0] -= overlap
        duration_range_list[i][1] += overlap
'''

# Shift boundary to make overlap.
for i in range(len(duration_range_list)):
    duration_range_list[i][0] -= overlap
    duration_range_list[i][1] += overlap

duration_range_tuple = tuple(duration_range_list)
duration_range_tuple = tuple(tuple(element) for element in duration_range_list)

# Set log.
fr_log_bufsize = 1 # 0 means unbuffered, 1 means line buffered.
fr_log_filename = 'fr' + namestr + '_' + size_start_str + '_' + size_end_str + '_' + size_step_str + '.log'
fr_log_path = rootDir +'log/'
fr_log_path_filename = fr_log_path + fr_log_filename
fr_log = open(fr_log_path_filename, 'w', fr_log_bufsize)
sys.stdout = fr_log

# Set search range.
searchDatetimeStart = datetime.strptime(time_head,'%Y-%m-%d %H:%M')
searchDatetimeEnd   = datetime.strptime(time_tail,'%Y-%m-%d %H:%M')

if searchDatetimeStart.day == searchDatetimeEnd.day:
    date_str = searchDatetimeStart.strftime('%d %B %Y %H:%M:%S') + '--' + searchDatetimeEnd.strftime('%H:%M:%S')
else:
    date_str = searchDatetimeStart.strftime('%d %B %Y %H:%M:%S') + '--' + searchDatetimeEnd.strftime('%d %B %Y %H:%M:%S')

print('\nDuration range is: {}'.format(duration_range_tuple))
print('Searching fluxrope in {}.'.format(date_str))

# If output folder for specific date does not exist, create it.
if not os.path.exists(rootDir + 'output/' + namestr[1:] + probe_str):
    os.makedirs(rootDir + 'output/' + namestr[1:] + probe_str)

# =================================== Read and Check data =======================================
# Read in data.
print('Reading data...')
GS_DataFrame = pd.read_csv(dataDir + 'data' + probe_str + namestr + '.csv', index_col=0)
dt = np.diff(pd.to_datetime(GS_DataFrame.index))[0] / pd.Timedelta('1s')
print(dt)

# Check data property.
print('Checking DataFrame keys... {}'.format(GS_DataFrame.keys()))
print('Checking DataFrame shape... {}'.format(GS_DataFrame.shape))
print('Data Time start: {}'.format(GS_DataFrame.index[0]))
print('Data Time end: {}'.format(GS_DataFrame.index[-1]))
print('Checking the number of NaNs in DataFrame...')
len_GS_DataFrame = len(GS_DataFrame)
for key in GS_DataFrame.keys():
    num_notNaN = GS_DataFrame[key].isnull().values.sum()
    percent_notNaN = 100.0 - num_notNaN * 100.0 / len_GS_DataFrame
    print('The number of NaNs in {} is {}, integrity is {}%'.format(key, num_notNaN, round(percent_notNaN, 2)))

'''
# DO NOT DELETE THIS COMMENT.
# Good for insertion, not good for range selection. if datetime = datetime(year, 1,1,2,59,59), then the returned index is from datetime(year, 1,1,2,59,0). However, if datetime = datetime(year, 1,1,2,59,0), then the returned index is from datetime(year, 1,1,2,58,0).
# Get the start and end DataFrame indices according to the start and end datetime.
index_start = GS_DataFrame.index.searchsorted(searchDatetimeStart)
index_end   = GS_DataFrame.index.searchsorted(searchDatetimeEnd) # index_end not include searchDatetimeEnd.
# Get the records between start and end time from DataFrame.
GS_DataFrame = GS_DataFrame.iloc[index_start:index_end+1] #.iloc works on location.
'''

# Get slice of data
selectedRange_mask = (pd.to_datetime(GS_DataFrame.index) >= searchDatetimeStart) & (pd.to_datetime(GS_DataFrame.index) <= searchDatetimeEnd)
GS_DataFrame  = GS_DataFrame.iloc[selectedRange_mask]
B_DataFrame   = GS_DataFrame[['Bx','By','Bz']] 
Vsw_DataFrame = GS_DataFrame[['Vx','Vy','Vz']]
Np_DataFrame  = GS_DataFrame['Np']
Tp_DataFrame  = GS_DataFrame['Tp']

# Multiprocessing
num_cpus = multiprocessing.cpu_count()
max_processes = num_cpus
print ('\nTotol CPU cores on this node = ', num_cpus)
pool = multiprocessing.Pool(processes=max_processes) # Create a multiprocessing pool with safe_lock.
results = []                                         # Create a list to save result.

num_FluxRope = 0
search_result_raw_true = {}
search_result_raw_false = {}
totalStartTime = datetime.now()

# Theta grid number: 90/9=10, d_theta=10(degree); 90/12=7.5, d_theta=7.5(degree)
n_theta_grid = 9
print('Grid size of theta & phi = {} & {}'.format(90/n_theta_grid, 180/n_theta_grid))

# First integer in tuple is minimum duration threshold, second integer in tuple is searching window width.
print('\nDuration range tuple is: {}'.format(duration_range_tuple))

# Apply GS detection in sliding window (window width).
for duration_range in duration_range_tuple:
    startTime = datetime.now()
    minDuration = duration_range[0]
    maxDuration = duration_range[1]
    print('Duration : {} ~ {} minutes.'.format(minDuration, maxDuration))
    
    '''
    # Choose a flexible savgol filter window width based on the length of minDuration.
    half_minDuration = minDuration//2
    half_maxDuration = maxDuration//2
    if (half_minDuration) % 2 == 0: # filter window must be odd.
        savgol_filter_window = half_minDuration + 1
    else:
        savgol_filter_window = half_minDuration
    print('savgol_filter_window = {}'.format(savgol_filter_window))
    '''
    
    # Flexible interpolation limit based on window length. 
    # The maximum gap tolerance is up to 30% of total points count.
    interp_limit = int(math.ceil(minDuration*3.0/10)) # interp_limit = 1
    # print('interp_limit = {}'.format(interp_limit))
    
    # Sliding window. If half_window=15, len(DataFrame)=60, range(15,45)=[15,...,44].
    for indexFluxRopeStart in range(len(B_DataFrame) - maxDuration): # minutes.
        indexFluxRopeEnd = indexFluxRopeStart + maxDuration - 1  # The end point is included, so -1.
        
        B_inWindow = B_DataFrame.iloc[indexFluxRopeStart : indexFluxRopeEnd+1] # End is not included.
        # If there is any NaN in this range, try to interpolate. For example, limit=3 means only interpolate the gap shorter than 4.
        if B_inWindow.isnull().values.sum():
            B_inWindow_copy = B_inWindow.copy(deep=True)
            B_inWindow_copy.interpolate(method='time', limit=interp_limit, inplace=True)
            if B_inWindow_copy.isnull().values.sum():
                print('Encounter NaN in B field data, skip this iteration.')
                continue
            else:
                B_inWindow = B_inWindow_copy

        Vsw_inWindow = Vsw_DataFrame.iloc[indexFluxRopeStart : indexFluxRopeEnd+1] # End is not included.
        if Vsw_inWindow.isnull().values.sum():
            Vsw_inWindow_copy = Vsw_inWindow.copy(deep=True) 
            Vsw_inWindow_copy.interpolate(method='time', limit=interp_limit, inplace=True)
            if Vsw_inWindow_copy.isnull().values.sum():
                print('Encounter NaN in Vsw data, skip this iteration.')
                continue
            else:
                Vsw_inWindow = Vsw_inWindow_copy

        Np_inWindow = Np_DataFrame.iloc[indexFluxRopeStart : indexFluxRopeEnd+1]
        if Np_inWindow.isnull().values.sum():
            Np_inWindow_copy = Np_inWindow.copy(deep=True)
            Np_inWindow_copy.interpolate(method='time', limit=interp_limit, inplace=True)
            if Np_inWindow_copy.isnull().values.sum():
                print('Encounter NaN in Np data, skip this iteration.')
                continue
            else:
                Np_inWindow = Np_inWindow_copy

        Tp_inWindow = Tp_DataFrame.iloc[indexFluxRopeStart : indexFluxRopeEnd+1]
        if Tp_inWindow.isnull().values.sum():
            Tp_inWindow_copy = Tp_inWindow.copy(deep=True)
            Tp_inWindow_copy.interpolate(method='time', limit=interp_limit, inplace=True)
            if Tp_inWindow_copy.isnull().values.sum():
                print('Encounter NaN in Tp data, skip this iteration.')
                continue
            else:
                Tp_inWindow = Tp_inWindow_copy
        
        # Calculate VHT in GSE frame. Calculating VHT takes very long time(0.02748s for 14 data points), so we use mean Vsw as VHT.
        VHT_inGSE = FR.findVHT(B_inWindow, Vsw_inWindow) # Very slow.
        # VHT_inGSE = np.array(Vsw_inWindow.mean())
        
        # Return value: timeRange, Residue, orientation
        result_temp = pool.apply_async(FR.searchFluxRopeInWindow, args=(B_inWindow, VHT_inGSE, n_theta_grid, minDuration, dt, flag_smoothA, Np_inWindow, Tp_inWindow, Vsw_inWindow))
        # print(result_temp.get()) # This statement will cause IO very slow.
        results.append(result_temp)
        # DO NOT unpack result here. It will block IO. Unpack in bulk.

    # Next we are going to save file We have to wait for all worker processes to finish.
    # Block main process to wait for worker processes to finish. This while loop will execute almost immediately when the innner for loop goes through. The inner for loop is non-blocked, so it finish in seconds.
    while len(pool._cache)!=0:
        print('{} - Waiting... There are {} worker processes in pool.'.format(time.ctime(), len(pool._cache)))
        time.sleep(5)
    print('{} - len(pool._cache) = {}'.format(time.ctime(), len(pool._cache)))
    print('{} - Duration range {}~{} minutes is completed!'.format(time.ctime(), minDuration, maxDuration))

    # Save result. One file per window size.
    results_true_tuple_list  = []
    results_false_tuple_list = []
    
    # Unpack results. Convert to tuple, and put into list.
    for one_result in results:
        results_tuple_temp = (one_result.get())
        if not np.isinf(results_tuple_temp[5]):                     # Check residue.
            if results_tuple_temp[4]:
                results_true_tuple_list.append(results_tuple_temp)  # Turn point on top.
            else:                                                       
                results_false_tuple_list.append(results_tuple_temp) # Turn point on bottom.

    # Save results to dictionary. One key per window size.
    key_temp = str(minDuration) + '~' + str(maxDuration)
    search_result_raw_true[key_temp]  = results_true_tuple_list
    search_result_raw_false[key_temp] = results_false_tuple_list
    results = [] # Empty container results

    #########
    pickle_filename = rootDir + 'output/' + namestr[1:] + probe_str + '/'
    pickle_filename_true  = pickle_filename + 'true_' + key_temp + 'min.p'
    pickle_filename_false = pickle_filename + 'false_' + key_temp + 'min.p'

    print('Save file to: {}'.format(pickle_filename_true))
    pickle.dump(results_true_tuple_list, open(pickle_filename_true, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
    print('Save file to: {}'.format(pickle_filename_false))
    pickle.dump(results_false_tuple_list, open(pickle_filename_false, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
    ##########

    endTime = datetime.now()
    time_spent_in_seconds = (endTime - startTime).total_seconds()
    print('Time spent on this window: {} seconds ({} hours).'.format(time_spent_in_seconds, time_spent_in_seconds/3600.0))

pool.close() # Close pool, prevent new worker process from joining.
pool.join()  # Block caller process until workder processes terminate.

totalEndTime = datetime.now()
time_spent_in_seconds = (totalEndTime - totalStartTime).total_seconds()
print('\n{}'.format(time.ctime()))
print('Number of CPU cores per node: {}.'.format(num_cpus))
print('Max number of worker processes in pool: {}.'.format(max_processes))
print('Total Time spent: {} seconds ({} hours).'.format(time_spent_in_seconds, time_spent_in_seconds/3600.0))

# Load pickle file.
result = pickle.load(open(pickle_filename_true, 'rb'))
for item in result:
    print(item)