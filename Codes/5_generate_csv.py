import os
import sys
import pickle
import numpy as np
import pandas as pd
from datetime import datetime

rootDir   = '/home/rharvey/Documents/GS-Analysis/'
inputDir  = 'SearchResult/selected_events/'
outputDir = rootDir + 'Events/'

def addTableRow(fileToWrite, data_list, background_color):
    fileToWrite.write('<tr> <!-- talbe row -->\n')
    for i in range(len(data_list)):
        if i == 1:
            trimed_time_str = data_list[i].replace('/', '')
            trimed_time_str = trimed_time_str.replace(' ', '')
            trimed_time_str = trimed_time_str.replace(':', '')
            trimed_time_str = trimed_time_str.replace('~', '_')
            URL = 'events/' + trimed_time_str + '.html'
            fileToWrite.write('<td align="center" valign="middle" bgcolor="' + background_color + '">' + '<a href="' + URL + '">' + data_list[i] + '</a>' + '</td>\n')
        else:
            fileToWrite.write('<td align="center" valign="middle" bgcolor="' + background_color + '">' + data_list[i] + '</td>\n')
    fileToWrite.write('</tr>\n')
#############################################################################################################

# Command line argument.
namestr   = sys.argv[1]
probe_str = sys.argv[2]

# Load DataFrame
inputFileName    = namestr[1:] + probe_str + '(|B|>0)_selected_events.p'
FR_detailed_info = pickle.load(open(rootDir + inputDir + inputFileName,'rb'))

# Save data to CSV file.
print('\nSaving results data to CSV file...')
columns_to_csv = ['startTime', 'endTime', 'duration','residue_fit', 'residue_diff', 'B_abs_mean','B_magnitude_max', 'Beta_mean', 'Beta_p_mean', 'Vsw_magnitude_mean', 'Tp_mean', 'theta_deg', 'phi_deg', 'Z_unitVector[0]', 'Z_unitVector[1]', 'Z_unitVector[2]','size_inAU','walenTest_slope','walenTest_slope_b4reverse','Np_mean', 'Mach_average','A_range','Bx_inFR_abs_mean', 'By_inFR_abs_mean', 'Bz_inFR_abs_mean','VHT_inGSE[0]', 'VHT_inGSE[1]', 'VHT_inGSE[2]','walenTest_r_value','walenTest_r_value_b4reverse','Cross_heli','Br_mean','Residue_energy','Jzmax','turnTime']
csv_header = ['Start Time', 'End Time', 'Duration','Residue_fit', 'Residue_diff', '<B>(nT)','Bmax(nT)', '<Beta>', '<Beta_p>', '<Vsw>(km/s)', '<Tp>(10^6K)', 'theta(deg)', 'phi(deg)', 'z_axis[0](in RTN)', 'z_axis[1](in RTN)', 'z_axis[2](in RTN)','size_inAU','walenTest_slope','walenTest_slope_b4reverse','Np_mean', 'Mach_average','A_range','Bx_inFR_abs_mean', 'By_inFR_abs_mean', 'Bz_inFR_abs_mean','VHT_inRTN[0]', 'VHT_inRTN[1]', 'VHT_inRTN[2]','walenTest_r_value','walenTest_r_value_b4reverse','Cross_heli','Br_mean','Residue_energy','Jzmax', 'turnTime']
FR_detailed_csv = FR_detailed_info[columns_to_csv]
FR_detailed_csv.columns = csv_header
FR_detailed_csv.to_csv(path_or_buf=outputDir + 'GS_events' + namestr + probe_str + '.csv')
print('\nGS_events' + namestr + probe_str + '.csv is saved!')
print(FR_detailed_csv[['Start Time', 'End Time', 'Duration', 'Residue_fit', 'Residue_diff', 'walenTest_slope', 'walenTest_r_value', 'Mach_average']])