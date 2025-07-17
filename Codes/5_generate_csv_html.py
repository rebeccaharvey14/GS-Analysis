import os
import sys
import pickle
import numpy as np
import pandas as pd
from datetime import datetime

rootDir   = '/home/rharvey/GS-Analsis/'
inputDir  = 'SearchResult/selected_events/'
outputDir = 'Events/'

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
year_str  = sys.argv[1]
year      = int(year_str)
namestr   = sys.argv[2]
probe_str = sys.argv[3]

# Load DataFrame
inputFileName = namestr[:-1] + probe_str + '(|B|>0)_selected_events.p'
FR_detailed_info = pickle.load(open(rootDir + inputDir + inputFileName,'rb'))

# Save data to CSV file.
print('\nSaving results data to CSV file...')
columns_to_csv = ['startTime', 'endTime', 'duration','residue_fit', 'residue_diff', 'B_abs_mean','B_magnitude_max', 'Beta_mean', 'Beta_p_mean', 'Vsw_magnitude_mean', 'Tp_mean', 'theta_deg', 'phi_deg', 'Z_unitVector[0]', 'Z_unitVector[1]', 'Z_unitVector[2]','size_inAU','walenTest_slope','walenTest_slope_b4reverse','Np_mean', 'Mach_average','A_range','Bx_inFR_abs_mean', 'By_inFR_abs_mean', 'Bz_inFR_abs_mean','VHT_inGSE[0]', 'VHT_inGSE[1]', 'VHT_inGSE[2]','walenTest_r_value','walenTest_r_value_b4reverse','Cross_heli','Br_mean','Residue_energy','Jzmax','turnTime']
csv_header = ['Start Time', 'End Time', 'Duration','Residue_fit', 'Residue_diff', '<B>(nT)','Bmax(nT)', '<Beta>', '<Beta_p>', '<Vsw>(km/s)', '<Tp>(10^6K)', 'theta(deg)', 'phi(deg)', 'z_axis[0](in RTN)', 'z_axis[1](in RTN)', 'z_axis[2](in RTN)','size_inAU','walenTest_slope','walenTest_slope_b4reverse','Np_mean', 'Mach_average','A_range','Bx_inFR_abs_mean', 'By_inFR_abs_mean', 'Bz_inFR_abs_mean','VHT_inRTN[0]', 'VHT_inRTN[1]', 'VHT_inRTN[2]','walenTest_r_value','walenTest_r_value_b4reverse','Cross_heli','Br_mean','Residue_energy','Jzmax', 'turnTime']
FR_detailed_csv = FR_detailed_info[columns_to_csv]
FR_detailed_csv.columns = csv_header
FR_detailed_csv.to_csv(path_or_buf=rootDir + outputDir + 'events' + namestr + probe_str + '.csv')
print('\nevents' + namestr + probe_str + '.csv is saved!')
print(FR_detailed_csv[['Start Time', 'End Time', 'Duration', 'Residue_fit', 'Residue_diff', 'walenTest_slope', 'walenTest_r_value', 'Mach_average']])

'''
# Constructing webpage.
print('\nCreating webpage...')
html_filename = 'year'+year_str+'.html'
html_webpage_file = open(rootDir + outputWebDir + html_filename,'w')
f = html_webpage_file # html heading.
f.write('<html>\n')   # python will convert \n to os.linesep

# link on top.
# f.write('<div style="width:100px; height:100px; position:fixed;">\n')
f.write('<a href="../index.html" style="font-family:arial; text-decoration:none; font-size:13"> HOME </a>\n')
f.write('<br />\n')
if year==year_start:
    f.write('<a href="../index.html" style="font-family:arial; text-decoration:none;"> &lsaquo;&ndash; </a> &nbsp;\n')
    f.write('<a href="../PSPE'+str(year+1)+'/E'+str(year+1)+'.html" style="font-family:arial; text-decoration:none;"> &ndash;&rsaquo; </a>\n')
elif year==year_end:
    f.write('<a href="../PSPE'+str(year-1)+'/E'+str(year-1)+'.html" style="font-family:arial; text-decoration:none;"> &lsaquo;&ndash; </a> &nbsp;\n')
    f.write('<a href="../index.html" style="font-family:arial; text-decoration:none;"> &ndash;&rsaquo; </a>\n')
else:
    f.write('<a href="../PSPE'+str(year-1)+'/E'+str(year-1)+'.html" style="font-family:arial; text-decoration:none;"> &lsaquo;&ndash; </a> &nbsp;\n')
    f.write('<a href="../PSPE'+str(year+1)+'/E'+str(year+1)+'.html" style="font-family:arial; text-decoration:none;"> &ndash;&rsaquo; </a>\n')
#f.write('</div>\n')
f.write('<br />\n')

# Webpage head.
f.write('<head>\n')
f.write('<title>Encounter' + year_str + '</title>\n')
f.write('</head>\n')

# Open webpage body.
f.write('<body>\n')

# Content title
f.write('<h2><center><font size="5">Small-scale Flux Rope Events in ' + year_str +'</font><center></h2>\n')
f.write('<hr width="1096px"/>\n')

# Open data table.
f.write('<table style="margin:1em auto; font-size:16; white-space:nowrap;">\n')                                                            # Table heading.
f.write('<tr> <!-- talbe row -->\n')                                                                                                       # Start table row.
f.write('<th scope="col" width="50" align="center" valign="middle" bgcolor="#6194C9"> No. </th>\n')                                        # Heading: 1) No.
f.write('<th scope="col" width="280" align="center" valign="middle" bgcolor="#6194C9"> Time Range </th>\n')                                # Heading: 2) Time range.
f.write('<th scope="col" width="70" align="center" valign="middle" bgcolor="#6194C9"> Duration </th>\n')                                   # Heading: 3) duration.
# f.write('<th scope="col" width="80" align="center" valign="middle" bgcolor="#6194C9"> Wait Time </th>\n')                                # Heading: 4) Wait Time.
f.write('<th scope="col" width="80" align="center" valign="middle" bgcolor="#6194C9"> Residue </th>\n')                                    # Heading: 5) Residue.
f.write('<th scope="col" width="80" align="center" valign="middle" bgcolor="#6194C9"> &lt;B&gt; (nT) </th>\n')                             # Heading: 6) Bmean.
f.write('<th scope="col" width="80" align="center" valign="middle" bgcolor="#6194C9"> B<sub>max</sub> (nT) </th>\n')                       # Heading: 7) Bmax.
f.write('<th scope="col" width="100" align="center" valign="middle" bgcolor="#6194C9"> &lt;&beta;&gt;, &lt;&beta;<sub>p</sub>&gt;</th>\n') # Heading: 8) <Beta>.
f.write('<th scope="col" width="100" align="center" valign="middle" bgcolor="#6194C9"> &lt;Vsw&gt; (km/s) </th>\n')                        # Heading: 9) <Vsw>.
f.write('<th scope="col" width="100" align="center" valign="middle" bgcolor="#6194C9"> &lt;Tp&gt; (10<sup>6</sup>K) </th>\n')              # Heading: 10) <Tp>.
f.write('<th scope="col" width="120" align="center" valign="middle" bgcolor="#6194C9"> Orientation </th>\n')                               # Heading: 11) Orientation.
f.write('</tr>\n')                                                                                                                         # End table row.

# '#D9E0EC' is darker.
color_list = ['#D9E0EC', '#EDF1F6']

for i in range(len(FR_detailed_info)):
    print('i = {}'.format(i))
    # Prepare table row list.
    # Extract one record.
    # DataFrame.iloc[[]] will keep original index.
    # DataFrame.iloc[] will create new index.
    oneRecord_DF = FR_detailed_info.iloc[[i]]
    
    # 1) No.
    number_str = str(oneRecord_DF.index[0] + 1)
    
    # 2) timeRange
    startTime = oneRecord_DF['startTime'].iloc[0]
    endTime = oneRecord_DF['endTime'].iloc[0]
    startTime_str = str(startTime.strftime('%Y/%m/%d %H:%M:%S'))
    endTime_str = str(endTime.strftime('%Y/%m/%d %H:%M:%S'))
    timeRange_str = startTime_str + ' ~ ' + endTime_str
    
    # 3) duration
    duration = oneRecord_DF['duration'].iloc[0]
    duration_str = str(int(duration/60.0))

    # 5) residue
    residue_fit = oneRecord_DF['residue_fit'].iloc[0]
    residue_fit = round(residue_fit,4)
    residue_fit_str = "{:.4f}".format(residue_fit)
    
    # 6) Bmean
    Bmean = oneRecord_DF['B_abs_mean'].iloc[0]
    Bmean_str = "{:.2f}".format(round(Bmean, 2))
    
    # 7) Bmax
    Bmax = oneRecord_DF['B_magnitude_max'].iloc[0]
    Bmax = round(Bmax, 2)
    Bmax_str = "{:.2f}".format(round(Bmax, 2))
    
    # 8) <Beta>
    Beta_mean = oneRecord_DF['Beta_mean'].iloc[0]
    Beta_p_mean = oneRecord_DF['Beta_p_mean'].iloc[0]
    Beta_mean = round(Beta_mean, 2)
    Beta_p_mean = round(Beta_p_mean, 2)
    Beta_mean_str = "{:.2f}".format(Beta_mean) + ', ' + "{:.2f}".format(Beta_p_mean)
    
    # 9) <Vsw>
    Vsw_mean = oneRecord_DF['Vsw_magnitude_mean'].iloc[0]
    Vsw_mean_str = str(int(Vsw_mean))
    
    # 10) <Tp>
    Tp_mean = oneRecord_DF['Tp_mean'].iloc[0]
    Tp_mean = round(Tp_mean, 2)
    Tp_mean_str = "{:.2f}".format(Tp_mean)
    
    # 11) orientation
    theta = oneRecord_DF['theta_deg'].iloc[0]
    phi = oneRecord_DF['phi_deg'].iloc[0]
    theta_str = str(int(theta))
    if np.isnan(phi):
        phi_str = str(9999)
    else:
        phi_str = str(round(phi))
    orientation_str = '&theta;=' + theta_str + ', ' + '&phi;=' + phi_str
    
    # 10) z_axis
    z_axis = oneRecord_DF['z_axis'].iloc[0]
    z_axis_str = str(z_axis)
    
    # Put into list
    data_list = [number_str, timeRange_str, duration_str, residue_fit_str, Bmean_str, Bmax_str, Beta_mean_str, Vsw_mean_str, Tp_mean_str, orientation_str]
    background_color = color_list[i%2]          # Set background color.
    addTableRow(f, data_list, background_color) # Add to table.

f.write('</table>\n') # Close table.
f.write('</body>\n')  # Close body.
f.write('</html>\n')  # Close html.
f.close()             # Close file.
'''