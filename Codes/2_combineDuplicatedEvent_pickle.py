# INPUT  : raw searching result, combined into one pickle file.
# OUTPUT : flux rope list without overlapped records.
# NOTE   : Modularize code. Use MyPythonPackage.fluxrope module.

import os
import sys
import pickle
import numpy as np
from datetime import datetime, timedelta
import pandas as pd
import fluxrope as FR
import warnings
warnings.filterwarnings(action="ignore", module="scipy", message="^internal gelsd")

# Terminal output format.
pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_columns', 100)
pd.set_option('display.width', 300)

rootDir   = '/home/rharvey/Documents/GS-Analysis/'
dataDir   = '/home/rharvey/data/'
outputDir = rootDir + 'SearchResult/no_overlap/'
inputDir  = 'SearchResult/raw_record_list_dict/'
#############################################################################################################

# Command line argument.
namestr   = sys.argv[1]
probe_str = sys.argv[2]

# Read in all data.
print('Reading data...')
DataFrame = pd.read_csv('/home/rharvey/data/data' + namestr + probe_str '.csv',index_col=0)

# Check the NaNs in each variable.
print('Checking DataFrame keys... {}'.format(DataFrame.keys()))
print('Checking DataFrame shape... {}'.format(DataFrame.shape))
print('Data Time start: {}'.format(DataFrame.index[0]))
print('Data Time end: {}'.format(DataFrame.index[-1]))
print('Checking the number of NaNs in DataFrame...')
len_DataFrame = len(DataFrame)
for key in Data_DataFrame.keys():
    num_notNaN = DataFrame[key].isnull().values.sum()
    percent_notNaN = 100.0 - num_notNaN * 100.0 / len_DataFrame
    print('The number of NaNs in {} is {}, integrity is {}%'.format(key, num_notNaN, round(percent_notNaN, 2)))

# Clean up overlapped records.
oneYear_dict_fileName = rootDir + inputDir + 'raw_dict' + namestr + probe_str + '.p'
# shockList_DF = pd.read_pickle(open('/home/rharvey/GS/GS_FluxRopeDetectionPackage/shockList/IPShock_ACE_or_WIND_or_Ulysses_1996_2016_DF.p', 'rb'))

search_result_no_overlap_DF = FR.clean_up_raw_result(DataFrame, oneYear_dict_fileName, min_residue_diff=0.2, min_residue_fit=0.2, walenTest_k_threshold=0.3, output_dir=outputDir, output_filename=name_str[:-1] + probe_str +'_no_overlap.p', isPrintIntermediateDF=False, isVerbose=True, isRemoveShock=False)






