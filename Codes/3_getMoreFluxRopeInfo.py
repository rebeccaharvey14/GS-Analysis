# INPUT  : no_overlap record.
# OUTPUT : flux rope records with detailed information.
# NOTE   : Modularize code. Use MyPythonPackage.fluxrope module.

import os
import sys
import pickle
import numpy as np
from datetime import datetime
import pandas as pd
import fluxrope as FR

pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 100)
pd.set_option('display.width', 200)

rootDir   = '/home/rharvey/Documents/GS-Analysis/'
dataDir   = '/home/rharvey/data/'
outputDir = rootDir + 'SearchResult/detailed_info/'
inputDir  = 'SearchResult/no_overlap/'
###############################################################################

# Command line argument.
namestr   = sys.argv[1]
probe_str = sys.argv[2]

# Read in all data.
print('Reading data...')
DataFrame = pd.read_csv('/home/rharvey/data/data' + namestr + probe_str '.csv',index_col=0)

# Get more flux rope information.
search_result_no_overlap_DF_fileName = rootDir + inputDir + namestr[:-1] + probe_str + '_no_overlap.p'
search_result_detail_info_DF = FR.get_more_flux_rope_info(DataFrame, search_result_no_overlap_DF_fileName, output_dir=outputDir, output_filename=namestr[:-1] + probe_str+'_detailed_info.p')