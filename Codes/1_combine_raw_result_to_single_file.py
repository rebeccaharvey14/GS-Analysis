import sys
import pickle
import pandas as pd
from datetime import datetime

# Terminal output format.
pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_columns', 100)
pd.set_option('display.width', 300)

rootDir   = '/home/rharvey/Documents/GS-Analysis/'
outputDir = '/SearchResult/raw_record_list_dict/'
#############################################################################################################

# Command line argument.
year_str  = sys.argv[1]
dt        = float(sys.argv[2]) # seconds
namestr   = sys.argv[3]
probe_str = sys.argv[4]

# Create an empty dictionary to store data.
one_dict = {'true':{}, 'timeRange':{'datetimeStart':datetime(year,1,1,0,0), 'datetimeEnd':datetime(year,12,31,23,59)}, 'false':{}}

for duration_str in ('10~16','13~24','19~33','31~45','41~55','53~66','62~76','74~88','84~98','96~109','105~119','117~131','126~172','170~216','212~258','256~302','298~344',):
                   # '342~388','384~430','428~491','487~550','548~611','607~670',):
    oneFileName_temp = rootDir + 'output/' + year_str + '/' + namestr[:-1] + '_true_' + duration_str + 'min.p'
    recordList_temp = pickle.load(open(oneFileName_temp, 'rb'))
    one_dict['true'][duration_str] = recordList_temp 

# Save raw record list
one_dict_fileName = rootDir + outputDir + 'raw_dict' + namestr + probe_str + '.p'
print(f'Saving raw record list to: {one_dict_fileName}')
pickle.dump(one_dict, open(one_dict_fileName, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
one_df = pd.DataFrame.from_dict(one_dict)
one_df.to_csv(one_dict_fileName[:-1] + 'csv')
exit()