import sys
import numpy as np
import pandas as pd
import pyspedas
from datetime import datetime, timedelta
from pytplot import get_data
####################################################################

time_head  = sys.argv[1] + ' ' + sys.argv[2]
time_tail  = sys.argv[3] + ' ' + sys.argv[4]
time_range = [time_head,time_tail]
namestr    = sys.argv[5]
probe      = sys.argv[6]
year       = sys.argv[7]

rootDir = '/home/rharvey/Documents/Research/GS-Analysis/'
dataFile = '/home/rharvey/data/' + 'data_MMS'+ probe + namestr + '.csv'

datetimeStart = datetime.strptime(time_head,'%Y-%m-%d %H:%M')
datetimeEnd   = datetime.strptime(time_tail,'%Y-%m-%d %H:%M')

fr_log_bufsize = 1 # 0 means unbuffered, 1 means line buffered.
fr_log_path_filename = rootDir + 'fr' + namestr + '.log'
fr_log = open(fr_log_path_filename, 'w', fr_log_bufsize)
sys.stdout = fr_log

######################### FGM measurements #########################
print('Getting FGM (magnetic field) data...')
pyspedas.mms.mms_load_fgm(probe=probe, trange=time_range, varformat='*_b_gse_srvy',time_clip=True)
[B_epoch, b_gse_vec] = get_data('mms'+probe+'_fgm_b_gse_srvy_l2') # nT

######################### FPI measurements #########################
print('Getting FPI (velocity, density, temperature) data...')
variables = pyspedas.mms.mms_load_fpi(probe=probe, trange=time_range, time_clip=True,center_measurement=True)
[V_epoch, u] = get_data('mms'+probe+'_dis_bulkv_gse_fast')             # km/s
[V_epoch, Np] = get_data('mms'+probe+'_dis_numberdensity_fast')        # cm^-3
[V_epoch, Temp] = get_data('mms'+probe+'_dis_temptensor_gse_fast')     # eV
if 'mms'+probe+'_des_numberdensity_fast' in variables:
    [V_epoch, Npe] = get_data('mms'+probe+'_des_numberdensity_fast')   # cm^-3
    [V_epoch, Temp_e] = get_data('mms'+probe+'_des_temptensor_gse_fast')# eV
    Te = np.empty(Temp_e.shape[0])
Vmag =  np.array([np.sqrt(u[i][0]**2 + u[i][1]**2 + u[i][2]**2) for i in range(V_epoch.size)])

# Proton and electron temperature
Tp = np.empty(Temp.shape[0])
for i in range(Temp.shape[0]):
    Tp[i] = np.trace(Temp[i])/3*11604.518    # eV --> K
    if 'mms'+probe+'_des_temptensor_gse_fast' in variables:
        Te[i] = np.trace(Temp_e[i])/3*11604.518  # eV --> K


####################################################################
# Trim data to specified range
####################################################################
print('Trimming data to specified time range...')
BTime = np.array([datetime.utcfromtimestamp(t) for t in B_epoch])
selected_index = [(BTime> datetimeStart) & (BTime < datetimeEnd)][0]
BTime = BTime[selected_index]
BGSE = b_gse_vec[selected_index]

VTime = np.array([datetime.utcfromtimestamp(t) for t in V_epoch])
selected_index = [(VTime > datetimeStart) & (VTime < datetimeEnd)][0]
VTime = VTime[selected_index]
VGSE = u[selected_index]
Np = Np[selected_index]
Tp = Tp[selected_index]
Vmag = Vmag[selected_index]
if 'mms'+probe+'_des_temptensor_gse_fast' in variables:
    Te = Te[selected_index]

# Time cadence
dt = datetime.utcfromtimestamp(V_epoch[1])-datetime.utcfromtimestamp(V_epoch[0])
dt = round(dt.seconds + dt.microseconds*1e-6,4)


####################################################################
# Put data into DataFrame
####################################################################
print('Putting BGSE, VGSE, Np, Tp (Te) into DataFrame...')
BGSE_DataFrame = pd.DataFrame(np.column_stack(BGSE.T), index = BTime, columns = ['Bx', 'By', 'Bz', 'Bmag']).resample(f'{dt}S').mean()
VGSE_DataFrame = pd.DataFrame(np.column_stack((VGSE,Vmag)), index = VTime, columns = ['Vx', 'Vy', 'Vz', 'Vmag']).resample(f'{dt}S').mean()
Np_DataFrame = pd.DataFrame(Np, index = VTime, columns = ['Np']).resample(f'{dt}S').mean()
Tp_DataFrame = pd.DataFrame(Tp, index = VTime, columns = ['Tp']).resample(f'{dt}S').mean()
if 'mms'+probe+'_des_temptensor_gse_fast' in variables:
    Te_DataFrame = pd.DataFrame(Te, index = VTime, columns = ['Te']).resample(f'{dt}S').mean()

# Alfven speed & magnitude, Thermal & Magnetic Pressure, Plasma Beta
print('Putting bGSE, beta into DataFrame...')
bGSE_DataFrame = (1e-9)*BGSE_DataFrame.div(np.sqrt(4*np.pi*Np_DataFrame['Np']*(1.67e-27)*(1e-1)),axis=0).rename(columns={'Bx':'bx','By':'by','Bz':'bz','Bmag':'bmag'})*(1e-3)
if 'mms'+probe+'_des_temptensor_gse_fast' in variables:
    p_DataFrame = Tp_DataFrame.add(Te_DataFrame['Te'],axis=0).mul(Np_DataFrame['Np']*(1.3806e-17),axis=0).rename(columns={'Tp':'p'})
else:
    p_DataFrame = Tp_DataFrame.mul(Np_DataFrame['Np']*(1.3806e-17),axis=0).rename(columns={'Tp':'p'})
Bp = ((BGSE[3]*(1e-9))**2)/2/(4*np.pi*(1e-7))
beta_DataFrame = p_DataFrame.div(((BGSE_DataFrame['Bmag']*(1e-9))**2)/2/(4*np.pi*(1e-7)),axis=0).rename(columns={'p':'beta'})


####################################################################
# Get spacecraft location
####################################################################
print('Getting ephemeris data...')
km_in_re = 6371.2
pyspedas.mms.mms_load_mec(probe=probe, trange=time_range, time_clip=True, varformat='*_r_gse')
[eph_epoch, pos] = get_data('mms' + probe + '_mec_r_gse') #km
x = pos[:,0]/6371.2
y = pos[:,1]/6371.2
z = pos[:,2]/6371.2
xyz_GSE = pd.DataFrame.from_dict(data = {'x_gse':x, 'y_gse':y, 'z_gse':z})
xyz_GSE.index = np.array([datetime.utcfromtimestamp(t) for t in eph_epoch])
xyz_GSE = xyz_GSE.resample(f'{dt}S').mean().interpolate(method='linear')


####################################################################
# Merge all DataFrames into one according to time index
####################################################################
print('Merging all DataFrames...')
DataFrame = pd.DataFrame(index=BGSE_DataFrame.index)
if 'mms'+probe+'_des_temptensor_gse_fast' in variables:
    DataFrame = pd.concat([DataFrame, BGSE_DataFrame, VGSE_DataFrame, Np_DataFrame, Tp_DataFrame, Te_DataFrame, bGSE_DataFrame, beta_DataFrame, xyz_GSE], axis=1)
else:
    DataFrame = pd.concat([DataFrame, BGSE_DataFrame, VGSE_DataFrame, Np_DataFrame, Tp_DataFrame,bGSE_DataFrame, beta_DataFrame, xyz_GSE], axis=1)

# Drop duplicated records and sort data by time index
print('Dropping duplicated records and sorting by time index...')
DataFrame.interpolate(method='linear',inplace=True)
DataFrame.drop_duplicates(keep='first', inplace=True)
DataFrame.sort_index(axis=0, ascending=True, inplace=True, kind='quicksort')
#DataFrame[~DataFrame.apply(lambda x: np.abs(x - x.mean()) / x.std() < 3).all(axis=1)] = np.nan

print('Checking the number of NaNs in DataFrame...')
len_DataFrame = len(DataFrame)
for key in DataFrame.keys():
    num_notNaN = DataFrame[key].isnull().values.sum()
    percent_notNaN = 100.0 - num_notNaN * 100.0 / len_DataFrame
    print('The number of NaNs in {} is {}, integrity is {}%'.format(key, num_notNaN, round(percent_notNaN, 2)))

DataFrame.dropna(inplace=True)


####################################################################
# Save merged DataFrame into csv file
####################################################################
print(f'Saving data file to: {dataFile}')
DataFrame.to_csv(dataFile)