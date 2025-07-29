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

rootDir  = '/home/rharvey/Documents/Research/GS-Analysis/'
dataFile = '/home/rharvey/data/' + 'data_THM'+ probe.upper() + namestr + '.csv'

datetimeStart = datetime.strptime(time_head,'%Y-%m-%d %H:%M')
datetimeEnd   = datetime.strptime(time_tail,'%Y-%m-%d %H:%M')

######################### FGM measurements #########################
print('Getting FGM (magnetic field) data...')
pyspedas.themis.fgm(probe=probe, trange=time_range, time_clip=True, varnames=['th'+probe+'_fgs_gse','th'+probe+'_fgs_btotal'])
[B_epoch, b_gse_vec, components] = get_data('th'+probe+'_fgs_gse') # nT
[B_epoch, Btot] = get_data('th'+probe+'_fgs_btotal') # nT
Bmag = np.array([np.sqrt(b_gse_vec[i][0]**2 + b_gse_vec[i][1]**2 + b_gse_vec[i][2]**2) for i in range(B_epoch.size)])

######################### ESA measurements #########################
print('Getting ESA (velocity, density, temperature) data...')
pyspedas.themis.esa(probe=probe, trange=time_range, time_clip=True, varnames=['th'+probe+'_peir_velocity_gse','th'+probe+'_peif_velocity_gse','th'+probe+'_peir_density','th'+probe+'_peer_density','th'+probe+'_peir_avgtemp','th'+probe+'_peer_avgtemp'])
[V_epoch, Np] = get_data('th'+probe+'_peir_density')                    # cm^-3
[V_epoch, Npe] = get_data('th'+probe+'_peer_density')                   # cm^-3
# Proton and electron temperature
[T_epoch, Tp] = get_data('th'+probe+'_peir_avgtemp')                    # eV
Tp *= 11604.518                                                         # eV --> K
[e_epoch, Te] = get_data('th'+probe+'_peer_avgtemp')                    # eV
Te *= 11604.518                                                         # eV --> K

[V_epoch, u, components] = get_data('th'+probe+'_peir_velocity_gse')    # km/s
if nan_checker(u[:,0]) > 0.05:
    [V_epoch, u, components] = get_data('th'+probe+'_peif_velocity_gse')# km/s
Vmag = np.array([np.sqrt(u[i][0]**2 + u[i][1]**2 + u[i][2]**2) for i in range(V_epoch.size)])
print(u)

####################################################################
# Trim data to specified range
####################################################################
print('Trimming data to specified time range...')
BTime = np.array([datetime.utcfromtimestamp(t) for t in B_epoch])
selected_index = [(BTime> datetimeStart) & (BTime < datetimeEnd)][0]
BTime = BTime[selected_index]
BGSE = b_gse_vec[selected_index]
Bmag = Bmag[selected_index]

TTime = np.array([datetime.utcfromtimestamp(t) for t in T_epoch])
selected_index = [(TTime > datetimeStart) & (TTime < datetimeEnd)][0]
Np = Np[selected_index]
Tp = Tp[selected_index]
if 'Te' in locals():
    eTime = np.array([datetime.utcfromtimestamp(t) for t in e_epoch])
    selected_index = [(eTime > datetimeStart) & (eTime < datetimeEnd)][0]
    Te = Te[selected_index]
    eTime = eTime[selected_index]

VTime = np.array([datetime.utcfromtimestamp(t) for t in V_epoch])
selected_index = [(VTime > datetimeStart) & (VTime < datetimeEnd)][0]
VTime = VTime[selected_index]
VGSE = u[selected_index]
Vmag = Vmag[selected_index]

# Time cadence
dt = datetime.utcfromtimestamp(B_epoch[1])-datetime.utcfromtimestamp(B_epoch[0])
dt = round(dt.seconds + dt.microseconds*1e-6,4)


####################################################################
# Put data into DataFrame
####################################################################
print('Putting BGSE, VGSE, Np, Tp (Te) into DataFrame...')
BGSE_DataFrame = pd.DataFrame(np.column_stack((BGSE,Bmag)), index = BTime, columns = ['Bx', 'By', 'Bz', 'Bmag']).resample(f'{dt}S').mean()
VGSE_DataFrame = pd.DataFrame(np.column_stack((VGSE,Vmag)), index = VTime, columns = ['Vx', 'Vy', 'Vz', 'Vmag']).resample(f'{dt}S').mean()
Np_DataFrame = pd.DataFrame(Np, index = TTime, columns = ['Np']).resample(f'{dt}S').mean()
Tp_DataFrame = pd.DataFrame(Tp, index = TTime, columns = ['Tp']).resample(f'{dt}S').mean()
Te_DataFrame = pd.DataFrame(Te, index = eTime, columns = ['Te']).resample(f'{dt}S').mean()

# Alfven speed & magnitude, Thermal & Magnetic Pressure, Plasma Beta
print('Putting bGSE, beta into DataFrame...')
bGSE_DataFrame = (1e-9)*BGSE_DataFrame.div(np.sqrt(4*np.pi*Np_DataFrame['Np']*(1.67e-27)*(1e-1)),axis=0).rename(columns={'Bx':'bx','By':'by','Bz':'bz','Bmag':'bmag'})*(1e-3)
if 'Te' in locals():
    p_DataFrame = Tp_DataFrame.add(Te_DataFrame['Te'],axis=0).mul(Np_DataFrame['Np']*(1.3806e-17),axis=0).rename(columns={'Tp':'p'})
else:
    p_DataFrame = Tp_DataFrame.mul(Np_DataFrame['Np']*(1.3806e-17),axis=0).rename(columns={'Tp':'p'})
Bp = ((Bmag*(1e-9))**2)/2/(4*np.pi*(1e-7))
beta_DataFrame = p_DataFrame.div(((BGSE_DataFrame['Bmag']*(1e-9))**2)/2/(4*np.pi*(1e-7)),axis=0).rename(columns={'p':'beta'})


####################################################################
# Get spacecraft location
####################################################################
print('Getting ephemeris data...')
pyspedas.themis.state(probe=probe, trange=time_range, time_clip=True, varnames=['th'+probe+'_pos_gse'])
[eph_epoch, pos, components] = get_data('th'+probe+'_pos_gse') #km
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
if 'Te' in locals():
    DataFrame = pd.concat([DataFrame, BGSE_DataFrame, VGSE_DataFrame, Np_DataFrame, Tp_DataFrame, Te_DataFrame, bGSE_DataFrame, beta_DataFrame, xyz_GSE], axis=1)
else:
    DataFrame = pd.concat([DataFrame, BGSE_DataFrame, VGSE_DataFrame, Np_DataFrame, Tp_DataFrame,bGSE_DataFrame, beta_DataFrame, xyz_GSE], axis=1)



# Drop duplicated records and sort data by time index
print('Dropping duplicated records and sorting by time index...')
DataFrame.interpolate(method='linear',inplace=True)
DataFrame.drop_duplicates(keep='first', inplace=True)
DataFrame.sort_index(axis=0, ascending=True, inplace=True, kind='quicksort')
#DataFrame = DataFrame.where(np.abs(DataFrame['Bx']) < 100).where(np.abs(DataFrame['By']) < 100).where(np.abs(DataFrame['Bz']) < 100)
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
