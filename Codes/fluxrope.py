from __future__ import division # Treat integer as float.
import os
import sys
import math
import time
import pickle
import astropy
import multiprocessing
import pandas as pd
import numpy as np
from numpy import linalg as la
import scipy as sp
from scipy import integrate,stats
from scipy.signal import savgol_filter # Savitzky-Golay filter
from datetime import datetime, timedelta
from collections import OrderedDict
import matplotlib
import matplotlib.dates as mdates
import matplotlib.cbook as cbook
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MaxNLocator, FuncFormatter
from matplotlib.lines import Line2D
from PIL import Image
import warnings
warnings.filterwarnings(action="ignore")
#############################################################################################################

# Physics constants.
global mu0              # (N/A^2) magnetic constant permeability of free space vacuum permeability
global m_proton         # Proton mass [kg]
global factor_deg2rad   # Convert degree to radians
global k_Boltzmann      # Boltzmann constant, in J/K.

mu0 = 4.0 * np.pi * 1e-7
m_proton = 1.6726219e-27 # kg
factor_deg2rad = np.pi/180.0 # radians
k_Boltzmann = 1.3806488e-23 # J/K.
#############################################################################################################

# Calculate the velocity of deHoffmann-Teller frame, VHT.
def findVHT(B_DF_inGSE, Vsw_DF_inGSE):
    N = len(B_DF_inGSE)
    B_square = np.square(B_DF_inGSE).sum(axis=1) # Take squre and sum row (axis=1 for row, axis=0 for column)
    KN = np.zeros((N,3,3))                       # np.zeros((layer, row, column)). Right most index change first.
    for n in range(N):
        for i in range(3):
            for j in range(3):
                if i == j:
                    KN[n,i,j] = B_square.iloc[n] - B_DF_inGSE.iloc[n][i] * B_DF_inGSE.iloc[n][j]
                else:
                    KN[n,i,j] = - B_DF_inGSE.iloc[n][i] * B_DF_inGSE.iloc[n][j]
    K = np.mean(KN, axis=0)                      # Take average on layer (axis=1 for row, axis=2 for column).
    KVN = np.zeros((N,3))                        # np.zeros((row, column)). Right most index change first.
    for n in range(N):
        KVN[n,:] = np.dot(KN[n,:], Vsw_DF_inGSE.iloc[n])
    
    # Average KVN over N to get KV.
    KV = np.mean(KVN, axis=0) # Take average on column.
    VHT = np.dot(np.linalg.inv(K), KV)
    return VHT
    
# Check VHT: find correlation coefficient between E = -v X B and EHT = -VHT X B
def check_VHT(VHT, V, B):
    VHT_array = np.array(VHT)
    V_array = np.array(V)
    B_array = np.array(B)               # EHT = B X VHT
    EHT = np.cross(B_array, VHT_array)  # E = B X v
    E = np.cross(B_array, V_array)
    EHT_1D = EHT.reshape(EHT.size)      # merge all component to 1-D array.
    E_1D = E.reshape(E.size)
    mask = ~np.isnan(EHT_1D) & ~np.isnan(E_1D)
    if mask.sum()>=5:
        # slope, intercept, r_value, p_value, std_err = stats.linregress(A,B)
        # scipy.stats.linregress(x, y=None). Put VA on X-axis, V_remaining on Y-axis.
        slope, intercept, r_value, p_value, std_err = stats.linregress(EHT_1D[mask], E_1D[mask])
        return slope, intercept, r_value
    else:
        return np.nan, np.nan, np.nan

#############################################################################################################

# Calculate the eignvectors and eigenvaluse of input matrix dataframe. This module is Python style.
def eigenMatrix(matrix_DataFrame, **kwargs):
    # Calculate the eigenvalues and eigenvectors of covariance matrix.
    eigenValue, eigenVector = la.eig(matrix_DataFrame)
    
    # Sort the eigenvalues and arrange eigenvectors by sorted eigenvalues.
    eigenValue_i = np.argsort(eigenValue) # covM_B_eigenValue_i is sorted index of covM_B_eigenValue
    lambda3 = eigenValue[eigenValue_i[0]] # lambda3, minimum variance
    lambda2 = eigenValue[eigenValue_i[1]] # lambda2, intermediate variance.
    lambda1 = eigenValue[eigenValue_i[2]] # lambda1, maximum variance.
    eigenVector3 = pd.DataFrame(eigenVector[:, eigenValue_i[0]], columns=['minVar(lambda3)']) # Eigenvector 3, along minimum variance
    eigenVector2 = pd.DataFrame(eigenVector[:, eigenValue_i[1]], columns=['interVar(lambda2)']) # Eigenvector 2, along intermediate variance.
    eigenVector1 = pd.DataFrame(eigenVector[:, eigenValue_i[2]], columns=['maxVar(lambda1)']) # Eigenvector 1, along maximum variance.
    
    if kwargs['formXYZ']==True:
        # Form an eigenMatrix with the columns: X = minimum variance direction, Y = Maximum variance direction, Z = intermediate variance direction.
        eigenMatrix = pd.concat([eigenVector3, eigenVector1, eigenVector2], axis=1)
        eigenValues = pd.DataFrame([lambda3, lambda1, lambda2], index=['X1(min)', 'X2(max)', 'X3(inter)'], columns=['eigenValue'])
    else:
        # Form a sorted eigenMatrix using three sorted eigenvectors. Columns are eigenvectors.
        eigenMatrix = pd.concat([eigenVector3, eigenVector2, eigenVector1], axis=1)
        eigenValues = pd.DataFrame([lambda3, lambda2, lambda1], index=['lambda3', 'lambda2', 'lambda1'], columns=['eigenValue'])
    
    eigenVectorMaxVar_lambda1 = (eigenVector[:, eigenValue_i[2]])
    eigenVectorInterVar_lambda2 = (eigenVector[:, eigenValue_i[1]])
    eigenVectorMinVar_lambda3 = (eigenVector[:, eigenValue_i[0]])
    # print('eigenVectorMaxVar_lambda1 = ', eigenVectorMaxVar_lambda1) # maxVar(lambda1)
    # print('eigenVectorInterVar_lambda2 = ', eigenVectorInterVar_lambda2) # interVar(lambda2)
    # print('eigenVectorMinVar_lambda3 = ', eigenVectorMinVar_lambda3) # minVar(lambda3)

    # return eigenValues, eigenMatrix
    return lambda1, lambda2, lambda3, eigenVectorMaxVar_lambda1, eigenVectorInterVar_lambda2, eigenVectorMinVar_lambda3

################################################################################################################

# Find X axis according to Z axis and V. The X axis is the projection of V on the plane perpendicular to Z axis.
# For this function, numba is slower than python.
def findXaxis(Z, V):
    Z = np.array(Z)
    V = np.array(V)

    # Both Z and V are unit vector representing the directions.
    z1 = Z[0]; z2 = Z[1]; z3 = Z[2]; v1 = V[0]; v2 = V[1]; v3 = V[2]
    # V, Z, and X must satisfy two conditions that lead to two equations with three unknow. We can solve for x1, x2, and x3, in which x1 is arbitrary. Let x1=1, then normalize X.
        # 1) They are co-planar : (Z cross V) dot X = 0
        # 2) Z perpendicular to X : Z dot X = 0
    x1 = 1.0 # Arbitray.
    x2 = -((x1*(v2*z1*z1 - v1*z1*z2 - v3*z2*z3 + v2*z3*z3))/(v2*z1*z2 - v1*z2*z2 + v3*z1*z3 - v1*z3*z3))
    x3 = -((x1*(v3*z1*z1 + v3*z2*z2 - v1*z1*z3 - v2*z2*z3))/(v2*z1*z2 - v1*z2*z2 + v3*z1*z3 - v1*z3*z3))
    
    # Normalization.
    X = np.array([float(x1), float(x2), float(x3)])
    X = X/(la.norm(X))
    if X.dot(V) < 0:
        X = - X
    return X

# Given two orthnormal vectors(Z and X), find the third vector(Y) to form right-hand side frame.
# For this function, numba is slower than python.
def formRightHandFrame(X, Z): # Z cross X = Y in right hand frame.
    X = np.array(X)
    Z = np.array(Z)
    Y = np.cross(Z, X)
    Y = Y/(la.norm(Y)) # Normalize.
    return Y

################################################################################################################

# Find how many turning points in an array.
def turningPoints(array):
    array = np.array(array)
    dx = np.diff(array)
    dx = dx[dx != 0] # if don't remove duplicate points, will miss the turning points with duplicate values.
    return np.sum(dx[1:] * dx[:-1] < 0)

################################################################################################################

# Usage: B_inFR = B_inGSE.dot(matrix_transToFluxRopeFrame)
def angle2matrix(theta_deg, phi_deg, VHT_inGSE):
    factor_deg2rad = np.pi/180.0 # Convert degree to rad.
    # Direction cosines:
        # x = rcos(alpha) = rsin(theta)cos(phi) => cos(alpha) = sin(theta)cos(phi)
        # y = rcos(beta)  = rsin(theta)sin(phi) => cos(beta)  = sin(theta)sin(phi)
        # z = rcos(gamma) = rcos(theta)         => cos(gamma) = cos(theta)
        
    # Use direction cosines to construct a unit vector.
    theta_rad = factor_deg2rad * theta_deg
    phi_rad   = factor_deg2rad * phi_deg
    
    Z_unitVector = np.array([np.sin(theta_rad)*np.cos(phi_rad), np.sin(theta_rad)*np.sin(phi_rad), np.cos(theta_rad)]) # Form new Z_unitVector according to direction cosines.
    X_unitVector = findXaxis(Z_unitVector, -VHT_inGSE)                                                                 # Find X axis from Z axis and -VHT.
    Y_unitVector = formRightHandFrame(X_unitVector, Z_unitVector)                                                      # Find the Y axis to form a right-handed coordinater with X and Z.
    
    # Project B_inGSE into FluxRope Frame.
    matrix_transToFluxRopeFrame = np.array([X_unitVector, Y_unitVector, Z_unitVector]).T
    return matrix_transToFluxRopeFrame

################################################################################################################

def directionVector2angle(V):
    Z = np.array([0,0,1])
    X = np.array([1,0,0])
    cos_theta = np.dot(V,Z)/la.norm(V)/la.norm(Z)
    V_cast2XY = np.array([V[0], V[1], 0])
    cos_phi = np.dot(V_cast2XY,X)/la.norm(V_cast2XY)/la.norm(X)
    theta_deg = np.arccos(np.clip(cos_theta, -1, 1))*180/np.pi
    phi_deg = np.arccos(np.clip(cos_phi, -1, 1))*180/np.pi
    if V[1]<0:
        phi_deg = 360 - phi_deg
    return (theta_deg, phi_deg)
    
################################################################################################################

def angle2directionVector(theta_deg, phi_deg):
    factor_deg2rad = np.pi/180.0
    theta_rad      = factor_deg2rad * theta_deg
    phi_rad        = factor_deg2rad * phi_deg
    
    # Form new Z_unitVector according to direction cosines.
    Z_unitVector = np.array([np.sin(theta_rad)*np.cos(phi_rad), np.sin(theta_rad)*np.sin(phi_rad), np.cos(theta_rad)])
    return Z_unitVector
    
################################################################################################################

def flip_direction(theta_deg, phi_deg):
    new_theta_deg = 180 - theta_deg
    new_phi_deg = phi_deg + 180
    if new_phi_deg >= 360:
        new_phi_deg -= 360
    return (new_theta_deg, new_phi_deg)

################################################################################################################

# Walen test. Find the correlation coefficient and slop between the remainning velocity and Alfven speed. This function return the component-by-component correlation coefficient and slope of the plasma velocities and the Alfven velocities.
def walenTest(VA, V_remaining):
    V_remaining = np.array(V_remaining)                     # V_remaining reshaped time series of solar wind velocity [km/s].
    VA = np.array(VA)                                       # VA is the reshaped time series of Alfven wave [km/s].
    mask = ~np.isnan(VA) & ~np.isnan(V_remaining)
    if mask.sum()>=5:
        # slope, intercept, r_value, p_value, std_err = stats.linregress(A,B)
        # scipy.stats.linregress(x, y=None). Put VA on X-axis, V_remaining on Y-axis.
        slope, intercept, r_value, p_value, std_err = stats.linregress(VA[mask], V_remaining[mask])
        return slope, intercept, r_value
    else:
        return np.nan, np.nan, np.nan

################################################################################################################

# Loop for all directions to calculate residue, return the smallest residue and corresponding direction.
def searchFluxRopeInWindow(B_DataFrame, VHT, n_theta_grid, minDuration, dt, flag_smoothA, Np_DataFrame, Tp_DataFrame, Vsw_DataFrame):
    print('{} - [{}~{} minutes] searching: ({} ~ {})'.format(time.ctime(), minDuration, len(B_DataFrame), B_DataFrame.index[0], B_DataFrame.index[-1]))
    
    # Initialization. Caution: the type of return value will be different if the initial data is updated. If updated, timeRange_temp will become to tuple, plotData_dict_temp will becomes to dict, et, al.
    time_start_temp     = np.nan
    time_end_temp       = np.nan
    time_turn_temp      = np.nan
    turnPointOnTop_temp = np.nan
    Residue_diff_temp   = np.inf
    Residue_fit_temp    = np.inf
    duration_temp       = np.nan
    theta_temp          = np.nan
    phi_temp            = np.nan
    Tp_in_temp          = np.nan
    Np_in_temp          = np.nan
    time_start, time_end, time_turn, turnPointOnTop, Residue_diff, Residue_fit, duration = getResidueForCurrentAxial(0, 0, minDuration, B_DataFrame, VHT, dt, flag_smoothA, Np_DataFrame, Tp_DataFrame, Vsw_DataFrame)
    print('For current orientation, the returned residue is {}'.format(Residue))
    print('For current orientation, the returned duration is {}'.format(duration))
    if  Residue_diff < Residue_diff_temp:
        time_start_temp     = time_start
        time_end_temp       = time_end
        time_turn_temp      = time_turn
        turnPointOnTop_temp = turnPointOnTop
        Residue_diff_temp   = Residue_diff
        Residue_fit_temp    = Residue_fit
        theta_temp          = 0
        phi_temp            = 0
        duration_temp       = duration
        Tp_in_temp          = Tp_in
        Np_in_temp          = Np_in
    
    # This step loops all theta and phi except for theta = 0.
    thetaArray = np.linspace(0, 90, n_theta_grid+1)
    thetaArray = thetaArray[1:]
    phiArray = np.linspace(0, 360, n_theta_grid*2+1)
    phiArray = phiArray[1:]
    for theta_deg in thetaArray: # Not include theta = 0.
        for phi_deg in phiArray: # Include phi = 0.
            time_start, time_end, time_turn, turnPointOnTop, Residue_diff, Residue_fit, duration = getResidueForCurrentAxial(theta_deg, phi_deg, minDuration, B_DataFrame, VHT, dt, flag_smoothA, Np_DataFrame, Tp_DataFrame, Vsw_DataFrame)
            if Residue_diff < Residue_diff_temp:
                time_start_temp     = time_start
                time_end_temp       = time_end
                time_turn_temp      = time_turn
                turnPointOnTop_temp = turnPointOnTop
                Residue_diff_temp   = Residue_diff
                Residue_fit_temp    = Residue_fit
                theta_temp          = theta_deg
                phi_temp            = phi_deg
                duration_temp       = duration
               Tp_in_temp           = Tp_in
               Np_in_temp           = Np_in

    print('Residue_diff = {}'.format(Residue_diff_temp))
    print('Residue_fit  = {}\n'.format(Residue_fit_temp))
    return time_start_temp, time_turn_temp, time_end_temp, duration_temp, turnPointOnTop_temp, Residue_diff_temp, Residue_fit_temp, (theta_temp, phi_temp), (round(VHT[0],5),round(VHT[1],5),round(VHT[2],5))

################################################################################################################

# Calculate the residue for given theta and phi.
def getResidueForCurrentAxial(theta_deg, phi_deg, minDuration, B_DataFrame, VHT, dt, flag_smoothA, Np_DataFrame, Tp_DataFrame, Vsw_DataFrame):

    # Initialize
    time_start = np.nan
    time_end = np.nan
    time_turn = np.nan
    Residue_diff = np.inf
    Residue_fit = np.inf
    duration = np.nan
    turnPointOnTop = np.nan
    
    # Loop for half polar angle (theta(0~90 degree)), and azimuthal angle (phi(0~360 degree)) for Z axis orientations.
    # Direction cosines:
        # x = rcos(alpha) = rsin(theta)cos(phi) => cos(alpha) = sin(theta)cos(phi)
        # y = rcos(beta)  = rsin(theta)sin(phi) => cos(beta)  = sin(theta)sin(phi)
        # z = rcos(gamma) = rcos(theta)         => cos(gamma) = cos(theta)
    
    # Using direction cosines to form a unit vector.
    theta_rad = factor_deg2rad * theta_deg
    phi_rad   = factor_deg2rad * phi_deg
    
    Z_unitVector = np.array([np.sin(theta_rad)*np.cos(phi_rad), np.sin(theta_rad)*np.sin(phi_rad), np.cos(theta_rad)]) # Form new Z_unitVector according to direction cosines.
    X_unitVector = findXaxis(Z_unitVector, -VHT)                                                                       # Find X axis from Z axis and -VHT.
    Y_unitVector = formRightHandFrame(X_unitVector, Z_unitVector)                                                      # Find the Y axis to form a right-handed coordinater with X and Z.

    # Project B_DataFrame & VHT into new trial Frame.
    transToTrialFrame = np.array([X_unitVector, Y_unitVector, Z_unitVector]).T
    B_inTrialframe_DataFrame = B_DataFrame.dot(transToTrialFrame)
    VHT_inTrialframe = VHT.dot(transToTrialFrame)

    ########################################
    Vsw_inTrialframe = Vsw_DataFrame.dot(transToTrialFrame)
    V_remaining = np.array(Vsw_inTrialframe - VHT_inTrialframe)
    P_massDensity = Np_DataFrame['Np'] * m_proton * 1e6 # [kg/m^3]
    len_P_massDensity = len(P_massDensity)
    P_massDensity_array = np.array(P_massDensity)
    P_massDensity_array = np.reshape(P_massDensity_array, (len_P_massDensity, 1))
    VA_inTrialframe = np.array(B_inTrialframe_DataFrame * 1e-9) / np.sqrt(mu0 * P_massDensity_array) / 1000.0
    Mach = np.sqrt(np.square(V_remaining).sum(axis=1))/np.sqrt(np.square(VA_inTrialframe).sum(axis=1))
    a = (Mach.mean())**2     
    B_norm_DF = pd.DataFrame(np.sqrt(np.square(B_DataFrame).sum(axis=1)),columns=['|B|'])

    # Calculate A(x,0) by integrating By. A(x,0) = Integrate[-By(s,0)ds, {s, 0, x}], where ds = -Vht dot unit(X) dt.
    # Don't forget to convert km/s to m/s, and convert nT to T.
    # By = B_inTrialframe_DataFrame[1]
    # A = integrate.cumtrapz(-B_inTrialframe_DataFrame[1]*1e-9, dx=ds, initial=0)
    ds = - VHT_inTrialframe[0] * 1000.0 * dt # Space increment along X axis. Convert km/s to m/s.
    A = integrate.cumtrapz(-(1-a)*B_inTrialframe_DataFrame[1]*1e-9, dx=ds, initial=0)

    # Calculate Pt(x,0). Pt(x,0)=p(x,0)+Bz^2/(2mu0)
    # By = B_inTrialframe_DataFrame[2]
    Pp1 = np.array(Np_DataFrame['Np']) * 1e6 * 1.3806488e-23 * 1e9 * np.array(Tp_DataFrame['Tp'])
    Pb1 = np.array((B_inTrialframe_DataFrame[2] * 1e-9)**2 / (2.0*mu0) * 1e9)   # 1e9 convert unit form pa to npa. 
    PB1 = np.array((B_norm_DF['|B|'] * 1e-9)**2 / (2.0*mu0) * 1e9)              # 1e9 convert unit form pa to npa. 
    Pt = ((1-a)**2)*Pb1 + (1-a)*Pp1 + (a*(1-a))*PB1
    #######################################

    # Check how many turning points in original data.
    num_A_turningPoints = turningPoints(A)
    print('num_A_turningPoints = {}'.format(num_A_turningPoints))
    
    # Smooth A series to find the number of the turning point in main trend.
    if flag_smoothA == True:
        # Because the small scale turning points are not important, we only use A_smoothed to find turning points. When fitting Pt with A, original A is used.
        # Firstly, downsample A to 20 points, then apply savgol_filter, then upsample to original data points number.
        index_A = range(len(A))
        index_downsample = np.linspace(index_A[0],index_A[-1], 20)
        A_downsample = np.interp(index_downsample, index_A, A)          # Downsample A to 20 points.
        A_downsample = savgol_filter(A_downsample, 7, 3)                # Apply savgol_filter. # 7 is smooth window size, 3 is polynomial order.
        A_upsample = np.interp(index_A, index_downsample, A_downsample) # Upsample to original data points amount.

        # The smoothed A is just upsampled A.
        A_smoothed = A_upsample
    else:
        A_smoothed = A
    

    # Check how many turning points in smoothed data.
    num_A_smoothed_turningPoints = turningPoints(A_smoothed)
    print('num_A_smoothed_turningPoints = {}'.format(num_A_smoothed_turningPoints))

    # num_A_smoothed_turningPoints==0 means the A value is not double folded. It's monotonous. Skip.
    # num_A_smoothed_turningPoints > 1 means the A valuse is 3 or higher folded. Skip.
    if (num_A_smoothed_turningPoints==0)|(num_A_smoothed_turningPoints>1):
        return time_start, time_end, time_turn, turnPointOnTop, Residue_diff, Residue_fit, duration # Skip the rest commands in current iteration.
    print('Theta={}, Phi={}. Double-folding feature detected!\n'.format(theta_deg, phi_deg))
    
    # Find the boundary of A.
    A_smoothed_start = A_smoothed[0]                            # The first value of A.
    A_smoothed_end = A_smoothed[-1]                             # The last value of A.
    A_smoothed_max_index = A_smoothed.argmax()                  # The index of max A, return the index of first max(A).
    A_smoothed_max = A_smoothed[A_smoothed_max_index]           # The maximum A.
    A_smoothed_min_index = A_smoothed.argmin()                  # The index of min A, return the index of first min(A).
    A_smoothed_min = A_smoothed[A_smoothed_min_index]           # The minimum A.

    if (A_smoothed_min == min(A_smoothed_start, A_smoothed_end))&(A_smoothed_max == max(A_smoothed_start, A_smoothed_end)):
        # This means the A value is not double folded. It's monotonous. Skip.
        # Sometimes num_A_smoothed_turningPoints == 0 does not work well. This is double check.
        return time_start, time_end, time_turn, turnPointOnTop, Residue_diff, Residue_fit, duration
    elif abs(A_smoothed_min - ((A_smoothed_start + A_smoothed_end)/2)) < abs(A_smoothed_max - ((A_smoothed_start + A_smoothed_end)/2)):
        # This means the turning point is on the right side.
        A_turnPoint_index = A_smoothed_max_index
        turnPointOnRight = True
    elif abs(A_smoothed_min - ((A_smoothed_start + A_smoothed_end)/2)) > abs(A_smoothed_max - ((A_smoothed_start + A_smoothed_end)/2)):
        # This means the turning point is on the left side.
        A_turnPoint_index = A_smoothed_min_index
        turnPointOnLeft = True

    # Split A into two subarray from turning point.
    A_sub1  = A[:A_turnPoint_index+1]
    Pt_sub1 = Pt[:A_turnPoint_index+1]  # Pick corresponding Pt according to index of A.
    A_sub2  = A[A_turnPoint_index:]
    Pt_sub2 = Pt[A_turnPoint_index:]    # Pick corresponding Pt according to index of A.

    # Get time stamps and split into two subarray from turning point.
    timeStamp      = B_inTrialframe_DataFrame.index
    timeStamp_sub1 = timeStamp[:A_turnPoint_index+1]
    timeStamp_sub2 = timeStamp[A_turnPoint_index:]

    # Keep the time of turn point and the value of Pt turn point.
    Pt_turnPoint        = Pt[A_turnPoint_index]
    timeStamp_turnPoint = timeStamp[A_turnPoint_index]
    
    # Put two branches into DataFrame.
    Pt_vs_A_sub1_DataFrame = pd.DataFrame({'Pt_sub1':np.array(Pt_sub1).T,'timeStamp_sub1':np.array(timeStamp_sub1).T}, index=A_sub1)
    Pt_vs_A_sub2_DataFrame = pd.DataFrame({'Pt_sub2':np.array(Pt_sub2).T,'timeStamp_sub2':np.array(timeStamp_sub2).T}, index=A_sub2)

    # Sort by A. A is index in Pt_vs_A_sub1_DataFrame.
    Pt_vs_A_sub1_DataFrame.sort_index(ascending=True, inplace=True, kind='quicksort')
    Pt_vs_A_sub2_DataFrame.sort_index(ascending=True, inplace=True, kind='quicksort')
    
    # Trim two branches to get same boundary A value.
    # Note that, triming is by A value, not by lenght. After trimming, two branches may have different lengths.
    A_sub1_boundary_left  = Pt_vs_A_sub1_DataFrame.index.min()
    A_sub1_boundary_right = Pt_vs_A_sub1_DataFrame.index.max()
    A_sub2_boundary_left  = Pt_vs_A_sub2_DataFrame.index.min()
    A_sub2_boundary_right = Pt_vs_A_sub2_DataFrame.index.max()

    A_boundary_left  = max(A_sub1_boundary_left, A_sub2_boundary_left)
    A_boundary_right = min(A_sub1_boundary_right, A_sub2_boundary_right)

    Pt_vs_A_sub1_trimmed_DataFrame = Pt_vs_A_sub1_DataFrame.iloc[Pt_vs_A_sub1_DataFrame.index.get_loc(A_boundary_left,method='nearest'):Pt_vs_A_sub1_DataFrame.index.get_loc(A_boundary_right,method='nearest')+1]
    Pt_vs_A_sub2_trimmed_DataFrame = Pt_vs_A_sub2_DataFrame.iloc[Pt_vs_A_sub2_DataFrame.index.get_loc(A_boundary_left,method='nearest'):Pt_vs_A_sub2_DataFrame.index.get_loc(A_boundary_right,method='nearest')+1]

    # Get the time range of trimmed A.
    timeStamp_start = min(Pt_vs_A_sub1_trimmed_DataFrame['timeStamp_sub1'].min(skipna=True), Pt_vs_A_sub2_trimmed_DataFrame['timeStamp_sub2'].min(skipna=True))
    timeStamp_end   = max(Pt_vs_A_sub1_trimmed_DataFrame['timeStamp_sub1'].max(skipna=True), Pt_vs_A_sub2_trimmed_DataFrame['timeStamp_sub2'].max(skipna=True))
    time_start      = int(timeStamp_start.strftime('%Y%m%d%H%M%S'))
    time_end        = int(timeStamp_end.strftime('%Y%m%d%H%M%S'))
    time_turn       = int(timeStamp_turnPoint.strftime('%Y%m%d%H%M%S'))
    duration        = int((timeStamp_end - timeStamp_start).total_seconds()/sys.argv[2])+1 #dt=sys.argv[2] for detection codes
    
    # Skip if shorter than minDuration.
    if duration < minDuration:
        return time_start, time_end, time_turn, turnPointOnTop, Residue_diff, Residue_fit, duration

    # Calculate two residues. respectively: Residue_fit and Residue_diff.
    # Preparing for calculating Residue_fit, the residue of all data sample w.r.t. fitted PtA curve.
    A_sub1_array  = np.array(Pt_vs_A_sub1_trimmed_DataFrame.index)
    A_sub2_array  = np.array(Pt_vs_A_sub2_trimmed_DataFrame.index)
    Pt_sub1_array = np.array(Pt_vs_A_sub1_trimmed_DataFrame['Pt_sub1'])
    Pt_sub2_array = np.array(Pt_vs_A_sub2_trimmed_DataFrame['Pt_sub2'])
    
    # Combine two trimmed branches.
    Pt_array = np.concatenate((Pt_sub1_array, Pt_sub2_array))
    A_array  = np.concatenate((A_sub1_array, A_sub2_array))
   
    # The order must be in accordance.
    sortedIndex     = np.argsort(A_array)
    A_sorted_array  = A_array[sortedIndex]
    Pt_sorted_array = Pt_array[sortedIndex]

    # Fit a polynomial function (3rd order) and use it to calculate residue.
    Pt_array_float = Pt_array.astype(np.float64)
    Pt_A_coeff     = np.polyfit(A_array, Pt_array_float, 3)
    Pt_A           = np.poly1d(Pt_A_coeff)

    # Preparing for calculating Residue_diff, the residue get by compare two branches.
    Pt_vs_A_trimmed_DataFrame = pd.concat([Pt_vs_A_sub1_trimmed_DataFrame, Pt_vs_A_sub2_trimmed_DataFrame], axis=1) # Merge two subset into one DataFrame.
    Pt_vs_A_trimmed_DataFrame.drop(['timeStamp_sub1', 'timeStamp_sub2'], axis=1, inplace=True)                      # axis=1 for column.

    # Interpolation. "TypeError: Cannot interpolate with all NaNs" can occur if the DataFrame contains columns of object dtype. Convert data to numeric type and check data type by print(Pt_vs_A_trimmed_DataFrame.dtypes).
    for one_column in Pt_vs_A_trimmed_DataFrame:
        Pt_vs_A_trimmed_DataFrame[one_column] = pd.to_numeric(Pt_vs_A_trimmed_DataFrame[one_column], errors='coerce')
    
    # Interpolate according to index A.
    Pt_vs_A_trimmed_DataFrame.interpolate(method='index', axis=0, inplace=True) # axis=0:fill column-by-column
    
    # Drop leading and trailing NaNs. The leading NaN won't be filled by linear interpolation, however, the trailing NaN will be filled by forward copy of the last non-NaN values. So, for leading NaN, just use pd.dropna, and for trailing NaN, remove the duplicated values.
    Pt_vs_A_trimmed_DataFrame.dropna(inplace=True)                                                                  # Drop leading NaNs.
    trailing_NaN_mask_DataFrame = (Pt_vs_A_trimmed_DataFrame.diff()!=0)                                             # Get duplicate bool mask.
    trailing_NaN_mask = np.array(trailing_NaN_mask_DataFrame['Pt_sub1'] & trailing_NaN_mask_DataFrame['Pt_sub2'])   # Drop trailing NaNs.
    Pt_vs_A_trimmed_DataFrame = Pt_vs_A_trimmed_DataFrame.iloc[trailing_NaN_mask]

    # Get Pt_max and Pt_min. They will be used to normalize Residue for both Residue_fit and Residue_diff.
    Pt_max = Pt_sorted_array.max()
    Pt_min = Pt_sorted_array.min()
    Pt_max_min_diff = abs(Pt_max - Pt_min)
    
    # Check if turn point is on top.
    turnPointOnTop = (Pt_turnPoint>(Pt_max-(Pt_max-Pt_min)*0.15)) 
    
    # Use two different defination to calculate Residues. Note that, the definition of Residue_diff is different with Hu's paper. We divided it by 2 two make it comparable with Residue_fit. The definition of Residue_fit is same as that in Hu2004.
    if Pt_max_min_diff == 0:
        Residue_diff = np.inf
        Residue_fit  = np.inf
    else:
        Residue_diff = 0.5 * np.sqrt((1.0/len(Pt_vs_A_trimmed_DataFrame))*((Pt_vs_A_trimmed_DataFrame['Pt_sub1'] - Pt_vs_A_trimmed_DataFrame['Pt_sub2']) ** 2).sum()) / Pt_max_min_diff
        Residue_fit  = np.sqrt((1.0/len(A_array))*((Pt_sorted_array - Pt_A(A_sorted_array)) ** 2).sum()) / Pt_max_min_diff
        Residue_diff = round(Residue_diff, 5)
        Residue_fit  = round(Residue_fit, 5)
    
    return time_start, time_end, time_turn, turnPointOnTop, Residue_diff, Residue_fit, duration
        
################################################################################################################  
        
# Calculate r_VHT (correlation coefficient of EHT and E).
def calculate_r_VHT(data_DF, fluxRopeList_DF):
    fluxRopeList_with_r_VHT_DF = fluxRopeList_DF.copy()                                                         # Copy fluxRopeList_DF.
    fluxRopeList_with_r_VHT_DF = fluxRopeList_with_r_VHT_DF.assign(r_VHT=[0.0]*len(fluxRopeList_with_r_VHT_DF)) # Add 'r_VHT' column to fluxRopeList_DF.
    recordLength = len(fluxRopeList_with_r_VHT_DF)                                                              # Loop each record to calculate cv.
    for index_temp, oneRecord_temp in fluxRopeList_with_r_VHT_DF.iterrows():
        print('\nindex_temp = {}/{}'.format(index_temp, recordLength))
        startTime_temp = oneRecord_temp['startTime']
        endTime_temp = oneRecord_temp['endTime']
        VHT_temp = oneRecord_temp['VHT']
        
        # Grab the data for one FR_record.
        selectedRange_mask  = (data_DF.index >= oneRecord_temp['startTime']) & (data_DF.index <= oneRecord_temp['endTime'])
        oneRecord_data_temp = data_DF.iloc[selectedRange_mask]
        oneRecord_Vsw_temp  = oneRecord_data_temp[['Vx', 'Vy', 'Vz']].copy()
        Vsw_temp            = np.sqrt(np.square(oneRecord_Vsw_temp).sum(axis=1))
        oneRecord_B_temp    = oneRecord_data_temp[['Bx', 'By', 'Bz']].copy()
        
        # Calculate r_VHT and add it to fluxRopeList_with_r_VHT_DF.
        slope, intercept, r_VHT_temp = check_VHT(VHT_temp, oneRecord_Vsw_temp, oneRecord_B_temp)
        fluxRopeList_with_r_VHT_DF.loc[index_temp, 'r_VHT'] = r_VHT_temp
        
        if 0:
            if (r_VHT_temp<0.98):
                plt.plot(np.array(Vsw_temp))
                plt.ylim([300,700])
                plt.show()
                
    print(fluxRopeList_with_r_VHT_DF)
    fluxRopeList_with_r_VHT_DF['r_VHT'].plot.hist(bins=100)
    plt.show()

################################################################################################################

# Clean up raw_result.
def clean_up_raw_result(data_DF, dataObject_or_dataPath, **kwargs):
    # Check input datatype
    if isinstance(dataObject_or_dataPath, dict):
        print('\nYour input is a dictionary data.')
        search_result_raw = dataObject_or_dataPath
    elif isinstance(dataObject_or_dataPath, str):
        print('\nYour input is a path. Load the dictionary data via this path.')
        search_result_raw = pd.read_pickle(open(dataObject_or_dataPath, 'rb'))
    else:
        print('\nPlease input the correct datatype!')
        return None

    # Set default values.
    dt = float(sys.argv[3]) # seconds
    turnTime_tolerance  = 5*dt
    min_residue_diff    = 0.2
    min_residue_fit     = 0.2

    # Set fitted curve quality parameters.
    max_tailPercentile = 0.3
    max_tailDiff       = 0.3
    max_PtFitStd       = 0.3
    Vsw_std_threshold  = 10000  # Max allowed standard deviation for solar wind speed.
    Vsw_diff_threshold = 10000 # Max allowed solar wind max-min difference.

    # Remove discontinuity.
    isRemoveShock = False
    
    # Walen test.
    walenTest_r_threshold = 0.8 # correlation coefficient.
    walenTest_k_threshold = 0.3 # slope.
   
    # Display control.
    isVerbose = False
    isPrintIntermediateDF = True
    output_filename = 'search_result_no_overlap' + namestr
    output_dir = os.getcwd()
    
    # If keyword is specified, overwrite the default value.
    print('\nSetting parameters:')
    if 'turnTime_tolerance' in kwargs:
        turnTime_tolerance = kwargs['turnTime_tolerance']
        print('turnTime_tolerance is set to {}.'.format(turnTime_tolerance))
    if 'min_residue_diff' in kwargs:
        min_residue_diff = kwargs['min_residue_diff']
        print('min_residue_diff is set to {}.'.format(min_residue_diff))
    if 'min_residue_fit' in kwargs:
        min_residue_fit = kwargs['min_residue_fit']
        print('min_residue_fit is set to {}.'.format(min_residue_fit))
    if 'max_tailPercentile' in kwargs:
        max_tailPercentile = kwargs['max_tailPercentile']
        print('max_tailPercentile is set to {}.'.format(max_tailPercentile))
    if 'max_tailDiff' in kwargs:
        max_tailDiff = kwargs['max_tailDiff']
        print('max_tailDiff is set to {}.'.format(max_tailDiff))
    if 'max_PtFitStd' in kwargs:
        max_PtFitStd = kwargs['max_PtFitStd']
        print('max_PtFitStd is set to {}.'.format(max_PtFitStd))
    if 'isRemoveShock' in kwargs:
        isRemoveShock = kwargs['isRemoveShock']
        print('isRemoveShock set to {}.'.format(isRemoveShock))
        if isRemoveShock:
            if 'shockList_DF' in kwargs:
                shockList_DF = kwargs['shockList_DF']
                print('shockList_DF is loaded.')
            else:
                print('isRemoveShock is True, but shockList_DF is not provided.')
                return None
            if 'spacecraftID' in kwargs:
                spacecraftID = kwargs['spacecraftID']
                print('spacecraftID is set to {}'.format(spacecraftID))
            else:
                print('isRemoveShock is True, but spacecraftID is not provided.')
                return None
    if 'Vsw_std_threshold' in kwargs:
        Vsw_std_threshold = kwargs['Vsw_std_threshold']
        print('Vsw_std_threshold is set to {}.'.format(Vsw_std_threshold))
    if 'Vsw_diff_threshold' in kwargs:
        Vsw_diff_threshold = kwargs['Vsw_diff_threshold']
        print('Vsw_diff_threshold is set to {}.'.format(Vsw_diff_threshold))
    if 'walenTest_r_threshold' in kwargs:
        walenTest_r_threshold = kwargs['walenTest_r_threshold']
        print('walenTest_r_threshold is set to {}.'.format(walenTest_r_threshold))
    if 'walenTest_k_threshold' in kwargs:
        walenTest_k_threshold = kwargs['walenTest_k_threshold']
        print('walenTest_k_threshold is set to {}.'.format(walenTest_k_threshold))
    if 'isVerbose' in kwargs:
        isVerbose = kwargs['isVerbose']
        print('isVerbose is set to {}.'.format(isVerbose))
    if 'isPrintIntermediateDF' in kwargs:
        isPrintIntermediateDF = kwargs['isPrintIntermediateDF']
        print('isPrintIntermediateDF is set to {}.'.format(isPrintIntermediateDF))
    if 'output_dir' in kwargs:
        output_dir = kwargs['output_dir']
        print('output_dir is set to {}.'.format(output_dir))
    if 'output_filename' in kwargs:
        output_filename = kwargs['output_filename']
        print('output_filename is set to {}.'.format(output_filename))
    if 'spacecraftID' in kwargs:
        spacecraftID = kwargs['spacecraftID']
        print('spacecraftID is {}.'.format(spacecraftID))
    if 'shockList' in kwargs:
        shockList = kwargs['shockList']
        print('shockList is loaded.')

    # Check keyword parameters
    print('\nDefault parameters:')
    print('turnTime_tolerance    = {} seconds'.format(turnTime_tolerance))
    print('min_residue_diff      = {}'.format(min_residue_diff))
    print('min_residue_fit       = {}'.format(min_residue_fit))
    print('max_tailPercentile    = {}'.format(max_tailPercentile))
    print('max_tailDiff          = {}'.format(max_tailDiff))
    print('max_PtFitStd          = {}'.format(max_PtFitStd))
    print('isRemoveShock         = {}'.format(isRemoveShock))
    print('Vsw_std_threshold     = {} km/s'.format(Vsw_std_threshold))
    print('Vsw_diff_threshold    = {} km/s'.format(Vsw_diff_threshold))
    print('walenTest_r_threshold = {}'.format(walenTest_r_threshold))
    print('walenTest_k_threshold = {}'.format(walenTest_k_threshold))
    print('isVerbose             = {}'.format(isVerbose))
    print('isPrintIntermediateDF = {}'.format(isPrintIntermediateDF))
    print('output_dir            = {}.'.format(output_dir))
    print('output_filename       = {}.'.format(output_filename))
    
    # Set terminal display format.
    if isPrintIntermediateDF:
        pd.set_option('display.max_rows', 1000)
        pd.set_option('display.max_columns', 500)
        pd.set_option('display.width', 500)

    year_str = output_filename[:4]
    residue_str = f'{min_residue_diff}_{min_residue_fit}'

    # Get duration list.
    duration_list = search_result_raw['true'].keys()
    window_size_list = []
    for item in duration_list:
        window_size = int(item.split('~')[1])
        window_size_list.append(window_size)
    
    # Sort the duration_list with window_size_list by argsort: the duration is in descending order.
    sorted_index_window_size_array = np.argsort(window_size_list)
    sorted_index_window_size_array = sorted_index_window_size_array[::-1]
    duration_array = np.array(list(duration_list))
    duration_list = list(duration_array)
    search_iteration = len(duration_list)

    # Get start and end time.
    datetimeStart = search_result_raw['timeRange']['datetimeStart']
    datetimeEnd = search_result_raw['timeRange']['datetimeEnd']

    # Create empty eventList_no_overlap DataFrame, the cleaned lists will be appended to it.
    eventList_DF_noOverlap = pd.DataFrame(columns=['startTime', 'turnTime', 'endTime', 'duration', 'residue_diff', 'residue_fit', 'theta_phi', 'VHT', 'reduced_residue'])
    # Create empty progressionList for each duration so that the progression of events can be recorded
    progressionList = pd.DataFrame(columns=['beforeSlots','residuals', 'same_turnTime', 'turnTime<5pts']) #, 'walenTest', 'fluctuations', 'turnTime<5min', 'fineFittingCurve', 'slots', 'noOverlap'

    for i_iteration in range(search_iteration):
        print('\n======================================================================')
        print('\niteration = {}/{}'.format(i_iteration+1, search_iteration))

        duration_str_temp = duration_list[i_iteration]
        print('Combining events with {} minutes duration range:'.format(duration_str_temp))
        eventList_temp = search_result_raw['true'][duration_str_temp]
        
        # 1) Check point: If event list not empty, put it into DataFrame. An empty list is itself considered false in true value testing.
        if not eventList_temp:
            print('\nEvent list eventList_temp is empty!')
            print('Go the the next iteration!')
            continue

        # Create headers.
        eventList_temp_Header = ['startTime', 'turnTime', 'endTime', 'duration', 'topTurn', 'residue_diff', 'residue_fit', 'theta_phi', 'VHT']
        eventList_DF_0_original_temp              = pd.DataFrame(eventList_temp, columns=eventList_temp_Header)
        eventList_DF_0_original_temp['startTime'] = pd.to_datetime(eventList_DF_0_original_temp['startTime'], format="%Y%m%d%H%M%S")
        eventList_DF_0_original_temp['turnTime']  = pd.to_datetime(eventList_DF_0_original_temp['turnTime'], format="%Y%m%d%H%M%S")
        eventList_DF_0_original_temp['endTime']   = pd.to_datetime(eventList_DF_0_original_temp['endTime'], format="%Y%m%d%H%M%S")
        
        # Find all records from eventList_DF_0_original_temp that fit the slots of slotList_DF_temp.
        print('\nFitting the events into available slots...')
        if isVerbose:
            print('Before fitting, totoal records is {}'.format(len(eventList_DF_0_original_temp)))
        eventList_DF_2_fineResidue_temp = eventList_DF_0_original_temp.copy()
        # Add keepFlag column to eventList_DF_1_fitSlot_temp.
        eventList_DF_2_fineResidue_temp = eventList_DF_2_fineResidue_temp.assign(keepFlag=[False]*len(eventList_DF_2_fineResidue_temp))

        # Remove the event with residue_diff > min_residue_diff and residue_fit > min_residue_fit.
        print('\nRemoving events with residue_diff > {} and residue_fit > {}...'.format(min_residue_diff, min_residue_fit))
        if isVerbose:
            print('Before Removing, total records is {}.'.format(len(eventList_DF_2_fineResidue_temp)))
        eventList_DF_2_fineResidue_temp = eventList_DF_2_fineResidue_temp[(eventList_DF_2_fineResidue_temp['residue_diff']<=min_residue_diff)&(eventList_DF_2_fineResidue_temp['residue_fit']<=min_residue_fit)]
        # Reset index
        eventList_DF_2_fineResidue_temp.reset_index(drop=True, inplace=True)
        if isVerbose:
            print('After Removing, total records is {}.'.format(len(eventList_DF_2_fineResidue_temp)))
        print('Done.')

        # 2) Check point: After removing bad residue, check if eventList_DF_2_fineResidue_temp is empty. If DataFrame eventList_DF_2_fineResidue_temp is empty. Skip the rest operations.
        if eventList_DF_2_fineResidue_temp.empty:
            print('\nDataFrame eventList_DF_2_fineResidue_temp is empty!')
            print('Go the the next iteration!')
            continue

        # Clean up the records with same turnTime.
        print('\nCombining events with same turnTime...')
        if isVerbose:
            print('Before combining, total records is {}.'.format(len(eventList_DF_2_fineResidue_temp)))
        
        eventList_DF_2_fineResidue_temp = eventList_DF_2_fineResidue_temp.sort_values(by='turnTime') # Sort by turnTime.
        index_max_duration_inGrouped = eventList_DF_2_fineResidue_temp.groupby(['turnTime'], sort=False)['duration'].transform(max) == eventList_DF_2_fineResidue_temp['duration'] # Group by turnTime.
        
        # Pick the event with max duration among the events sharing same turnPoint.
        eventList_DF_3_combinedByTurnTime_temp = eventList_DF_2_fineResidue_temp[index_max_duration_inGrouped]
       
        # Some events with the same duration, same turnTime, but different residue_diff exist.
        index_min_Residue_diff_inGrouped = eventList_DF_3_combinedByTurnTime_temp.groupby(['turnTime'], sort=False)['residue_diff'].transform(min) == eventList_DF_3_combinedByTurnTime_temp['residue_diff'] # Group by turnTime.
        
        # Pick the event with min residue_diff among the events sharing same turnPoint.
        eventList_DF_3_combinedByTurnTime_temp = eventList_DF_3_combinedByTurnTime_temp[index_min_Residue_diff_inGrouped]
        eventList_DF_3_combinedByTurnTime_temp.reset_index(drop=True, inplace=True) # Reset index
        if isVerbose:
            print('After combining, total records is {}.'.format(len(eventList_DF_3_combinedByTurnTime_temp)))
        print('Done.')

        # No need to check whether eventList_DF_3_combinedByTurnTime_temp is empty: if eventList_DF_1_fitSlot_temp is not empty, eventList_DF_3_combinedByTurnTime_temp cannot be empty.
        print('\nRemoving events failed in walen test...')
        # Add an walen test result column: .assign() always returns a copy of the data, leaving the original DataFrame untouched.
        eventList_DF_4_passWalenTest_temp = eventList_DF_3_combinedByTurnTime_temp.copy()
        eventList_DF_4_passWalenTest_temp = eventList_DF_4_passWalenTest_temp.assign(r = len(eventList_DF_4_passWalenTest_temp)*[np.nan])
        eventList_DF_4_passWalenTest_temp = eventList_DF_4_passWalenTest_temp.assign(k = len(eventList_DF_4_passWalenTest_temp)*[np.nan])
        eventList_DF_4_passWalenTest_temp = eventList_DF_4_passWalenTest_temp.assign(MA = len(eventList_DF_4_passWalenTest_temp)*[np.nan])
        
        # Cacluate walen test r (correlation coefficient) value.
        eventList_DF_4_passWalenTest_temp.reset_index(drop=True, inplace=True)
        len_eventList_DF_4_before = len(eventList_DF_4_passWalenTest_temp)
        for index, FR_record in eventList_DF_4_passWalenTest_temp.iterrows():
            if isVerbose:
                print('Walen test: checking duration {} minutes, {}/{}...'.format(duration_str_temp, index+1, len_eventList_DF_4_before))
            theta_deg, phi_deg = FR_record['theta_phi']
            VHT_inGSE = np.array(FR_record['VHT'])
            
            # Grab the data for one fluxrope candidate.
            selectedRange_mask = (data_DF.index >= FR_record['startTime']) & (data_DF.index <= FR_record['endTime'])
            FR_record_data = data_DF.iloc[selectedRange_mask]

            # Interpolate FR_record_data if there is any NaNs in any column.
            FR_record_data_interpolated = FR_record_data.copy(deep=True)
            FR_record_data_interpolated.interpolate(method='time', limit=None, inplace=True)
            FR_record_data_interpolated.bfill(inplace=True)
            FR_record_data_interpolated.ffill(inplace=True)
            
            # Apply Walen test on the result(in FR frame).
            B_inGSE = FR_record_data_interpolated[['Bx', 'By', 'Bz']]
            Vsw_inGSE = FR_record_data_interpolated[['Vx', 'Vy', 'Vz']]
            
            # Project B_inGSE, VHT_inGSE, and Vsw_inFR into FR Frame.
            matrix_transToFluxRopeFrame = angle2matrix(theta_deg, phi_deg, np.array(VHT_inGSE))
            B_inFR   = B_inGSE.dot(matrix_transToFluxRopeFrame)
            VHT_inFR = VHT_inGSE.dot(matrix_transToFluxRopeFrame)
            Vsw_inFR = Vsw_inGSE.dot(matrix_transToFluxRopeFrame)
           
            # Proton mass density. Original Np is in #/cc ( cc = cubic centimeter). Multiply by 1e6 to convert cc to m^3.
            P_massDensity       = FR_record_data['Np'] * m_proton * 1e6 # kg/m^3.
            len_P_massDensity   = len(P_massDensity)
            P_massDensity_array = np.array(P_massDensity)
            P_massDensity_array = np.reshape(P_massDensity_array, (len_P_massDensity, 1))
            
            # Alfven speed. Multiply by 1e-9 to convert nT to T. Divided by 1000.0 to convert m/s to km/s.
            VA              = np.array(B_inFR * 1e-9) / np.sqrt(mu0 * P_massDensity_array) / 1000.0
            VA_1D           = np.reshape(VA, VA.size)
            V_remaining     = np.array(Vsw_inFR - VHT_inFR)
            V_remaining_1D  = np.reshape(V_remaining, V_remaining.size) # VA_inFR
            Mach_average    = np.array((np.sqrt(np.square(V_remaining).sum(axis=1))/np.sqrt(np.square(VA).sum(axis=1))).mean())
            
            # Call Walen test function: First row is x component, second row is y component, third row is z component.
            walenTest_slope, walenTest_intercept, walenTest_r_value = walenTest(VA_1D, V_remaining_1D)
            eventList_DF_4_passWalenTest_temp.loc[index, 'r']  = round(walenTest_r_value, 4) # r.
            eventList_DF_4_passWalenTest_temp.loc[index, 'k']  = round(walenTest_slope, 4) # k.
            eventList_DF_4_passWalenTest_temp.loc[index, 'MA'] = Mach_average

        # Remove the records with |k| > 0.3; except for when |k| > 0.3 & |r| >= 0.8 & MA <= 0.9
        # Keep the records with |k|<=0.3 & |r| >= 0.8 & MA =< 0.9
        print(f'\nKeeping records with |k| > {walenTest_k_threshold} & |r| >= {walenTest_r_threshold} & MA <= 0.9')
        eventList_DF_4_passWalenTest_temp_FRFF = eventList_DF_4_passWalenTest_temp[( (abs(eventList_DF_4_passWalenTest_temp['k']) > walenTest_k_threshold) & ((abs(eventList_DF_4_passWalenTest_temp['r'])>= walenTest_r_threshold) & (eventList_DF_4_passWalenTest_temp['MA'] <= 0.900)))]

        print(f'Keeping records with |k| <= {walenTest_k_threshold}')
        eventList_DF_4_passWalenTest_temp = pd.concat([eventList_DF_4_passWalenTest_temp[(abs(eventList_DF_4_passWalenTest_temp['k']) <= walenTest_k_threshold)],eventList_DF_4_passWalenTest_temp_FRFF], axis=0)
        eventList_DF_4_passWalenTest_temp = eventList_DF_4_passWalenTest_temp.sort_values(by='startTime')
        eventList_DF_4_passWalenTest_temp.reset_index(drop=True, inplace=True)
        len_eventList_DF_4_after = len(eventList_DF_4_passWalenTest_temp)
        if isVerbose:
            print('Before Walen test, total records is {}.'.format(len_eventList_DF_4_before))
            print('After Walen test, total records is {}.'.format(len_eventList_DF_4_after))
        print('Done.')

        # 4) Check point: After Walen test, check if eventList_DF_4_passWalenTest_temp is empty.
        if eventList_DF_4_passWalenTest_temp.empty:
            print('\nDataFrame eventList_DF_4_passWalenTest_temp is empty!')
            print('Go the the next iteration!')
            continue

        ################################################################################################
        # Until now, we still may have duplicated events with same residue_diff but different residue_fit. Keep this in mind when perform furture operations.
        
        # Clean events with less than turnTime_tolerance minutes turnTime difference.
        eventList_DF_6_cleanedTurnTime_temp = eventList_DF_4_passWalenTest_temp.copy() # Default is deep copy.
        
        # Add difference of turn time column. Combine close turnTime.
        eventList_DF_6_cleanedTurnTime_temp = eventList_DF_6_cleanedTurnTime_temp.assign(diff=eventList_DF_6_cleanedTurnTime_temp['turnTime'].diff())
        eventList_DF_6_cleanedTurnTime_temp = eventList_DF_6_cleanedTurnTime_temp.assign(keepFlag=[True]*len(eventList_DF_6_cleanedTurnTime_temp))
        index_column_keepFlag = eventList_DF_6_cleanedTurnTime_temp.columns.get_loc('keepFlag')
        print('\nCombining events with less than {} seconds turnTime difference...'.format(turnTime_tolerance))
        if isVerbose:
            print('Before combining, total records is {}.'.format(len(eventList_DF_6_cleanedTurnTime_temp)))
        i_index = 1
        while(i_index < len(eventList_DF_6_cleanedTurnTime_temp)):
            if(eventList_DF_6_cleanedTurnTime_temp['diff'].iloc[i_index] <= timedelta(seconds=turnTime_tolerance)):
                cluster_begin_temp = i_index - 1
                while((eventList_DF_6_cleanedTurnTime_temp['diff'].iloc[i_index] <= timedelta(seconds=turnTime_tolerance)) ):
                    i_index += 1
                    if (i_index > len(eventList_DF_6_cleanedTurnTime_temp)-1):
                        break
                cluster_end_temp = i_index - 1
                
                # Get minimum residue_diff: .iloc[a:b]=[a,b), .loc[a:b]=[a,b]
                min_residue_diff_index_temp = eventList_DF_6_cleanedTurnTime_temp['residue_diff'].iloc[cluster_begin_temp:cluster_end_temp+1].idxmin()
                
                # Set record with min_residue as true, others as false.
                eventList_DF_6_cleanedTurnTime_temp.iloc[cluster_begin_temp:cluster_end_temp+1, index_column_keepFlag] = False
                eventList_DF_6_cleanedTurnTime_temp.loc[min_residue_diff_index_temp, 'keepFlag'] = True
            else:
                # For .iloc, can only use integer to specify row and column. No column string allowed.
                eventList_DF_6_cleanedTurnTime_temp.iloc[i_index, index_column_keepFlag] = True
                i_index += 1

        # Remove the records labeled as false.
        eventList_DF_6_cleanedTurnTime_temp = eventList_DF_6_cleanedTurnTime_temp[eventList_DF_6_cleanedTurnTime_temp['keepFlag']==True]
        eventList_DF_6_cleanedTurnTime_temp.reset_index(drop=True, inplace=True)
        eventList_DF_6_cleanedTurnTime_temp = eventList_DF_6_cleanedTurnTime_temp.drop('diff', axis=1) # Drop diff.
        eventList_DF_6_cleanedTurnTime_temp.reset_index(drop=True, inplace=True) # Reset index
        if isVerbose:
            print('After combining, total records is {}.'.format(len(eventList_DF_6_cleanedTurnTime_temp)))
        print('Done.')
        eventList_DF_7_fineFittingCurve_temp = eventList_DF_6_cleanedTurnTime_temp.copy()
        
        ###############################################################################################################
        
        # Sort eventList_DF_7_fineFittingCurve_temp by endTime.
        eventList_DF_7_fineFittingCurve_temp = eventList_DF_7_fineFittingCurve_temp.sort_values(by='endTime')
        eventList_DF_7_fineFittingCurve_copy_temp = eventList_DF_7_fineFittingCurve_temp.copy()
        
        # Use the interval scheduling greedy algorithm to remove the overlapes.
        print('\nRemoving overlapped events...')
        eventList_DF_noOverlap = eventList_DF_noOverlap._append(eventList_DF_7_fineFittingCurve_copy_temp, ignore_index=True)
        eventList_DF_noOverlap.sort_values(by='endTime', inplace=True)
        print('Done.')
        
        print('\nAppending records to eventList_DF_noOverlap...')
        print('Done.')
        progressionList = pd.concat([progressionList, pd.DataFrame({'beforeSlots':[len(eventList_DF_0_original_temp)], 'residuals':[len(eventList_DF_2_fineResidue_temp)], 'same_turnTime':[len(eventList_DF_3_combinedByTurnTime_temp)], 'turnTime<5pts':[len(eventList_DF_6_cleanedTurnTime_temp)]},index=[duration_str_temp])])

        if isPrintIntermediateDF:
            print('\neventList_DF_0_original_temp:')
            print(eventList_DF_0_original_temp)
            print('\neventList_DF_2_fineResidue_temp:')
            print(eventList_DF_2_fineResidue_temp)
            print('\neventList_DF_3_combinedByTurnTime_temp:')
            print(eventList_DF_3_combinedByTurnTime_temp)
            print('\neventList_DF_6_cleanedTurnTime_temp:')
            print(eventList_DF_6_cleanedTurnTime_temp)
            print('\neventList_DF_7_fineFittingCurve_temp:')
            print(eventList_DF_7_fineFittingCurve_temp)

    # Reset index.
    eventList_DF_noOverlap.reset_index(drop=True, inplace=True)

    if isPrintIntermediateDF:
        print('\neventList_DF_noOverlap:')
        print(eventList_DF_noOverlap)
        
    # Move on to cleaning up entire list
    print('\n======================================================================')
    print('\nCleaning up entire list')

    # Clean up the records with same turnTime [for all windows].
    print('\nCombining events with same turnTime...')
    if isVerbose:
        print('Before combining, total records is {}.'.format(len(eventList_DF_noOverlap)))
    
    eventList_DF_noOverlap = eventList_DF_noOverlap.sort_values(by='turnTime') # Sort by turnTime.
    index_max_duration_inGrouped = eventList_DF_noOverlap.groupby(['turnTime'], sort=False)['duration'].transform(max) == eventList_DF_noOverlap['duration'] # Group by turnTime.
    
    # Pick the event with max duration among the events sharing same turnPoint.
    eventList_DF_noOverlap = eventList_DF_noOverlap[index_max_duration_inGrouped]
    index_min_Residue_diff_inGrouped =eventList_DF_noOverlap.groupby(['turnTime'], sort=False)['residue_diff'].transform(min) == eventList_DF_noOverlap['residue_diff'] # Group by turnTime
    
    # Pick the event with min residue_diff among the events sharing same turnPoint.
    eventList_DF_noOverlap = eventList_DF_noOverlap[index_min_Residue_diff_inGrouped]
    eventList_DF_noOverlap.reset_index(drop=True, inplace=True) # Reset index
    if isVerbose:
        print('After combining, total records is {}.'.format(len(eventList_DF_noOverlap)))
    print('Done.')

    # Clean events with less than turnTime_tolerance minutes turnTime difference [for all windows].
    eventList_DF_noOverlap_temp = eventList_DF_noOverlap.copy() # Default is deep copy.
    eventList_DF_noOverlap_temp = eventList_DF_noOverlap_temp.assign(diff=eventList_DF_noOverlap_temp['turnTime'].diff())
    eventList_DF_noOverlap_temp = eventList_DF_noOverlap_temp.assign(keepFlag=[True]*len(eventList_DF_noOverlap_temp))
    index_column_keepFlag = eventList_DF_noOverlap_temp.columns.get_loc('keepFlag')
    print('\nCombining events with less than {} seconds turnTime difference...'.format(turnTime_tolerance))
    if isVerbose:
        print('Before combining, total records is {}.'.format(len(eventList_DF_noOverlap_temp)))
    i_index = 1
    while(i_index < len(eventList_DF_noOverlap_temp)):
        if(eventList_DF_noOverlap_temp['diff'].iloc[i_index] <= timedelta(seconds=turnTime_tolerance)):
            cluster_begin_temp = i_index - 1
            while((eventList_DF_noOverlap_temp['diff'].iloc[i_index] <= timedelta(seconds=turnTime_tolerance)) ):
                i_index += 1
                if (i_index > len(eventList_DF_noOverlap_temp)-1):
                    break
            cluster_end_temp = i_index - 1
            
            # Get minimum residue_diff: .iloc[a:b]=[a,b), .loc[a:b]=[a,b]
            min_residue_diff_index_temp = eventList_DF_noOverlap_temp['residue_diff'].iloc[cluster_begin_temp:cluster_end_temp+1].idxmin()
            
            # Set record with min_residue as true, others as false.
            eventList_DF_noOverlap_temp.iloc[cluster_begin_temp:cluster_end_temp+1, index_column_keepFlag] = False
            eventList_DF_noOverlap_temp.loc[min_residue_diff_index_temp, 'keepFlag'] = True
        else:
            # For .iloc, can only use integer to specify row and column. No column string allowed.
            eventList_DF_noOverlap_temp.iloc[i_index, index_column_keepFlag] = True
            i_index += 1

    # Remove the records labeled as false.
    eventList_DF_noOverlap_temp = eventList_DF_noOverlap_temp[eventList_DF_noOverlap_temp['keepFlag']==True]
    eventList_DF_noOverlap_temp.reset_index(drop=True, inplace=True)
    eventList_DF_noOverlap_temp = eventList_DF_noOverlap_temp.drop('diff', axis=1)
    eventList_DF_noOverlap_temp.reset_index(drop=True, inplace=True)
    if isVerbose:
        print('After combining, total records is {}.'.format(len(eventList_DF_noOverlap_temp)))
    print('Done.')
    eventList_DF_noOverlap = eventList_DF_noOverlap_temp.copy()
    eventList_DF_noOverlap = eventList_DF_noOverlap.sort_values(by='startTime')

    # Eliminate events that overlap
    print('\nRemoving overlapped events...')
    if isVerbose:
        print('Before removing, total records is {}.'.format(len(eventList_DF_noOverlap)))
    overlap_length = len(eventList_DF_noOverlap) +1
    while overlap_length > len(eventList_DF_noOverlap):
        # Create array of intervals
        start_array = np.array([pd.Timestamp(elem) for elem in eventList_DF_noOverlap['startTime'].values])
        end_array = np.array([pd.Timestamp(elem) for elem in eventList_DF_noOverlap['endTime'].values])
        interval_array = pd.arrays.IntervalArray.from_arrays(start_array,end_array,closed='both')

        # Retrieving column index from column name.
        eventList_DF_noOverlap = eventList_DF_noOverlap.assign(keepFlag=[True]*len(eventList_DF_noOverlap))
        idx_column_keepFlag = eventList_DF_noOverlap.columns.get_loc('keepFlag')

        # Keep event with highest duration
        for i_index in range(len(eventList_DF_noOverlap)):
            interval = interval_array[i_index]
            idx_overlap = np.where(interval_array.overlaps(interval)==True)[0]
            max_duration_idx = eventList_DF_noOverlap.loc[idx_overlap]['duration'].astype(float).idxmax()
            eventList_DF_noOverlap.iloc[idx_overlap, idx_column_keepFlag] = False
            eventList_DF_noOverlap.loc[max_duration_idx, 'keepFlag'] = True

        eventList_DF_noOverlap = eventList_DF_noOverlap[eventList_DF_noOverlap['keepFlag']==True]
        eventList_DF_noOverlap = eventList_DF_noOverlap.drop(columns=['keepFlag'])
        eventList_DF_noOverlap.reset_index(drop=True,inplace=True)

        start_array = np.array([pd.Timestamp(elem) for elem in eventList_DF_noOverlap['startTime'].values])
        end_array = np.array([pd.Timestamp(elem) for elem in eventList_DF_noOverlap['endTime'].values])
        interval_array = pd.arrays.IntervalArray.from_arrays(start_array,end_array,closed='both')
        overlap_length = np.array([interval_array.overlaps(interval_array[i]).sum() for i in range(len(interval_array))]).sum()
        if isVerbose:
            print('After removing overlapped events, total records is {}.'.format(len(eventList_DF_noOverlap)))
        print('Done.')

    # Save DataFrame to output file.
    print('\nSaving eventList_DF_noOverlap to pickle file...')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    eventList_DF_noOverlap.to_csv(output_dir + '/' + output_filename + '.csv')
    progressionList.to_csv('/home/rharvey/GS/GS_SearchResult/progression/' + year_str + '_progression' + residue_str + '.csv')
    print('Done.')

    return eventList_DF_noOverlap

#########################################################################

# Calculate more information of given flux rope.
def get_more_flux_rope_info(data_DF, dataObject_or_dataPath, **kwargs):
    # Check input datatype: dataObject_or_dataPath should be the data or path of no overlapped eventlist.
    if isinstance(dataObject_or_dataPath, pd.DataFrame):
        print('\nYour input is a DataFrame data.')
        search_result_no_overlap_DF = dataObject_or_dataPath
    elif isinstance(dataObject_or_dataPath, str):
        print('\nYour input is a path. Load the dictionary data via this path.')
        search_result_no_overlap_DF = pd.read_pickle(open(dataObject_or_dataPath, 'rb'))
    else:
        print('\nPlease input the correct datatype! The input data must be a DataFrame or the path of a DataFrame!')
        return None
    
    output_filename = 'search_result_detailed_info' + namestr
    output_dir = os.getcwd()
    isVerbose = False
    
    print('\nDefault parameters:')
    print('output_dir            = {}.'.format(output_dir))
    print('output_filename       = {}.'.format(output_filename))
    
    # If keyword is specified, overwrite the default value.
    if 'output_dir' in kwargs:
        output_dir = kwargs['output_dir']
        print('output_dir is set to {}.'.format(output_dir))
    if 'output_filename' in kwargs:
        output_filename = kwargs['output_filename']
        print('output_filename is set to {}.'.format(output_filename))
    if 'isVerbose' in kwargs:
        isVerbose = kwargs['isVerbose']
        print('isVerbose is set to {}'.format())
    
    # Create an empty dataframe.
    eventList_DF_detailedInfo = pd.DataFrame(columns=['startTime', 'turnTime', 'endTime', \
        'duration', 'size_inAU','residue_diff', 'residue_fit', 'theta_deg', 'phi_deg', 'A_range', \
        'Pt_coeff', 'Path_length', 'VHT_inGSE[0]', 'VHT_inGSE[1]', 'VHT_inGSE[2]', \
        'X_unitVector[0]', 'X_unitVector[1]', 'X_unitVector[2]', 'Y_unitVector[0]', \
        'Y_unitVector[1]', 'Y_unitVector[2]', 'Z_unitVector[0]', 'Z_unitVector[1]', \
        'Z_unitVector[2]', 'walenTest_slope', 'walenTest_intercept', 'walenTest_r_value', \
        'walenTest_slope_b4reverse','walenTest_intercept_b4reverse', 'walenTest_r_value_b4reverse', \
        'B_abs_mean', 'Bx_abs_mean', 'By_abs_mean', 'Bz_abs_mean', 'B_std', \
        'Bx_std', 'By_std', 'Bz_std', 'Bx_inFR_abs_mean', 'By_inFR_abs_mean', 'Bz_inFR_abs_mean', \
        'Bx_inFR_std', 'By_inFR_std', 'Bz_inFR_std', 'B_magnitude_max', 'Vsw_magnitude_mean', \
        'Tp_mean', 'Np_mean', 'Te_mean', 'Beta_mean', 'Beta_p_mean', \
        'lambda1', 'lambda2', 'lambda3', 'eigenVectorMaxVar_lambda1[0]', 'eigenVectorMaxVar_lambda1[1]', \
        'eigenVectorMaxVar_lambda1[2]', 'eigenVectorInterVar_lambda2[0]', 'eigenVectorInterVar_lambda2[1]', \
        'eigenVectorInterVar_lambda2[2]', 'eigenVectorMinVar_lambda3[0]', 'eigenVectorMinVar_lambda3[1]', \
        'eigenVectorMinVar_lambda3[2]', 'Mach_average', \
        'Vswlambda1', 'Vswlambda2', 'Vswlambda3', 'VsweigenVectorMaxVar_lambda1[0]', \
        'VsweigenVectorMaxVar_lambda1[1]', 'VsweigenVectorMaxVar_lambda1[2]', 'VsweigenVectorInterVar_lambda2[0]', \
        'VsweigenVectorInterVar_lambda2[1]', 'VsweigenVectorInterVar_lambda2[2]', 'VsweigenVectorMinVar_lambda3[0]', \
        'VsweigenVectorMinVar_lambda3[1]', 'VsweigenVectorMinVar_lambda3[2]','Vrlambda1', 'Vrlambda2', 'Vrlambda3', \
        'VreigenVectorMaxVar_lambda1[0]', 'VreigenVectorMaxVar_lambda1[1]', 'VreigenVectorMaxVar_lambda1[2]', \
        'VreigenVectorInterVar_lambda2[0]', 'VreigenVectorInterVar_lambda2[1]', 'VreigenVectorInterVar_lambda2[2]', \
        'VreigenVectorMinVar_lambda3[0]', 'VreigenVectorMinVar_lambda3[1]', 'VreigenVectorMinVar_lambda3[2]',\
        'RD_mean', 'NaNp_mean', 'Alpha2Proton_ratio_mean', \
        'Jzmax','VA_mean','Cross_heli','Br_mean','Bt_mean','Bn_mean','Residue_energy','Am'])

    for index_FR in range(len(search_result_no_overlap_DF)):
        print('\nCalculating detailed information of flux ropes: {}/{}...'.format(index_FR+1, len(search_result_no_overlap_DF)))
        oneEvent            = search_result_no_overlap_DF.iloc[index_FR]
        startTime           = oneEvent['startTime']
        turnTime            = oneEvent['turnTime']
        endTime             = oneEvent['endTime']
        # duration          = oneEvent['duration']
        duration            = (endTime - startTime).total_seconds()/60
        residue_diff        = oneEvent['residue_diff']
        residue_fit         = oneEvent['residue_fit']
        theta_deg, phi_deg  = oneEvent['theta_phi']
        VHT_inGSE           = np.array(oneEvent['VHT'])
        
        if isVerbose:
            print('startTime = {}'.format(startTime))
            print('turnTime = {}'.format(turnTime))
            print('endTime = {}'.format(endTime))
            print('(theta_deg, phi_deg) = ({},{})'.format(theta_deg, phi_deg))
            print('residue_diff = {}'.format(residue_diff))
            print('residue_fit = {}'.format(residue_fit))
                    
        # Grab data in specific range of fluxrope candidate.
        selectedRange_mask = (data_DF.index >= startTime) & (data_DF.index <= endTime)
        data_oneFR_DF = data_DF.iloc[selectedRange_mask]
        if selectedRange_mask.sum() <=1:
                continue
        dt = float(sys.argv[3]) # seconds.
        
        # Get slice of data.
        B_inGSE   = data_oneFR_DF.loc[:,['Bx', 'By', 'Bz']] # [nT]
        Vsw_inGSE = data_oneFR_DF.loc[:,['Vx', 'Vy', 'Vz']] # [km/s]
        Np        = data_oneFR_DF.loc[:,['Np']] # [cm^-3]
        Tp        = data_oneFR_DF.loc[:,['Tp']] # [K]

        # RD = data_oneFR_DF.ix[:,['RD']]
        # Na = data_oneFR_DF.ix[:,['N_alpha']]
        # Alpha2Proton_ratio = data_oneFR_DF.loc[:,['Alpha2Proton_ratio']]
        # Alpha2Proton_ratio =  data_oneFR_DF.ix[:,['Alpha2Proton_ratio']]

        if 'Te' in data_oneFR_DF.keys():
            Te = data_oneFR_DF.loc[:,['Te']] # [K]
        
        # If there is any NaN in B_inGSE, try to interpolate.
        if B_inGSE.isnull().values.sum():
            if isVerbose:
                print('Found NaNs, interpolate B.')
            B_inGSE_copy = B_inGSE.copy()
            # limit=3 means only interpolate the gap shorter than 4.
            B_inGSE_copy.interpolate(method='time', limit=None, inplace=True)
            # interpolate won't fill leading NaNs, so we use backward fill.
            B_inGSE_copy.bfill(inplace=True)
            B_inGSE_copy.ffill(inplace=True)
            if B_inGSE_copy.isnull().values.sum():
                print('Too many NaNs in B. Skip this record. If this situation happens, please check.')
                detailed_info_dict = None
                continue
            else:
                B_inGSE = B_inGSE_copy

        # If there is any NaN in Vsw_inGSE, try to interpolate.
        if Vsw_inGSE.isnull().values.sum():
            if isVerbose:
                print('Found NaNs, interpolate Vsw.')
            Vsw_inGSE_copy = Vsw_inGSE.copy()
            Vsw_inGSE_copy.interpolate(method='time', limit=None, inplace=True)
            Vsw_inGSE_copy.bfill(inplace=True)
            Vsw_inGSE_copy.ffill(inplace=True)
            if Vsw_inGSE_copy.isnull().values.sum():
                print('Too many NaNs in Vsw. Skip this record. If this situation happens, please check.')
                detailed_info_dict = None
                continue
            else:
                Vsw_inGSE = Vsw_inGSE_copy
                
        # If there is any NaN in Np, try to interpolate.
        if Np.isnull().values.sum():
            if isVerbose:
                print('Found NaNs, interpolate Np.')
            Np_copy = Np.copy()
            Np_copy.interpolate(method='time', limit=None, inplace=True)
            Np_copy.bfill(inplace=True)
            Np_copy.ffill(inplace=True)
            if Np_copy.isnull().values.sum():
                print('Too many NaNs in Vsw. Skip this record. If this situation happens, please check.')
                detailed_info_dict = None
                continue
            else:
                Np = Np_copy

        # Use direction cosines to construct a unit vector.
            # x = rcos(alpha) = rsin(theta)cos(phi) => cos(alpha) = sin(theta)cos(phi)
            # y = rcos(beta)  = rsin(theta)sin(phi) => cos(beta)  = sin(theta)sin(phi)
            # z = rcos(gamma) = rcos(theta)         => cos(gamma) = cos(theta)
        theta_rad = factor_deg2rad * theta_deg
        phi_rad   = factor_deg2rad * phi_deg

        # Form new Z_unitVector according to direction cosines.
        Z_unitVector = np.array([np.sin(theta_rad)*np.cos(phi_rad), np.sin(theta_rad)*np.sin(phi_rad), np.cos(theta_rad)])
        X_unitVector = findXaxis(Z_unitVector, -VHT_inGSE)            # Find X axis from Z axis and -VHT.
        Y_unitVector = formRightHandFrame(X_unitVector, Z_unitVector) # Find the Y axis to form a right-handed coordinater with X and Z.

        # Project B_inGSE & Vsw_inFR into FR Frame.
        matrix_transToFluxRopeFrame = np.array([X_unitVector, Y_unitVector, Z_unitVector]).T
        B_inFR   = B_inGSE.dot(matrix_transToFluxRopeFrame)
        VHT_inFR = VHT_inGSE.dot(matrix_transToFluxRopeFrame).
        Vsw_inFR = Vsw_inGSE.dot(matrix_transToFluxRopeFrame)
        
        # Apply Walen test on the result (in optimal frame).
        # Proton mass density. Original Np is in #/cc ( cc = cubic centimeter). Multiply by 1e6 to convert cc to m^3.
        P_massDensity = Np * m_proton * 1e6 # kg/m^3.
        len_P_massDensity = len(P_massDensity)
        P_massDensity_array = np.array(P_massDensity)
        P_massDensity_array = np.reshape(P_massDensity_array, (len_P_massDensity, 1))

        # Alfven speed. Multiply by 1e-9 to convert nT to T. Divided by 1000.0 to convert m/s to km/s.
        VA_inFR         = np.array(B_inFR * 1e-9) / np.sqrt(mu0 * P_massDensity_array) / 1000.0
        VA_inFR_1D      = np.reshape(VA_inFR, VA_inFR.size)
        V_remaining     = np.array(Vsw_inFR - VHT_inFR)
        V_remaining_1D  = np.reshape(V_remaining, V_remaining.size)
 
        walenTest_slope_b4reverse, walenTest_intercept_b4reverse, walenTest_r_value_b4reverse = walenTest(VA_inFR_1D, V_remaining_1D)        

        # Check if Bz has negative values, if does, flip Z-axis direction.
        num_Bz_lt0 = (B_inFR[2]<0).sum()
        num_Bz_gt0 = (B_inFR[2]>0).sum()
        # If the negative Bz is more than positive Bz, filp.
        if (num_Bz_lt0 > num_Bz_gt0):
            # Reverse the direction of Z-axis.
            print('Reverse the direction of Z-axis!')
            Z_unitVector = -Z_unitVector

            # Recalculate theta and phi with new Z_unitVector.
            theta_deg, phi_deg = directionVector2angle(Z_unitVector) 
            X_unitVector = findXaxis(Z_unitVector, -VHT_inGSE)                                  # Refind X axis frome Z axis and -Vsw.
            Y_unitVector = formRightHandFrame(X_unitVector, Z_unitVector)                        # Refind the Y axis to form a right-handed coordinater with X and Z.  
            matrix_transToFluxRopeFrame = np.array([X_unitVector, Y_unitVector, Z_unitVector]).T # Reproject B_inGSE_DataFrame into FR frame.
            B_inFR = B_inGSE.dot(matrix_transToFluxRopeFrame)
        
        # Project VHT_inGSE & Vsw_inFR into FR Frame.
        VHT_inFR = VHT_inGSE.dot(matrix_transToFluxRopeFrame)
        Vsw_inFR = Vsw_inGSE.dot(matrix_transToFluxRopeFrame)
        
        # Alfven speed. Multiply by 1e-9 to convert nT to T. Divided by 1000.0 to convert m/s to km/s.
        # VA_inGSE = np.array(B_inGSE * 1e-9) / np.sqrt(mu0 * P_massDensity_array) / 1000.0
        VA_inFR    = np.array(B_inFR * 1e-9) / np.sqrt(mu0 * P_massDensity_array) / 1000.0
        VA_inFR_1D = np.reshape(VA_inFR, VA_inFR.size)

        size = - VHT_inFR[0] * 1000.0 * duration * 60.0 # Space increment along X axis. Convert km/s to m/s.
        AU = 149597870700
        size_inAU = size/AU

        V_remaining = np.array(Vsw_inFR - VHT_inFR)
        V_remaining_1D = np.reshape(V_remaining, V_remaining.size)
        
        # Call Walen test function: first row is x component, second row is y component, third row is z component.
        walenTest_slope, walenTest_intercept, walenTest_r_value = walenTest(VA_inFR_1D, V_remaining_1D)
        
        # Calculate the covariance matrix of Magnetic field.
        covM_B_inGSE = B_inGSE.cov()
        
        # Calculate the eigenvalues and eigenvectors of convariance matrix of B field.
        lambda1, lambda2, lambda3, eigenVectorMaxVar_lambda1, eigenVectorInterVar_lambda2, eigenVectorMinVar_lambda3 = eigenMatrix(covM_B_inGSE, formXYZ=True)
        covM_Vsw_inGSE = Vsw_inGSE.cov()
        Vswlambda1, Vswlambda2, Vswlambda3, VsweigenVectorMaxVar_lambda1, VsweigenVectorInterVar_lambda2, VsweigenVectorMinVar_lambda3 = eigenMatrix(covM_Vsw_inGSE, formXYZ=True)

        Vr_inGSE = Vsw_inGSE - VHT_inGSE
        covM_Vr_inGSE = Vr_inGSE.cov()
        Vrlambda1, Vrlambda2, Vrlambda3, VreigenVectorMaxVar_lambda1, VreigenVectorInterVar_lambda2, VreigenVectorMinVar_lambda3 = eigenMatrix(covM_Vr_inGSE, formXYZ=True)

        '''
        Project B_DataFrame & VHt_inFR onto new frame (MVB frame).The dot product of two dataframe requires the columns and indices are same, so we convert to np.array:
            B_inMVB = B_inGSE.dot(np.array(eigenVectors_covM_B_inGSE))
            VHT_inMVB = VHT_inGSE.dot(np.array(eigenVectors_covM_B_inGSE))
        '''

        Mach = np.sqrt(np.square(V_remaining).sum(axis=1))/np.sqrt(np.square(VA_inFR).sum(axis=1))
        a = Mach.mean()**2
        B_norm_DF = pd.DataFrame(np.sqrt(np.square(B_inGSE).sum(axis=1)),columns=['|B|'])

        # Calculate A(x,0) by integrating By. A(x,0) = Integrate[-By(s,0)ds, {s, 0, x}], where ds = -Vht dot unit(X) dt.
        ds = - VHT_inFR[0] * 1000.0 * dt
        A = integrate.cumulative_trapezoid(-(1-a)*B_inFR[1]*1e-9, dx=ds, initial=0)
        
        # Calculate Pt(x,0). Pt(x,0)=p(x,0)+Bz^2/(2mu0). By = B_inFR[2]        
        Pp1 = np.array(Np['Np']) * 1e6 * 1.3806488e-23 * 1e9 * np.array(Tp['Tp'])
        Pb1 = np.array((B_inFR[2] * 1e-9)**2 / (2.0*mu0) * 1e9) # 1e9 convert unit form pa to npa.
        PB1 = np.array((B_norm_DF['|B|'] * 1e-9)**2 / (2.0*mu0) * 1e9) # 1e9 convert unit form pa to npa.         
        Pt  = ((1-a)**2)*Pb1 + (1-a)*Pp1 + (a*(1-a))*PB1
        # Pt = np.array((B_inFR[2] * 1e-9)**2 / (2.0*mu0) * 1e9) # 1e9 convert unit form pa to npa.
        
        # Find the index of turnPoint.
        index_turnTime = B_inFR.index.get_indexer([turnTime],method='nearest')[0]
        # Split A and Pt into two branches.
        A_sub1  = A[:index_turnTime+1]
        A_sub2  = A[index_turnTime:]
        Pt_sub1 = Pt[:index_turnTime+1]
        Pt_sub2 = Pt[index_turnTime:]
        
        z         = np.polyfit(A, Pt, 3)
        Func_Pt_A = np.poly1d(z)
        Func_Jz   = np.polyder(Func_Pt_A) 
        Pt_coeff  = list(z)
        
        A_phy   = A / (1-a)
        A_range = [min(A_phy), max(A_phy)]
        A_norm  = A_range/max(A_phy)
        if abs(min(A_phy)) > abs(max(A_phy)):
            Am = min(A_phy)
        else:
            Am = max(A_phy)
        Path_length = ds * duration # The lenght of spacecraft trajectory across the flux rope.
        Jz = Func_Jz(A_range)
        Jzmax = np.nanmax(Jz)

        '''
        plt.plot(A_sub1, Pt_sub1, 'ro-', A_sub2, Pt_sub2, 'bo-', np.sort(A), Func_Pt_A(np.sort(A)),'g--')
        plt.title('residue_diff = {},  residue_fit = {}'.format(residue_diff, residue_fit))
        plt.show()
        '''
        
        # Get B statistical properties.
        B_magnitude_max = B_norm_DF['|B|'].max(skipna=True)
        B_inGSE = pd.concat([B_inGSE, B_norm_DF], axis=1)
        B_std_Series = B_inGSE.std(axis=0,skipna=True,numeric_only=True)
        B_abs_mean_Series = B_inGSE.abs().mean(axis=0,skipna=True,numeric_only=True)
        
        B_mean_Series = B_inGSE.mean(axis=0,skipna=True,numeric_only=True)
        Br_mean = round(B_mean_Series[0],4)
        Bt_mean = round(B_mean_Series[1],4)
        Bn_mean = round(B_mean_Series[2],4)

        B_abs_mean  = round(B_abs_mean_Series['|B|'],4)
        Bx_abs_mean = round(B_abs_mean_Series[0],4)
        By_abs_mean = round(B_abs_mean_Series[1],4)
        Bz_abs_mean = round(B_abs_mean_Series[2],4)
        B_std  = round(B_std_Series['|B|'],4)
        Bx_std = round(B_std_Series[0],4)
        By_std = round(B_std_Series[1],4)
        Bz_std = round(B_std_Series[2],4)
        
        # B_inFR.
        B_inFR_std_Series = B_inFR.std(axis=0,skipna=True,numeric_only=True)
        B_inFR_abs_mean_Series = B_inFR.abs().mean(axis=0,skipna=True,numeric_only=True)
        Bx_inFR_abs_mean = round(B_inFR_abs_mean_Series[0],4)
        By_inFR_abs_mean = round(B_inFR_abs_mean_Series[1],4)
        Bz_inFR_abs_mean = round(B_inFR_abs_mean_Series[2],4)
        Bx_inFR_std = round(B_inFR_std_Series[0],4)
        By_inFR_std = round(B_inFR_std_Series[1],4)
        Bz_inFR_std = round(B_inFR_std_Series[2],4)
        
        Mach_average = np.array((np.sqrt(np.square(V_remaining).sum(axis=1))/np.sqrt(np.square(VA_inFR).sum(axis=1))).mean())
        Vsw_norm_DF  = pd.DataFrame(np.sqrt(np.square(Vsw_inGSE).sum(axis=1)),columns=['|Vsw|'])
        Vsw_magnitude_mean = Vsw_norm_DF['|Vsw|'].mean(skipna=True)    
        VA_mean = np.mean(np.sqrt(np.square(VA_inFR).sum(axis=1)))
        RD_mean = np.mean(np.array(RD['RD']))
        NaNp_mean = np.mean(np.array(Na['N_alpha'])/np.array(Np['Np']))
        
        # Get Plasma Beta statistical properties.
        Tp_mean = np.mean(np.ma.masked_invalid(np.array(Tp['Tp'])))/1e6 # Convert unit to 10^6K.
        Tp_mean = round(Tp_mean, 6)
        Np_mean = float(Np.mean(skipna=True, numeric_only=True)) # [#/cc]
        Alpha2Proton_ratio_mean = float(Alpha2Proton_ratio.mean(skipna=True, numeric_only=True))
        Alpha2Proton_ratio_mean = np.mean(np.ma.masked_invalid(np.array(Alpha2Proton_ratio['Alpha2Proton_ratio'])))
        
        # Calculate Te_mean. Convert unit to 10^6K.
        if 'Te' in data_oneFR_DF.keys():
            Te_mean = float(Te.mean(skipna=True, numeric_only=True))/1e6
        else:
            Te_mean = None
        
        # Get Other statistical properties.
        Bx_inFR_VA = np.array(B_inFR[0] * 1e-9)/np.sqrt(mu0 * np.array(P_massDensity['Np'])) / 1000.0 # km/s
        By_inFR_VA = np.array(B_inFR[1] * 1e-9)/np.sqrt(mu0 * np.array(P_massDensity['Np'])) / 1000.0 # km/s
        Bz_inFR_VA = np.array(B_inFR[2] * 1e-9)/np.sqrt(mu0 * np.array(P_massDensity['Np'])) / 1000.0 # km/s
        Cross_heli=2*(V_remaining[:,0]*Bx_inFR_VA+V_remaining[:,1]*By_inFR_VA+V_remaining[:,2]*Bz_inFR_VA).mean()/((V_remaining[:,0]**2+V_remaining[:,1]**2+V_remaining[:,2]**2).mean()+(Bx_inFR_VA**2+By_inFR_VA**2+Bz_inFR_VA**2).mean())
        Residue_energy=((V_remaining[:,0]**2+V_remaining[:,1]**2+V_remaining[:,2]**2).mean()-(Bx_inFR_VA**2+By_inFR_VA**2+Bz_inFR_VA**2).mean())/((V_remaining[:,0]**2+V_remaining[:,1]**2+V_remaining[:,2]**2).mean()+(Bx_inFR_VA**2+By_inFR_VA**2+Bz_inFR_VA**2).mean())
        
        # Calculate plasma Dynamic Pressure PD.
        Pp = np.array(Np['Np']) * 1e6 * k_Boltzmann * np.array(Tp['Tp']) # Proton pressure.
        if 'Te' in data_oneFR_DF.keys():
            Pe = np.array(Np['Np']) * 1e6 * k_Boltzmann * np.array(Te['Te']) # Electron pressure.
            PD = Pp + Pe # Total dynamic pressure.
        else:
            PD = Pp

        # Calculate plasma Magnetic pressure PB.
        PB = (np.array(B_norm_DF['|B|'])*1e-9)**2/(2*mu0)
        
        # Calculate plasma Beta = PD/PB
        Beta        = PD/PB
        Beta_mean   = np.mean(np.ma.masked_invalid(Beta)) # Exclude nan and inf.
        Beta_p      = Pp/PB
        Beta_p_mean = np.mean(np.ma.masked_invalid(Beta_p))

        #detailed_info_dict = {'startTime':startTime, 'turnTime':turnTime, 'endTime':endTime, 'duration':duration, 'residue_diff':residue_diff, 'residue_fit':residue_fit, 'theta_deg':theta_deg, 'phi_deg':phi_deg, 'A_range':A_range, 'Pt_coeff':Pt_coeff, 'Path_length':Path_length, 'VHT_inGSE[0]':VHT_inGSE[0], 'VHT_inGSE[1]':VHT_inGSE[1], 'VHT_inGSE[2]':VHT_inGSE[2], 'X_unitVector[0]':X_unitVector[0], 'X_unitVector[1]':X_unitVector[1], 'X_unitVector[2]':X_unitVector[2], 'Y_unitVector[0]':Y_unitVector[0], 'Y_unitVector[1]':Y_unitVector[1], 'Y_unitVector[2]':Y_unitVector[2], 'Z_unitVector[0]':Z_unitVector[0], 'Z_unitVector[1]':Z_unitVector[1], 'Z_unitVector[2]':Z_unitVector[2], 'walenTest_slope':walenTest_slope, 'walenTest_intercept':walenTest_intercept, 'walenTest_r_value':walenTest_r_value,  'B_abs_mean':B_abs_mean, 'Bx_abs_mean':Bx_abs_mean, 'By_abs_mean':By_abs_mean, 'Bz_abs_mean':Bz_abs_mean, 'B_std':B_std, 'Bx_std':Bx_std, 'By_std':By_std, 'Bz_std':Bz_std, 'Bx_inFR_abs_mean':Bx_inFR_abs_mean, 'By_inFR_abs_mean':By_inFR_abs_mean, 'Bz_inFR_abs_mean':Bz_inFR_abs_mean, 'Bx_inFR_std':Bx_inFR_std, 'By_inFR_std':By_inFR_std, 'Bz_inFR_std':Bz_inFR_std, 'B_magnitude_max':B_magnitude_max, 'Vsw_magnitude_mean':Vsw_magnitude_mean, 'Tp_mean':Tp_mean, 'Np_mean':Np_mean, 'Te_mean':Te_mean, 'Beta_mean':Beta_mean, 'Beta_p_mean':Beta_p_mean, 'lambda1':lambda1, 'lambda2':lambda2, 'lambda3':lambda3, 'eigenVectorMaxVar_lambda1[0]':eigenVectorMaxVar_lambda1[0], 'eigenVectorMaxVar_lambda1[1]':eigenVectorMaxVar_lambda1[1], 'eigenVectorMaxVar_lambda1[2]':eigenVectorMaxVar_lambda1[2], 'eigenVectorInterVar_lambda2[0]':eigenVectorInterVar_lambda2[0], 'eigenVectorInterVar_lambda2[1]':eigenVectorInterVar_lambda2[1], 'eigenVectorInterVar_lambda2[2]':eigenVectorInterVar_lambda2[2], 'eigenVectorMinVar_lambda3[0]':eigenVectorMinVar_lambda3[0], 'eigenVectorMinVar_lambda3[1]':eigenVectorMinVar_lambda3[1], 'eigenVectorMinVar_lambda3[2]':eigenVectorMinVar_lambda3[2], 'Vswlambda1':Vswlambda1, 'Vswlambda2':Vswlambda2, 'Vswlambda3':Vswlambda3, 'VsweigenVectorMaxVar_lambda1[0]':VsweigenVectorMaxVar_lambda1[0], 'VsweigenVectorMaxVar_lambda1[1]':VsweigenVectorMaxVar_lambda1[1], 'VsweigenVectorMaxVar_lambda1[2]':VsweigenVectorMaxVar_lambda1[2], 'VsweigenVectorInterVar_lambda2[0]':VsweigenVectorInterVar_lambda2[0], 'VsweigenVectorInterVar_lambda2[1]':VsweigenVectorInterVar_lambda2[1], 'VsweigenVectorInterVar_lambda2[2]':VsweigenVectorInterVar_lambda2[2], 'VsweigenVectorMinVar_lambda3[0]':VsweigenVectorMinVar_lambda3[0], 'VsweigenVectorMinVar_lambda3[1]':VsweigenVectorMinVar_lambda3[1], 'eigenVectorMinVar_lambda3[2]':eigenVectorMinVar_lambda3[2], 'Vrlambda1':Vrlambda1, 'Vrlambda2':Vrlambda2, 'Vrlambda3':Vrlambda3, 'VreigenVectorMaxVar_lambda1[0]':VreigenVectorMaxVar_lambda1[0], 'VreigenVectorMaxVar_lambda1[1]':VreigenVectorMaxVar_lambda1[1], 'VreigenVectorMaxVar_lambda1[2]':VreigenVectorMaxVar_lambda1[2], 'VreigenVectorInterVar_lambda2[0]':VreigenVectorInterVar_lambda2[0], 'VreigenVectorInterVar_lambda2[1]':VreigenVectorInterVar_lambda2[1], 'VreigenVectorInterVar_lambda2[2]':VreigenVectorInterVar_lambda2[2], 'VreigenVectorMinVar_lambda3[0]':VreigenVectorMinVar_lambda3[0], 'VreigenVectorMinVar_lambda3[1]':VreigenVectorMinVar_lambda3[1], 'VreigenVectorMinVar_lambda3[2]':VreigenVectorMinVar_lambda3[2], 'Mach_average':Mach_average, 'Average_Mach':Average_Mach, 'Alpha2Proton_ratio_mean': Alpha2Proton_ratio_mean, 'Alpha2Proton_ratio_max':Alpha2Proton_ratio_max, 'Alpha2Proton_ratio_min':Alpha2Proton_ratio_min}
        detailed_info_dict = {'startTime':startTime, 'turnTime':turnTime, 'endTime':endTime, \
        'duration':duration, 'size_inAU':size_inAU,'residue_diff':residue_diff, 'residue_fit':residue_fit, 'theta_deg':theta_deg, \
        'phi_deg':phi_deg, 'A_range':A_range, 'Pt_coeff':Pt_coeff, 'Path_length':Path_length, 'VHT_inGSE[0]':VHT_inGSE[0], \
        'VHT_inGSE[1]':VHT_inGSE[1], 'VHT_inGSE[2]':VHT_inGSE[2], 'X_unitVector[0]':X_unitVector[0], \
        'X_unitVector[1]':X_unitVector[1], 'X_unitVector[2]':X_unitVector[2], 'Y_unitVector[0]':Y_unitVector[0], \
        'Y_unitVector[1]':Y_unitVector[1], 'Y_unitVector[2]':Y_unitVector[2], 'Z_unitVector[0]':Z_unitVector[0], \
        'Z_unitVector[1]':Z_unitVector[1], 'Z_unitVector[2]':Z_unitVector[2], 'walenTest_slope':walenTest_slope, \
        'walenTest_intercept':walenTest_intercept, 'walenTest_r_value':walenTest_r_value,  'B_abs_mean':B_abs_mean, \
        'walenTest_slope_b4reverse':walenTest_slope_b4reverse, 'walenTest_intercept_b4reverse':walenTest_intercept_b4reverse, 'walenTest_r_value_b4reverse':walenTest_r_value_b4reverse, \
        'Bx_abs_mean':Bx_abs_mean, 'By_abs_mean':By_abs_mean, 'Bz_abs_mean':Bz_abs_mean, 'B_std':B_std, 'Bx_std':Bx_std, \
        'By_std':By_std, 'Bz_std':Bz_std, 'Bx_inFR_abs_mean':Bx_inFR_abs_mean, 'By_inFR_abs_mean':By_inFR_abs_mean, \
        'Bz_inFR_abs_mean':Bz_inFR_abs_mean, 'Bx_inFR_std':Bx_inFR_std, 'By_inFR_std':By_inFR_std, 'Bz_inFR_std':Bz_inFR_std, \
        'B_magnitude_max':B_magnitude_max, 'Vsw_magnitude_mean':Vsw_magnitude_mean, 'Tp_mean':Tp_mean, 'Np_mean':Np_mean, \
        'Te_mean':Te_mean, 'Beta_mean':Beta_mean, 'Beta_p_mean':Beta_p_mean, 'lambda1':lambda1, 'lambda2':lambda2, \
        'lambda3':lambda3, 'eigenVectorMaxVar_lambda1[0]':eigenVectorMaxVar_lambda1[0], \
        'eigenVectorMaxVar_lambda1[1]':eigenVectorMaxVar_lambda1[1], 'eigenVectorMaxVar_lambda1[2]':eigenVectorMaxVar_lambda1[2], \
        'eigenVectorInterVar_lambda2[0]':eigenVectorInterVar_lambda2[0], 'eigenVectorInterVar_lambda2[1]':eigenVectorInterVar_lambda2[1], \
        'eigenVectorInterVar_lambda2[2]':eigenVectorInterVar_lambda2[2], 'eigenVectorMinVar_lambda3[0]':eigenVectorMinVar_lambda3[0], \
        'eigenVectorMinVar_lambda3[1]':eigenVectorMinVar_lambda3[1], 'eigenVectorMinVar_lambda3[2]':eigenVectorMinVar_lambda3[2], \
        'Mach_average':Mach_average, 'Vswlambda1':Vswlambda1, \
        'Vswlambda2':Vswlambda2, 'Vswlambda3':Vswlambda3, 'VsweigenVectorMaxVar_lambda1[0]':VsweigenVectorMaxVar_lambda1[0], \
        'VsweigenVectorMaxVar_lambda1[1]':VsweigenVectorMaxVar_lambda1[1], 'VsweigenVectorMaxVar_lambda1[2]':VsweigenVectorMaxVar_lambda1[2], \
        'VsweigenVectorInterVar_lambda2[0]':VsweigenVectorInterVar_lambda2[0], 'VsweigenVectorInterVar_lambda2[1]':VsweigenVectorInterVar_lambda2[1], \
        'VsweigenVectorInterVar_lambda2[2]':VsweigenVectorInterVar_lambda2[2], 'VsweigenVectorMinVar_lambda3[0]':VsweigenVectorMinVar_lambda3[0], \
        'VsweigenVectorMinVar_lambda3[1]':VsweigenVectorMinVar_lambda3[1], 'VsweigenVectorMinVar_lambda3[2]':VsweigenVectorMinVar_lambda3[2],\
        'Vrlambda1':Vrlambda1, 'Vrlambda2':Vrlambda2, 'Vrlambda3':Vrlambda3, 'VreigenVectorMaxVar_lambda1[0]':VreigenVectorMaxVar_lambda1[0], \
        'VreigenVectorMaxVar_lambda1[1]':VreigenVectorMaxVar_lambda1[1], 'VreigenVectorMaxVar_lambda1[2]':VreigenVectorMaxVar_lambda1[2], \
        'VreigenVectorInterVar_lambda2[0]':VreigenVectorInterVar_lambda2[0], 'VreigenVectorInterVar_lambda2[1]':VreigenVectorInterVar_lambda2[1], \
        'VreigenVectorInterVar_lambda2[2]':VreigenVectorInterVar_lambda2[2], 'VreigenVectorMinVar_lambda3[0]':VreigenVectorMinVar_lambda3[0], \
        'VreigenVectorMinVar_lambda3[1]':VreigenVectorMinVar_lambda3[1], 'VreigenVectorMinVar_lambda3[2]':VreigenVectorMinVar_lambda3[2],\
        'RD_mean':RD_mean,'NaNp_mean':NaNp_mean, 'Alpha2Proton_ratio_mean':Alpha2Proton_ratio_mean, \
        'Jzmax':Jzmax,'VA_mean':VA_mean,'Cross_heli':Cross_heli,\
        'Br_mean':Br_mean,'Bt_mean':Bt_mean,'Bn_mean':Bn_mean,'Residue_energy':Residue_energy,'Am':Am}

        # Append detailed_info_dict to FR_detailed_info_DF.
        if not (detailed_info_dict is None):
            eventList_DF_detailedInfo = eventList_DF_detailedInfo._append(detailed_info_dict, ignore_index=True)

    # Save DataFrame to output file.
    print('\nSaving eventList_DF_detailedInfo to output file...')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    eventList_DF_detailedInfo.to_pickle(output_dir + '/' + output_filename + '.p')
    print('Done.')
    return eventList_DF_detailedInfo


###############################################################################

# Stack two images verticlly.
def vstack_images(file1, file2):
    """Merge two images into one, displayed side by side
    :param file1: path to first image file
    :param file2: path to second image file
    :return: the merged Image object
    """
    image1 = Image.open(file1)
    image2 = Image.open(file2)

    (width1, height1) = image1.size
    (width2, height2) = image2.size

    result_width = max(width1, width2)
    result_height = height1 + height2

    result = Image.new('RGB', (result_width, result_height), color=(255,255,255))
    result.paste(im=image1, box=(0, 0))
    result.paste(im=image2, box=(0, height1))
    return result