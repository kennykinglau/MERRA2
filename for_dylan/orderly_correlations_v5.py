#! /usr/bin/env python

import netCDF4 as nc
import numpy as np
import scipy.integrate as grate
import scipy.interpolate as polate
import scipy.stats as sts
import datetime as dt
import os
import argparse
from matplotlib import pyplot as plt
import pickle as pkl
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.layers.experimental import preprocessing
#Some of the machine learning code is taken from a TensorFlow tutorial here: https://www.tensorflow.org/tutorials/keras/regression

parser = argparse.ArgumentParser(epilog = 'The program assumes we have a circle centered around the South Pole.')
monthChoices = [1,2,3,4,5,6,7,8,9,10,11,12]
dayChoices = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31]
testChoices = ['Spearman','linregress','ks']
SPTchoices = ['Q','U','QUratio','Q-U']
DIM_DICT = {'CPT':'[J m^-2]', 'KE':'[J m^-2]', 'MASS':'[kg m^-2]', 'THV':'[K]', 'TOX':'[kg m^-2]', 'TQI':'[kg m^-2]', 'TQL':'[kg m^-2]', 'TQV':'[kg m]', \
    'CLOUD':'[1]', 'DELP':'[Pa]', 'EPV':'[K m^2 kg^-1 s^-1]', 'H':'[m]', 'O3':'[kg kg^-1]', 'OMEGA':'[Pa s^-1]', 'PHIS':'[m^2 s^-2]', 'PL':'[Pa]', 'PS':'[Pa]', \
    'QI':'[kg kg^-1]', 'QL':'[kg kg^-1]', 'QV':'[kg kg^-1]', 'RH':'[1]', 'SLP':'[Pa]', 'T':'[K]', 'U':'[m s^-1]', 'V':'[m s^-1]'}

LAT_RES = 0.5
LON_RES = 0.625
R = 6357000 #radius of Earth at the poles in meters
ZERO_MEAN = False
filters = ['None','<100m','100-1000m','>1000m']

NUM_OF_TICKS = 5

def plot_loss(history):
  plt.plot(history.history['loss'], label='loss')
  plt.plot(history.history['val_loss'], label='val_loss')
  plt.xlabel('Epoch')
  plt.ylabel('Error')
  plt.legend()
  plt.grid(True)

def compute_1d_psd(input_array, sample_spacing = 1.0): #This function is defined here intead of imported from demonstrate_psd.py because I figured it's better if this program can run without having to remember to move over another program as well.
    frequencies = np.fft.fftshift(np.fft.fftfreq(input_array.shape[0], d=sample_spacing))
    fourier_transform = np.fft.fftshift(np.fft.fft(input_array))
    power_spectrum = np.abs(fourier_transform)**2
    return frequencies,power_spectrum

def convert2time(dateString,START_DATE,END_DATE):
    #Finds the 3-hour interval that the time is in and then returns the hour number corresponding to the beginning of the 3-hour interval
    #Returns -1 if the time is after END_DATE so those data can easily be thrown out with the data before START_DATE
    #datetime.fromisoformat is not available in Python versions older than 3.7, so I did it myself 
    #observationTime = dt.fromisoformat(dateString)
    dateStrings = dateString.split(sep='T')
    dayString = dateStrings[0]
    timeString = dateStrings[1]
    
    dayStrings = dayString.split(sep='-')
    year = int(dayStrings[0])
    month = int(dayStrings[1])
    day = int(dayStrings[2])
    
    timeStrings = timeString.split(sep=':')
    hour = int(timeStrings[0])
    minute = int(timeStrings[1])
    second = int(round(float(timeStrings[2])))
    
    observationTime = dt.datetime(year,month,day,hour,minute,second)
    
    #Give the start and end-dates times so they can be compared to datetime objects
    START_DATE, END_DATE = dt.datetime(year=START_DATE.year,month=START_DATE.month,day=START_DATE.day,hour=0,minute=0,second=0), dt.datetime(year=END_DATE.year,month=END_DATE.month,day=END_DATE.day,hour=0,minute=0,second=0)
    
    t = (observationTime - START_DATE) / dt.timedelta(hours=1)
    #I think date objects without a time will have time = 0, so I'm going to add a day to END_DATE when checking if we've gone past it
    END_DATE = END_DATE + dt.timedelta(days=1)
    if (observationTime - END_DATE) > dt.timedelta(hours=0):
        t = -1
    return t

def get_merra_data(MERRA_VAR,DIM,FOLDER_NAME,FINAL_T_INTERVAL=-1, #FINAL_T_INTERVAL only has an effect if going from 1-hour intervals to 3-hour intervals
    FILE_BASE_INST2 = 'MERRA2_400.inst1_2d_asm_Nx.', FILE_END_INST2 = '.SUB.nc',FILE_BASE_TAVG2 = 'MERRA2_400.tavg1_2d_slv_Nx.',FILE_END_TAVG2 = '.SUB.nc',
    FILE_BASE_INST3 = 'MERRA2_400.inst3_3d_asm_Nv.', FILE_END_INST3 = '.SUB.nc',FILE_BASE_TAVG3 = 'MERRA2_400.tavg3_3d_asm_Nv.',FILE_END_TAVG3 = '.SUB.nc',
    MIN_FREQ = 0,MAX_FREQ=float('inf'),DELTA = True,time_type='tavg',
    NUM_OF_SAMPS=300, PRINT_MESS = False,
    START_YEAR = 2019, START_MONTH=3,START_DAY=22,
    END_YEAR = 2019, END_MONTH = 12, END_DAY=31):
    
    #This creates the folder where everything will go
    dir_path = os.path.dirname(os.path.realpath(__file__))
    folderPath = os.path.join(dir_path,FOLDER_NAME)
    if not os.path.isdir(folderPath):
        os.mkdir(folderPath)
    
    messages = []
    if DIM==3:
        T_INTERVAL = 3
    if DIM==2:
        T_INTERVAL = 1
    NUM_OF_TIMES = 24 // T_INTERVAL
    #Generate the dates
    START_DATE = dt.date(year=START_YEAR,month=START_MONTH,day=START_DAY)
    END_DATE = dt.date(year=END_YEAR,month=END_MONTH,day=END_DAY)
    date = START_DATE
    interval = dt.timedelta(days=1)
    dates = []
    while date <= END_DATE:
        dates.append(date.isoformat().replace('-',''))
        date = date + interval
    
    message = "Loading MERRA2 data..."
    if PRINT_MESS: print(message)
    messages.append(message)
    #Load the data
    setsOfData_inst = {}
    setsOfData_tavg = {}
    for date in dates:
        if DIM==3:
            setsOfData_inst[date] = nc.Dataset(FILE_BASE_INST3+date+FILE_END_INST3)
            setsOfData_tavg[date] = nc.Dataset(FILE_BASE_TAVG3+date+FILE_END_TAVG3)
        elif DIM==2:
            setsOfData_inst[date] = nc.Dataset(FILE_BASE_INST2+date+FILE_END_INST2)
            setsOfData_tavg[date] = nc.Dataset(FILE_BASE_TAVG2+date+FILE_END_TAVG2)
    
    message = "Loaded "+str(len(setsOfData_tavg))+ " days of data."
    if PRINT_MESS: print(message)
    messages.append(message)
    
    merra_data = np.zeros(NUM_OF_TIMES*len(dates))
    i = 0
    for date in dates:
        for t in range(NUM_OF_TIMES):
            if DIM==2:
                if DELTA:
                    merra_data[i] = setsOfData_tavg[date][MERRA_VAR][t,0,0] - setsOfData_inst[date][MERRA_VAR][t,0,0]
                elif time_type=='tavg':
                    merra_data[i] = setsOfData_tavg[date][MERRA_VAR][t,0,0]
                elif time_type=='inst':
                    merra_data[i] = setsOfData_inst[date][MERRA_VAR][t,0,0]
            elif DIM==3:
                #Load the data and then, because the pressure levels may move up and down, interpolate the difference on a common set of sample heights
                if DELTA:
                    SP_array_inst = setsOfData_inst[date][MERRA_VAR][t,:,0,0]
                    heights_inst = setsOfData_inst[date]['H'][t,:,0,0]
                    
                    SP_array_tavg = setsOfData_tavg[date][MERRA_VAR][t,:,0,0]
                    heights_tavg = setsOfData_tavg[date]['H'][t,:,0,0]
                    
                    #Remove high altitudes where they are both zero
                    nonzero_inst = SP_array_inst!=0
                    nonzero_heights_inst = heights_inst[nonzero_inst]
                    
                    nonzero_tavg = SP_array_tavg!=0
                    nonzero_heights_tavg =  heights_tavg[nonzero_tavg]
                    
                    message = "Number of nonzero instantaneous values: %i, time-averaged: %i" % (nonzero_heights_inst.shape[0],nonzero_heights_tavg.shape[0])
                    if PRINT_MESS: print(message)
                    messages.append(message)
                    
                    interpolater_inst = polate.interp1d(heights_inst,SP_array_inst)
                    interpolater_tavg = polate.interp1d(heights_tavg,SP_array_tavg)
                    
                    minHeight = min(np.amin(heights_inst),np.amin(heights_tavg))
                    
                    
                    if nonzero_heights_inst.shape[0]==0 and nonzero_heights_tavg.shape[0]==0:
                        lowerBound = max(np.amin(heights_tavg),np.amin(heights_inst))
                        upperBound = min(np.amax(heights_tavg),np.amax(heights_inst))
                    elif nonzero_heights_inst.shape[0]==0 and nonzero_heights_tavg.shape[0] != 0:
                        lowerBound = max(np.amin(nonzero_heights_tavg),np.amin(heights_inst))
                        upperBound = min(np.amax(nonzero_heights_tavg),np.amax(heights_inst))
                    elif nonzero_heights_inst.shape[0]!=0 and nonzero_heights_tavg.shape[0] == 0:
                        lowerBound = max(np.amin(nonzero_heights_inst),np.amin(heights_tavg))
                        upperBound = min(np.amax(nonzero_heights_inst),np.amax(heights_tavg))
                    else:
                        lowerBound = max(np.amin(heights_inst),np.amin(heights_tavg),min(np.amin(nonzero_heights_inst),np.amin(nonzero_heights_tavg))) #The lowest height where at least one is nonzero and both are defined for interpolation
                        upperBound = min(np.amax(heights_inst),np.amax(heights_tavg),max(np.amax(nonzero_heights_inst),np.amax(nonzero_heights_tavg))) #The highest height where at least one is nonzero and both are defined for interpolation
                    intervalLength = upperBound - lowerBound
                    
                    height_samples,height_interval = np.linspace(lowerBound,upperBound,NUM_OF_SAMPS,retstep = True)
                    
                    message = "The power spectrum is on %i heights going from %.4f m to %.4f m with spacing %.4f m." % (NUM_OF_SAMPS,lowerBound,upperBound,height_interval)
                    if PRINT_MESS: print(message)
                    messages.append(message)
                    
                    SP_array = interpolater_tavg(height_samples) - interpolater_inst(height_samples)
                    if ZERO_MEAN: SP_array = SP_array - np.mean(SP_array)
                
                else:
                    if time_type=='tavg':
                        SP_array = setsOfData_tavg[date][MERRA_VAR][t,:,0,0]
                        heights = setsOfData_tavg[date]['H'][t,:,0,0]
                    if time_type=='inst':
                        SP_array = setsOfData_inst[date][MERRA_VAR][t,:,0,0]
                        heights = setsOfData_inst[date]['H'][t,:,0,0]
                    
                    #Remove high altitudes beyond which it is always zero
                    
                    nonzero = SP_array!=0
                    nonzero_heights =  heights[nonzero]
                    
                    message = "Number of nonzero values: %i" % nonzero_heights.shape[0]
                    if PRINT_MESS: print(message)
                    messages.append(message)
                    
                    interpolater = polate.interp1d(heights,SP_array)
                    
                    minHeight = np.amin(heights)
                    
                    if nonzero_heights.shape[0]==0:
                        lowerBound = np.amin(heights)
                        upperBound = np.amax(heights)
                    elif nonzero_heights.shape[0] != 0:
                        lowerBound = np.amin(nonzero_heights)
                        upperBound = np.amax(nonzero_heights)

                    intervalLength = upperBound - lowerBound
                    
                    height_samples,height_interval = np.linspace(lowerBound,upperBound,NUM_OF_SAMPS,retstep = True)
                    
                    message = "The power spectrum is on %i heights going from %.4f m to %.4f m with spacing %.4f m." % (NUM_OF_SAMPS,lowerBound,upperBound,height_interval)
                    if PRINT_MESS: print(message)
                    messages.append(message)
                    
                    SP_array = interpolater(height_samples)
                    if ZERO_MEAN: SP_array = SP_array - np.mean(SP_array)
                
                sample_freqs,power = compute_1d_psd(SP_array,sample_spacing = height_interval)
                
                if MIN_FREQ != 0:
                    MIN_INV_FREQ = 1/MIN_FREQ
                else:
                    MIN_INV_FREQ = float('inf')
                if MAX_FREQ != float('inf'):
                    MAX_INV_FREQ = 1/MAX_FREQ
                else:
                    MAX_INV_FREQ = 0
                
                if MIN_INV_FREQ==float('inf') and MAX_INV_FREQ==0:

                    psd_integral = grate.trapz(power,sample_freqs)/(NUM_OF_SAMPS**2)*intervalLength #I don't know why you need to divide by (NUM_OF_SAMPS**2) and multiply by the interval length, but that made it work
        
                elif MIN_INV_FREQ==float('inf'):
        
                    valid_neg_freqs = (MAX_INV_FREQ<=np.abs(sample_freqs))&(sample_freqs<0)
                    valid_pos_freqs = (MAX_INV_FREQ<=np.abs(sample_freqs))&(sample_freqs>0)
                    
                    filtered_neg_power = power[valid_neg_freqs]
                    filtered_pos_power = power[valid_pos_freqs]
                    filtered_sample_neg_freqs = sample_freqs[valid_neg_freqs]
                    filtered_sample_pos_freqs = sample_freqs[valid_pos_freqs]
                    
                    
                    psd_integral = (grate.trapz(filtered_neg_power,filtered_sample_neg_freqs) + grate.trapz(filtered_pos_power,filtered_sample_pos_freqs))/(NUM_OF_SAMPS**2)*intervalLength
                
                elif MAX_INV_FREQ==0:
        
                    valid_freqs = np.abs(sample_freqs)<=MIN_INV_FREQ
        
                    filtered_power = power[valid_freqs]
                    filtered_sample_freqs = sample_freqs[valid_freqs]
                    
        
                    psd_integral = grate.trapz(filtered_power,filtered_sample_freqs)/(NUM_OF_SAMPS**2)*intervalLength
        
                else:
                    
                    valid_neg_freqs = (MIN_INV_FREQ>=np.abs(sample_freqs))&(np.abs(sample_freqs)>=MAX_INV_FREQ)&(sample_freqs<0)
                    valid_pos_freqs = (MIN_INV_FREQ>=np.abs(sample_freqs))&(np.abs(sample_freqs)>=MAX_INV_FREQ)&(sample_freqs>0)
                    #Now we need to add in the edge values that just barely weren't in the interval, since otherwise part of the power spectrum isn't included in any filter
                    firstNegValid = np.amin(np.where(valid_neg_freqs==True))
                    lastNegValid = np.amax(np.where(valid_neg_freqs==True))
                    
                    firstPosValid = np.amin(np.where(valid_pos_freqs==True))
                    lastPosValid = np.amax(np.where(valid_pos_freqs==True))
                    
                    if firstNegValid>0:
                        valid_neg_freqs[firstNegValid-1] = True
                    valid_neg_freqs[lastNegValid+1] = True
                    
                    valid_pos_freqs[firstPosValid-1] = True
                    if lastPosValid < sample_freqs.shape[0] - 1:
                        valid_pos_freqs[lastPosValid+1] = True
                    
                    filtered_neg_power = power[valid_neg_freqs]
                    filtered_pos_power = power[valid_pos_freqs]
                    filtered_sample_neg_freqs = sample_freqs[valid_neg_freqs]
                    filtered_sample_pos_freqs = sample_freqs[valid_pos_freqs]
                    
        
                    psd_integral = (grate.trapz(filtered_neg_power,filtered_sample_neg_freqs) + grate.trapz(filtered_pos_power,filtered_sample_pos_freqs))/(NUM_OF_SAMPS**2)*intervalLength
                
                merra_data[i] = psd_integral
            
            message = "Processed day %s, t=%i." % (date,t)
            if PRINT_MESS: print(message)
            messages.append(message)
            
            i += 1
    
    if T_INTERVAL==1 and FINAL_T_INTERVAL==3:
        merra_data = np.mean(merra_data.reshape(-1,3),axis=1)
    
    if MIN_FREQ==0 and MAX_FREQ==float('inf'):
        rangeString = ''
    if MIN_FREQ==0 and MAX_FREQ != float('inf'):
        rangeString = '<%im_' % MAX_FREQ
    if MIN_FREQ!=0 and MAX_FREQ==float('inf'):
        rangeString = '>%im_' % MIN_FREQ
    if MIN_FREQ != 0 and MAX_FREQ != float('inf'):
        rangeString = '%im-%im_' % (MIN_FREQ,MAX_FREQ)
    
    xLabel = MERRA_VAR
    if DIM==2 and DELTA:
        xLabel = "Delta_"+xLabel
    if DIM==2 and not DELTA:
        xLabel = time_type+'_'+xLabel
    if DIM==3 and DELTA:
        xLabel = "Variance_of_Delta_"+xLabel
    if DIM==3 and not DELTA:
        xLabel = "Variance_of_"+time_type+'_'+xLabel
    
    fileName = '%s%s_%s-%s_processing_results.txt' % (rangeString,xLabel,START_DATE.strftime('%Y%m%d'),END_DATE.strftime('%Y%m%d'))
    
    txtfile = open(os.path.join(folderPath,fileName),'w')
    messages.reverse()
    txtfile.write('\n'.join(messages))
    txtfile.close()
    
    return merra_data,rangeString,xLabel

def get_spt_data(SPT_VAR,T_INTERVAL,FOLDER_NAME,
    SPT_FREQ = 220, LOWEST_ELL = 110, HIGHEST_ELL = 250,
    RANDOM = False, PRINT_MESS = False,
    START_YEAR = 2019, START_MONTH=3,START_DAY=22,
    END_YEAR = 2019, END_MONTH = 12, END_DAY=31):
    
    #This creates the folder where everything will go
    dir_path = os.path.dirname(os.path.realpath(__file__))
    folderPath = os.path.join(dir_path,FOLDER_NAME)
    if not os.path.isdir(folderPath):
        os.mkdir(folderPath)
    
    messages = []
    NUM_OF_TIMES = 24 // T_INTERVAL
    #Generate the dates
    START_DATE = dt.date(year=START_YEAR,month=START_MONTH,day=START_DAY)
    END_DATE = dt.date(year=END_YEAR,month=END_MONTH,day=END_DAY)
    date = START_DATE
    interval = dt.timedelta(days=1)
    dates = []
    while date <= END_DATE:
        dates.append(date.isoformat().replace('-',''))
        date = date + interval
    
    message = "Loading SPT data..."
    if PRINT_MESS: print(message)
    messages.append(message)
    
    pklFileName = 'map_PSDs_index_isoformat_dates.pkl'
    infile = open(pklFileName,'rb')
    raw_data = pkl.load(infile)
    infile.close()
    
    ell_bins = raw_data.pop('ell_bins')
    SPT_time_strings = np.array(list(raw_data.keys()))
    
    Q = np.zeros(SPT_time_strings.shape[0])
    U = np.zeros(SPT_time_strings.shape[0])
    T = np.zeros(SPT_time_strings.shape[0])
    
    high_enough = ell_bins>=LOWEST_ELL
    low_enough = ell_bins<=HIGHEST_ELL
    bins_used = high_enough & low_enough
    
    for i,date in enumerate(raw_data):
        Qpower = raw_data[date]['Q'][SPT_FREQ]
        Q[i] = 10*np.sum(Qpower[bins_used])
        Upower = raw_data[date]['U'][SPT_FREQ]
        U[i] = 10*np.sum(Upower[bins_used])
        Tpower = raw_data[date]['T'][SPT_FREQ]
        T[i] = 10*np.sum(Tpower[bins_used])
    
    date_sort_indices = np.argsort(SPT_time_strings)
    SPT_time_strings = SPT_time_strings[date_sort_indices]
    Q = Q[date_sort_indices]
    U = U[date_sort_indices]
    T = T[date_sort_indices]
    
    if SPT_VAR=='T':
        SPT_data = T
    elif SPT_VAR=='Q':
        SPT_data = Q
    elif SPT_VAR=='U':
        SPT_data = U
    elif SPT_VAR=='QUratio':
        SPT_data = Q/U
    elif SPT_VAR=='Q-U':
        SPT_data = Q-U
    
    SPT_time_numbers = np.zeros(SPT_time_strings.shape[0])
    for i in range(SPT_time_strings.shape[0]): #I used a for loop instead of vectorizing the function because the function has the additional argumnets of START_DATE and END_DATE
        SPT_time_numbers[i] = convert2time(SPT_time_strings[i],START_DATE,END_DATE)
    
    #Remove data that is outside of the interval [START_DATE,END_DATE]
    data_to_keep = SPT_time_numbers>=0
    
    SPT_data = SPT_data[data_to_keep]
    SPT_time_numbers = SPT_time_numbers[data_to_keep]
    
    #Create an array with the SPT data at the same indices as the MERRA data will end up
    
    SPT_time_indices = (SPT_time_numbers // T_INTERVAL).astype(int)
    
    SPT_for_correlation = np.ma.array(np.zeros(NUM_OF_TIMES*len(dates)))
    SPT_for_correlation.mask = True
    for observationIndex,timeIndex in enumerate(SPT_time_indices):
        if SPT_for_correlation.mask[timeIndex]: #This makes it so that if there are several observations in a time interval, only the first goes into SPT_for_correlation
            SPT_for_correlation[timeIndex] = SPT_data[observationIndex]
    
    rangeString = '%sGHz_ell=%s-%s_' % (SPT_FREQ,LOWEST_ELL,HIGHEST_ELL)
    txtfile = open(os.path.join(folderPath,'%s%s_%s-%s_processing_results.txt' % (rangeString,SPT_VAR,START_DATE.strftime('%Y%m%d'),END_DATE.strftime('%Y%m%d'))),'w')
    messages.reverse()
    txtfile.write('\n'.join(messages))
    txtfile.close()
    
    yLabel = SPT_VAR
    return SPT_for_correlation,SPT_data,SPT_time_numbers,rangeString,yLabel

def process_data(x_for_correlation,x_data,x_time_numbers,y_for_correlation,y_data,y_time_numbers,
    xLabel,yLabel,STAT_TEST,FOLDER_NAME,
    colorLabel = 'None',colorRangeString='',color_data = False, #colorRangeString should include the percentile cutoff used.
    condition = 'None',conditionLabel = '',
    USE_RANK=False, INVY = False, ABSX = False, ABSY = False,LOGX = False, LOGY = False,
    xRangeString='',yRangeString='',
    RANDOM = False, PRINT_MESS = False,
    START_YEAR = 2019, START_MONTH=3,START_DAY=22,
    END_YEAR = 2019, END_MONTH = 12, END_DAY=31):
    
    if isinstance(condition,str):
        condition = np.ones(x_for_correlation.size,dtype=bool)
    
    messages = []
    START_DATE = dt.date(year=START_YEAR,month=START_MONTH,day=START_DAY)
    END_DATE = dt.date(year=END_YEAR,month=END_MONTH,day=END_DAY)
    date = START_DATE
    interval = dt.timedelta(days=1)
    dates = []
    while date <= END_DATE:
        dates.append(date.isoformat().replace('-',''))
        date = date + interval
    
    #This creates the folder where everything will go
    dir_path = os.path.dirname(os.path.realpath(__file__))
    folderPath = os.path.join(dir_path,FOLDER_NAME)
    if not os.path.isdir(folderPath):
        os.mkdir(folderPath)
    #Do statistical test and make scatter and time series graphs
    if USE_RANK:
        yLabel = 'rank(%s)' % yLabel
        xLabel = 'rank(%s)' % xLabel
        
    else:
        yLabel = '(%s)^-1' % yLabel if INVY else yLabel
        
        yLabel = 'abs(%s)' % yLabel if ABSY else yLabel
        yLabel = 'log(%s)' % yLabel if LOGY else yLabel
        
        if RANDOM:
            xLabel = 'random_var'
        xLabel = 'abs(%s)' % xLabel if ABSX else xLabel
        xLabel = 'log(%s)' % xLabel if LOGX else xLabel
        
        if INVY:
            y_for_correlation = 1/y_for_correlation
            y_data = 1/y_data
        if ABSY:
            y_for_correlation - np.abs(y_for_correlation)
            y_data = np.abs(y_data)
        
        if LOGY:
            y_for_correlation = np.log10(y_for_correlation)
            y_data = np.log10(y_data)
        
        if ABSX:
            x_for_correlation = np.abs(x_for_correlation)
            x_data = np.abs(x_data)
        
        if LOGX:
            x_for_correlation = np.log10(x_for_correlation)
            x_data = np.log10(x_data)

    #first get rid of the values which are masked because there wasn't an SPT reading from that interval, as well as any nan or infinite values
    #First deal with masked values
    invalidValues = np.zeros(x_for_correlation.size,dtype=bool)
    if isinstance(condition,np.ma.MaskedArray):
        invalidValues = invalidValues | condition.mask
    if isinstance(y_for_correlation,np.ma.MaskedArray):
        invalidValues = invalidValues | y_for_correlation.mask
    if isinstance(x_for_correlation,np.ma.MaskedArray):
        invalidValues = invalidValues | x_for_correlation.mask
    if isinstance(color_data,np.ma.MaskedArray):
        invalidValues = invalidValues | color_data
    #If there's a condition, take that into account
    if not isinstance(condition,str):
        invalidValues = invalidValues | ~condition
    #If there's color data, take that into account
    if not isinstance(color_data,bool):
        invalidValues = invalidValues | np.isnan(color_data) | np.isinf(color_data)
    invalidValues = invalidValues | np.isnan(y_for_correlation) | np.isinf(y_for_correlation) | np.isnan(x_for_correlation) | np.isinf(x_for_correlation)
    x_for_correlation = x_for_correlation[~invalidValues]
    y_for_correlation = y_for_correlation[~invalidValues]
    if not isinstance(color_data,bool):
        color_data = color_data[~invalidValues]
    
    if USE_RANK:
        y_for_correlation = sts.rankdata(y_for_correlation)
        x_for_correlation = sts.rankdata(x_for_correlation)

    if STAT_TEST=='linregress':

        regression = sts.linregress(x_for_correlation,y_for_correlation)
        slope = regression.slope
        intercept = regression.intercept
        r_sq = regression.rvalue**2
        pvalue = regression.pvalue
        message = "r^2: %f, slope: %f, intercept: %f, p-value: %f" % (r_sq,slope,intercept,pvalue)
        if PRINT_MESS: print(message)
        messages.append(message)
        graphTitle = "r^2: %.4f, p-value: %.4f" % (r_sq,pvalue)
        statNumber = r_sq
        
    elif STAT_TEST=='Spearman':
        
        rho, pvalue = sts.spearmanr(x_for_correlation,y_for_correlation)
        message = "rho: %f, p-value: %f" % (rho,pvalue)
        if PRINT_MESS: print(message)
        messages.append(message)
        graphTitle = "rho: %.4f, p-value: %.4f" % (rho,pvalue)
        statNumber = rho
    
    elif STAT_TEST=='ks':
        
        #Split the MERRA data into greater than average and less than average for the two sets of data
        aboveMedian = x_for_correlation>np.median(x_for_correlation)
        ks_stat,pvalue = sts.ks_2samp(y_for_correlation[aboveMedian],y_for_correlation[~aboveMedian])
        message = "KS stat: %f, p-value: %f" % (ks_stat,pvalue)
        if PRINT_MESS: print(message)
        messages.append(message)
        graphTitle = "KS stat: %.4f, p-value: %.4f" % (ks_stat,pvalue)
        statNumber = ks_stat
    
    fig, ax = plt.subplots()
    if not type(color_data) is bool:
        graphing_order = np.argsort(color_data)
        ordered_x = x_for_correlation[graphing_order]
        ordered_y = y_for_correlation[graphing_order]
        ordered_c = color_data[graphing_order]
        pos = ax.scatter(ordered_x,ordered_y,c=ordered_c,label='Data points')
        #fig.colorbar(pos,ax)
    else:
        ax.plot(x_for_correlation,y_for_correlation,'bo',label='Data points')
    
    if STAT_TEST=='linregress':
        ax.plot(x_for_correlation,intercept + slope*x_for_correlation,'r-',label = 'Best fit')
    if STAT_TEST=='ks':
        ax.vlines(np.median(x_for_correlation),0,1,transform = ax.get_xaxis_transform(),colors='r',label='KS divider')
    
    ax.set_title(graphTitle)
    ax.set_xlabel(xLabel)
    ax.set_ylabel(yLabel)
    
    ax.legend()
    ax.grid()
    fig.tight_layout()
    
    if not type(color_data) is bool:
        imageName = '%s%s%s_vs_%s%s_%s_%s-%s_color=%s%s_scatter.png' % (conditionLabel,yRangeString,yLabel,xRangeString,xLabel,STAT_TEST,START_DATE.strftime('%Y%m%d'),END_DATE.strftime('%Y%m%d'),colorRangeString,colorLabel)
    else:
        imageName = '%s%s%s_vs_%s%s_%s_%s-%s_scatter.png' % (conditionLabel,yRangeString,yLabel,xRangeString,xLabel,STAT_TEST,START_DATE.strftime('%Y%m%d'),END_DATE.strftime('%Y%m%d'))

    imagePath = os.path.join(folderPath,imageName)
    fig.savefig(imagePath)
    plt.close(fig)
    
    #Graph functions of time
    #Now get rid of invalid values for the time series data
    invalidValues = np.isnan(x_data) | np.isinf(x_data)
    
    x_data = x_data[~invalidValues]
    x_time_numbers = x_time_numbers[~invalidValues]
    
    invalidValues = np.isnan(y_data) | np.isinf(y_data)
    y_data = y_data[~invalidValues]
    y_time_numbers = y_time_numbers[~invalidValues]
    
    if USE_RANK:
        x_data = sts.rankdata(x_data)
        y_data = sts.rankdata(y_data)
    
    fig, ax1 = plt.subplots()
    
    ax1.plot(x_time_numbers,x_data,'b-',label=xLabel)
    ax1.plot([],[],'r-',label=yLabel)
    ax1.set_xlabel('Time [hrs]')
    ax1.set_ylabel(xLabel)
    ax1.set_title("%s and %s vs. Time" % (xLabel,yLabel))
    ax1.legend()
    
    ax2 = ax1.twinx()
    ax2.plot(y_time_numbers,y_data,'r-',label=yLabel)
    ax2.set_ylabel(yLabel)
    
    fig.tight_layout()
    
    imageName = '%s%s_vs_%s%s_%s-%s_time_series.png' % (yRangeString,yLabel,xRangeString,xLabel,START_DATE.strftime('%Y%m%d'),END_DATE.strftime('%Y%m%d'))
    imagePath = os.path.join(folderPath,imageName)
    fig.savefig(imagePath)
    plt.close(fig)
    
    txtfile = open(os.path.join(folderPath,'%s%s_vs_%s%s_%s_%s-%s_results.txt' % (yRangeString,yLabel,xRangeString,xLabel,STAT_TEST,START_DATE.strftime('%Y%m%d'),END_DATE.strftime('%Y%m%d'))),'w')
    messages.reverse()
    txtfile.write('\n'.join(messages))
    txtfile.close()
    
    return yRangeString,yLabel,xRangeString,xLabel,statNumber,pvalue

def make_model(xLabels,xRangeStrings,xDataTuple,yLabel,yRangeString,yData,FOLDER_NAME,
    PRINT_MESS = False,hiddenLayers=[],epoch_num = 200,
    USE_RANK = False,ABSX=False,LOGX=False,LOGY=False):
    
    #Because there is an element of randomness in the models, it includes the time in the file names so multiple attempts don't overwrite one another
    timeStamp = dt.datetime.now().isoformat()
    messages = []
    dir_path = os.path.dirname(os.path.realpath(__file__))
    folderPath = os.path.join(dir_path,FOLDER_NAME)
    if not os.path.isdir(folderPath):
        os.mkdir(folderPath)
    
    features = np.column_stack(xDataTuple)
    labels = yData
    if USE_RANK:
        features = sts.rankdata(features)
        labels = sts.rankdata(labels)
        for i in range(len(xLabels)):
            xLabels[i] = 'rank(%s)' % xLabels[i]
        yLabel = 'rank(%s)' % yLabel
    else:
        if ABSX:
            features = np.abs(features)
            for i in range(len(xLabels)):
                xLabels[i] = 'abs('+xLabels[i]+')'
        if LOGX: 
            features = np.log10(features)
            for i in range(len(merraVars)):
                merraVars[i] = 'log('+merraVars[i]+')'
        if LOGY:
            labels = np.log10(labels)
            SPT_VAR = 'log('+SPT_VAR+')'
        
    #Eliminate invalid values
    invalidIndices = np.any(np.isnan(features),axis=1) | np.any(np.isinf(features),axis=1) | np.isnan(labels) | np.isinf(labels)
    
    features = features[~invalidIndices]
    labels = labels[~invalidIndices]
    
    #Assign some of the data for testing instead of training
    test_indices = np.sort(np.random.choice(labels.shape[0],int(labels.shape[0]/5),replace=False))

    train_features = features[~test_indices]
    train_labels = labels[~test_indices]
    test_features = features[test_indices]
    test_labels = labels[test_indices]
    
    #Create normalizer
    normalizer = preprocessing.Normalization()
    normalizer.adapt(train_features)
    
    #Bring together layers
    layerList = [normalizer]
    for neuronNumber in hiddenLayers:
        layerList.append(layers.Dense(neuronNumber,activation='relu'))
    layerList.append(layers.Dense(1))
    
    #Build and compile model
    model = keras.Sequential(layerList)
    model.compile(loss='mean_absolute_error', optimizer=tf.keras.optimizers.Adam(0.001))
    
    #Fit model
    history = model.fit(
        train_features, train_labels, 
        epochs=epoch_num,
        verbose=1, #Logs results
        # Calculate validation results on 20% of the training data
        validation_split = 0.2)
    
    #Save model
    for i in range(len(xLabels)):
        xLabels[i] = xRangeStrings[i]+xLabels[i]
    independentVariableNames = ','.join(xLabels)
    if len(hiddenLayers)==0:
        architectureString = "linear"
    else:
        architectureString = "hidden_"+str(hiddenLayers)
    
    modelFileName = '%s%svs%s_%s_model_%s' % (yRangeString,yLabel,independentVariableNames,architectureString,timeStamp)
    model.save(os.path.join(folderPath,modelFileName))
    
    #Graph history
    fig, ax = plt.subplots()
    ax.plot(history.history['loss'],label='loss')
    ax.plot(history.history['val_loss'],label='validation loss')
    ax.set_xlabel('Epoch')
    ax.set_ylabel('Error')
    ax.legend()
    ax.grid()
    
    imageName = '%s%svs%s_%s_history_%s.png' % (yRangeString,yLabel,independentVariableNames,architectureString,timeStamp)
    imagePath = os.path.join(folderPath,imageName)
    fig.savefig(imagePath)
    plt.close(fig)
    
    #Test on training data and graph results
    prediction = model.predict(test_features)
    fig, ax = plt.subplots()
    ax.scatter(test_labels,prediction)
    ax.set_xlabel("True Value")
    ax.set_ylabel("Model Prediction")
    
    #Prediction is has shape (n,1) so flatten it
    prediction = prediction.flatten()
    
    #Adjust axes
    lowerLimit = min(np.amin(test_labels),np.amin(prediction))
    upperLimit = max(np.amax(test_labels),np.amax(prediction))
    
    ax.set_xlim(lowerLimit,upperLimit)
    ax.set_ylim(lowerLimit,upperLimit)
    
    regression = sts.linregress(test_labels,prediction)
    slope = regression.slope
    intercept = regression.intercept
    r_sq = regression.rvalue**2
    pvalue = regression.pvalue
    
    message = "r^2: %f, slope: %f, intercept: %f, p-value: %f" % (r_sq,slope,intercept,pvalue)
    if PRINT_MESS: print(message)
    messages.append(message)
    
    graphTitle = "r^2: %.4f" % r_sq
    ax.set_title(graphTitle)
    imageName = '%s%svs%s_%s_test_scatter_%s.png' % (LyRangeString,yLabel,independentVariableNames,architectureString,timeStamp)
    imagePath = os.path.join(folderPath,imageName)
    fig.savefig(imagePath)
    plt.close(fig)
    
    txtfile = open(os.path.join(folderPath,'%s%svs%s_%s__results_%s.txt' % (yRangeString,yLabel,independentVariableNames,architectureString,timeStamp)),'w')
    messages.reverse()
    txtfile.write('\n'.join(messages))
    txtfile.close()
    
    #Graph the model's function if there's just one independent variable
    if len(xDataTuple)==1:
        x = np.linspace(np.amin(features),np.amax(features),num=200)
        y = model.predict(x)
        fig,ax = plt.subplots()
        ax.scatter(features,labels,color='b',label='data')
        ax.plot(x,y,color='r',label='prediction')
        ax.set_xlabel(independentVariableNames)
        ax.set_ylabel(SPT_VAR)
        ax.legend()
        ax.grid()
        
        imageName = '%s%svs%s_%s__function_%s.png' % (yRangeString,yLabel,independentVariableNames,architectureString,timeStamp)
        imagePath = os.path.join(folderPath,imageName)
        fig.savefig(imagePath)
        plt.close(fig)
    
    return r_sq
if __name__ == '__main__':
    
    parser.add_argument('-fbi2','--file_base_instantaneous_2d',default='MERRA2_400.inst1_2d_asm_Nx.')
    parser.add_argument('-fei2','--file_end_instantaneous_2d',default='.SUB.nc')
    parser.add_argument('-fbt2','--file_base_time_averaged_2d',default='MERRA2_400.tavg1_2d_slv_Nx.')
    parser.add_argument('-fet2','--file_end_time_averaged_2d',default='.SUB.nc')
    
    parser.add_argument('-fbi3','--file_base_instantaneous_3d',default='MERRA2_400.inst3_3d_asm_Nv.')
    parser.add_argument('-fei3','--file_end_instantaneous_3d',default='.SUB.nc')
    parser.add_argument('-fbt3','--file_base_time_averaged_3d',default='MERRA2_400.tavg3_3d_asm_Nv.')
    parser.add_argument('-fet3','--file_end_time_averaged_3d',default='.SUB.nc')
    
    parser.add_argument('-fn','--folder_name',help='The name of the folder where the program\'s output is stored')
    
    parser.add_argument('-mv','--merra_variable',help='The name of the MERRA-2 variable you wish to use')
    parser.add_argument('-cv','--color_variable',default="None")
    parser.add_argument('-cd','--color_dimension',type=int,default=2)
    parser.add_argument('-cl','--color_delta',action='store_true')
    parser.add_argument('-sv','--spt_variable',help='The SPT noise variable you want',choices=SPTchoices)
    parser.add_argument('-d','--dimensions',help='Either the 2- or 3-D MERRA-2 Data',type=int,choices=[2,3])
    parser.add_argument('-ff','--frequency_filter',choices = filters,default='None')
    parser.add_argument('-sf','--spt_frequency',type=int,choices=[90,150,220],default=220)
    parser.add_argument('-le','--lowest_ell',type=float,default=110)
    parser.add_argument('-he','--highest_ell',type=float,default=250)
    
    parser.add_argument('-r','--random',action='store_true')
    parser.add_argument('-mi','--make_images',action='store_true')
    parser.add_argument('-ns','--number_of_samples',type=int,default=300)
    parser.add_argument('-pm','--print_messages',action='store_true')
    
    parser.add_argument('-st','--statistical_test',choices=testChoices)
    parser.add_argument('-delta','--delta',action='store_true')
    parser.add_argument('-invy','--inverse_y',action='store_true')
    parser.add_argument('-logx','--logarithm_x',action='store_true')
    parser.add_argument('-logy','--logarithm_y',action='store_true')
    
    parser.add_argument('-sy','--start_year',type=int,default=2019)
    parser.add_argument('-sm','--start_month',type=int,choices=monthChoices,default=3)
    parser.add_argument('-sd','--start_day',type=int,choices=dayChoices,default=22)
    
    parser.add_argument('-ey','--end_year',type=int,default=2019)
    parser.add_argument('-em','--end_month',type=int,choices=monthChoices,default=12)
    parser.add_argument('-ed','--end_day',type=int,choices=dayChoices,default=31)
    
    args = parser.parse_args()
    
    FILE_BASE_INST2 = args.file_base_instantaneous_2d
    FILE_END_INST2 = args.file_end_instantaneous_2d
    FILE_BASE_TAVG2 = args.file_base_time_averaged_2d
    FILE_END_TAVG2 = args.file_end_time_averaged_2d
    
    FILE_BASE_INST3 = args.file_base_instantaneous_3d
    FILE_END_INST3 = args.file_end_instantaneous_3d
    FILE_BASE_TAVG3 = args.file_base_time_averaged_3d
    FILE_END_TAVG3 = args.file_end_time_averaged_3d
    
    START_YEAR = args.start_year
    START_MONTH = args.start_month
    START_DAY = args.start_day
    END_YEAR = args.end_year
    END_MONTH = args.end_month
    END_DAY = args.end_day
    
    FOLDER_NAME = args.folder_name #The name of the folder where the images are stored, which will be created if it doesn't already exist
    
    MERRA_VAR = args.merra_variable #Which meteorological quantity we're checking
    COLOR_VAR = args.color_variable

    RANDOM = args.random
    SPT_VAR = args.spt_variable
    COLOR_DIM = args.color_dimension
    DIM = args.dimensions
    
    FREQ_FILTER = args.frequency_filter
    SPT_FREQ = args.spt_frequency
    LOWEST_ELL = args.lowest_ell
    HIGHEST_ELL = args.highest_ell
    
    MAKE_IMAGES = args.make_images
    NUM_OF_SAMPS = args.number_of_samples
    PRINT_MESS = args.print_messages
    
    DELTA = args.delta
    COLOR_DELTA = args.color_delta
    STAT_TEST = args.statistical_test
    INVY = args.inverse_y
    LOGX = args.logarithm_x
    LOGY = args.logarithm_y
    
    color_data = False
    if not COLOR_VAR=='None':
        color_data = get_merra_data(COLOR_VAR,COLOR_DIM,FOLDER_NAME,
        FILE_BASE_INST2, FILE_END_INST2, FILE_BASE_TAVG2, FILE_END_TAVG2,
        FILE_BASE_INST3, FILE_END_INST3, FILE_BASE_TAVG3, FILE_END_TAVG3,
        FREQ_FILTER, COLOR_DELTA,
        MAKE_IMAGES, NUM_OF_SAMPS, PRINT_MESS,
        START_YEAR, START_MONTH, START_DAY,
        END_YEAR, END_MONTH, END_DAY)
        #Just mark the top few %
        color_data = (color_data>np.percentile(color_data,80)).astype(int)
    
    merra_data = get_merra_data(MERRA_VAR,DIM,FOLDER_NAME,
        FILE_BASE_INST2, FILE_END_INST2, FILE_BASE_TAVG2, FILE_END_TAVG2,
        FILE_BASE_INST3, FILE_END_INST3, FILE_BASE_TAVG3, FILE_END_TAVG3,
        FREQ_FILTER, DELTA,
        MAKE_IMAGES, NUM_OF_SAMPS, PRINT_MESS,
        START_YEAR, START_MONTH, START_DAY,
        END_YEAR, END_MONTH, END_DAY)

    SPT_for_correlation,SPT_data,SPT_time_numbers = get_spt_data(SPT_VAR,DIM,FOLDER_NAME,
        SPT_FREQ, LOWEST_ELL, HIGHEST_ELL,
        RANDOM, PRINT_MESS,
        START_YEAR, START_MONTH, START_DAY,
        END_YEAR, END_MONTH, END_DAY)

    process_data(merra_data,SPT_for_correlation,SPT_data,SPT_time_numbers,
        SPT_VAR,MERRA_VAR,DIM,STAT_TEST,FOLDER_NAME,
        COLOR_VAR,color_data, INVY, LOGX, LOGY, DELTA,
        FREQ_FILTER,SPT_FREQ,LOWEST_ELL,HIGHEST_ELL,
        RANDOM, MAKE_IMAGES, NUM_OF_SAMPS, PRINT_MESS,
        START_YEAR, START_MONTH,START_DAY,
        END_YEAR, END_MONTH, END_DAY)