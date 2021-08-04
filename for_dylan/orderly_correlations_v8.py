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
from PIL import Image
import pickle as pkl
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.layers.experimental import preprocessing
#Some of the machine learning code is taken from a TensorFlow tutorial here: https://www.tensorflow.org/tutorials/keras/regression

ZERO_MEAN = False

def plot_loss(history): #to see the progress of the ML training
  plt.plot(history.history['loss'], label='loss')
  plt.plot(history.history['val_loss'], label='val_loss')
  plt.xlabel('Epoch')
  plt.ylabel('Error')
  plt.legend()
  plt.grid(True)

def compute_1d_psd(input_array, sample_spacing = 1.0): #This function is defined here intead of imported from demonstrate_psd.py because I figured it's better if this program can run without having to remember to move over another program as well.
    #Does what it says on the tin
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
    MIN_FREQ = 0,MAX_FREQ=float('inf'),DELTA = True,time_type='tavg', RETURN_3D_PROFILES = False,
    NUM_OF_SAMPS=300, PRINT_MESS = False,
    START_YEAR = 2019, START_MONTH=3,START_DAY=22,
    END_YEAR = 2019, END_MONTH = 12, END_DAY=31):
    
    #This creates the folder where everything will go
    dir_path = os.path.dirname(os.path.realpath(__file__))
    folderPath = os.path.join(dir_path,FOLDER_NAME)
    if not os.path.isdir(folderPath):
        os.mkdir(folderPath)
    
    messages = []
    #3D MERRA-2 data is on 3-hour time intervals, whereas 2D MERRA-2 data is on 1-hour time intervals.
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
    if RETURN_3D_PROFILES:
        profiles_3d = []
        heights_3d = []
        rangeString = ''
    i = 0
    for date in dates:
        for t in range(NUM_OF_TIMES):
            #Gets the appropriate dataset for the dimension and whether it's the instantaneous data, the time-averaged data, or their difference
            if DIM==2:
                if MERRA_VAR == 'WIND2M':
                    data_tavg = np.sqrt(np.power(setsOfData_tavg[date]['U2M'][t,0,0],2)+np.power(setsOfData_tavg[date]['V2M'][t,0,0],2))
                    data_inst = np.sqrt(np.power(setsOfData_inst[date]['U2M'][t,0,0],2)+np.power(setsOfData_inst[date]['V2M'][t,0,0],2))
                elif MERRA_VAR == 'WIND10M':
                    data_tavg = np.sqrt(np.power(setsOfData_tavg[date]['U10M'][t,0,0],2)+np.power(setsOfData_tavg[date]['V10M'][t,0,0],2))
                    data_inst = np.sqrt(np.power(setsOfData_inst[date]['U10M'][t,0,0],2)+np.power(setsOfData_inst[date]['V10M'][t,0,0],2))
                elif MERRA_VAR == 'WIND50M':
                    data_tavg = np.sqrt(np.power(setsOfData_tavg[date]['U50M'][t,0,0],2)+np.power(setsOfData_tavg[date]['V50M'][t,0,0],2))
                    data_inst = np.sqrt(np.power(setsOfData_inst[date]['U50M'][t,0,0],2)+np.power(setsOfData_inst[date]['V50M'][t,0,0],2))
                else:
                    data_tavg = setsOfData_tavg[date][MERRA_VAR][t,0,0]
                    data_inst = setsOfData_inst[date][MERRA_VAR][t,0,0]
                if DELTA:
                    merra_data[i] = data_tavg - data_inst
                elif time_type=='tavg':
                    merra_data[i] = data_tavg
                elif time_type=='inst':
                    merra_data[i] = data_inst
            elif DIM==3:
                #Load the data and then, because the pressure levels may move up and down, interpolate the difference on a common set of sample heights
                if DELTA:
                    if MERRA_VAR == 'TOTAL_WIND':
                        SP_array_inst = np.sqrt(np.power(setsOfData_inst[date]['U'][t,:,0,0],2) + np.power(setsOfData_inst[date]['V'][t,:,0,0],2))
                        SP_array_tavg = np.sqrt(np.power(setsOfData_tavg[date]['U'][t,:,0,0],2) + np.power(setsOfData_tavg[date]['V'][t,:,0,0],2))
                    else:
                        SP_array_inst = setsOfData_inst[date][MERRA_VAR][t,:,0,0]
                        SP_array_tavg = setsOfData_tavg[date][MERRA_VAR][t,:,0,0]

                    heights_inst = setsOfData_inst[date]['H'][t,:,0,0] #heights in meters corresponding to the pressure levels of the MERRA-2 grid
                    heights_tavg = setsOfData_tavg[date]['H'][t,:,0,0]
                    
                    #Remove high altitudes where they are both zero
                    nonzero_inst = SP_array_inst!=0
                    nonzero_heights_inst = heights_inst[nonzero_inst]
                    
                    nonzero_tavg = SP_array_tavg!=0
                    nonzero_heights_tavg =  heights_tavg[nonzero_tavg]
                    
                    message = "Number of nonzero instantaneous values: %i, time-averaged: %i" % (nonzero_heights_inst.shape[0],nonzero_heights_tavg.shape[0])
                    if PRINT_MESS: print(message)
                    messages.append(message)
                    
                    interpolater_inst = polate.interp1d(heights_inst,SP_array_inst,fill_value='extrapolate')
                    interpolater_tavg = polate.interp1d(heights_tavg,SP_array_tavg,fill_value='extrapolate')
                    
                    minHeight = min(np.amin(heights_inst),np.amin(heights_tavg))
                    
                    #Deal with edge cases where sometime one of the arrays is all zeros
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
                    
                    #Generate the even samples for interpolation, which is then used for the 1D power spectrum
                    height_samples,height_interval = np.linspace(lowerBound,upperBound,NUM_OF_SAMPS,retstep = True)
                    
                    message = "The power spectrum is on %i heights going from %.4f m to %.4f m with spacing %.4f m." % (NUM_OF_SAMPS,lowerBound,upperBound,height_interval)
                    if PRINT_MESS: print(message)
                    messages.append(message)
                    
                    SP_array = interpolater_tavg(height_samples) - interpolater_inst(height_samples)
                    if ZERO_MEAN: SP_array = SP_array - np.mean(SP_array)
                
                else:
                    #This works similarly to the above, but without taking the difference of tavg and inst
                    if time_type=='tavg':
                        if MERRA_VAR == 'TOTAL_WIND':
                            SP_array = np.sqrt(np.power(setsOfData_tavg[date]['U'][t,:,0,0],2) + np.power(setsOfData_tavg[date]['V'][t,:,0,0],2))
                        else:
                            SP_array = setsOfData_tavg[date][MERRA_VAR][t,:,0,0]
                        heights = setsOfData_tavg[date]['H'][t,:,0,0]
                    if time_type=='inst':
                        if MERRA_VAR == 'TOTAL_WIND':
                            SP_array = np.sqrt(np.power(setsOfData_inst[date]['U'][t,:,0,0],2) + np.power(setsOfData_inst[date]['V'][t,:,0,0],2))
                        else:
                            SP_array = setsOfData_inst[date][MERRA_VAR][t,:,0,0]
                        heights = setsOfData_inst[date]['H'][t,:,0,0]
                    
                    #Remove high altitudes beyond which it is always zero
                    
                    nonzero = SP_array!=0
                    nonzero_heights =  heights[nonzero]
                    
                    message = "Number of nonzero values: %i" % nonzero_heights.shape[0]
                    if PRINT_MESS: print(message)
                    messages.append(message)
                    
                    interpolater = polate.interp1d(heights,SP_array,fill_value='extrapolate')
                    
                    minHeight = np.amin(heights)
                    
                    if nonzero_heights.size<=1:
                        lowerBound = np.amin(heights)
                        upperBound = np.amax(heights)
                    else:
                        lowerBound = np.amin(nonzero_heights)
                        upperBound = np.amax(nonzero_heights)
                    
                    intervalLength = upperBound - lowerBound
                    
                    height_samples,height_interval = np.linspace(lowerBound,upperBound,NUM_OF_SAMPS,retstep = True)
                    
                    message = "The power spectrum is on %i heights going from %.4f m to %.4f m with spacing %.4f m." % (NUM_OF_SAMPS,lowerBound,upperBound,height_interval)
                    if PRINT_MESS: print(message)
                    messages.append(message)
                    
                    SP_array = interpolater(height_samples)
                    if ZERO_MEAN: SP_array = SP_array - np.mean(SP_array)
                
                if RETURN_3D_PROFILES:
                    profiles_3d.append(SP_array)
                    heights_3d.append(height_samples)
                else:
                    sample_freqs,power = compute_1d_psd(SP_array,sample_spacing = height_interval)
                    
                    #Convert the ranges in meters to inverse meters, and then only integrate over the inverse distances which are within the range
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
    
                        if np.where(valid_neg_freqs==True)[0].size != 0:
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
                        else:
                            psd_integral = 0
                    
                    merra_data[i] = psd_integral
            
            message = "Processed day %s, t=%i." % (date,t)
            if PRINT_MESS: print(message)
            messages.append(message)
            
            i += 1
    
    
    if not RETURN_3D_PROFILES:
        if T_INTERVAL==1 and FINAL_T_INTERVAL==3:
            merra_data = np.mean(merra_data.reshape(-1,3),axis=1) #Averages out 1-hour data in order to compare it with 3-hour data
        
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

    if RETURN_3D_PROFILES:
        xLabel = MERRA_VAR
        if DELTA:
            xLabel = 'Delta_'+xLabel
        else:
            xLabel = time_type+'_'+xLabel

    fileName = '%s%s_%s-%s_processing_results.txt' % (rangeString,xLabel,START_DATE.strftime('%Y%m%d'),END_DATE.strftime('%Y%m%d'))
    
    txtfile = open(os.path.join(folderPath,fileName),'w')
    messages.reverse()
    txtfile.write('\n'.join(messages))
    txtfile.close()
    if RETURN_3D_PROFILES:
        return profiles_3d,heights_3d
    else:
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
    
    #Opens the pickle file and processes the arrays
    pklFileName = 'map_PSDs_index_isoformat_dates.pkl'
    infile = open(pklFileName,'rb')
    raw_data = pkl.load(infile)
    infile.close()
    
    ell_bins = raw_data.pop('ell_bins')
    SPT_time_strings = np.array(list(raw_data.keys()))
    
    Q = np.zeros(SPT_time_strings.shape[0])
    U = np.zeros(SPT_time_strings.shape[0])
    T = np.zeros(SPT_time_strings.shape[0])
    
    #Integrates only over the specified ell-ranges
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
    for i in range(SPT_time_strings.shape[0]): #I used a for loop instead of vectorizing the function because the function has the additional arguments of START_DATE and END_DATE
        SPT_time_numbers[i] = convert2time(SPT_time_strings[i],START_DATE,END_DATE)
    
    #Remove data that is outside of the interval [START_DATE,END_DATE]
    data_to_keep = SPT_time_numbers>=0 #Dates which are after END_DATE map to -1, so they're also cut  by this procedure
    
    SPT_data = SPT_data[data_to_keep]
    SPT_time_numbers = SPT_time_numbers[data_to_keep]
    
    #Create an array with the SPT data at the same indices as the MERRA data will end up
    
    SPT_time_indices = (SPT_time_numbers // T_INTERVAL).astype(int)
    
    SPT_for_correlation = np.ma.array(np.zeros(NUM_OF_TIMES*len(dates)))
    SPT_for_correlation.mask = True
    #Only unmasks if there is indeed an SPT observation in that time bin.
    #This will average all observations falling within a particular bin
    for timeIndex in SPT_time_indices:
        in_the_bin = (SPT_time_indices==timeIndex)
        SPT_for_correlation[timeIndex] = np.mean(SPT_data[in_the_bin])
    
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
    RANDOM = False, PRINT_MESS = False, TIME_STREAM = False, ImageQuality=95, #1 to 95
    ImageType = 'png',
    START_YEAR = 2019, START_MONTH=3,START_DAY=22,
    END_YEAR = 2019, END_MONTH = 12, END_DAY=31):
    
    if isinstance(condition,str):
        condition = np.ones(x_for_correlation.size,dtype=bool)
        #If condition is 'None' then this makes it so that the condition is true everywhere, so that every data point is included
    
    messages = []
    #Generate the dates
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
    
    #Do any functions on the x and y data that are specified, with rank and log not mixing, cause rank(log(x)) = rank(x) and log(rank(x)) is sort of meaningless
    if ABSY:
        y_for_correlation = np.abs(y_for_correlation)
        y_data = np.abs(y_data)
        yLabel = 'abs(%s)' % yLabel
    if ABSX:
        x_for_correlation = np.abs(x_for_correlation)
        x_data = np.abs(x_data)
        xLabel = 'abs(%s)' % xLabel if ABSX else xLabel
    if USE_RANK:
        yLabel = 'rank(%s)' % yLabel
        xLabel = 'rank(%s)' % xLabel
        
    else:
        yLabel = '(%s)^-1' % yLabel if INVY else yLabel
        yLabel = 'log(%s)' % yLabel if LOGY else yLabel
        
        if RANDOM:
            xLabel = 'random_var'
        xLabel = 'log(%s)' % xLabel if LOGX else xLabel
        
        if INVY:
            y_for_correlation = 1/y_for_correlation
            y_data = 1/y_data
        
        if LOGY:
            y_for_correlation = np.log10(y_for_correlation)
            y_data = np.log10(y_data)
        
        if LOGX:
            x_for_correlation = np.log10(x_for_correlation)
            x_data = np.log10(x_data)

    #first get rid of the values which are masked because there wasn't an SPT reading from that interval, as well as any nan or infinite values

    invalidValues = np.zeros(x_for_correlation.size,dtype=bool)
    #First deal with masked values
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
    invalidValues = invalidValues | np.isnan(y_for_correlation) | np.isinf(y_for_correlation) | np.isnan(x_for_correlation) | np.isinf(x_for_correlation) #Get rid of nan or inf values
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
    if not type(color_data) is bool: #Only if color_data is provided as an array
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
        imageName = '%s%s%s_vs_%s%s_%s_%s-%s_color=%s%s_scatter.%s' % (conditionLabel,yRangeString,yLabel,xRangeString,xLabel,STAT_TEST,START_DATE.strftime('%Y%m%d'),END_DATE.strftime('%Y%m%d'),colorRangeString,colorLabel,ImageType)
    else:
        imageName = '%s%s%s_vs_%s%s_%s_%s-%s_scatter.%s' % (conditionLabel,yRangeString,yLabel,xRangeString,xLabel,STAT_TEST,START_DATE.strftime('%Y%m%d'),END_DATE.strftime('%Y%m%d'),ImageType)

    imagePath = os.path.join(folderPath,imageName)
    fig.savefig(imagePath,quality=ImageQuality)
    plt.close(fig)
    
    #Graph functions of time
    if TIME_STREAM:
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
        
        imageName = '%s%s_vs_%s%s_%s-%s_time_series.%s' % (yRangeString,yLabel,xRangeString,xLabel,START_DATE.strftime('%Y%m%d'),END_DATE.strftime('%Y%m%d'),ImageType)
        imagePath = os.path.join(folderPath,imageName)
        fig.savefig(imagePath,quality=ImageQuality)
        plt.close(fig)
    
    txtfile = open(os.path.join(folderPath,'%s%s_vs_%s%s_%s_%s-%s_results.txt' % (yRangeString,yLabel,xRangeString,xLabel,STAT_TEST,START_DATE.strftime('%Y%m%d'),END_DATE.strftime('%Y%m%d'))),'w')
    messages.reverse()
    txtfile.write('\n'.join(messages))
    txtfile.close()
    
    return yRangeString,yLabel,xRangeString,xLabel,statNumber,pvalue

def make_model(xLabels,xRangeStrings,xDataList,yLabel,yRangeString,yData,FOLDER_NAME,
    PRINT_MESS = False,hiddenLayers=[],epoch_num = 100, ImageQuality=95, #1 to 95
    ImageType = 'png',ABSX=False,ABSY=False,LOGX=False,LOGY=False,runID=''):
    
    #Because there is an element of randomness in the models, it includes a distinct runID so multiple attempts don't overwrite one another
    messages = []
    dir_path = os.path.dirname(os.path.realpath(__file__))
    folderPath = os.path.join(dir_path,FOLDER_NAME)
    if not os.path.isdir(folderPath):
        os.mkdir(folderPath)
    
    if ABSX:
        for i in range(len(xLabels)):
            xDataList[i] = np.abs(xDataList[i])
            xLabels[i] = 'abs('+xLabels[i]+')'
    if ABSY:
        labels = np.abs(labels)
        yLabel = 'abs(%s)' % yLabel
    
    if LOGX:
        for i in range(len(xLabels)):
            xDataList[i] = np.log10(xDataList[i])
            xLabels[i] = 'log('+xLabels[i]+')'
    if LOGY:
        yData = np.log10(yData)
        yLabel = 'log('+yLabel+')'
    if np.ma.is_masked(yData):
        invalidValues = np.isnan(yData) | np.isinf(yData) | yData.mask
    else:
        invalidValues = np.isnan(yData) | np.isinf(yData)
    for xData in xDataList:
        invalidValues = invalidValues | np.isnan(xData) | np.isinf(xData)
    
    yData = yData[~invalidValues]
    for i in range(len(xDataList)):
        print(xData.size)
        xDataList[i] = xDataList[i][~invalidValues]
        print(xData.size)
    
    #It tries to predict yData based on xDataTuple
    features = np.column_stack(tuple(xDataList))
    labels = yData
    
    print(labels.size)
    print(features.size)
    #Assign some of the data for testing instead of training
    pool = np.array([0,0,0,0,1],dtype=bool)
    test_indices = np.random.choice(pool,size=labels.size,replace=True)
    #test_indices = np.sort(np.random.choice(labels.shape[0],int(labels.shape[0]/5),replace=False))

    train_features = features[~test_indices]
    train_labels = labels[~test_indices]
    test_features = features[test_indices]
    test_labels = labels[test_indices]
    
    print(train_features.shape,train_labels.shape)
    print(test_features.shape,test_labels.shape)
    
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
    for i in range(len(xDataList)):
        xLabels[i] = xRangeStrings[i]+xLabels[i]
    independentVariableNames = ','.join(xLabels)
    if len(hiddenLayers)==0:
        architectureString = "linear"
    else:
        architectureString = "hidden_"+str(hiddenLayers)
    
    modelFileName = '%s%svs%s_%s_model_%s' % (yRangeString,yLabel,independentVariableNames,architectureString,runID)
    model.save(os.path.join(folderPath,modelFileName))
    
    #Graph history
    fig, ax = plt.subplots()
    ax.plot(history.history['loss'],label='loss')
    ax.plot(history.history['val_loss'],label='validation loss')
    ax.set_xlabel('Epoch')
    ax.set_ylabel('Error')
    ax.legend()
    ax.grid()
    
    imageName = '%s%svs%s_%s_history_%s.%s' % (yRangeString,yLabel,independentVariableNames,architectureString,runID,ImageType)
    imagePath = os.path.join(folderPath,imageName)
    fig.savefig(imagePath,quality=ImageQuality)
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
    
    rho,unusedp = sts.spearmanr(test_labels,prediction)
    
    message = "rho: %f; r^2: %f, slope: %f, intercept: %f, p-value: %f" % (rho,r_sq,slope,intercept,pvalue)
    if PRINT_MESS: print(message)
    messages.append(message)
    
    graphTitle = "rho: %.4f, r^2: %.4f" % (rho,r_sq)
    ax.set_title(graphTitle)
    imageName = '%s%svs%s_%s_test_scatter_%s.%s' % (yRangeString,yLabel,independentVariableNames,architectureString,runID,ImageType)
    imagePath = os.path.join(folderPath,imageName)
    fig.savefig(imagePath,quality=ImageQuality)
    plt.close(fig)
    
    txtfile = open(os.path.join(folderPath,'%s%svs%s_%s__results_%s.txt' % (yRangeString,yLabel,independentVariableNames,architectureString,runID)),'w')
    messages.reverse()
    txtfile.write('\n'.join(messages))
    txtfile.close()
    
    #Graph the model's function if there's just one independent variable
    if len(xDataList)==1:
        x = np.linspace(np.amin(features),np.amax(features),num=200)
        y = model.predict(x)
        fig,ax = plt.subplots()
        ax.scatter(features,labels,color='b',label='data')
        ax.plot(x,y,color='r',label='prediction')
        ax.set_xlabel(independentVariableNames)
        ax.set_ylabel(yLabel)
        ax.legend()
        ax.grid()
        
        imageName = '%s%svs%s_%s__function_%s.%s' % (yRangeString,yLabel,independentVariableNames,architectureString,runID,ImageType)
        imagePath = os.path.join(folderPath,imageName)
        fig.savefig(imagePath,quality=ImageQuality)
        plt.close(fig)
    
    return yRangeString,yLabel,independentVariableNames,architectureString,rho,r_sq

if __name__ == '__main__':
    SPT_for_correlation,SPT_data,SPT_time_numbers,yRangeString,yLabel = get_spt_data('QUratio',1,'MLtest2')
    merra_data,xRangeString,xLabel = get_merra_data('TQI',2,'MLtest2')
    make_model([xLabel],[xRangeString],[merra_data],yLabel,yRangeString,SPT_for_correlation,'MLtest2', hiddenLayers = [16],
    PRINT_MESS = True,LOGX=True,LOGY=True,runID='1A')
