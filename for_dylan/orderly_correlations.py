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

parser = argparse.ArgumentParser(description = 'Checks the integration procedure',epilog = 'The program assumes we have a circle centered around the South Pole.')
monthChoices = [1,2,3,4,5,6,7,8,9,10,11,12]
dayChoices = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31]
testChoices = ['Spearman','linregress','ks']
SPTchoices = ['QUratio']
DIM_DICT = {'CPT':'[J m^-2]', 'KE':'[J m^-2]', 'MASS':'[kg m^-2]', 'THV':'[K]', 'TOX':'[kg m^-2]', 'TQI':'[kg m^-2]', 'TQL':'[kg m^-2]', 'TQV':'[kg m]', \
    'CLOUD':'[1]', 'DELP':'[Pa]', 'EPV':'[K m^2 kg^-1 s^-1]', 'H':'[m]', 'O3':'[kg kg^-1]', 'OMEGA':'[Pa s^-1]', 'PHIS':'[m^2 s^-2]', 'PL':'[Pa]', 'PS':'[Pa]', \
    'QI':'[kg kg^-1]', 'QL':'[kg kg^-1]', 'QV':'[kg kg^-1]', 'RH':'[1]', 'SLP':'[Pa]', 'T':'[K]', 'U':'[m s^-1]', 'V':'[m s^-1]'}

messages = []
LAT_RES = 0.5
LON_RES = 0.625
R = 6357000 #radius of Earth at the poles in meters
ZERO_MEAN = False
filters = ['None','<100m','100-1000m','>1000m']

NUM_OF_TICKS = 5

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

def make_images(height_samples,SP_array,sample_freqs,power,FILE_BASE_INST,VAR,DIM_DICT,date,t,messages,folderPath,NUM_OF_TICKS = 5,PRINT_MESS = True,test=False):
    testString = ".test" if test else ".%s.t=%i" % (date,t)
    imageName = "%s%s%s.png" % (FILE_BASE_INST,VAR,testString)
    message = "Making image %s..." % imageName
    if PRINT_MESS: print(message)
    messages.append(message)
    
    fig, (press_graph,oned_power) = plt.subplots(1,2)
    fig.suptitle("Power Spectra of %s on %s at t=%i" % (VAR,date,t),y=1.15)
    fig.subplots_adjust(wspace=0.5)
    
    #Altitude plot
    press_graph.plot(height_samples,SP_array,'b-')
    tick_positions = np.around(np.linspace(np.amin(height_samples),np.amax(height_samples),NUM_OF_TICKS)).astype(int)
    press_graph.set_xticks(tick_positions)
    
    pressTitle = "Delta %s %s vs. Height" % (VAR,DIM_DICT[VAR])
    press_graph.set_title(pressTitle,y=1.15)
    press_graph.set_xlabel("Height [m]")
    press_graph.set_ylabel("Delta %s %s" % (VAR,DIM_DICT[VAR]))
    press_graph.grid()
    
    #1-dimensional power spectra
    #Convert from inverse hPa to regular hPa
    
    oned_power.plot(sample_freqs,power,'b-')
    tick_positions = np.around(np.linspace(np.amin(sample_freqs),np.amax(sample_freqs),NUM_OF_TICKS),decimals=3)
    oned_power.set_xticks(tick_positions)
    
    freqTitle = "Power Spectrum of Delta %s" % VAR
    oned_power.set_title(freqTitle,y=1.15)
    oned_power.set_xlabel("Frequency [m^-1]")
    oned_power.set_ylabel("Power")
    oned_power.set_ybound(lower=0,upper=np.amax(power))
    oned_power.grid()
    
    fig.tight_layout()
    imagePath = os.path.join(folderPath,imageName)
    fig.savefig(imagePath)
    plt.close(fig)
    
    return messages

parser.add_argument('-fbi2','--file_base_instantaneous_2d',default='MERRA2_400.inst1_2d_asm_Nx.')
parser.add_argument('-fei2','--file_end_instantaneous_2d',default='.SUB.nc')
parser.add_argument('-fbt2','--file_base_time_averaged_2d',default='MERRA2_400.tavg1_2d_slv_Nx.')
parser.add_argument('-fet2','--file_end_time_averaged_2d',default='.SUB.nc')

parser.add_argument('-fbi3','--file_base_instantaneous_3d',default='MERRA2_400.inst3_3d_asm_Nv.')
parser.add_argument('-fei3','--file_end_instantaneous_3d',default='.SUB.nc')
parser.add_argument('-fbt3','--file_base_time_averaged_3d',default='MERRA2_400.tavg3_3d_asm_Nv.')
parser.add_argument('-fet3','--file_end_time_averaged_3d',default='.SUB.nc')
parser.add_argument('-sfn','--spt_file_name',default='ratios_dates.txt')

parser.add_argument('-fn','--folder_name',help='The name of the folder where the program\'s output is stored')

parser.add_argument('-mv','--merra_variable',help='The name of the MERRA-2 variable you wish to use')
parser.add_argument('-sv','--spt_variable',help='The SPT noise variable you want',choices=SPTchoices)
parser.add_argument('-d','--dimensions',help='Either the 2- or 3-D MERRA-2 Data',type=int,choices=[2,3])
parser.add_argument('-ff','--frequency_filter',choices = filters,default='None')
parser.add_argument('-r','--random',action='store_true')

parser.add_argument('-mi','--make_images',action='store_true')
parser.add_argument('-ns','--number_of_samples',type=int,default=300)
parser.add_argument('-pm','--print_messages',action='store_true')

parser.add_argument('-st','--statistical_test',choices=testChoices)
parser.add_argument('-invy','--inverse_y',action='store_true')
parser.add_argument('-logx','--logarithm_x',action='store_true')
parser.add_argument('-logy','--logarithm_y',action='store_true')

parser.add_argument('-sy','--start_year',type=int,default=2019)
parser.add_argument('-sm','--start_month',type=int,choices=monthChoices,default=3)
parser.add_argument('-sd','--start_day',type=int,choices=dayChoices,default=22)

parser.add_argument('-ey','--end_year',type=int,default=2019)
parser.add_argument('-em','--end_month',type=int,choices=monthChoices,default=8)
parser.add_argument('-ed','--end_day',type=int,choices=dayChoices,default=11)

args = parser.parse_args()

FILE_BASE_INST2 = args.file_base_instantaneous_2d
FILE_END_INST2 = args.file_end_instantaneous_2d
FILE_BASE_TAVG2 = args.file_base_time_averaged_2d
FILE_END_TAVG2 = args.file_end_time_averaged_2d

FILE_BASE_INST3 = args.file_base_instantaneous_3d
FILE_END_INST3 = args.file_end_instantaneous_3d
FILE_BASE_TAVG3 = args.file_base_time_averaged_3d
FILE_END_TAVG3 = args.file_end_time_averaged_3d

SPT_FILE_NAME = args.spt_file_name

START_DATE = dt.date(year=args.start_year,month=args.start_month,day=args.start_day)
END_DATE = dt.date(year=args.end_year,month=args.end_month,day=args.end_day)
FOLDER_NAME = args.folder_name #The name of the folder where the images are stored, which will be created if it doesn't already exist

MERRA_VAR = args.merra_variable #Which meteorological quantity we're checking
RANDOM = args.random
SPT_VAR = args.spt_variable
DIM = args.dimensions
if DIM==3:
    T_INTERVAL = 3
if DIM==2:
    T_INTERVAL = 1
NUM_OF_TIMES = 24 // T_INTERVAL
FREQ_FILTER = args.frequency_filter

MAKE_IMAGES = args.make_images
NUM_OF_SAMPS = args.number_of_samples
PRINT_MESS = args.print_messages

STAT_TEST = args.statistical_test
INVY = args.inverse_y
LOGX = args.logarithm_x
LOGY = args.logarithm_y

#Generate the dates
date = START_DATE
interval = dt.timedelta(days=1)
dates = []
while date <= END_DATE:
    dates.append(date.isoformat().replace('-',''))
    date = date + interval

message = "Loading SPT data..."
if PRINT_MESS: print(message)
messages.append(message)

if SPT_VAR=='QUratio':
    SPT = np.loadtxt(SPT_FILE_NAME,dtype=str)
    SPT_data = SPT[:,0]
    SPT_time_strings = SPT[:,1]

SPT_data = SPT_data.astype(float)

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


#This creates the folder where everything will go

dir_path = os.path.dirname(os.path.realpath(__file__))
folderPath = os.path.join(dir_path,FOLDER_NAME)
if not os.path.isdir(folderPath):
    os.mkdir(folderPath)

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
            merra_data[i] = abs(setsOfData_tavg[date][MERRA_VAR][t,0,0] - setsOfData_inst[date][MERRA_VAR][t,0,0])
        elif DIM==3:
            #Load the data and then, because the pressure levels may move up and down, interpolate the difference on a common set of sample heights
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
            
            sample_freqs,power = compute_1d_psd(SP_array,sample_spacing = height_interval)
            
            if MAKE_IMAGES:
                messages = make_images(height_samples,SP_array,sample_freqs,power,FILE_BASE_INST,VAR,DIM_DICT,date,t,messages,folderPath,NUM_OF_TICKS,PRINT_MESS,test=False)
    
            if FREQ_FILTER=='None':
    
                psd_integral = grate.trapz(power,sample_freqs)/(NUM_OF_SAMPS**2)*intervalLength #I don't know why you need to divide by (NUM_OF_SAMPS**2) and multiply by the interval length, but that made it work
    
            if FREQ_FILTER=='<100m':
    
                valid_neg_freqs = ((1/100)<=np.abs(sample_freqs))&(sample_freqs<0)
                valid_pos_freqs = ((1/100)<=np.abs(sample_freqs))&(sample_freqs>0)
                
                filtered_neg_power = power[valid_neg_freqs]
                filtered_pos_power = power[valid_pos_freqs]
                filtered_sample_neg_freqs = sample_freqs[valid_neg_freqs]
                filtered_sample_pos_freqs = sample_freqs[valid_pos_freqs]
                
                
                psd_integral = (grate.trapz(filtered_neg_power,filtered_sample_neg_freqs) + grate.trapz(filtered_pos_power,filtered_sample_pos_freqs))/(NUM_OF_SAMPS**2)*intervalLength
    
            if FREQ_FILTER=='100-1000m':
                
                valid_neg_freqs = ((1/100)>=np.abs(sample_freqs))&(np.abs(sample_freqs)>=(1/1000))&(sample_freqs<0)
                valid_pos_freqs = ((1/100)>=np.abs(sample_freqs))&(np.abs(sample_freqs)>=(1/1000))&(sample_freqs>0)
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
    
            if FREQ_FILTER=='>1000m':
    
                valid_freqs = np.abs(sample_freqs)<=(1/1000)
    
                filtered_power = power[valid_freqs]
                filtered_sample_freqs = sample_freqs[valid_freqs]
                
    
                psd_integral = grate.trapz(filtered_power,filtered_sample_freqs)/(NUM_OF_SAMPS**2)*intervalLength
            
            merra_data[i] = psd_integral
        
        if RANDOM:
            merra_data[i] = np.random.normal()
        
        message = "Processed day %s, t=%i." % (date,t)
        if PRINT_MESS: print(message)
        messages.append(message)
        
        i += 1

#Do statistical test and make scatter and time series graphs
yLabel = SPT_VAR
yLabel = '(%s)^-1' % yLabel if INVY else yLabel

yLabel = 'log(%s)' % yLabel if LOGY else yLabel

if DIM==2:
    xLabel = "Delta_"+MERRA_VAR
elif DIM==3:
    xLabel = "Variance_of_Delta_"+MERRA_VAR
    if FREQ_FILTER != 'None':
        xLabel = FREQ_FILTER+" "+xLabel
if RANDOM:
    xLabel = 'random_var'
xLabel = 'log(%s)' % xLabel if LOGX else xLabel

if INVY:
    SPT_for_correlation = 1/SPT_for_correlation
    SPT_data = 1/SPT_data

if LOGY:
    SPT_for_correlation = np.log10(SPT_for_correlation)
    SPT_data = np.log10(SPT_data)

if LOGX:
    merra_data = np.log10(merra_data)

merra_for_correlation = merra_data.copy()
merra_for_timeseries = merra_data.copy()
#first get rid of the values which are masked because there wasn't an SPT reading from that interval, as well as any nan or infinite values
invalidValues = SPT_for_correlation.mask | np.isnan(SPT_for_correlation) | np.isinf(SPT_for_correlation) | np.isnan(merra_for_correlation) | np.isinf(merra_for_correlation)

merra_for_correlation = merra_for_correlation[~invalidValues]
SPT_for_correlation = SPT_for_correlation[~invalidValues]
    
if STAT_TEST=='linregress':

    regression = sts.linregress(merra_for_correlation,SPT_for_correlation)
    slope = regression.slope
    intercept = regression.intercept
    r_sq = regression.rvalue**2
    pvalue = regression.pvalue
    message = "r^2: %f, slope: %f, intercept: %f, p-value: %f" % (r_sq,slope,intercept,pvalue)
    if PRINT_MESS: print(message)
    messages.append(message)
    graphTitle = "r^2: %.4f, p-value: %.4f" % (r_sq,pvalue)
    
    
elif STAT_TEST=='Spearman':
    
    rho, pvalue = sts.spearmanr(merra_for_correlation,SPT_for_correlation)
    message = "rho: %f, p-value: %f" % (rho,pvalue)
    if PRINT_MESS: print(message)
    messages.append(message)
    graphTitle = "rho: %.4f, p-value: %.4f" % (rho,pvalue)

elif STAT_TEST=='ks':
    
    #Split the MERRA data into greater than average and less than average for the two sets of data
    aboveMedian = merra_for_correlation>np.median(merra_for_correlation)
    ks_stat,pvalue = sts.ks_2samp(SPT_for_correlation[aboveMedian],SPT_for_correlation[~aboveMedian])
    message = "KS stat: %f, p-value: %f" % (ks_stat,pvalue)
    if PRINT_MESS: print(message)
    messages.append(message)
    graphTitle = "KS stat: %.4f, p-value: %.4f" % (ks_stat,pvalue)

fig, ax = plt.subplots()
ax.plot(merra_for_correlation,SPT_for_correlation,'bo',label='Data points')

if STAT_TEST=='linregress':
    ax.plot(merra_for_correlation,intercept + slope*merra_for_correlation,'r-',label = 'Best fit')
if STAT_TEST=='ks':
    ax.vlines(np.median(merra_for_correlation),0,1,transform = ax.get_xaxis_transform(),colors='r',label='KS divider')

ax.set_title(graphTitle)
ax.set_xlabel(xLabel)
ax.set_ylabel(yLabel)

ax.legend()
ax.grid()
fig.tight_layout()

imageName = '%s.vs.%s.%s.%s-%s.scatter.png' % (yLabel,xLabel,STAT_TEST,START_DATE.strftime('%Y%m%d'),END_DATE.strftime('%Y%m%d'))
imagePath = os.path.join(folderPath,imageName)
fig.savefig(imagePath)
plt.close(fig)

#Graph functions of time
#Now get rid of invalid values for the time series data
invalidValues = np.isnan(SPT_data) | np.isinf(SPT_data)
SPT_data = SPT_data[~invalidValues]
SPT_time_numbers = SPT_time_numbers[~invalidValues]

invalidValues = np.isnan(merra_for_timeseries) | np.isinf(merra_for_timeseries)

times = T_INTERVAL*np.arange(NUM_OF_TIMES*len(dates))
times = times[~invalidValues]
merra_for_timeseries = merra_for_timeseries[~invalidValues]

fig, ax1 = plt.subplots()

ax1.plot(times,merra_for_timeseries,'b-',label=xLabel)
ax1.plot([],[],'r-',label=yLabel)
ax1.set_xlabel('Time [hrs]')
ax1.set_ylabel(xLabel)
ax1.set_title("%s and %s vs. Time" % (xLabel,yLabel))
ax1.legend()

ax2 = ax1.twinx()
ax2.plot(SPT_time_numbers,SPT_data,'r-',label=yLabel)
ax2.set_ylabel(yLabel)

fig.tight_layout()

imageName = '%s.vs.%s.%s-%s.time_series.png' % (yLabel,xLabel,START_DATE.strftime('%Y%m%d'),END_DATE.strftime('%Y%m%d'))
imagePath = os.path.join(folderPath,imageName)
fig.savefig(imagePath)
plt.close(fig)

txtfile = open(os.path.join(folderPath,'%s.vs.%s.%s.%s-%s.results.txt' % (yLabel,xLabel,STAT_TEST,START_DATE.strftime('%Y%m%d'),END_DATE.strftime('%Y%m%d'))),'w')
messages.reverse()
txtfile.write('\n'.join(messages))
txtfile.close()