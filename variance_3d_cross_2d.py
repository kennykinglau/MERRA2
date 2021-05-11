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
DIM_DICT = {'CPT':'[J m^-2]', 'KE':'[J m^-2]', 'MASS':'[kg m^-2]', 'THV':'[K]', 'TOX':'[kg m^-2]', 'TQI':'[kg m^-2]', 'TQL':'[kg m^-2]', 'TQV':'[kg m]', \
    'CLOUD':'[1]', 'DELP':'[Pa]', 'EPV':'[K m^2 kg^-1 s^-1]', 'H':'[m]', 'O3':'[kg kg^-1]', 'OMEGA':'[Pa s^-1]', 'PHIS':'[m^2 s^-2]', 'PL':'[Pa]', 'PS':'[Pa]', \
    'QI':'[kg kg^-1]', 'QL':'[kg kg^-1]', 'QV':'[kg kg^-1]', 'RH':'[1]', 'SLP':'[Pa]', 'T':'[K]', 'U':'[m s^-1]', 'V':'[m s^-1]'}

messages = []
LAT_RES = 0.5
LON_RES = 0.625
R = 6357000 #radius of Earth at the poles in meters

NUM_OF_TICKS = 5
numOfLats = 11 #TODO: Make it smarter about this
numOfLons = 576 #TODO: Make it smarter about this
colorMap = 'Blues'
zeroMean = False #TODO: Put in the argparser?
filters = ['None','<100m','100-1000m','>1000m']

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
    
    #1-dimensional power spectru
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

parser.add_argument('-fbi','--file_base_instantaneous')
parser.add_argument('-fei','--file_end_instantaneous')
parser.add_argument('-fbt','--file_base_time_averaged')
parser.add_argument('-fet','--file_end_time_averaged')
parser.add_argument('-fb2d','--file_base_2d')
parser.add_argument('-fe2d','--file_end_2d')

parser.add_argument('-fn','--folder_name',help='The name of the folder where the program\'s output is stored')

parser.add_argument('-v3','--variable_3d',help='The name of the 3D MERRA-2 variable you wish to use')
parser.add_argument('-v2','--variable_2d')
parser.add_argument('-mi','--make_images',action='store_true')
parser.add_argument('-ns','--number_of_samples',type=int,default=300)
parser.add_argument('-pm','--print_messages',action='store_true')

parser.add_argument('-sy','--start_year',type=int)
parser.add_argument('-sm','--start_month',type=int,choices=monthChoices)
parser.add_argument('-sd','--start_day',type=int,choices=dayChoices)

parser.add_argument('-ey','--end_year',type=int)
parser.add_argument('-em','--end_month',type=int,choices=monthChoices)
parser.add_argument('-ed','--end_day',type=int,choices=dayChoices)
parser.add_argument('-ff','--frequency_filter',choices = filters,default='None')

args = parser.parse_args()

FILE_BASE_INST = args.file_base_instantaneous
FILE_END_INST = args.file_end_instantaneous
FILE_BASE_TAVG = args.file_base_time_averaged
FILE_END_TAVG = args.file_end_time_averaged
FILE_BASE_2D = args.file_base_2d
FILE_END_2D = args.file_end_2d

START_DATE = dt.date(year=args.start_year,month=args.start_month,day=args.start_day)
END_DATE = dt.date(year=args.end_year,month=args.end_month,day=args.end_day)
FOLDER_NAME = args.folder_name #The name of the folder where the images are stored, which will be created if it doesn't already exist

VAR_3D = args.variable_3d
VAR_2D = args.variable_2d
MAKE_IMAGES = args.make_images
NUM_OF_SAMPS = args.number_of_samples
PRINT_MESS = args.print_messages
FREQ_FILTER = args.frequency_filter

#Generate the dates
date = START_DATE
interval = dt.timedelta(days=1)
dates = []
while date <= END_DATE:
    dates.append(date.isoformat().replace('-',''))
    date = date + interval

#This creates the folder where everything will go

if MAKE_IMAGES:
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
setsOfData_2d = {}
for date in dates:
    setsOfData_inst[date] = nc.Dataset(FILE_BASE_INST+date+FILE_END_INST)
    setsOfData_tavg[date] = nc.Dataset(FILE_BASE_TAVG+date+FILE_END_TAVG)
    setsOfData_2d[date] = nc.Dataset(FILE_BASE_2D+date+FILE_END_2D)

message = "Loaded "+str(len(setsOfData_tavg))+ " days of data."
if PRINT_MESS: print(message)
messages.append(message)

variances = np.zeros(len(dates)*8)
values_2d = np.zeros(len(dates)*8)
i = 0
for date in dates:
    for t in range(8):
        #Load the data and then, because the pressure levels may move up and down, interpolate the difference on a common set of sample heights
        SP_array_inst = setsOfData_inst[date][VAR_3D][t,:,0,0]
        heights_inst = setsOfData_inst[date]['H'][t,:,0,0]
        SP_array_tavg = setsOfData_tavg[date][VAR_3D][t,:,0,0]
        heights_tavg = setsOfData_tavg[date]['H'][t,:,0,0]
        #Load the 2D instantaneous value at the beginning of the interval
        values_2d[i] = setsOfData_2d[date][VAR_2D][3*t,0,0]
        
        message = "Loaded day %s, t=%i." % (date,t)
        print(message)
        messages.append(message)
        
        #Remove high altitudes where they are both zero
        nonzero_inst = SP_array_inst!=0
        nonzero_heights_inst = heights_inst[nonzero_inst]
        
        message = "Number of nonzero instantaneous values: %i" % nonzero_heights_inst.shape[0]
        print(message)
        messages.append(message)
        
        nonzero_tavg = SP_array_tavg!=0
        nonzero_heights_tavg =  heights_tavg[nonzero_tavg]
        
        message = "Number of nonzero time-averaged values: %i" % nonzero_heights_tavg.shape[0]
        print(message)
        messages.append(message)
        
        interpolater_inst = polate.interp1d(heights_inst,SP_array_inst)
        interpolater_tavg = polate.interp1d(heights_tavg,SP_array_tavg)
        
        minHeight = min(np.amin(heights_inst),np.amin(heights_tavg))
        
        message = "The minimum height in the data for t=%i is %f" % (t,minHeight)
        print(message)
        messages.append(message)
        
        
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
        print(message)
        messages.append(message)
        SP_array = interpolater_tavg(height_samples) - interpolater_inst(height_samples)
        
        #Remove high altitudes where the difference is zero
        #nonZero = SP_array!=0
        #nonZeroHeights = height_samples[nonZero]
        
        #Generate new sample heights only in the range where the difference was nonzero
        #height_samples, height_interval = np.linspace(np.amin(nonZeroHeights),np.amax(nonZeroHeights),NUM_OF_SAMPS,retstep = True)
        #intervalLength = np.amax(nonZeroHeights)-np.amin(nonZeroHeights)
        #SP_array = interpolater_tavg(height_samples) - interpolater_inst(height_samples)
        if zeroMean: SP_array = SP_array - np.mean(SP_array)

        variance = np.var(SP_array)
        variances[i] = variance
        
        sample_freqs,power = compute_1d_psd(SP_array,sample_spacing = height_interval)
        
        if MAKE_IMAGES:
            messages = make_images(height_samples,SP_array,sample_freqs,power,FILE_BASE_INST,VAR_3D,DIM_DICT,date,t,messages,folderPath,NUM_OF_TICKS,PRINT_MESS,test=False)

        if FREQ_FILTER=='None':

            psd_integral = grate.trapz(power,sample_freqs)/(NUM_OF_SAMPS**2)*intervalLength #I don't know why you need to divide by (NUM_OF_SAMPS**2) and multiply by the interval length, but that made it work
            total = psd_integral

        if FREQ_FILTER=='<100m':

            valid_neg_freqs = ((1/100)<=np.abs(sample_freqs))&(sample_freqs<0)
            valid_pos_freqs = ((1/100)<=np.abs(sample_freqs))&(sample_freqs>0)
            
            filtered_neg_power = power[valid_neg_freqs]
            filtered_pos_power = power[valid_pos_freqs]
            filtered_sample_neg_freqs = sample_freqs[valid_neg_freqs]
            filtered_sample_pos_freqs = sample_freqs[valid_pos_freqs]
            
            
            psd_integral = (grate.trapz(filtered_neg_power,filtered_sample_neg_freqs) + grate.trapz(filtered_pos_power,filtered_sample_pos_freqs))/(NUM_OF_SAMPS**2)*intervalLength
            firstTerm = psd_integral

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
            secondTerm = psd_integral

        if FREQ_FILTER=='>1000m':

            valid_freqs = np.abs(sample_freqs)<=(1/1000)

            filtered_power = power[valid_freqs]
            filtered_sample_freqs = sample_freqs[valid_freqs]
            

            psd_integral = grate.trapz(filtered_power,filtered_sample_freqs)/(NUM_OF_SAMPS**2)*intervalLength
            thirdTerm = psd_integral
            

        variances[i] = psd_integral
        i += 1
            
#Do linear regression

regression = sts.linregress(values_2d,variances)

slope = regression.slope
intercept = regression.intercept
r_sq = regression.rvalue**2

fig, ax = plt.subplots()
ax.plot(values_2d,variances,'bo',label='Data points')
ax.plot(values_2d,intercept + slope*values_2d,'r-',label = 'Best fit')

message = "r^2: "+str(r_sq)
print(message)
messages.append(message)

message = "Slope: "+str(slope)
print(message)
messages.append(message)

message = "Intercept: "+str(intercept)
print(message)
messages.append(message)


ax.set_title("r^2=%.3f, slope=%.3f, intercept = %.3f" % (r_sq,slope,intercept))
ax.set_ylabel('%s Filtered Variance' % FREQ_FILTER)
ax.set_xlabel("%s %s" % (VAR_2D,DIM_DICT[VAR_2D]))

ax.legend()
ax.grid()
fig.tight_layout()

imageName = '0%s.%s,%s.%svs%s.%s-%s.scatter.png' % (FREQ_FILTER,FILE_BASE_INST,FILE_BASE_TAVG,VAR_3D,VAR_2D,START_DATE.strftime('%Y%m%d'),END_DATE.strftime('%Y%m%d'))
imagePath = os.path.join(folderPath,imageName)
fig.savefig(imagePath)
plt.close(fig)

#Graph functions of time
#Now get rid of invalid values for the time series data

times = 3*np.arange(8*len(dates))

fig, ax1 = plt.subplots()
varianceLabel = '%s Filtered variance' % FREQ_FILTER
twoDLabel = '%s %s' % (VAR_2D,DIM_DICT[VAR_2D])
ax1.plot(times,variances,'b-',label=varianceLabel)
ax1.plot([],[],'r-',label=twoDLabel)
ax1.set_xlabel('Time [hrs]')
ax1.set_ylabel(varianceLabel)
ax1.set_title("%s Filtered Variance and %s %s vs. Time" % (FREQ_FILTER,VAR_2D,DIM_DICT[VAR_2D]))
ax1.legend()

ax2 = ax1.twinx()
ax2.plot(times,values_2d,'r-',label=twoDLabel)
ax2.set_ylabel(twoDLabel)

fig.tight_layout()

imageName = '0%s.%s,%s.%svs%s.%s-%s.time_series.png' % (FREQ_FILTER,FILE_BASE_INST,FILE_BASE_TAVG,VAR_3D,VAR_2D,START_DATE.strftime('%Y%m%d'),END_DATE.strftime('%Y%m%d'))
imagePath = os.path.join(folderPath,imageName)
fig.savefig(imagePath)
plt.close(fig)

txtfile = open(os.path.join(folderPath,'0%s.%s,%s.%svs%s.%s-%s.results.txt' % (FREQ_FILTER,FILE_BASE_INST,FILE_BASE_TAVG,VAR_3D,VAR_2D,START_DATE.strftime('%Y%m%d'),END_DATE.strftime('%Y%m%d'))),'w')
txtfile.write('\n'.join(messages))
txtfile.close()