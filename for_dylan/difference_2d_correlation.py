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

parser.add_argument('-fn','--folder_name',help='The name of the folder where the program\'s output is stored')

parser.add_argument('-v','--variable',help='The name of the 2D MERRA-2 variable you wish to use')
parser.add_argument('-mi','--make_images',action='store_true')
parser.add_argument('-ns','--number_of_samples',type=int,default=300)
parser.add_argument('-pm','--print_messages',action='store_true')
parser.add_argument('-in','--inverse',action='store_true')
parser.add_argument('-log','--logarithm',action='store_true')
parser.add_argument('-ttt','--ten_to_the',action='store_true')
parser.add_argument('-wf','--weird_function',action='store_true')

parser.add_argument('-sy','--start_year',type=int)
parser.add_argument('-sm','--start_month',type=int,choices=monthChoices)
parser.add_argument('-sd','--start_day',type=int,choices=dayChoices)

parser.add_argument('-ey','--end_year',type=int)
parser.add_argument('-em','--end_month',type=int,choices=monthChoices)
parser.add_argument('-ed','--end_day',type=int,choices=dayChoices)

args = parser.parse_args()

FILE_BASE_INST = args.file_base_instantaneous
FILE_END_INST = args.file_end_instantaneous
FILE_BASE_TAVG = args.file_base_time_averaged
FILE_END_TAVG = args.file_end_time_averaged
SPT_FILE_NAME = 'ratios_dates.txt' #Could put this in the argument parser

START_DATE = dt.date(year=args.start_year,month=args.start_month,day=args.start_day)
END_DATE = dt.date(year=args.end_year,month=args.end_month,day=args.end_day)
FOLDER_NAME = args.folder_name #The name of the folder where the images are stored, which will be created if it doesn't already exist

VAR = args.variable #Which meteorological quantity we're checking
MAKE_IMAGES = args.make_images
NUM_OF_SAMPS = args.number_of_samples
PRINT_MESS = args.print_messages
TAKE_INVERSE = args.inverse
TAKE_LOG = args.logarithm
TAKE_TTT = args.ten_to_the
TAKE_WEIRD = args.weird_function

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

SPT_time_indices = (SPT_time_numbers // 1).astype(int)

SPT_for_correlation = np.ma.array(np.zeros(24*len(dates)))
SPT_for_correlation.mask = True
for observationIndex,timeIndex in enumerate(SPT_time_indices):
    if SPT_for_correlation.mask[timeIndex]: #This makes it so that if there are several observations in a time interval, only the first goes into SPT_for_correlation
        SPT_for_correlation[timeIndex] = SPT_data[observationIndex]


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
for date in dates:
    setsOfData_inst[date] = nc.Dataset(FILE_BASE_INST+date+FILE_END_INST)
    setsOfData_tavg[date] = nc.Dataset(FILE_BASE_TAVG+date+FILE_END_TAVG)

message = "Loaded "+str(len(setsOfData_tavg))+ " days of data."
if PRINT_MESS: print(message)
messages.append(message)

differences = np.zeros(24*len(dates))
i = 0
for date in dates:
    for t in range(24):
        
        differences[i] = abs(setsOfData_tavg[date][VAR][t,0,0] - setsOfData_inst[date][VAR][t,0,0])
        
        message = "Loaded day %s, t=%i." % (date,t)
        print(message)
        messages.append(message)
        
        i += 1

#Do linear regression
ratio = '1/(Q/U ratio)' if TAKE_INVERSE else 'Q/U ratio'
ratio = '10^('+ratio+')' if TAKE_TTT else ratio
modeString = '_ten_to_the' if TAKE_TTT else ''
modeString = modeString+'_inverse' if TAKE_INVERSE else modeString
modeString = modeString+'_log10' if TAKE_LOG else modeString
if TAKE_WEIRD:
    ratio = 'log_10(|-1+(1/(Q/U))|)'
    modeString = modeString+'_weird'
description = 'log_10(tavg-inst)' if TAKE_LOG else 'tavg-inst'

if TAKE_WEIRD:
    SPT_for_correlation = np.log10(np.abs(-1+np.divide(1,SPT_for_correlation)))
    SPT_data = np.log10(np.abs(-1+np.divide(1,SPT_data)))
else:
    if TAKE_INVERSE:
        SPT_for_correlation = 1/SPT_for_correlation
        SPT_data = 1/SPT_data
    
    if TAKE_TTT:
        SPT_for_correlation = np.power(10,SPT_for_correlation)
        SPT_data = np.power(10,SPT_data)

#first get rid of the values which are masked because there wasn't an SPT reading from that interval, as well as any nan or infinite alues
invalidValues = SPT_for_correlation.mask | np.isnan(SPT_for_correlation) | np.isinf(SPT_for_correlation)

differences_for_correlation = differences[~invalidValues]
SPT_for_correlation = SPT_for_correlation[~invalidValues]

if TAKE_LOG:
    differences_for_correlation = np.log10(differences_for_correlation)
    invalidValues = np.isnan(differences_for_correlation) | np.isinf(differences_for_correlation)
    differences_for_correlation = differences_for_correlation[~invalidValues]
    SPT_for_correlation = SPT_for_correlation[~invalidValues]

regression = sts.linregress(differences_for_correlation,SPT_for_correlation)

slope = regression.slope
intercept = regression.intercept
r_sq = regression.rvalue**2

fig, ax = plt.subplots()
ax.plot(differences_for_correlation,SPT_for_correlation,'bo',label='Data points')
ax.plot(differences_for_correlation,intercept + slope*differences_for_correlation,'r-',label = 'Best fit')

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
ax.set_xlabel('%s for %s' % (description,VAR))
ax.set_ylabel(ratio)

ax.legend()
ax.grid()
fig.tight_layout()

imageName = '0%s,%s.%s.%s-%s.scatter%s.png' % (FILE_BASE_INST,FILE_BASE_TAVG,VAR,START_DATE.strftime('%Y%m%d'),END_DATE.strftime('%Y%m%d'),modeString)
imagePath = os.path.join(folderPath,imageName)
fig.savefig(imagePath)
plt.close(fig)

#Graph functions of time
#Now get rid of invalid values for the time series data
invalidValues = np.isnan(SPT_data) | np.isinf(SPT_data)
SPT_data = SPT_data[~invalidValues]
SPT_time_numbers = SPT_time_numbers[~invalidValues]

times = np.arange(24*len(dates))
if TAKE_LOG:
    differences = np.log10(differences)
    invalidValues = np.isnan(differences) | np.isinf(differences)
    times = times[~invalidValues]
    differences = differences[~invalidValues]

fig, ax1 = plt.subplots()
varianceLabel = '%s for %s' % (description,VAR)
ratioLabel = ratio
ax1.plot(times,differences,'b-',label=varianceLabel)
ax1.plot([],[],'r-',label=ratioLabel)
ax1.set_xlabel('Time [hrs]')
ax1.set_ylabel(varianceLabel)
ax1.set_title("%s %s and %s vs. Time" % (VAR,description,ratio))
ax1.legend()

ax2 = ax1.twinx()
ax2.plot(SPT_time_numbers,SPT_data,'r-',label=ratioLabel)
ax2.set_ylabel(ratioLabel)

fig.tight_layout()

imageName = '0%s,%s.%s.%s-%s.time_series%s.png' % (FILE_BASE_INST,FILE_BASE_TAVG,VAR,START_DATE.strftime('%Y%m%d'),END_DATE.strftime('%Y%m%d'),modeString)
imagePath = os.path.join(folderPath,imageName)
fig.savefig(imagePath)
plt.close(fig)

txtfile = open(os.path.join(folderPath,'0%s,%s.%s.%s-%s.results%s.txt' % (FILE_BASE_INST,FILE_BASE_TAVG,VAR,START_DATE.strftime('%Y%m%d'),END_DATE.strftime('%Y%m%d'),modeString)),'w')
txtfile.write('\n'.join(messages))
txtfile.close()