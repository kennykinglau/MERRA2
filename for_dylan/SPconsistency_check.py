import netCDF4 as nc
import numpy as np
import scipy.integrate as grate
import scipy.stats as sts
import datetime
import os
import argparse
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser(description = 'Checks a dataset for consistency at the South Pole',epilog = 'The program assumes we have a circle centered around the South Pole.')
monthChoices = [1,2,3,4,5,6,7,8,9,10,11,12]
dayChoices = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31]
DIM_DICT = {'CPT':'[J m^-2]', 'KE':'[J m^-2]', 'MASS':'[kg m^-2]', 'THV':'[K]', 'TOX':'[kg m^-2]', 'TQI':'[kg m^-2]', 'TQL':'[kg m^-2]', 'TQV':'[kg m]', \
    'CLOUD':'[1]', 'DELP':'[Pa]', 'EPV':'[K m^2 kg^-1 s^-1]', 'H':'[m]', 'O3':'[kg kg^-1]', 'OMEGA':'[Pa s^-1]', 'PHIS':'[m^2 s^-2]', 'PL':'[Pa]', 'PS':'[Pa]', \
    'QI':'[kg kg^-1]', 'QL':'[kg kg^-1]', 'QV':'[kg kg^-1]', 'RH':'[1]', 'SLP':'[Pa]', 'T':'[K]', 'U':'[m s^-1]', 'V':'[m s^-1]'}
LAT_RES = 0.5
LON_RES = 0.625
messages = []

parser.add_argument('-fb','--file_base')
parser.add_argument('-fe','--file_end')
parser.add_argument('-fn','--folder_name',help='The name of the folder where the program\'s output is stored')

parser.add_argument('-v','--variable',help='The name of the MERRA-2 variable you wish to check')
parser.add_argument('-mi','--make_images',action='store_true')
parser.add_argument('-pm','--print_messages',action='store_true')

parser.add_argument('-sy','--start_year',type=int)
parser.add_argument('-sm','--start_month',type=int,choices=monthChoices)
parser.add_argument('-sd','--start_day',type=int,choices=dayChoices)

parser.add_argument('-ey','--end_year',type=int)
parser.add_argument('-em','--end_month',type=int,choices=monthChoices)
parser.add_argument('-ed','--end_day',type=int,choices=dayChoices)

args = parser.parse_args()

FILE_BASE = args.file_base
FILE_END = args.file_end
START_DATE = datetime.date(year=args.start_year,month=args.start_month,day=args.start_day)
END_DATE = datetime.date(year=args.end_year,month=args.end_month,day=args.end_day)
FOLDER_NAME = args.folder_name #The name of the folder where the images are stored, which will be created if it doesn't already exist
VAR = args.variable #Which meteorological quantity we're checking
MAKE_IMAGES = args.make_images
PRINT_MESS = args.print_messages

#This creates the folder where everything will go

if MAKE_IMAGES:
    dir_path = os.path.dirname(os.path.realpath(__file__))
    folderPath = os.path.join(dir_path,FOLDER_NAME)
if not os.path.isdir(folderPath):
    os.mkdir(folderPath)

#Generate the dates
date = START_DATE
interval = datetime.timedelta(days=1)
dates = []
while date <= END_DATE:
    dates.append(date.isoformat().replace('-',''))
    date = date + interval

message = "Loading data..."
if PRINT_MESS: print(message)
messages.append(message)

#Load the data
setsOfData = {}
for date in dates:
    setsOfData[date] = nc.Dataset(FILE_BASE+date+FILE_END)
message = "Loaded "+str(len(setsOfData))+ " days of data."
if PRINT_MESS: print(message)
messages.append(message)
data_merra = setsOfData[dates[0]] #load a dataset to determine the dimensions

#Determine the dimensions of the data
if len(data_merra[VAR].shape)==4:
    IS_3D = True
elif len(data_merra[VAR].shape)==3:
    IS_3D = False

numOfTimes = len(data_merra[VAR])
#Generate time strings for later use
timeStrings = {}
for t in range(numOfTimes):
    if numOfTimes >= 10 and t<10:
        timeStrings[t] = "0"+str(t)
    else:
        timeStrings[t] = str(t)

if IS_3D:
    numOfLevs = len(data_merra[VAR][0])
    numOfLats = len(data_merra[VAR][0][0])
    numOfLons = len(data_merra[VAR][0][0][0])
else:
    numOfLats = len(data_merra[VAR][0])
    numOfLons = len(data_merra[VAR][0][0])

message = "Checking consistency..."
if PRINT_MESS: print(message)
messages.append(message)

for date in dates:
    
    data_merra = setsOfData[date]
    
    for t in range(numOfTimes):
        timeString = timeStrings[t]
        ice_array = data_merra[VAR][t]
        
        if IS_3D:
            for lev in range(numOfLevs):
                imageName = '%s%s.SPinconsistency.%s.t=%s.lev=%s.png' % (FILE_BASE,VAR,date,timeString,str(lev))
                SP_array = ice_array[lev,0,:]
                
                if np.amax(SP_array)==np.amin(SP_array):
                    message = "Consistency confirmed for %s on %s at t=%s, lev=%s" % (VAR,date,timeString,str(lev))
                    if PRINT_MESS: print(message)
                    messages.append(message)
                else:
                    message = "Inconsistency found for %s on %s at t=%s" % (VAR,date,timeString)
                    if PRINT_MESS: print(message)
                    messages.append(message)
                    
                    if MAKE_IMAGES:
                        fig,ax = plt.subplots()
                        lon = -180+LON_RES*np.arange(numOfLons)
                        ax.plot(lon,SP_array)
                        
                        ax.set_title('%s vs. Longitude, %s, t=%s' % (VAR,date,timeString))
                        ax.set_xlabel("Longitude")
                        ax.set_ylabel('%s %s' % (VAR,DIM_DICT[VAR]))
                        fig.savefig(os.path.join(folderPath,imageName))
                        plt.close()
        else:
            imageName = '%s%s.SPinconsistency.%s.t=%s.png' % (FILE_BASE,VAR,date,timeString)
            SP_array = ice_array[0,:]
            
            if np.amax(SP_array)==np.amin(SP_array):
                message = "Consistency confirmed for %s on %s at t=%s" % (VAR,date,timeString)
                if PRINT_MESS: print(message)
                messages.append(message)
            else:
                message = "Inconsistency found for %s on %s at t=%s" % (VAR,date,timeString)
                if PRINT_MESS: print(message)
                messages.append(message)
                
                if MAKE_IMAGES:
                    fig,ax = plt.subplots()
                    lon = -180+LON_RES*np.arange(numOfLons)
                    ax.plot(lon,SP_array)
                    
                    ax.set_title('%s vs. Longitude, %s, t=%s' % (VAR,date,timeString))
                    ax.set_xlabel("Longitude")
                    ax.set_ylabel('%s %s' % (VAR,DIM_DICT[VAR]))
                    fig.savefig(os.path.join(folderPath,imageName))
                    plt.close()

txtfile = open(os.path.join(folderPath,'0%s.%s.%s-%s.results.txt' % (FILE_BASE,VAR,START_DATE.strftime('%Y%m%d'),END_DATE.strftime('%Y%m%d'))),'w')
txtfile.write('\n'.join(messages))
txtfile.close()