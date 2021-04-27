#! /usr/bin/env python

import netCDF4 as nc
import numpy as np
import scipy.integrate as grate
import scipy.stats as sts
import datetime
import os
import argparse
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser(description = 'Checks the integration procedure',epilog = 'The program assumes we have a circle centered around the South Pole.')
monthChoices = [1,2,3,4,5,6,7,8,9,10,11,12]
dayChoices = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31]
DIM_DICT = {'CPT':'[J m^-2]', 'KE':'[J m^-2]', 'MASS':'[kg m^-2]', 'THV':'[K]', 'TOX':'[kg m^-2]', 'TQI':'[kg m^-2]', 'TQL':'[kg m^-2]', 'TQV':'[kg m]', \
    'CLOUD':'[1]', 'DELP':'[Pa]', 'EPV':'[K m^2 kg^-1 s^-1]', 'H':'[m]', 'O3':'[kg kg^-1]', 'OMEGA':'[Pa s^-1]', 'PHIS':'[m^2 s^-2]', 'PL':'[Pa]', 'PS':'[Pa]', \
    'QI':'[kg kg^-1]', 'QL':'[kg kg^-1]', 'QV':'[kg kg^-1]', 'RH':'[1]', 'SLP':'[Pa]', 'T':'[K]', 'U':'[m s^-1]', 'V':'[m s^-1]'}
PRESS_LEVELS = np.array([0.0100, 0.0200, 0.0327, 0.0476, 0.0660, 0.0893, 0.1197, 0.1595, 0.2113, 0.2785, 0.3650, 0.4758, \
0.6168, 0.7951, 1.0194, 1.3005, 1.6508, 2.0850, 2.6202, 3.2764, 4.0766, 5.0468, 6.2168, 7.6198, \
9.2929, 11.2769, 13.6434, 16.4571, 19.7916, 23.7304, 28.3678, 33.8100, 40.1754, 47.6439, 56.3879, 66.6034, \
78.5123, 92.3657, 108.663, 127.837, 150.393, 176.930, 208.152, 244.875, 288.083, 337.500, 375.000, 412.500, \
450.000, 487.500, 525.000, 562.500, 600.000, 637.500, 675.000, 700.000, 725.000, 750.000, 775.000, 800.000, \
820.000, 835.000, 850.000, 865.000, 880.000, 895.000, 910.000, 925.000, 940.000, 955.000, 975.000, 985.000])
messages = []
modes = ['TRAP']
long_mode_names = {'TRAP':'Trapezoidal','SIMP':'Simpson\'s Rule'}
LAT_RES = 0.5
LON_RES = 0.625

def vertically_integrate(input_array,pressure_levels,mode):
    #I defined a function to take care of the masked values issue
    g = 9.80665 #m/s^2

    #Get rid of masked values
    if np.ma.isMaskedArray(input_array):
        masks = input_array.mask
        input_array = input_array[~masks]
        pressure_levels = pressure_levels[~masks]

    if mode=='TRAP': return grate.trapz(input_array,100*pressure_levels)/g
    elif mode=='SIMP': return grate.simps(input_array,100*pressure_levels,even='avg')/g

parser.add_argument('-fb2d','--file_base_2d')
parser.add_argument('-fe2d','--file_end_2d')
parser.add_argument('-fb3d','--file_base_3d')
parser.add_argument('-fe3d','--file_end_3d')
parser.add_argument('-fn','--folder_name',help='The name of the folder where the program\'s output is stored')

parser.add_argument('-v2d','--variable_2d',help='The name of the 2D MERRA-2 variable you wish to check')
parser.add_argument('-v3d','--variable_3d',help='The name of the 3D MERRA-2 variable you wish to check')
parser.add_argument('-mi','--make_images',action='store_true')
parser.add_argument('-pm','--print_messages',action='store_true')
parser.add_argument('-i','--instantaneous',action='store_true')

parser.add_argument('-sy','--start_year',type=int)
parser.add_argument('-sm','--start_month',type=int,choices=monthChoices)
parser.add_argument('-sd','--start_day',type=int,choices=dayChoices)

parser.add_argument('-ey','--end_year',type=int)
parser.add_argument('-em','--end_month',type=int,choices=monthChoices)
parser.add_argument('-ed','--end_day',type=int,choices=dayChoices)

args = parser.parse_args()

FILE_BASE_2D = args.file_base_2d
FILE_END_2D = args.file_end_2d
FILE_BASE_3D = args.file_base_3d
FILE_END_3D = args.file_end_3d
START_DATE = datetime.date(year=args.start_year,month=args.start_month,day=args.start_day)
END_DATE = datetime.date(year=args.end_year,month=args.end_month,day=args.end_day)
FOLDER_NAME = args.folder_name #The name of the folder where the images are stored, which will be created if it doesn't already exist
VAR_2D = args.variable_2d #Which meteorological quantity we're checking
VAR_3D = args.variable_3d #Which meteorological quantity we're checking
MAKE_IMAGES = args.make_images
PRINT_MESS = args.print_messages
INST = args.instantaneous

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
setsOfData_2d = {}
setsOfData_3d = {}
for date in dates:
    setsOfData_2d[date] = nc.Dataset(FILE_BASE_2D+date+FILE_END_2D)
    setsOfData_3d[date] = nc.Dataset(FILE_BASE_3D+date+FILE_END_3D)
message = "Loaded "+str(len(setsOfData_3d))+ " days of data."
if PRINT_MESS: print(message)
messages.append(message)

r_sqs,slopes,intercepts = {},{},{}
for mode in modes:
    r_sqs[mode],slopes[mode],intercepts[mode] = np.zeros(len(dates)*8),np.zeros(len(dates)*8),np.zeros(len(dates)*8)

i = 0
for date in dates:
    data_merra_3d = setsOfData_3d[date][VAR_3D]
    data_merra_2d = setsOfData_2d[date][VAR_2D]
    #Now because there are 24 times in the 2d data but only 8 in the 3d, some of the 2d data is eliminated or averaged
    if INST:
        data_merra_2d = data_merra_2d[::3,:,:]
    else:
        #TODO: There ought to be a more efficient way to do this with numpy arrays, but I can't figure it out
        averaged_array = np.zeros((8,data_merra_2d.shape[1],data_merra_2d.shape[2]))
        for j in range(8):
            averaged_array[j] = (data_merra_2d[3*j]+data_merra_2d[3*j+1]+data_merra_2d[3*j+2])/3
        data_merra_2d = averaged_array
    for t in range(8):
        if MAKE_IMAGES:
            imageName = "%svs.%s,%s.%s.t=%i.png" % (FILE_BASE_2D,FILE_BASE_3D,VAR_2D,date,t)
            message = "Making image %s..." % imageName
            if PRINT_MESS: print(message)
            messages.append(message)
            if len(modes) > 1: fig, axes = plt.subplots(nrows = len(modes),sharey=True)
            else: fig, ax = plt.subplots()
        array_3d = data_merra_3d[t]
        array_3d = np.swapaxes(np.swapaxes(array_3d,0,1),1,2) #This changes it from zyx to yxz to facilitate the integration
        array_2d = data_merra_2d[t]
        
        for k,mode in enumerate(modes):
            integrated_array = np.apply_along_axis(vertically_integrate,2,array_3d,pressure_levels = PRESS_LEVELS,mode = mode)
            
            #Remove South Pole
            integrated_array = integrated_array[1:,:]
            array_2d = array_2d[1:,:]
            
            #Flatten the arrays to make the regression simpler
            integrated_array = integrated_array.flatten()
            array_2d = array_2d.flatten()
            
            regression = sts.linregress(array_2d,integrated_array)
            slope = regression.slope
            slopes[mode][i] = slope
            
            intercept = regression.intercept
            intercepts[mode][i] = intercept
            
            r_sq = regression.rvalue**2
            r_sqs[mode][i] = r_sq
            
            message = "On %s at t=%s, we had r^2=%.3f, slope=%.3f, and intercept=%.3f." % (date,t,r_sq,slope,intercept)
            if PRINT_MESS: print(message)
            messages.append(message)
            
            if MAKE_IMAGES and len(modes) > 1:
                axes[k].plot(array_2d,integrated_array,'bo',label='Data points')
                axes[k].plot(array_2d,array_2d,'k-',label="Agreement")
                axes[k].plot(array_2d,intercept + slope*array_2d,'r-',label = 'Best fit')
                axes[k].set_title("%s:r^2=%.3f, slope=%.3f, intercept = %.3f" % (long_mode_names[mode],r_sq,slope,intercept))
                axes[k].set_xlabel('%s %s' % (VAR_2D,DIM_DICT[VAR_2D]))
                axes[k].set_ylabel('Integrated %s %s' % (VAR_3D,DIM_DICT[VAR_3D]))
                axes[k].legend()
                axes[k].grid()
            elif MAKE_IMAGES:
                ax.plot(array_2d,integrated_array,'bo',label='Data points')
                ax.plot(array_2d,array_2d,'k-',label="Agreement")
                ax.plot(array_2d,intercept + slope*array_2d,'r-',label = 'Best fit')
                ax.set_title("%s:r^2=%.3f, slope=%.3f, intercept = %.3f" % (long_mode_names[mode],r_sq,slope,intercept))
                ax.set_xlabel('%s %s' % (VAR_2D,DIM_DICT[VAR_2D]))
                ax.set_ylabel('Integrated %s %s' % (VAR_3D,DIM_DICT[VAR_3D]))
                ax.legend()
                ax.grid()
        if MAKE_IMAGES:
            fig.subplots_adjust(hspace = 0.3)
            imagePath = os.path.join(folderPath,imageName)
            fig.savefig(imagePath)
            plt.close(fig)
        i += 1
for mode in modes:
    if MAKE_IMAGES:
        imageName = '%svs.%s,%s.r_sq_histogram.%s.png' % (FILE_BASE_2D,FILE_BASE_3D,VAR_2D,mode)
        fig, ax = plt.subplots()
        ax.hist(r_sqs[mode])
        ax.set_title("r^2 Histogram")
        ax.set_xlabel('r^2')
        ax.set_ylabel('Number of times')
        ax.grid()
        imagePath = os.path.join(folderPath,imageName)
        fig.savefig(imagePath)
        plt.close(fig)
        
        imageName = '%svs.%s,%s.slope_histogram.%s.png' % (FILE_BASE_2D,FILE_BASE_3D,VAR_2D,mode)
        fig, ax = plt.subplots()
        ax.hist(slopes[mode])
        ax.set_title("Slope Histogram")
        ax.set_xlabel('slope')
        ax.set_ylabel('Number of times')
        ax.grid()
        imagePath = os.path.join(folderPath,imageName)
        fig.savefig(imagePath)
        plt.close(fig)
        
        imageName = '%svs.%s,%s.intercept_histogram.%s.png' % (FILE_BASE_2D,FILE_BASE_3D,VAR_2D,mode)
        fig, ax = plt.subplots()
        ax.hist(intercepts[mode])
        ax.set_title("Intercept Histogram")
        ax.set_xlabel('intercept')
        ax.set_ylabel('Number of times')
        ax.grid()
        imagePath = os.path.join(folderPath,imageName)
        fig.savefig(imagePath)
        plt.close(fig)
        
    avg_r_sq = np.mean(r_sqs[mode])
    avg_slope = np.mean(slopes[mode])
    avg_intercept = np.mean(intercepts[mode])
    messages.append("With mode %s, average r^2 is %.3f, average slope is %.3f, and average intercept is %.3f." % (mode,avg_r_sq,avg_slope,avg_intercept))
    
txtfile = open(os.path.join(folderPath,'0%s,%s.%s,%s.%s-%s.results.txt' % (FILE_BASE_3D,FILE_BASE_2D,VAR_3D,VAR_2D,START_DATE.strftime('%Y%m%d'),END_DATE.strftime('%Y%m%d'))),'w')
txtfile.write('\n'.join(messages))
txtfile.close()