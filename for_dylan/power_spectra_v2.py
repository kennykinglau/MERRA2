#! /usr/bin/env python

import netCDF4 as nc
import numpy as np
import scipy.integrate as grate
import scipy.interpolate as polate
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
LAT_RES = 0.5
LON_RES = 0.625
R = 6357000 #radius of Earth at the poles in meters
PRESSURE_INTERVAL = 10 #I could make this an argument, but it's so long already
NUM_OF_TICKS = 5
numOfLats = 11 #TODO: Make it smarter about this
numOfLons = 576 #TODO: Make it smarter about this
colorMap = 'Blues'

def compute_1d_psd(input_array, sample_spacing = 1.0): #This function is defined here intead of imported from demonstrate_psd.py because I figured it's better if this program can run without having to remember to move over another program as well.
    frequencies = np.fft.fftshift(np.fft.fftfreq(input_array.shape[0], d=sample_spacing))
    fourier_transform = np.fft.fftshift(np.fft.fft(input_array))
    power_spectrum = np.abs(fourier_transform)**2
    return frequencies,power_spectrum

def compute_2d_psd(input_array, x_spacing=1.0, y_spacing=1.0):
    freq_x = np.fft.fftshift(np.fft.fftfreq(input_array.shape[0], d=x_spacing))
    freq_y = np.fft.fftshift(np.fft.fftfreq(input_array.shape[1], d=y_spacing))


    fourier_transform = np.fft.fftshift(np.fft.fft2(input_array))
    power_spectrum = np.abs(fourier_transform)**2
    return freq_x,freq_y,power_spectrum

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

def latlon_2_cart(lat,lon):
    R = 6357000 #radius of Earth at the poles in meters
    x,y = R*np.cos((90-lon)*np.pi/180)*np.cos(lat*np.pi/180),R*np.sin((90-lon)*np.pi/180)*np.cos(lat*np.pi/180)
    return x,y

def custom_inverse(a): #When converting from frequency to inverse frequency, just sends zero to zero instead of having something weird happen
    if a != 0:
        return 1/a
    else:
        return 0
custom_inverse_vectorized = np.vectorize(custom_inverse)


parser.add_argument('-fb2d','--file_base_2d')
parser.add_argument('-fe2d','--file_end_2d')
parser.add_argument('-fb3d','--file_base_3d')
parser.add_argument('-fe3d','--file_end_3d')
parser.add_argument('-fn','--folder_name',help='The name of the folder where the program\'s output is stored')

parser.add_argument('-v2d','--variable_2d',help='The name of the 2D MERRA-2 variable you wish to check')
parser.add_argument('-v3d','--variable_3d',help='The name of the 3D MERRA-2 variable you wish to check')
parser.add_argument('-mi','--make_images',action='store_true')
parser.add_argument('-sl','--side_length',type=int,default=100)
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
SIDE_LENGTH = args.side_length
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

#Calculate the Cartesian coordinates of the MERRA2 grid points
data_merra_2d = setsOfData_2d[dates[0]][VAR_2D][0]
lat = np.linspace(-90,-85,num=data_merra_2d.shape[0])
lon = np.linspace(-180,180,num=data_merra_2d.shape[1])
llon,llat = np.meshgrid(lon,lat)
merra_xx,merra_yy = latlon_2_cart(llat,llon)
#Flatten so we can use griddata later
merra_xx = merra_xx.flatten()
merra_yy = merra_yy.flatten()

#Create the Cartesian grid to graph on

maxLat = LAT_RES*(numOfLats-1)
circleRadius = 0.99*R*np.cos((-90+maxLat)*np.pi/180) #using spherical coordinates; I think multiplying by 0.99 prevents it from going out of bounds for the interpolation?
conversionRatio = np.sqrt(2)*circleRadius/SIDE_LENGTH #meters per pixel of the Cartesian grid
x = (-1*circleRadius/np.sqrt(2))+conversionRatio*np.arange(SIDE_LENGTH)
y = (-1*circleRadius/np.sqrt(2))+conversionRatio*np.arange(SIDE_LENGTH)
graphing_xx,graphing_yy = np.meshgrid(x,y)

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
        
        #Interpolate 3D variable at South Pole to have even spacing for FFT
        SP_array = data_merra_3d[t,:,0,0]
        f = polate.interp1d(PRESS_LEVELS,SP_array)
        pressure_samples = np.arange(0.0100,985.000,10) #I think the pressures on the grid are not actual, but refer to what they would be if the pressure at the surface were 1000 hPa
        
        #Interpolate the 2D variable
        array_2d = data_merra_2d[t]
        array_2d = array_2d.flatten() #For griddata, we first need to flatten
        graphing_array = polate.griddata(points=(merra_xx,merra_yy),values=array_2d,xi=(graphing_xx,graphing_yy))
        
        if MAKE_IMAGES:
            imageName = "%s,%s%s,%s.%s.t=%i.png" % (FILE_BASE_2D,FILE_BASE_3D,VAR_2D,VAR_3D,date,t)
            message = "Making image %s..." % imageName
            if PRINT_MESS: print(message)
            messages.append(message)
            fig, ((press_graph,oned_power),(twod_graph,twod_power)) = plt.subplots(2,2)
            fig.suptitle("Power Spectra of %s and %s on %s at t=%i" % (VAR_3D,VAR_2D,date,t),y=1.15)
            fig.subplots_adjust(hspace=0.6,wspace=0.5)
            
            #Pressure plot
            press_graph.plot(pressure_samples,f(pressure_samples),'b-')
            tick_positions = np.linspace(np.amin(pressure_samples),np.amax(pressure_samples),NUM_OF_TICKS)
            press_graph.set_xticks(tick_positions)
            
            pressTitle = "%s %s vs. Pressure" % (VAR_3D,DIM_DICT[VAR_3D])
            press_graph.set_title(pressTitle,y=1.15)
            press_graph.set_xlabel("Pressure [hPa]")
            press_graph.set_ylabel("%s %s" % (VAR_3D,DIM_DICT[VAR_3D]))
            press_graph.grid()
            
            #1-dimensional power spectrum
            psd = compute_1d_psd(f(pressure_samples),sample_spacing = PRESSURE_INTERVAL)
            powers = psd[1]
            sample_freqs = psd[0]
            #Convert from inverse hPa to regular hPa
            
            sample_distances = custom_inverse_vectorized(sample_freqs)
            oned_power.plot(sample_distances,powers,'b-')
            tick_positions = np.linspace(np.amin(sample_distances),np.amax(sample_distances),NUM_OF_TICKS)
            oned_power.set_xticks(tick_positions)
            
            freqTitle = "Power Spectrum of %s" % VAR_3D
            oned_power.set_title(freqTitle,y=1.15)
            oned_power.set_xlabel("Inverse Frequency [hPa]")
            oned_power.set_ylabel("Power")
            oned_power.set_ybound(lower=0,upper=np.amax(powers))
            oned_power.grid()
            
            twod_graph.imshow(graphing_array,origin="lower",cmap=colorMap)
            twodgraphTitle = "%s %s" % (VAR_2D,DIM_DICT[VAR_2D])
            twod_graph.set_title(twodgraphTitle)
            
            #Axis labels in kilometers
            twod_graph.set_xlabel("km in the plane")
            twod_graph.set_ylabel("km in the plane")
            
            tick_positions = np.arange(0,SIDE_LENGTH,SIDE_LENGTH/NUM_OF_TICKS)
            tick_labels = np.around(tick_positions*conversionRatio/1000).astype(int)
            twod_graph.set_xticks(tick_positions)
            twod_graph.set_xticklabels(tick_labels)
            twod_graph.set_yticks(tick_positions)
            twod_graph.set_yticklabels(tick_labels)
            twod_graph.grid()
            
            twoDpsd = compute_2d_psd(graphing_array,x_spacing = conversionRatio/1000,y_spacing=conversionRatio/1000)
            x_frequencies = twoDpsd[0]
            y_frequencies = twoDpsd[1]
            
            powers_2d = twoDpsd[2]
            xx_frequencies, yy_frequencies = np.meshgrid(x_frequencies,y_frequencies)
            #Convert from inverse km to regular km
            xx_distances,yy_distances = custom_inverse_vectorized(xx_frequencies),custom_inverse_vectorized(yy_frequencies)
            
            #Now we need to sort the arrays so it graphs nicely; some of the following code is from https://github.com/numpy/numpy/issues/4724
            xx_sort_indices = np.argsort(xx_distances,axis=1)
            i = np.arange(len(xx_sort_indices))[:, np.newaxis]
            
            xx_distances = xx_distances[i,xx_sort_indices]
            yy_distances = yy_distances[i,xx_sort_indices]
            powers_2d = powers_2d[i,xx_sort_indices]

            yy_sort_indices = np.argsort(yy_distances,axis=0) #axis=0 yields the correct indices
            j = np.arange(len(yy_sort_indices))[np.newaxis,:]
            
            xx_distances = xx_distances[yy_sort_indices,j]
            yy_distances = yy_distances[yy_sort_indices,j]
            powers_2d = powers_2d[yy_sort_indices,j]
            
            twod_power.imshow(np.transpose(powers_2d),origin='lower',cmap=colorMap,aspect = 'auto')
            twod_freqs_title = "Power Spectrum of %s" % VAR_2D
            twod_power.set_title(twod_freqs_title,y=1.15)
            
            twod_power.set_xlabel("Inverse x-frequency [km]")
            twod_power.set_ylabel("Inverse y-frequency [km]")
            x_tick_positions = np.arange(0,len(x_frequencies),int(len(x_frequencies)/NUM_OF_TICKS))
            y_tick_positions = np.arange(0,len(y_frequencies),int(len(y_frequencies)/NUM_OF_TICKS))

            x_tick_labels = np.around(xx_distances[0][x_tick_positions],1)
            j = np.arange(len(y_tick_positions))[:, np.newaxis]
            y_tick_labels = np.around(yy_distances[y_tick_positions,j][0],1)

            twod_power.set_xticks(x_tick_positions)
            twod_power.set_yticks(y_tick_positions)
            twod_power.set_xticklabels(x_tick_labels)
            twod_power.set_yticklabels(y_tick_labels)
            twod_power.grid()
            
            imagePath = os.path.join(folderPath,imageName)
            fig.savefig(imagePath)
            plt.close(fig)
            
#Make test graphics
imageName = "%s,%s%s,%s.test.png" % (FILE_BASE_2D,FILE_BASE_3D,VAR_2D,VAR_3D)
message = "Making image %s..." % imageName
if PRINT_MESS: print(message)
messages.append(message)
fig, ((press_graph,oned_power),(twod_graph,twod_power)) = plt.subplots(2,2)
fig.suptitle("Power Spectra of %s and %s on %s at t=%i" % (VAR_3D,VAR_2D,date,t),y=1.15)
fig.subplots_adjust(hspace=0.6,wspace=0.5)

test_data = np.sin(pressure_samples/40)+2*np.cos(pressure_samples/80)
graphing_array = np.sin(graphing_xx/200)+2*np.cos(graphing_yy/400)

#Pressure plot
test_data = np.sin(pressure_samples/40)+2*np.cos(pressure_samples/80)
press_graph.plot(pressure_samples,test_data,'b-')
tick_positions = np.linspace(np.amin(pressure_samples),np.amax(pressure_samples),NUM_OF_TICKS)
press_graph.set_xticks(tick_positions)

pressTitle = "Test vs. Pressure"
press_graph.set_title(pressTitle,y=1.15)
press_graph.set_xlabel("Pressure [hPa]")
press_graph.set_ylabel("Test")
press_graph.grid()

#1-dimensional power spectrum
psd = compute_1d_psd(test_data,sample_spacing = PRESSURE_INTERVAL)
powers = psd[1]
sample_freqs = psd[0]
#Convert from inverse hPa to regular hPa

sample_distances = custom_inverse_vectorized(sample_freqs)
oned_power.plot(sample_distances,powers,'b-')
tick_positions = np.linspace(np.amin(sample_distances),np.amax(sample_distances),NUM_OF_TICKS)
oned_power.set_xticks(tick_positions)

freqTitle = "Power Spectrum of Test"
oned_power.set_title(freqTitle,y=1.15)
oned_power.set_xlabel("Inverse Frequency [hPa]")
oned_power.set_ylabel("Power")
oned_power.set_ybound(lower=0,upper=np.amax(powers))
oned_power.grid()

twod_graph.imshow(graphing_array,origin="lower",cmap=colorMap)
twodgraphTitle = "Test"
twod_graph.set_title(twodgraphTitle)

#Axis labels in kilometers
twod_graph.set_xlabel("km in the plane")
twod_graph.set_ylabel("km in the plane")

tick_positions = np.arange(0,SIDE_LENGTH,SIDE_LENGTH/NUM_OF_TICKS)
tick_labels = np.around(tick_positions*conversionRatio/1000).astype(int)
twod_graph.set_xticks(tick_positions)
twod_graph.set_xticklabels(tick_labels)
twod_graph.set_yticks(tick_positions)
twod_graph.set_yticklabels(tick_labels)
twod_graph.grid()

twoDpsd = compute_2d_psd(graphing_array,x_spacing = conversionRatio/1000,y_spacing=conversionRatio/1000)
x_frequencies = twoDpsd[0]
y_frequencies = twoDpsd[1]

powers_2d = twoDpsd[2]
xx_frequencies, yy_frequencies = np.meshgrid(x_frequencies,y_frequencies)
#Convert from inverse km to regular km
xx_distances,yy_distances = custom_inverse_vectorized(xx_frequencies),custom_inverse_vectorized(yy_frequencies)

#Now we need to sort the arrays so it graphs nicely; some of the following code is from https://github.com/numpy/numpy/issues/4724
xx_sort_indices = np.argsort(xx_distances,axis=1)
i = np.arange(len(xx_sort_indices))[:, np.newaxis]

xx_distances = xx_distances[i,xx_sort_indices]
yy_distances = yy_distances[i,xx_sort_indices]
powers_2d = powers_2d[i,xx_sort_indices]

yy_sort_indices = np.argsort(yy_distances,axis=0) #axis=0 yields the correct indices
j = np.arange(len(yy_sort_indices))[np.newaxis,:]

xx_distances = xx_distances[yy_sort_indices,j]
yy_distances = yy_distances[yy_sort_indices,j]
powers_2d = powers_2d[yy_sort_indices,j]

twod_power.imshow(np.transpose(powers_2d),origin='lower',cmap=colorMap,aspect = 'auto')
twod_freqs_title = "Power Spectrum of Test"
twod_power.set_title(twod_freqs_title,y=1.15)

twod_power.set_xlabel("Inverse x-frequency [km]")
twod_power.set_ylabel("Inverse y-frequency [km]")
x_tick_positions = np.arange(0,len(x_frequencies),int(len(x_frequencies)/NUM_OF_TICKS))
y_tick_positions = np.arange(0,len(y_frequencies),int(len(y_frequencies)/NUM_OF_TICKS))

x_tick_labels = np.around(xx_distances[0][x_tick_positions],1)
j = np.arange(len(y_tick_positions))[:, np.newaxis]
y_tick_labels = np.around(yy_distances[y_tick_positions,j][0],1)

twod_power.set_xticks(x_tick_positions)
twod_power.set_yticks(y_tick_positions)
twod_power.set_xticklabels(x_tick_labels)
twod_power.set_yticklabels(y_tick_labels)
twod_power.grid()

imagePath = os.path.join(folderPath,imageName)
fig.savefig(imagePath)
plt.close(fig)

txtfile = open(os.path.join(folderPath,'0%s,%s.%s,%s.%s-%s.results.txt' % (FILE_BASE_3D,FILE_BASE_2D,VAR_3D,VAR_2D,START_DATE.strftime('%Y%m%d'),END_DATE.strftime('%Y%m%d'))),'w')
txtfile.write('\n'.join(messages))
txtfile.close()