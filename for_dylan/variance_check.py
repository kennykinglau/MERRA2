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

parser.add_argument('-fbi','--file_base_instantaneous')
parser.add_argument('-fei','--file_end_instantaneous')
parser.add_argument('-fbt','--file_base_time_averaged')
parser.add_argument('-fet','--file_end_time_averaged')

parser.add_argument('-fn','--folder_name',help='The name of the folder where the program\'s output is stored')

parser.add_argument('-v','--variable',help='The name of the 3D MERRA-2 variable you wish to use')
parser.add_argument('-mi','--make_images',action='store_true')
parser.add_argument('-ns','--number_of_samples',type=int,default=5000)
parser.add_argument('-pm','--print_messages',action='store_true')

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

START_DATE = datetime.date(year=args.start_year,month=args.start_month,day=args.start_day)
END_DATE = datetime.date(year=args.end_year,month=args.end_month,day=args.end_day)
FOLDER_NAME = args.folder_name #The name of the folder where the images are stored, which will be created if it doesn't already exist

VAR = args.variable #Which meteorological quantity we're checking
MAKE_IMAGES = args.make_images
NUM_OF_SAMPS = args.number_of_samples
PRINT_MESS = args.print_messages

#This creates the folder where everything will go

if MAKE_IMAGES:
    dir_path = os.path.dirname(os.path.realpath(__file__))
    folderPath = os.path.join(dir_path,FOLDER_NAME)
if not os.path.isdir(folderPath):
    os.mkdir(folderPath)

NUM_OF_TIMES = 1000
variances = np.zeros(NUM_OF_TIMES)
integrals = np.zeros(NUM_OF_TIMES)

centers_sq = np.zeros(NUM_OF_TIMES)
differences = np.zeros(NUM_OF_TIMES)
for i in range(NUM_OF_TIMES):
    if (i+1)%50==0:
        percent = 100*(i+1)/NUM_OF_TIMES
        print('%i%% done' % percent)
    scale = np.random.uniform(0.1,10)
    NUM_OF_SAMPS = int(np.random.uniform(100,10000))
    intervalLength = np.random.uniform(1,20)
    start = np.random.uniform(-20,20)
    stop = np.random.uniform(21,61)
    center = np.random.uniform(-5,5)
    
    height_samples, height_interval = np.linspace(start,stop,NUM_OF_SAMPS,retstep = True)
    intervalLength = stop-start
    SP_array = np.random.normal(loc = center,scale = scale,size = NUM_OF_SAMPS)
    SP_array = SP_array - np.mean(SP_array)
    variance = np.var(SP_array)
    
    sample_freqs,power = compute_1d_psd(SP_array,sample_spacing = height_interval)
    psd_integral = grate.trapz(power,sample_freqs)/(NUM_OF_SAMPS**2)*intervalLength #I don't know why you need to divide by (NUM_OF_SAMPS**2) and multiply by the interval length, but that made it work
    
    variances[i] = variance
    integrals[i] = psd_integral
    
    centers_sq[i] = center**2
    differences[i] = psd_integral-variance

regression = sts.linregress(variances,integrals)

slope = regression.slope
intercept = regression.intercept
r_sq = regression.rvalue**2

fig, ax = plt.subplots()
ax.plot(variances,integrals,'bo',label='Data points')
ax.plot(variances,variances,'k-',label="Agreement")
ax.plot(variances,intercept + slope*variances,'r-',label = 'Best fit')

print("r^2: "+str(r_sq))
print("Slope: "+str(slope))
print("Intercept: "+str(intercept))


ax.set_title("r^2=%.3f, slope=%.3f, intercept = %.3f" % (r_sq,slope,intercept))
ax.set_xlabel('Variance')
ax.set_ylabel('PSD Integral')

ax.legend()
ax.grid()

imageName = 'variance_vs_integral_corrected.png'
imagePath = os.path.join(folderPath,imageName)
fig.savefig(imagePath)
plt.close(fig)