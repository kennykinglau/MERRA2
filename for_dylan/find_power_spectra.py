#! /usr/bin/env python

import netCDF4 as nc
import numpy as np
import math
import datetime
import os
#import demonstrate_psd as psd
from matplotlib import pyplot as plt

LAT_RES = 0.5 #latitude resolution of the MERRA grid
LON_RES = 0.625 #longitude resolution
PRESS_LEVELS = np.array([1000,975,950,925,900,875,850,825,800,775,750,725,700,650,600,550,500,450,400,350,300,250,200,150,100,70,50,40,30,20,10,7,5,4,3,2,1,0.7,0.5,0.4,0.3,0.1])
#The above are in hPa
is3D = True
FILE_BASE = 'MERRA2_400.tavg3_3d_cld_Np.'
startDate = datetime.date(year=2019,month=3,day=1)
endDate = datetime.date(year=2019,month=3,day=7)
FILE_END = '.SUB.nc'
FOLDER_NAME = 'Power_Spectra_Adjusted_Spectrum_Graphs' #The name of the folder where the images are stored, which will be created if it doesn't already exist
param = 'QI' #Which meteorological quantity we're graphing

colorMap = 'Blues'
adjustPowerSpectra = True
makeImages = True

#This program assumes we have a circle cented at the South Pole

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

def vertically_integrate(input_array,pressure_levels):
	#Using trapezoidal sums
	g = 9.81 #m/s^2
	#Get rid of masked values
	if np.ma.isMaskedArray(input_array):
		masks = input_array.mask
		input_array = input_array[~masks]
		pressure_levels = pressure_levels[~masks]
	
	ans = 0
	for i in range(len(input_array)-1):
		avg = (input_array[i]+input_array[i+1])/2
		ans += avg*abs(pressure_levels[i+1]-pressure_levels[i])*100 #Multiplying by 100 converts from hPa to Pa, which are needed for the units to work out.
	return ans/g

print("Loading data...")
date = startDate #This generates the list of dates in the proper format from the start and end dates
interval = datetime.timedelta(days=1)
dates = []
while date <= endDate:
	dates.append(date.isoformat().replace('-',''))
	date = date + interval

setsOfData = {}
for date in dates:
	setsOfData[date] = nc.Dataset(FILE_BASE+date+FILE_END)
print("Loaded "+str(len(setsOfData))+ " days of data.")
data_merra = setsOfData[dates[0]] #load a dataset to determine the dimensions

numOfTimes = data_merra[param].shape[0]
numOfLats = data_merra[param].shape[2]
numOfLons = data_merra[param].shape[3]

print("Preparing to make graphs...")
maxVal = 0 #This finds the maximum of the whole dataset to scale the axes
maxIntVal = 0
maxPwr = 0
maxIntPwr = 0
pressure_interval = 50
SP_arrays = {}
integrated_arrays = {}
#This code finds the maxima to scale the axes and also stores the relevant arrays in a dictionary of dictionaries
for date in dates:
	data_merra = setsOfData[date]
	SP_arrays[date] = {}
	integrated_arrays[date] = {}
	for t in range(numOfTimes):
		ice_array = data_merra[param][t]
		#Because QI appears to be consistent among the different South Pole longitudes, I'll just use the first one.
		SP_array = ice_array[:,0,0]
		#Remove some of the points so that the spacing is at even pressure intervals, and remove pressures below 200 hPa, at which nothing significant seems to occur
		SP_array = SP_array[13:23]
		SP_arrays[date][t] = SP_array
		
		ice_array = np.swapaxes(np.swapaxes(ice_array,0,1),1,2)
		integrated_array = np.apply_along_axis(vertically_integrate,2,ice_array,pressure_levels = PRESS_LEVELS)
		integrated_arrays[date][t] = integrated_array
		
		thisMaxVal = np.amax(SP_array)
		thisMaxIntVal = np.amax(integrated_array)
		thisMaxPwr = np.amax(compute_1d_psd(SP_array,sample_spacing=pressure_interval)[1])
		thisMaxIntPwr = np.amax(compute_2d_psd(integrated_array,x_spacing=LON_RES,y_spacing=LAT_RES)[2])
		if thisMaxVal > maxVal:
			maxVal = thisMaxVal
		if thisMaxIntVal > maxIntVal:
			maxIntVal = thisMaxIntVal
		if thisMaxPwr > maxPwr:
			maxPwr = thisMaxPwr
		if thisMaxIntPwr > maxIntPwr:
			maxIntPwr = thisMaxIntPwr
		
#This creates the folder where everything will go
if makeImages or makeMovie:
	dir_path = os.path.dirname(os.path.realpath(__file__))
	folderPath = os.path.join(dir_path,FOLDER_NAME)
	if not os.path.isdir(folderPath):
		os.mkdir(folderPath)

#Generates an array with the pressure levels used
first_pressure = 650 #hPa
numOfPressures = SP_arrays[dates[0]][0].shape[0]
pressure_interval = 50 #hPa
pressures = np.arange(numOfPressures)
pressures = first_pressure - pressure_interval*pressures

#This prints out some tests of the power spectrum
numbs = np.linspace(-2*np.pi,2*np.pi,100)
test_array = np.sin(2*np.pi*numbs) + 0.5*np.cos(3*np.pi*numbs)

y = np.arange(integrated_arrays[dates[0]][0].shape[0])
x = np.arange(integrated_arrays[dates[0]][0].shape[1])
xx,yy = np.meshgrid(x,y)
twod_test_array = np.sin(xx)+np.sin(3*xx)-2*np.cos(2*yy)

fig, ((press_graph,oned_power),(twod_graph,twod_power)) = plt.subplots(2,2)
fig.suptitle("Test")
imageName = "test.png"
press_graph.plot(numbs,test_array,'b-')

psd = compute_1d_psd(test_array,sample_spacing = 4*np.pi/100)
bar_width = (np.amax(psd[0])-np.amin(psd[0]))/(2*psd[0].shape[0])
oned_power.bar(psd[0],psd[1],width=bar_width)

twod_graph.imshow(twod_test_array,origin='lower',cmap=colorMap,aspect='auto')

twoDpsd = compute_2d_psd(twod_test_array,x_spacing = 1.0,y_spacing=1.0)
twod_power.imshow(np.transpose(twoDpsd[2]),origin='lower',cmap=colorMap,aspect = 'auto')

if makeImages:
	imagePath = os.path.join(folderPath,imageName)
	fig.savefig(imagePath)
	plt.close(fig)
#This actually makes the images
for date in dates: 
	data_merra = setsOfData[date]
	for t in range(numOfTimes):
		if numOfTimes>=10 and t<10:
			timeString = "0"+str(t)
		else:
			timeString = str(t)
		imageName = FILE_BASE+param+'.'+date+'.t='+timeString+'.png'
		print("Making image "+imageName+"...")
		SP_array = SP_arrays[date][t]
		integrated_array = integrated_arrays[date][t]
		
		fig, ((press_graph,oned_power),(twod_graph,twod_power)) = plt.subplots(2,2)
		figTitle = date+" at t="+timeString
		fig.suptitle(figTitle,y=1.0)
		fig.subplots_adjust(wspace=0.4,hspace=0.6)
		
		press_graph.plot(pressures,SP_array,'b-')
		pressTitle = param+" vs. Pressure"
		press_graph.set_title(pressTitle,y=1.15)
		press_graph.set_xlabel("Pressure [hPa]")
		press_graph.set_ylabel(param)
		press_graph.set_ybound(upper=maxVal)
		
		psd = compute_1d_psd(SP_array,sample_spacing = pressure_interval)
		bar_width = (np.amax(psd[0])-np.amin(psd[0]))/(2*psd[0].shape[0])
		oned_power.bar(psd[0],psd[1],width=bar_width)
		freqTitle = "Power Spectrum of "+param
		oned_power.set_title(freqTitle,y=1.15)
		oned_power.set_xlabel("Frequency [(hPa)^-1]")
		oned_power.set_ylabel("Power")
		if adjustPowerSpectra:
			oned_power.set_ybound(lower=0,upper=maxPwr)
		else:
			oned_power.set_ybound(lower=0,upper=np.amax(psd[1]))
		twod_graph.imshow(integrated_array,origin='lower',cmap=colorMap,vmax=maxIntVal,aspect='auto')
		twodTitle = "Vertical integration of "+param
		twod_graph.set_title(twodTitle,y=1.15)
		
		twod_graph.set_xlabel("Longitude")
		twod_graph.set_ylabel("Latitude")
		x_tick_positions = np.arange(0,numOfLons,int(numOfLons/4))
		y_tick_positions = np.arange(0,numOfLats,int(numOfLats/4))
		x_tick_labels = -180+LON_RES*x_tick_positions
		y_tick_labels = -90+LAT_RES*y_tick_positions
		twod_graph.set_xticks(x_tick_positions)
		twod_graph.set_yticks(y_tick_positions)
		twod_graph.set_xticklabels(x_tick_labels)
		twod_graph.set_yticklabels(y_tick_labels)
		#fig.colorbar(mappable=im1,cax=twod_graph)
		#For some reason the colorbar is covering up the whole heatmap
		
		twoDpsd = compute_2d_psd(integrated_array,x_spacing = LON_RES,y_spacing=LAT_RES)
		if adjustPowerSpectra:
			twod_power.imshow(np.transpose(twoDpsd[2]),origin='lower',cmap=colorMap,vmax = maxIntPwr,aspect = 'auto')
		else:
			twod_power.imshow(np.transpose(twoDpsd[2]),origin='lower',cmap=colorMap,aspect = 'auto')
		twod_freqs_title = "Power Spectrum of Integrated "+param
		twod_power.set_title(twod_freqs_title,y=1.15)
		
		twod_power.set_xlabel("x-frequency [(longitude)^-1]")
		twod_power.set_ylabel("y-frequency [(latitude)^-1]")
		x_tick_positions = np.arange(0,len(twoDpsd[0]),int(len(twoDpsd[0])/4))
		y_tick_positions = np.arange(0,len(twoDpsd[1]),int(len(twoDpsd[1])/4))
		x_tick_labels = np.around(twoDpsd[0][x_tick_positions],2)
		y_tick_labels = np.around(twoDpsd[1][y_tick_positions],2)
		twod_power.set_xticks(x_tick_positions)
		twod_power.set_yticks(y_tick_positions)
		twod_power.set_xticklabels(x_tick_labels)
		twod_power.set_yticklabels(y_tick_labels)
		#fig.colorbar(mappable=im2,cax=twod_power)
		#For some reason the colorbar is covering up the whole heatmap
		
		if makeImages:
			imagePath = os.path.join(folderPath,imageName)
			fig.savefig(imagePath)
		plt.close(fig)
