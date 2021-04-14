import netCDF4 as nc
import numpy as np
import math
import datetime
import os
from matplotlib import pyplot as plt

R = 6357000 #radius of Earth at the poles in meters
LAT_RES = 0.5 #latitude resolution of the MERRA grid
LON_RES = 0.625 #longitude resolution

is3D = True
FILE_BASE = 'MERRA2_400.tavg3_3d_cld_Np.'
startDate = datetime.date(year=2019,month=3,day=1)
endDate = datetime.date(year=2019,month=3,day=1)
FILE_END = '.SUB.nc'
FOLDER_NAME = 'WedTest3D' #The name of the folder where the images are stored, which will be created if it doesn't already exist
param = 'QI' #Which meteorological quantity we're graphing
elevation = 18 #Only needed for 3D data, 10=750hPa,16=500hPa,21=250hPa

sideLength = 400 #number of pixels used on the Cartesian grid; should be an integer multiple of numOfTicks
numOfTicks = 5
colorMap = 'Blues'
makeMovie = False #Not yet implemented
makeImages = True
checkSouthPoleConsistency = True

#This program assumes we have a circle cented at the South Pole
#TODO clarify variable names for when they're actually latitude and longitude and when they're indices for the MERRA2 numpy array

def convert2Cartesian(lat,lon): #using spherical coordinates
	R = 6357000 #radius of Earth at the poles in meters
	x,y = R*np.cos((90-lon)*np.pi/180)*np.cos(lat*np.pi/180),R*np.sin((90-lon)*np.pi/180)*np.cos(lat*np.pi/180)
	return x,y

def convert2LatLon(x,y):
	R = 6357000
	z = -np.sqrt(R**2-x**2-y**2)
	#First find theta, then convert to longitude
	theta = np.arctan2(y,x)*180/np.pi
	if theta > -90:
		lon = 90-theta
	elif theta <= -90:
		lon = -270-theta
	if np.sqrt(x**2+y**2) != 0:
		lat = np.arctan(z/np.sqrt(x**2+y**2))*180/np.pi
	else:
		lat = -90
	return lat,lon

#The following code establishes the correspondence between the Cartesian grid used for graphing and the spherical coordinate grid of the original datase
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

numOfTimes = len(data_merra[param])
if is3D:
	numOfLats = len(data_merra[param][0][elevation])
	numOfLons = len(data_merra[param][0][elevation][0])
else:
	numOfLats = len(data_merra[param][0])
	numOfLons = len(data_merra[param][0][0])

print("Creating Cartesian grid to graph on...")
cartesianGrid = np.zeros((sideLength,sideLength,2)) #at each point it'll have the lat and lon indices of the nearest point on the MERRA2 grid
maxLat = LAT_RES*(numOfLats-1)
circleRadius = R*np.cos((-90+maxLat)*np.pi/180) #using spherical coordinates
conversionRatio = np.sqrt(2)*circleRadius/sideLength #meters per pixel of the Cartesian grid

merra2_in_Cartesian = np.zeros((numOfLats,numOfLons,2)) #At each latitude, longitude point it gives the x,y coordinates
#TODO: is there a more efficient way to do this using numpy arrays?
for lat in range(numOfLats):
	for lon in range(numOfLons):
		x,y = convert2Cartesian(-90+LAT_RES*lat,-180+LON_RES*lon)
		merra2_in_Cartesian[lat][lon][0] = x
		merra2_in_Cartesian[lat][lon][1] = y

for i in range(sideLength):
	for j in range(sideLength):
		#Finds the four closest MERRA-2 grid points and finds the closest among them
		currentx,currenty = conversionRatio*(j-(sideLength/2)),conversionRatio*(i-(sideLength/2)) #Not sure if one choice of which is i and which is j is superior to another, but this choice seems to lead to better outcomes
		currentlat,currentlon = convert2LatLon(currentx,currenty)

		floorlat,floorlon = math.floor(currentlat/LAT_RES)*LAT_RES,math.floor(currentlon/LON_RES)*LON_RES
		minDistance = math.inf
		lon,lat = round((floorlon+180)/(LON_RES)),round((floorlat+90)/(LAT_RES))

		if lat >= numOfLats-1: #Prevents us from going outside of the circle
			lat -= 1
		
		if lon != numOfLons-1: #This avoids going over the bounds of the array containing the lattitudes and longitudes
			for k in range(2):
				for ell in range(2):
					x,y = merra2_in_Cartesian[lat+k][lon+ell][0],merra2_in_Cartesian[lat+k][lon+ell][1]
					currentDistance = np.sqrt((currentx-x)**2+(currenty-y)**2)
					if currentDistance < minDistance:
						minDistance = currentDistance
						minLat,minLon = lat+k,lon+ell
			cartesianGrid[i][j][0] = minLat
			cartesianGrid[i][j][1] = minLon
		else:
			edges = (numOfLons-1,0)
			for k in range(2):
				for lon in edges:
					x,y = merra2_in_Cartesian[lat+k][lon][0],merra2_in_Cartesian[lat+k][lon][1]
					currentDistance = np.sqrt((currentx-x)**2+(currenty-y)**2)
					if currentDistance < minDistance:
						minDistance = currentDistance
						minlat,minLon = lat+k,lon
			cartesianGrid[i][j][0] = minLat
			cartesianGrid[i][j][1] = minLon

ice_graph = np.zeros((sideLength,sideLength))

print("Determining consistent scale for color map...")
maxIce = 0 #This finds the maximum of the whole dataset to scale the color maps consistently
lastProgressUpdate = 0
for date in dates:
	data_merra = setsOfData[date]
	for t in range(numOfTimes):
		if is3D:
			ice_array = data_merra[param][t][elevation][:][:]
		else:
			ice_array = data_merra[param][t][:][:]
		for i in range(sideLength):
			for j in range(sideLength):
				lat = int(round(cartesianGrid[i][j][0]))
				lon = int(round(cartesianGrid[i][j][1]))
				if ice_array[lat][lon] > maxIce:
					maxIce = ice_array[lat][lon]
	print("Checked "+str(dates.index(date)+1)+" out of "+str(len(dates))+" days.")
print("The maximum value for "+param+" in the box is "+str(maxIce)+".")

if makeImages and makeMovie:
	print("Making images and movie...")
elif makeImages:
	print("Making images...")
elif makeMovie:
	print("Making movie...")
#This creates the folder where everything will go

if makeImages or makeMovie:
	dir_path = os.path.dirname(os.path.realpath(__file__))
	folderPath = os.path.join(dir_path,FOLDER_NAME)
	if not os.path.isdir(folderPath):
		os.mkdir(folderPath)

#Adds longitude and latitude graphs for scale
lat_graph = np.zeros((sideLength,sideLength))
lon_graph = np.zeros((sideLength,sideLength))
for i in range(sideLength):
	for j in range(sideLength):
		lat = int(round(cartesianGrid[i][j][0]))
		lon = int(round(cartesianGrid[i][j][1]))
		lat_graph[i][j] = -90+LAT_RES*lat
		lon_graph[i][j] = -180+LON_RES*lon

plt.imshow(lat_graph,origin="lower")
plt.title('Latitude in the region being graphed')
plt.colorbar()
plt.savefig(os.path.join(folderPath,'latitude.png'))

plt.clf()
plt.imshow(lon_graph,origin="lower")
plt.title('Longitude in the region being graphed')
plt.colorbar()
plt.savefig(os.path.join(folderPath,'longitude.png'))
plt.clf()
#TODO Get rid of number labelling of axes

haveASouthPoleGraph = False
for date in dates:
	data_merra = setsOfData[date]
	for t in range(1): #Change back!
		if numOfTimes>=10 and t<10:
			timeString = "0"+str(t)
		else:
			timeString = str(t)
		if is3D:
			imageName = FILE_BASE+param+'.'+date+'.elev='+str(elevation)+'.t='+timeString+'.png'
			ice_array = data_merra[param][t][elevation][:][:]
		else:
			imageName = FILE_BASE+param+'.'+date+'.t='+timeString+'.png'
			ice_array = data_merra[param][t][:][:]
		
		print("Making image "+imageName+"...")
		
		if checkSouthPoleConsistency and not haveASouthPoleGraph:
			print("Checking South Pole Consistency...")
			southPolePoints = ice_array[0]
			upper = np.amax(southPolePoints)
			lower = np.amin(southPolePoints)
			if upper == lower:
				print("Consistency confirmed.")
			else:
				print("The value at the South Pole varies from "+str(lower)+" to "+str(upper)+".")
				firstCirclePoints = ice_array[1]
				secondCirclePoints = ice_array[2]
				x = np.arange(0,numOfLons)
				y1 = southPolePoints[x]
				y2 = firstCirclePoints[x]
				y3 = secondCirclePoints[x]
				x = -180+LON_RES*x
				plt.plot(x,y1,"r-")
				plt.plot(x,y2,"g-")
				plt.plot(x,y3,"b-")
				plt.title("Red is the South Pole, green is at 89.5 S, blue is at 89 S")
				plt.xlabel("Longitude")
				plt.ylabel(param)
				plt.savefig(os.path.join(folderPath,'SouthPoleGraph.png'))
				haveASouthPoleGraph = True
				plt.clf()
		
		for i in range(sideLength):
			for j in range(sideLength):
				lat = int(round(cartesianGrid[i][j][0]))
				lon = int(round(cartesianGrid[i][j][1]))
				ice_graph[i][j] = ice_array[lat][lon]
		plt.imshow(ice_graph,origin="lower",vmin=0,vmax=maxIce,cmap=colorMap)
		if date==dates[0] and t==0: #So only one colorbar gets added
			plt.colorbar()
		if is3D:
			plt.title(param+' on '+date+', elevation='+str(elevation)+', time='+str(t))
		else:
			plt.title(param+' on '+date+', time='+str(t))
		
		#Axis labels in kilometers
		plt.xlabel("km, measured in the plane onto which the sphere is being projected")
		plt.ylabel("km, also measured in the plane")
		
		tick_positions = np.arange(0,sideLength,sideLength/numOfTicks)
		tick_labels = np.around(tick_positions*conversionRatio/1000).astype(int)
		plt.xticks(tick_positions,tick_labels)
		plt.yticks(tick_positions,tick_labels)
		 
		if makeImages:
			imagePath = os.path.join(folderPath,imageName)
			plt.savefig(imagePath)
		#Maybe add movie-making capability?
