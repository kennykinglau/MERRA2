from collections import Counter
import pickle
import glob
import datetime
import time
import os
import merra2Player as m2p
#from pylab import *

currentSite = "ChajnantorCerro"
m = m2p.merra2Player(site = currentSite)
#for i in range(1996,2015):
#    print i, m.findFailedSpectra(str(i))

spData = {}
cpData = {}
keys = ['t','pwv', 'BK100','BK150','BK220','tipper850']
for i in keys:
    spData[i] = []
    cpData[i] = []

spData = Counter(spData)
cpData = Counter(cpData)

#South Pole
for i in glob.glob("/n/home10/anwang16/merra2_work/merra2_products/%s/farmfiles/pickles/%s*"%(currentSite,currentSite)):
    spData = spData+Counter(pickle.load(open(i,'rb')))
spData = dict(spData)
#print len(spData['t'])
#print len(spData['pwv'])
#print len(spData['tipper850'])

#Chajnantor
#for i in glob.glob("/n/home10/anwang16/merra2_work/merra2_products/ChajnantorPlateau/farmfiles/pickles/ChajnantorPlateau*"):
#    cpData = cpData+Counter(pickle.load(open(i,'rb')))
#cpData = dict(cpData)


#Finding missing data
delta = datetime.timedelta(days = 1)
start = datetime.datetime(1996,1,1)
end = datetime.datetime(2017,6,1)
f = open("Missing%sFiles.txt"%currentSite,'w')
date = start
while date < end:
    if not (date in spData['t']):
        f.write(date.strftime("%Y-%m-%d")+"\n")
    date += delta

#Deleting empty MERRA2 files and re-downloading them
#date = start
#while date < end:
#    datestr = date.strftime("%Y%m%d")
#    try:
#        merrafile = glob.glob("/n/home10/anwang16/merra2_work/merra2_products/ChajnantorPlateau/merra2_raw_data/*3d*%s*"%datestr)[0]
#        if os.stat(merrafile).st_size < 10000:
#            os.remove(merrafile)
#            m.retrieve_merra2_data_for_date(datestr)
#            time.sleep(10)
#        date += delta
#    except:
#        m.retrieve_merra2_data_for_date(datestr)
#        time.sleep(10)
#        date += delta
    


#Reading APEX PWV data
#apexData = []
#for i in glob.glob("/n/home10/anwang16/merra2_work/merra2_products/ChajnantorPlateau/Apex_pwv/*2010-01*radiometer*"):
#    apexData.append(genfromtxt(i))
#apexData = [x for x in apexData if len(shape(x))==2]
#apexData = concatenate(apexData, axis = 0)
#apexData = [x for x in apexData if x[1] < 10]
#apexData = zip(*apexData)
#apexData[0] = [datetime.datetime.fromtimestamp(x) for x in apexData[0]]

#Reading and Averaging Tipper Data
#datetime, Tsky, tau, Tatm = m.readTipper(start = "20100101",end = "20100131")
#Tsky_m, Tsky_s = m.averageTipper(cpData['t'],datetime,Tsky)

#scatter(spData['t'],spData['pwv'], color = 'blue')
#scatter(cpData['t'], [x/1000 for x in cpData['pwv']], color = 'red')
#plot(cpData['t'], cpData['pwv']
#ylabel("PWV (mm)")
#xlabel("Date")
#title("PWV over Chajnantor Plateau predicted by MERRA2")

#scatter(apexData[0][::60],apexData[1][::60],color = 'blue')
#show()

#print shape(cpData['t'])
#print shape(cpData['tipper850'])
#print shape(Tsky_m)
#print cpData['t']
#print Tsky_m
#scatter(cpData['t'], cpData['tipper850'], color = 'red')
#scatter(cpData['t'], Tsky_m,color = 'blue')
#show()

#scatter(cpData['tipper850'],Tsky_m)
#show()

#print len(apexData[0])
#print len(cpData['t'])
#print apexData[0][::180]
#print len(apexData[0][::180])

#scatter([x/1000 for x in cpData['pwv'][:162]],apexData[1][::180])
#show()
