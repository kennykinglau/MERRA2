#Frequency Step Testimport merra2Player as m2pimport timeimport osfrom pylab import *#for arbitrary purposes, we choose 2010-07-01 00:00:00 as our datetime to run am on repeatedly. #dfList = [100, 200, 500, 1000, 2000, 5000, 8000, 10000, 20000, 50000, 100000, 200000, 500000]
dfList = [100, 1000, 2000, 10000]data = []m = m2p.merra2Player()
m.defineSite(sitestr = "SouthPole")m.debug = True
date = '20150601'amcFile = 'SouthPole_%s_000000_gndData1.amc'%date
outFile = 'SouthPole_%s_000000_gndData1.out'%date
errFile = 'SouthPole_%s_000000_gndData1.err'%datef850,band850f,band850rj = m.readBandpass(expt = 'tipper', band = '850')
f220,band220f,band220rj = m.readBandpass(band='220')
band220rj = band220f  # hack because rj band is missing
f150,band150f,band150rj = m.readBandpass(band='150')
f100,band100f,band100rj = m.readBandpass(band='100')
for step in dfList:    os.remove(m.amcDir+outFile)    #print(os.listdir(m.amcDir))    os.remove(m.amcDir+errFile)    start = time.time()    a, b, c, d = m.run_am(amcFile, f0 = 0.0, f1 = 1200.0, df = step)
    #fs, tb, trj, pwv
    plot(a,b, color = 'red')
    plot(a,c, color = 'blue')
    title("Brightness temperature spectrum from am: %s MHz step size"%step)
    show()
    runtime = time.time() - start    print runtime    tsky850 = m.integBand(f850,band850rj,b,c)
    tsky220 = m.integBand(f220,band220rj,b,c)
    tsky150 = m.integBand(f150,band150rj,b,c)
    tsky100 = m.integBand(f100,band100rj,b,c)    data += [[ step, runtime, tsky850, tsky220, tsky150, tsky100, d ]]transposed = zip(*data)print(transposed)print()print()print(transposed[0])print(transposed[1])print(transposed[2])print(transposed[3])plot(transposed[0], transposed[1])show()

scatter(transposed[0],100*(transposed[3]-transposed[3][0])/transposed[3][0])axhline(y=0, xmin=0, xmax=transposed[0][-1], hold=None)ylim([-20,20])xscale('log')ylabel('percent deviation from Tsky at 100MHz step size')xlabel('frequency step size (MHz)')title('percent divergence of Tsky at South Pole, integrated over 220 GHz band: 2015-06-01 00:00:00')show()

scatter(transposed[0],100*(transposed[4]-transposed[4][0])/transposed[4][0])axhline(y=0, xmin=0, xmax=transposed[0][-1], hold=None)ylim([-20,20])xscale('log')ylabel('percent deviation from Tsky at 100MHz step size')xlabel('frequency step size (MHz)')title('percent divergence of Tsky at South Pole, 150 GHz band: 2015-06-01 00:00:00')show()

scatter(transposed[0],100*(transposed[5]-transposed[5][0])/transposed[5][0])axhline(y=0, xmin=0, xmax=transposed[0][-1], hold=None)ylim([-20,20])xscale('log')ylabel('percent deviation from Tsky at 100MHz step size')xlabel('frequency step size (MHz)')title('percent divergence of Tsky at South Pole, 100 GHz band: 2015-06-01 00:00:00')show()

scatter(transposed[0], transposed[1])plot(transposed[0], transposed[1])xscale('log')yscale('log')ylabel('time (s)')xlabel('frequency step size (MHz)')title('time to run_am')show()


