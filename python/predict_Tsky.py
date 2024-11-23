#! /usr/bin/env python

import argparse
import datetime
import time
import os

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-l",
        dest="location",
        default='None',
        type=str,
        help="-l. Comma-separated string to define site. Can be either -l [pre-selected site], -l NewSite,shortname,lat,lon,alt. Pre-selected sites are: SouthPole, ChajnantorPlateau, CerroChajnantor, MaunaKea, Summit, Qubic. Default: None",
    )

    parser.add_argument(
        "-b",
        dest="bandpass",
        default="None",
        type=str,
        help="-b. Bandpass option. If you have a custom passband that can be modeled as a tophat, please input -b bandname, band_center(GHz), bandwidth(GHz). Default: None",
    )

    parser.add_argument(
        "-d",
        dest="dateopt",
        default="20150123",
        type=str,
        help="-d. Can be either: 'start, end' or 'start', or 'year'.  Start and end can be either YYYYMMDD or \"YYYYMMDDTHH:MM:SS\". Year: YYYY. Default: 20160123",
    )
    # TODO: Detect when the daterange goes out of the range of MERRA2, and switch to using GFS.
    parser.add_argument(
        "-p",
        dest="plotFig",
        action="store_true",
        default=False,
        help="-p option will plot a figure, default: False",
    )

    parser.add_argument(
        "-g",
        dest="groundData",
        default=2,
        type=int,
        help="-g. Type of surface level data. 1: historical measured data when available. 2(default): MERRA2 surface data. 3: Extrapolation from MERRA2 3d data.",
    )

    parser.add_argument(
        "-c",
        dest="cloud",
        action="store_true",
        default=False,
        help="-c. adds ice and liquid water layers in the am profile. Default (false) is for am profile to only contain water-vapor data.",
    )

    parser.add_argument(
        "-s",
        dest="saveResults",
        action="store_false",
        default=True,
        help="-s. Default=True, Forces the output of a run to be saved in a csv and pickle file",
    )

    options = parser.parse_args()

    start_datetime = datetime.datetime.now()
    start_time_process = time.process_time_ns()
    start_time_total = time.perf_counter_ns()
    import merra2Player as m2p

    m = m2p.merra2Player()

    # Option parsing
    bandopt = {}
    dateopt = {}
    siteopt = {}

    bd = options.bandpass.split(',')
    if len(bd) == 1:
        bandopt['name'] = bd[0]
    elif len(bd) == 3:
        bandopt['name'] = 'custom'
        bandopt['v_cen'] = float(bd[1])
        bandopt['bw'] = float(bd[2])
    else:
        print("ERROR in bandpass selection")
        exit()

    st = options.location.split(',')
    if len(st) == 1:
        siteopt['type'] = 'presel'
        siteopt['name'] = st[0]
    elif len(st) == 5:
        siteopt['type'] = 'custom'
        siteopt['name'] = st[0]
        siteopt['sname'] = st[1]
        siteopt['lat'] = float(st[2])
        siteopt['lon'] = float(st[3])
        siteopt['alt'] = float(st[4])
    else:
        print("ERROR in site selection")
        exit()

    if options.cloud == False:
        siteopt['cldtype'] = 'vapor'
    else:
        siteopt['cldtype'] = 'icewater'
    siteopt['groundData'] = options.groundData

    dt = options.dateopt.split(',')
    if len(dt) == 1:
        if len(dt[0]) == 4:
            dateopt['type'] = 'datelist'
            dateopt['start'] = dt[0] + '0101'
            dateopt['end'] = dt[0] + '1231'
            dateopt['single'] = False
        else:
            dateopt['type'] = 'datelist'
            dateopt['start'] = dt[0]
            dateopt['end'] = dt[0]
            dateopt['single'] = True
    elif len(dt) == 2:
        dateopt['type'] = 'datelist'
        dateopt['start'] = dt[0]
        dateopt['end'] = dt[1]
        dateopt['single'] = False
    else:
        print("ERROR in dateopt selection")
        exit()

    # site has to be defined first as self.site is used throughout the class
    sitopt = m.defineSite(site=siteopt)
    print(dateopt, siteopt, bandopt)
    # Links and directory names have to be defined first as they are used throughout
    m.checkLinks()
    pre_run_time_process = time.process_time_ns()
    pre_run_time_total = time.perf_counter_ns()
    Tsky_merra = m.runMERRA(dateopt=dateopt, bandopt=bandopt)
    if siteopt['type'] == 'presel':
        Tsky_merra = m.runTipper(Tsky_merra, dateopt)
    post_run_time_process = time.process_time_ns()
    post_run_time_total = time.perf_counter_ns()

    opts = {}
    opts['web'] = False
    opts['now'] = start_datetime
    opts['bandopt'] = bandopt
    opts['siteopt'] = siteopt
    opts['dateopt'] = dateopt

    if dateopt['single']:
        m.printSingle(dateopt['start'], Tsky_merra)
    else:
        if options.saveResults:
            csvFile = m.save2csv(Tsky_merra, opts)
            pickleFile = m.save2pickle(Tsky_merra, opts)
        else:
            print("Results not saved to any file")

    if options.plotFig:
        tskyPlot = m.save2plot(Tsky_merra, opts)
        os.system('eog %s' % tskyPlot)

    end_time_total = time.perf_counter_ns()
    end_time_process = time.process_time_ns()
    end_datetime = datetime.datetime.now()
    print(
        "Total time seconds (datetime, not very precise): %.2f" % (end_datetime - start_datetime).total_seconds()
    )
    print(
        "Total time seconds (perf_counter, precise): %.5f\n"
        "  Initialization: %.5f\n"
        "  Running (generating amcFiles): %.5f\n"
        "  Saving (generating final): %.5f" % (
            (end_time_total - start_time_total) * 1e-9,
            (pre_run_time_total - start_time_total) * 1e-9,
            (post_run_time_total - pre_run_time_total) * 1e-9,
            (end_time_total - post_run_time_total) * 1e-9
        )
    )
    print(
        "Process time seconds: %.5f\n"
        "  Initialization: %.5f\n"
        "  Running (generating amcFiles): %.5f\n"
        "  Saving (generating final): %.5f" % (
            (end_time_process - start_time_process) * 1e-9,
            (pre_run_time_process - start_time_process) * 1e-9,
            (post_run_time_process - pre_run_time_process) * 1e-9,
            (end_time_process - post_run_time_process) * 1e-9
        )
    )


    # TODO add unit tests
    # predict_Tsky.py -l SouthPole -d 20170801
    # predict_Tsky.py -l SouthPole -d 20170801, 20170804
    # predict_Tsky.py -l SouthPole -d 2017
    # predict_Tsky.py -l SouthPole -d 20170501 04:05:06

    # predict_Tsky.py -d 20170801 04:05:06 -l "SMT, smt, 32.7016, -109.891, 3185"
