#! /usr/bin/env python

import argparse
import datetime
import os
import merra2Player as m2p
import netCDF4


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--location",
        default='SouthPole',
        help="Sites are: SouthPole, ChajnantorPlateau, CerroChajnantor, MaunaKea, Summit, Qubic. Default: SouthPole",
    )

    parser.add_argument(
        "--date",
        default=["20190101", "20201130"],
        nargs=2,
        help="Date range in the format YYYYMMDD, YYYYMMDD",
    )

    # parser.add_argument(
    #    "-g",
    #    dest="groundData",
    #    default=2,
    #    type=int,
    #    help="-g. Type of surface level data. 1: historical measured data when available. 2(default): MERRA2 surface data. 3: Extrapolation from MERRA2 3d data.",
    # )

    # parser.add_argument(
    #    "-c",
    #    dest="cloud",
    #    action="store_true",
    #    default=False,
    #    help="-c. adds ice and liquid water layers in the am profile. Default (false) is for am profile to only contain water-vapor data.",
    # )

    # parser.add_argument(
    #    "-s",
    #    dest="saveResults",
    #    action="store_false",
    #    default=True,
    #    help="-s. Default=True, Forces the output of a run to be saved in a csv and pickle file",
    # )

    args = parser.parse_args()

    import merra2Player as m2p

    m = m2p.merra2Player()

    # This should load the cloud data (?)
    siteopt = {}
    siteopt['type'] = 'presel'
    siteopt['name'] = args.location
    siteopt['cldtype'] = 'icewater'
    siteopt['groundData'] = 2
    siteopt = m.defineSite(site=siteopt)

    # Simplify the datelist interface
    dateopt = {}
    dateopt['type'] = 'datelist'
    dateopt['start'] = args.date[0]
    dateopt['end'] = args.date[1]
    dateopt['single'] = False

    # Links and directory names have to be defined first as they are used throughout
    print(dateopt, siteopt)
    m.checkLinks()

    # Now let's load the MERRA2 data for the correct date range
    datelist = m.retrieve_merra2_data_for_dateRange(dateopt['start'], dateopt['end'])

    exit()

    merradata = {}

    for date in datelist:

        # Let's load the netDCF4 filenames
        url, merraFilename_m = m.get_url_for_date(date.strftime('%Y%m%d'), 'multiLevel')
        url, merraFilename_s = m.get_url_for_date(
            date.strftime('%Y%m%d'), 'singleLevel'
        )
        del url

        print(merraFilename_m)
        print(merraFilename_s)

        data_s = netCDF4.Dataset(os.path.join(m.merraDir, merraFilename_s))
        data_m = netCDF4.Dataset(os.path.join(m.merraDir, merraFilename_m))

        print(data_s)
        print(data_m)

        print(data_s.variables.keys())

        # This seems to be the right number to pull from the MERRA2 data
        print(data_s.variables['TQI'][:, 1, 1])

        exit()

    # Tsky_merra = m.runMERRA(dateopt=dateopt, bandopt=bandopt)
    exit()

    # if siteopt['type'] == 'presel':
    #    Tsky_merra = m.runTipper(Tsky_merra, dateopt)

    # opts = {}
    # opts['web'] = False
    # opts['now'] = now
    # opts['bandopt'] = bandopt
    # opts['siteopt'] = siteopt
    # opts['dateopt'] = dateopt

    # if dateopt['single']:
    #    m.printSingle(dateopt['start'], Tsky_merra)
    # else:
    #    if options.saveResults:
    #        csvFile = m.save2csv(Tsky_merra, opts)
    #        pickleFile = m.save2pickle(Tsky_merra, opts)
    #    else:
    #        print("Results not saved to any file")

    # if options.plotFig:
    #    tskyPlot = m.save2plot(Tsky_merra, opts)
    #    os.system('eog %s' % tskyPlot)

    # print(
    #    "Total time: %.2f seconds" % ((datetime.datetime.now() - now).total_seconds())
    # )

    # TODO add unit tests
    # predict_Tsky.py -l SouthPole -d 20170801
    # predict_Tsky.py -l SouthPole -d 20170801, 20170804
    # predict_Tsky.py -l SouthPole -d 2017
    # predict_Tsky.py -l SouthPole -d 20170501 04:05:06

    # predict_Tsky.py -d 20170801 04:05:06 -l "SMT, smt, 32.7016, -109.891, 3185"
