import numpy
import obspy
import pytz
from datetime import datetime, timedelta


def data2miniseedGeneral(files='', newfolder='', newname='', delimiter='', verbose=True, precision=9,
                         save=False, network='9E', station='ZPT', location='RN'):
    file = numpy.loadtxt(open(files, 'rt'), delimiter=delimiter)
    stats = obspy.core.trace.Stats()
    stats.network = network
    stats.station = station
    stats.location = location

    instrument_codes = 'MMM'  # disp uncertainty
    orientation_codes = 'ENZ'

    # makes empty stream to append to
    newstream = obspy.core.stream.Stream()

    # find record sections with constant sampling rate
    delta = round(numpy.median(numpy.diff(file[:, 1])), precision)
    sections = [-1]
    for i in range(file.shape[0] - 1):
        test = abs(file[i + 1, 1] - file[i, 1])
        if test > delta * 1.2 or test < delta * .8:
            sections.append(i)
    sections.append(i)

    for i in range(len(sections) - 1):
        begin = sections[i] + 1
        end = sections[i + 1]
        stats.npts = len(file[begin:end, 1])
        if stats.npts > 1:
            mjd = file[begin, 0]
            sec_of_day = file[begin, 1] - 18  # leap seconds
            utc_tz = pytz.timezone('UTC')
            unix_epoch = datetime(1970, 1, 1, tzinfo=utc_tz)
            delta_t = timedelta(days=mjd - 40587, seconds=sec_of_day)
            unix_time = (unix_epoch + delta_t).timestamp()

            stats.starttime = obspy.UTCDateTime(unix_time)
            stats.delta = round(numpy.median(numpy.diff(file[begin:end, 1])), precision)

            for columnindex, data in enumerate(file.transpose()[2:5]):
                trace = obspy.core.trace.Trace(data=data[begin:end],
                                               header=stats)
                trace.stats.channel = 'H' + instrument_codes[columnindex] + orientation_codes[columnindex]
                if verbose:
                    print(trace)
                newstream += trace

                # newstream._cleanup()
    newstream.merge(fill_value='interpolate')  # , interpolation_samples=)
    if save:
        newstream.write('%s/%s.mseed' % (newfolder, newname), format='MSEED')
    return newstream