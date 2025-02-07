import obspy
from obspy import UTCDateTime, read
from obspy.clients.fdsn import Client
import pandas
import numpy
from numpy import sin, cos, pi
from obspy.clients.fdsn.header import FDSNException

'''

# Get subset of stations
client = Client('IRIS') # this can be NCEDC, SCEDC, EIDA, IRISDMC etc...
# website to iris map of stations and network
# https://ds.iris.edu/gmap/#network=PY&planet=earth

tstart=UTCDateTime('2024-04-02T23:58:12.0')
tend=UTCDateTime('2024-04-03T02:58:12.0')

datetime = '%s.%s.%s' % (tstart.year, tstart.month, tstart.day)

inventory=client.get_stations(network ='PY', level='channel',
				starttime=tstart, endtime=tend) # you can many more parameters here (e.g. a radius etc…)
inv=inventory.select(station="PFOIX", location = '*', channel='*')
print(inv)

'''


def get_PFO_data(root='folderpath', start = '2021-08-23T00:00:00', end = '2021-08-24T00:00:00', blueSeis=True, strongmotion=True, weakmotion=True):
    """
    Fetches seismic waveform data from IRIS and saves it in MiniSEED format for the network PY. PY is the Piñon flat observatory in central California.

    Parameters:
    root (str): Directory to save the waveform files.
    start (str): Start time in ISO format ('YYYY-MM-DDTHH:MM:SS').
    end (str): End time in ISO format ('YYYY-MM-DDTHH:MM:SS').
    blueSeis (bool): If True, fetches BlueSeis data (HJ* channels).
    strongmotion (bool): If True, fetches strong motion data (HH* channels).
    weakmotion (bool): If True, fetches weak motion data (BH* channels).
    """
    # Create a client instance for IRIS data center
    client = Client('IRIS')

    # Fetch blueSeis-3A data from BSPF station if enabled
    if blueSeis:
        try:
            waveforms = client.get_waveforms(network = 'PY', station = 'BSPF', location = '', channel='HJ*', starttime=start, endtime=end, attach_response=True)
            print(waveforms)# Print retrieved waveforms

            for trace in waveforms:
                starttime=trace.stats.starttime.timetuple()
                trace.write("%s/%s.D.%4d.%03d" %(root, trace.id,starttime.tm_year,starttime.tm_yday),format='MSEED',reclen=512)
        except:
            print('failed to download BSPF waveforms for time: %s'  %start)

    # Fetch strongmotion data from PFOIX station if enabled
    if strongmotion:
        try:
            waveforms = client.get_waveforms(network = 'PY', station = 'PFOIX', location = '', channel='HH*', starttime=start, endtime=end, attach_response=True)
            print(waveforms)# Print retrieved waveforms

            for trace in waveforms:
                starttime=trace.stats.starttime.timetuple()
                trace.write("%s/%s.D.%4d.%03d" %(root, trace.id,starttime.tm_year,starttime.tm_yday),format='MSEED',reclen=512)
        except:
            print('failed to download PFOIX waveforms for time: %s'  %start)

    # Fetch weakmotion data from PFOIX station if enabled
    if weakmotion:
        try:
            waveforms = client.get_waveforms(network = 'PY', station = 'PFOIX', location = '', channel='BH*', starttime=start, endtime=end, attach_response=True)
            print(waveforms)# Print retrieved waveforms

            for trace in waveforms:
                starttime=trace.stats.starttime.timetuple()
                trace.write("%s/%s.D.%4d.%03d" %(root, trace.id,starttime.tm_year,starttime.tm_yday),format='MSEED',reclen=512)
        except:
            print('failed to download PFOIX waveforms for time: %s'  %start)


root_data = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/PFO_data'
root_event = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data'


df_events = pandas.read_csv("%s/event_times.txt" % root_event, parse_dates=['eventtime', 'tstart', 'tend', 'location'])

for tstart, tend, location in zip(df_events.tstart, df_events.tend, df_events.location):
    print(tstart, tend, location)

    if location in ['Atiquipa Peru2']:
        # there is a data gap for the blueSeis-3A during the Atiquipa EQ.
        # "Atiquipa Peru2" allows to download the full day, making the datagap visible.
        continue

    try:
        get_PFO_data(root=root_data, start = UTCDateTime(tstart), end = UTCDateTime(tend), blueSeis=False, strongmotion=True, weakmotion=False)
    except:
        print('failed to download data for location %s' % location)