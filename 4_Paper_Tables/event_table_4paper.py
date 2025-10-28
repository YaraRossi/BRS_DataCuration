
import pandas as pd
from obspy.geodetics.base import gps2dist_azimuth


path = '' # path to event_times_LIV.txt file
# the columns in this file are:
# Eventtime,   tstart,       tend,        eventstart,            eventend,          location, latitude, longitude, mag, name, country, type
# which is:
# origintime, data start, data end, event start at station, event end at station, etc...
events = f'{path}event_times_LIV.txt'
data = pd.read_csv(events)

## get distances for events

dist_han, dist_liv = [],[]
for _, data_event in data.iterrows():
    source_latitude = float(data_event[' latitude'])
    source_longitude = float(data_event[' longitude'])
    # HAN
    station_latitude = 46.455#47.5
    station_longitude = -119.407 #-120.7
    dist_han.append(int(gps2dist_azimuth(source_latitude, source_longitude, station_latitude, station_longitude)[0]/1000))
    # LIV
    station_latitude = 30.562775410180993
    station_longitude = -90.77469251757329
    dist_liv.append(int(gps2dist_azimuth(source_latitude, source_longitude, station_latitude, station_longitude)[0]/1000))


table = data[[' mag']].copy()
#table['name'] = '\\makecell{' + data[' name'].astype(str) + ',\\\\' + data[' country'].astype(str) + '}'
table['name'] = data[' name'].astype(str) + ',' + data[' country'].astype(str)
order = [' latitude', ' longitude', 'Eventtime']
table[order] = data[order]
table['dist_han'] = dist_han
table['dist_liv'] = dist_liv
table['type'] = data[' type']

table.drop(table[table['type'].str.contains('tornado')].index, inplace=True)
table.drop(table[table['type'].str.contains('eruption')].index, inplace=True)
print(table)

# sort table according to magnitude
table = table.sort_values(by=' mag', ascending=False)

output_path = f'{path}event_table_latex.txt'
table.to_csv(output_path, index=False)

with open(f'{path}event_table_latex.txt', 'w') as f:
    for _, row in table.iterrows():
        row_values = [str(val) for val in row]
        line = ' & '.join(row_values) + ' \\\\'
        f.write(line + '\n')
print(table)