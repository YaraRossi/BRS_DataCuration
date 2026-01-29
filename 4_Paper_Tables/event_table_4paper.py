
import pandas as pd
from obspy.geodetics.base import gps2dist_azimuth
import os


path = '' #'/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/' # path to event_times_LIV.txt file
root_data = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/miniseed_LIGO'
# the columns in this file are:
# Eventtime,   tstart,       tend,        eventstart,            eventend,          location, latitude, longitude, mag, name, country, type
# which is:
# origintime, data start, data end, event start at station, event end at station, etc...
liv_file = f'{path}event_times_LIV.txt'
han_file = f'{path}event_times_HAN.txt'
data_liv = pd.read_csv(liv_file)
data_han = pd.read_csv(han_file)

data = pd.concat([data_liv, data_han[~data_han['location'].isin(data_liv['location'])]], ignore_index=True)


# add note for stations that are not good quality:
list_of_errors = {
    '2024-04-02': r'\makecell{LIV.IX.DJY\\abnormal resonance}',
    '2024-12-05': r'\makecell{clipping of\\seismometers}',
    '2025-02-08': r'\makecell{clipping of\\seismometers}',
    '2024-05-28': r'\makecell{HAN.VA.DHX\\high noise level}',
    '2024-09-23': r'\makecell{HAN.VA.DHX\\high noise level}',
    '2025-01-20': r'\makecell{HAN.VA.DHX\\high noise level}',
    '2024-07-31': r'\makecell{HAN.VA.DHX\\high noise level}',
    '2023-08-20': r'\makecell{HAN.VA.DHX\\high noise level}',
    '2024-12-17': r'\makecell{HAN.VA.DHX\\high noise level}',
    '2021-03-04': r'\makecell{LIV.EY.DJX\\unusable data}',
    '2020-06-23': r'\makecell{LIV.EY.DJX\\unusable data}',
    '2020-10-19': r'\makecell{LIV.EY.DJX\\unusable data}',
    '2026-01-16': r'\makecell{all LIV BRS\\unusable data}',
    '2022-09-19': r'\makecell{all LIV BRS\\unusable data}',
    '2025-12-08': r'only 1 sps data'
}
# LIV.EY.DJX for some events, i.e. 2021-03-04, 2020-06-23, 2020-10-19

# list of events for Andreas
events_of_interest = ['2024-04-02', '2024-06-28','2024-12-17','2025-07-29','2024-09-15',]

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

# which data is available:
availability = []
all_folders = [f for f in os.listdir(root_data) if os.path.isdir(os.path.join(root_data, f))]
for event_name in data['location']:
    # look for folder with this in its name. if name has _L then append (Han, LIV)
    # if han in its name append HAN and if liv in its name append LIV.
    # if there are two folders then also append HAN, LIV
    matching_folders = [f for f in all_folders if event_name in f]
    if len(matching_folders) > 1:
        has_B = True
        has_han = False
        has_liv = False
    else:
        has_han = any('han' in f for f in matching_folders)
        has_liv = any('liv' in f for f in matching_folders)
        has_B = False


    if has_B:
        availability.append('HAN, LIV')
    elif has_han:
        availability.append('HAN')
    elif has_liv:
        availability.append('LIV')
    else:
        availability.append('HAN, LIV')  # no data found


table['availability'] = availability
# add error column (empty if no issue)
table['errors'] = data['location'].map(list_of_errors).fillna('')


table.drop(table[table['type'].str.contains('tornado', na=False)].index, inplace=True)
table.drop(table[table['type'].str.contains('eruption', na=False)].index, inplace=True)
print(table)

# sort table according to magnitude
table = table.sort_values(by=' mag', ascending=False)

output_path = f'{path}event_table_latex.txt'
table.to_csv(output_path, index=False)
output_path_local = f'event_table_latex.txt'
table.to_csv(output_path_local, index=False)

with open(f'{path}event_table_latex.txt', 'w') as f:
    for _, row in table.iterrows():
        row_values = [str(val) for val in row]
        line = ' & '.join(row_values) + ' \\\\'
        f.write(line + '\n')
print(table)