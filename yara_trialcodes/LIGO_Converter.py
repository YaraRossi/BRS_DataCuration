import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
from obspy import Trace, Stream, UTCDateTime

# Directory containing the .csv files
#directory = "./"


# Check if the script received an argument
if len(sys.argv) != 2:
    print("Usage: python process_files.py <directory_path>")
    sys.exit(1)

# Get the directory path from the command-line arguments
directory = sys.argv[1]

# Check if the provided path is a valid directory
if not os.path.isdir(directory):
    print(f"The provided path '{directory}' is not a valid directory.")
    sys.exit(1)

# Process files in the directory
#for filename in os.listdir(directory):
#    if filename.endswith(".csv"):
#

# Loop through all files in the directory
for filename in os.listdir(directory):
    # Only process files ending in .csv
    if filename.endswith(".csv"):
        # Split the filename into parts
        parts = filename.rstrip(".csv").split('_')
        
            # Ensure the filename has the expected number of parts
        if len(parts) >= 8:
            # Assign the parts to variables
            site = parts[0]          # e.g., HAN
            loc = parts[1]           # e.g., VA
            chan = parts[2]          # e.g., DH2
            
            # Parse the datetime components
            year = parts[5]          # e.g., 2024
            month = parts[3]         # e.g., Sep
            day = parts[4]           # e.g., 23
            time = parts[6]          # e.g., 03-17-33-000
            
            # Convert month abbreviation to a zero-padded numeric month
            from datetime import datetime
            month_num = datetime.strptime(month, "%b").month  # Convert 'Sep' to 09
            
            # Split the time into components and reformat
            time_components = time.split('-')
            hour, minute, second, microsecond = time_components[0], time_components[1], time_components[2], time_components[3]
            
            # Build the formatted date_string
            date_string = f"{year}/{month_num:02}/{day} {hour}:{minute}:{second}.{microsecond}"
            
            samprate = parts[7]      # e.g., 256
            
            net = "LI"               #Unless from source other than LIGO!!


        else:
            print(f"Filename {filename} does not have the expected format.")
        df = pd.read_csv(directory + filename, usecols=[0, 1])  
        df.columns = ['ordinate', 'Rot']
        trace_out = Trace(data=df['Rot'].values, header={'starttime': date_string, \
                                             'sampling_rate': samprate, 'network': net, \
                                             'station': site, 'channel': chan, 'location':loc})
        # Create an ObsPy Stream object and add the Trace objects
        #stream = Stream(traces=[trace_x, trace_y, trace_temp])
        stream = Stream(traces=[trace_out])

# Plot the data in the stream object

# this Kludge seems to make the script pause after drawing a plot of a channel until
# the user dismisses it ...
        x=stream.plot(color='black', size=(2000,1200), linewidth=0.2, number_of_ticks=10)
        
        # Print the values to the screen
        #print(f"filename")
        print(f"{filename}")
        print(f"Site: {site}")
        print(f"Location: {loc}")
        print(f"Channel: {chan}")
        print(f"Date String: {date_string}")
        print(f"Sample Rate: {samprate}")
        print("-" * 40)  # Separator for clarity
        
        print(stream[0].stats)
        
# output miniseed files
        stream.write(directory + filename + '.ms',format='MSEED')
