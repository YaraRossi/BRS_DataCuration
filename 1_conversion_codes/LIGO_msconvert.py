import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from obspy import Trace, Stream, UTCDateTime
import sys
from datetime import datetime, timedelta

# Constants for conversion
GPS_EPOCH_OFFSET = 315964800  # GPS epoch starts at 1980-01-06T00:00:00Z
GPS_UTC_OFFSET = 18  # Leap second difference as of 2024

# Convert GPS time to UTC (Unix epoch)
def gps_to_utc(gps_timestamp):
    unix_timestamp = gps_timestamp + (GPS_EPOCH_OFFSET - GPS_UTC_OFFSET)
    return datetime.utcfromtimestamp(unix_timestamp)  # Convert to UTC datetime

# Ensure correct usage
if len(sys.argv) != 2:
    print("Usage: python script.py <directory_path>")
    sys.exit(1)

directory = sys.argv[1]

# Check if the directory exists
if not os.path.isdir(directory):
    print(f"Error: The provided path '{directory}' is not a valid directory.")
    sys.exit(1)

# Process each CSV file in the directory
for filename in os.listdir(directory):
    if filename.endswith(".csv"):
        file_path = os.path.join(directory, filename)

        # Extract metadata from filename
        parts = filename.rstrip(".csv").split('_')
        if len(parts) != 4:
            print(f"Skipping file {filename}: Incorrect filename format.")
            continue

        NET, STA, LOC, CHAN = parts  # Assign metadata

        # Read CSV file
        df = pd.read_csv(file_path, header=None, names=["gps_timestamp", "amplitude"])

        # Ensure there are at least two rows to compute sampling rate
        if len(df) < 2:
            print(f"Skipping file {filename}: Not enough data points.")
            continue

        # Convert first GPS timestamp to UTC for MiniSEED metadata
        start_utc = gps_to_utc(df["gps_timestamp"].iloc[0])
        starttime = UTCDateTime(start_utc)

        # Compute sampling rate (1 / time difference between first two timestamps)
        sampling_rate = 1 / (df["gps_timestamp"].iloc[1] - df["gps_timestamp"].iloc[0])

        # Convert GPS timestamps to human-readable format for plotting
        df["gps_datetime"] = df["gps_timestamp"].apply(gps_to_utc)

        # Print header variables
        print(f"Processing: {filename}")
        print(f"Network: {NET}")
        print(f"Station: {STA}")
        print(f"Location: {LOC}")
        print(f"Channel: {CHAN}")
        print(f"Start Time (UTC): {starttime}")
        print(f"Sampling Rate: {sampling_rate:.6f} Hz")
        print("-" * 40)

        # Create an ObsPy Trace and Stream
        trace = Trace()
        trace.data = df["amplitude"].values.astype(np.float32)  # Convert data to float32
        trace.stats.network = NET
        trace.stats.station = STA
        trace.stats.location = LOC
        trace.stats.channel = CHAN
        trace.stats.sampling_rate = sampling_rate
        trace.stats.starttime = starttime

        stream = Stream(traces=[trace])

        # Write to MiniSEED file
        mseed_filename = file_path.replace(".csv", ".ms")
        stream.write(mseed_filename, format="MSEED")
        print(f"Saved MiniSEED: {mseed_filename}")

        # Plot the time series using GPS time
#       plt.figure(figsize=(10, 5))
        plt.figure(figsize=(15, 7.5))
#       plt.plot(df["gps_datetime"], df["amplitude"], label="Ground Motion", color='black', size=(2000,1200), linewidth=0.2, number_of_ticks=10)
        plt.plot(df["gps_datetime"], df["amplitude"], label="Ground Motion", color='red',linewidth='0.2')
        plt.xlabel("GPS Time (YYYY-MM-DD HH:MM:SS)")
        plt.ylabel("Amplitude")
#       plt.title(f"Seismogram: {NET}.{STA}.{LOC}.{CHAN} (GPS Time)")
        plt.title(f"Seismogram: {NET}.{STA}.{LOC}.{CHAN} {starttime}" )
        plt.xticks(rotation=45)
        plt.legend()
        plt.grid()
        plt.show()

print("Processing complete.")
