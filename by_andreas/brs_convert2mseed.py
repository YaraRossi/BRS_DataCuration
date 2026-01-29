#!/usr/bin/env python3
"""
Convert BRS data to MSEED format.

This script reads BRS data files and converts them to MSEED format,
with options for writing to an SDS file system or a single MSEED file.
"""

import numpy as np
import os
import matplotlib.pyplot as plt
import obspy as obs
import argparse
from datetime import datetime
import sys

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Convert BRS data to MSEED format.')
    
    # Input arguments
    parser.add_argument('--input', type=str, required=True,
                        help='Path to input data directory or file')
    parser.add_argument('--event', type=str, default=None,
                        help='Event subdirectory (e.g., "2024-06-28/")')
    
    # Output arguments
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--output_sds', type=str, 
                       help='Output directory for SDS file system')
    group.add_argument('--output_file', type=str,
                       help='Output file path for MSEED file')
    
    # Processing arguments
    parser.add_argument('--sps', type=float, default=None,
                        help='Target sampling rate in Hz for downsampling')

    # Additional options
    parser.add_argument('--plot', action='store_true',
                        help='Plot the data after processing')
    parser.add_argument('--verbose', action='store_true',
                        help='Print detailed information during processing')
    
    args = parser.parse_args()
    
    return args

def read_data(input_path, event=None, verbose=False):
    """
    Read data from input path.
    
    Parameters:
    -----------
    input_path : str
        Path to input data directory or file
    event : str, optional
        Event subdirectory
    verbose : bool
        Print detailed information
        
    Returns:
    --------
    obspy.Stream
        Stream containing the data
    """
    st = obs.Stream()
    
    # Determine the full path
    if event:
        full_path = os.path.join(input_path, event)
    else:
        full_path = input_path
    
    # Check if it's a directory or file
    if os.path.isdir(full_path):
        if verbose:
            print(f"Reading files from directory: {full_path}")
        
        for f in os.listdir(full_path):
            file_path = os.path.join(full_path, f)
            try:
                if verbose:
                    print(f"Reading file: {file_path}")
                st += obs.read(file_path, "MSEED")
            except Exception as e:
                print(f"Failed to read {file_path}")
                print(e)
    else:
        # It's a single file
        if verbose:
            print(f"Reading file: {full_path}")
        try:
            st = obs.read(full_path, "MSEED")
        except Exception as e:
            print(f"Failed to read {full_path}")
            print(e)
            sys.exit(1)
    
    if verbose:
        print(f"Read {len(st)} traces")
    
    return st

def process_stream(st, sps=None, verbose=False):
    """
    Process the stream with optional downsampling and response removal.
    
    Parameters:
    -----------
    st : obspy.Stream
        Stream to process
    sps : float, optional
        Target sampling rate in Hz for downsampling
    remove_response : bool
        Whether to remove instrument response
    response_file : str, optional
        Path to response file
    verbose : bool
        Print detailed information
        
    Returns:
    --------
    obspy.Stream
        Processed stream
    """
    # Make a copy to avoid modifying the original
    st_proc = st.copy()
    
    # Downsample if requested
    if sps is not None:
        if verbose:
            print(f"Downsampling to {sps} Hz...")
        
        for tr in st_proc:
            current_sps = tr.stats.sampling_rate
            if current_sps > sps:
                try:
                    # Calculate decimation factor
                    # factor = int(current_sps / sps)
                    # if factor > 1:
                    #     tr.decimate(factor)
                    #     if verbose:
                    #         print(f"Decimated {tr.id} from {current_sps} Hz to {tr.stats.sampling_rate} Hz")

                    # Use resample instead of decimate
                    tr = tr.filter("lowpass", freq=sps//2, zerophase=True, corners=2)
                    tr.resample(sps, no_filter=True)
                    if verbose:
                        print(f"Resampled {tr.id} from {current_sps} Hz to {tr.stats.sampling_rate} Hz")
                    
                    # Adjust channel code
                    if int(sps) == 100 or int(sps) == 200:
                        tr.stats.channel = "H" + tr.stats.channel[1:]
                    elif int(sps) == 20 or int(sps) == 40:
                        tr.stats.channel = "B" + tr.stats.channel[1:]
                    elif int(sps) == 1 or int(sps) == 5:
                        tr.stats.channel = "L" + tr.stats.channel[1:]
                    else:
                        print(f"Unknown sampling rate: {sps}. Channel code not adjusted.")
                    
                except Exception as e:
                    print(f"Failed to resample {tr.id}: {e}")
    
    return st_proc

def write_to_sds(st, output_dir, verbose=False):
    """
    Write stream to SDS file structure.
    
    Parameters:
    -----------
    st : obspy.Stream
        Stream to write
    output_dir : str
        Output directory for SDS structure
    verbose : bool
        Print detailed information
    """
    if verbose:
        print(f"Writing to SDS structure in {output_dir}...")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Write each trace to SDS structure
    for tr in st:
        # Get trace info
        net = tr.stats.network
        sta = tr.stats.station
        loc = tr.stats.location if tr.stats.location else ""
        cha = tr.stats.channel
        year = tr.stats.starttime.year
        day = tr.stats.starttime.julday
        
        # Create SDS directory structure
        sds_dir = os.path.join(output_dir, str(year), net, sta, cha + "." + loc)
        os.makedirs(sds_dir, exist_ok=True)
        
        # Create filename
        filename = f"{net}.{sta}.{loc}.{cha}.D.{year}.{day:03d}"
        filepath = os.path.join(sds_dir, filename)
        
        # Write trace
        tr_stream = obs.Stream(traces=[tr])
        tr_stream.write(filepath, format="MSEED")
        
        if verbose:
            print(f"Wrote {filepath}")

def main():
    """Main function."""
    # Parse arguments
    args = parse_arguments()
    
    # Read data
    st = read_data(args.input, args.event, args.verbose)
    
    if args.verbose:
        print(st.__str__(extended=True))
    
    # Process data
    st_proc = process_stream(st, args.sps, args.verbose)
    
    if args.verbose:
        print(st_proc.__str__(extended=True))
    
    # Write output
    if args.output_sds:
        write_to_sds(st_proc, args.output_sds, args.verbose)
    else:
        if args.verbose:
            print(f"Writing to MSEED file: {args.output_file}")
        st_proc.write(args.output_file, format="MSEED")
    
    # Plot if requested
    if args.plot:
        st_proc.plot(equal_scale=False)
        plt.show()

    if args.verbose:
        print("Done!")

if __name__ == "__main__":
    main()