#!/usr/bin/env python3
"""
Prepare the Data for Analysis

This script processes waveform data for all dates in a given folder.
It loads rotation and translation data, rotates to ZNE, and outputs miniseed files.
"""

import os
import sys
import argparse
import re
import obspy as obs
from pathlib import Path
from datetime import datetime

from sixdegrees.sixdegrees import sixdegrees
from obspy.signal.rotate import rotate2zne


def extract_date_from_folder(folder_name):
    """
    Extract date from folder name using various date patterns.
    
    Parameters:
    -----------
    folder_name : str
        Folder name that may contain a date
        
    Returns:
    --------
    str or None
        Date string in YYYY-MM-DD format, or None if no date found
    """
    # Common date patterns to try
    patterns = [
        r'(\d{4}-\d{2}-\d{2})',           # YYYY-MM-DD
        r'(\d{4}_\d{2}_\d{2})',           # YYYY_MM_DD
        r'(\d{4}\.\d{2}\.\d{2})',         # YYYY.MM.DD
        r'(\d{8})',                        # YYYYMMDD
        r'(\d{4}/\d{2}/\d{2})',           # YYYY/MM/DD
    ]
    
    for pattern in patterns:
        match = re.search(pattern, folder_name)
        if match:
            date_str = match.group(1)
            # Normalize to YYYY-MM-DD format
            if len(date_str) == 8 and date_str.isdigit():
                # YYYYMMDD format
                date_str = f"{date_str[:4]}-{date_str[4:6]}-{date_str[6:8]}"
            elif '_' in date_str:
                # YYYY_MM_DD format
                date_str = date_str.replace('_', '-')
            elif '.' in date_str:
                # YYYY.MM.DD format
                date_str = date_str.replace('.', '-')
            elif '/' in date_str:
                # YYYY/MM/DD format
                date_str = date_str.replace('/', '-')
            
            # Validate the date
            try:
                datetime.strptime(date_str, '%Y-%m-%d')
                return date_str
            except ValueError:
                continue
    
    return None


def extract_suffix_from_folder(folder_name, date):
    """
    Return an explicit event suffix from the folder name, if present.

    Matches patterns like 2023-02-06_01, 20230206_02, or trailing _01 on the name.
    """
    date_compact = date.replace("-", "")
    patterns = [
        rf"{re.escape(date)}_(\d{{2}})(?:_|$)",
        rf"{date_compact}_(\d{{2}})(?:_|$)",
        rf"_(\d{{2}})$",
    ]
    for pattern in patterns:
        match = re.search(pattern, folder_name)
        if match:
            return f"_{match.group(1)}"
    return ""


def assign_output_suffixes(date_folders):
    """
    Assign output suffixes for event directories.

    One folder per date -> no suffix (event_20240628).
    Multiple folders with the same date -> _01, _02, ... in sorted folder order,
    unless the folder name already contains an explicit _NN suffix.
    """
    from collections import defaultdict

    by_date = defaultdict(list)
    for folder_name, date in date_folders:
        by_date[date].append(folder_name)

    assigned = []
    for date in sorted(by_date):
        folders = sorted(by_date[date])
        if len(folders) == 1:
            suffix = extract_suffix_from_folder(folders[0], date)
            assigned.append((folders[0], date, suffix))
            continue

        used = set()
        auto_index = 1
        for folder in folders:
            suffix = extract_suffix_from_folder(folder, date)
            if suffix and suffix not in used:
                used.add(suffix)
                assigned.append((folder, date, suffix))
                continue
            while f"_{auto_index:02d}" in used:
                auto_index += 1
            suffix = f"_{auto_index:02d}"
            used.add(suffix)
            assigned.append((folder, date, suffix))
            auto_index += 1

    return assigned


def event_output_dir(output_base, date, suffix):
    """Build output directory name, e.g. event_20230206_01."""
    return os.path.join(output_base, f"event_{date.replace('-', '')}{suffix}")


def scan_data_path(path_to_data_base):
    """
    Scan the data path and extract dates from folder names.
    
    Parameters:
    -----------
    path_to_data_base : str
        Base path containing date folders
        
    Returns:
    --------
    list of tuples
        List of (folder_name, date_string) tuples
    """
    date_folders = []
    
    if not os.path.exists(path_to_data_base):
        return date_folders
    
    print(f"Scanning for date folders in: {path_to_data_base}")
    
    for item in os.listdir(path_to_data_base):
        item_path = os.path.join(path_to_data_base, item)
        if os.path.isdir(item_path):
            date = extract_date_from_folder(item)
            if date:
                date_folders.append((item, date))
                print(f"  Found: {item} -> {date}")
            else:
                print(f"  Skipped (no date found): {item}")
    
    return date_folders


def process_date(path_to_data, date, stationxml, rot_seed, tra_seed, path_to_output, code):
    """
    Process data for a single date.
    
    Parameters:
    -----------
    path_to_data : str
        Path to the data folder for this date
    date : str
        Date string (e.g., "2020-07-22")
    stationxml : str
        Path to StationXML file
    rot_seed : list
        List of rotation seed names
    tra_seed : list
        List of translation seed names
    path_to_output : str
        Path to output directory
    code : str
        Channel code as in .xml metadata file
    """
    print(f"\n{'='*60}")
    print(f"Processing date: {date}")
    print(f"{'='*60}")
    
    # Load Metadata
    print("Loading metadata...")
    inv = obs.read_inventory(stationxml)
    
    # Load & Combine Waveform Data
    print("Loading rotation data...")
    rot_st = obs.Stream()
    for seed in rot_seed:
        filename = seed.replace(".", "_")+".ms"
        data_to_load = os.path.join(path_to_data, filename)
        try:
            rot_st += obs.read(data_to_load)
        except Exception as e:
            if os.path.isfile(data_to_load):
                print(f"Error while loading {data_to_load}: {e}")
            else:
                print(f"File does not exist: {data_to_load}")
                return False
    
    if len(rot_st) == 0:
        print(f"Warning: No rotation data loaded for {date}")
        return False
    
    # Remove instrument response
    rot_st = rot_st.remove_sensitivity(inventory=inv)
    
    # Add dummy BJZ to rot_st
    dummy = rot_st.select(station="LIV", channel="*JY").copy()
    if len(dummy) > 0:
        dummy[0].data = dummy[0].data * 0
        dummy[0].stats.channel = dummy[0].stats.channel.replace("JY", "JZ")
        rot_st += dummy
    
    # Load translation data
    print("Loading translation data...")
    tra_st = obs.Stream()
    for seed in tra_seed:
        filename = seed.replace(".", "_")+".ms"
        data_to_load = os.path.join(path_to_data, filename)
        try:
            tra_st += obs.read(data_to_load)
        except Exception as e:
            if os.path.isfile(data_to_load):
                print(f"Error while loading {data_to_load}: {e}")
            else:
                print(f"File does not exist: {data_to_load}")
    
    if len(tra_st) == 0:
        print(f"Warning: No translation data loaded for {date}")
        return False
    
    # Remove instrument response
    tra_st = tra_st.remove_response(inventory=inv, output="VEL", water_level=10)
    
    # Rotate Data to ZNE
    print("Rotating data to ZNE...")
    angles = {}
    for _seed in rot_seed + tra_seed:
        angles[_seed] = inv.get_orientation(_seed)
    
    
    # Rotate rotation data
    try:
        rot_jz = rot_st.select(station="LIV", channel=f"*JZ")
        rot_jy = rot_st.select(station="LIV", channel=f"*JY")
        rot_jx = rot_st.select(station="LIV", channel=f"*JX")
        
        if len(rot_jz) > 0 and len(rot_jy) > 0 and len(rot_jx) > 0:
            Z, N, E = rotate2zne(
                rot_jz[0].data, 0, 90,
                rot_jy[0].data, angles[f"LI.LIV.IX.{code}JY"]['azimuth'], angles[f"LI.LIV.IX.{code}JY"]['dip'],
                rot_jx[0].data, angles[f"LI.LIV.IY.{code}JX"]['azimuth'], angles[f"LI.LIV.IY.{code}JX"]['dip'],
            )
            
            # ZYX > ZNE
            rot_jy[0].data = N
            rot_jy[0].stats.channel = code+"JN"
            rot_jx[0].data = E
            rot_jx[0].stats.channel = code+"JE"

    except Exception as e:
        print(f"Error rotating rotation data: {e}")
        return False
    
    # Rotate translation data
    try:
        tra_hz = tra_st.select(station="LIV", channel="*HZ")
        tra_hy = tra_st.select(station="LIV", channel="*HY")
        tra_hx = tra_st.select(station="LIV", channel="*HX")
        
        if len(tra_hz) > 0 and len(tra_hy) > 0 and len(tra_hx) > 0:
            Z, N, E = rotate2zne(
                tra_hz[0].data, angles[f"LI.LIV.IX.{code}HZ"]['azimuth'], angles[f"LI.LIV.IX.{code}HZ"]['dip'],
                tra_hy[0].data, angles[f"LI.LIV.IX.{code}HY"]['azimuth'], angles[f"LI.LIV.IX.{code}HY"]['dip'],
                tra_hx[0].data, angles[f"LI.LIV.IX.{code}HX"]['azimuth'], angles[f"LI.LIV.IX.{code}HX"]['dip'],
            )
            
            # ZYX > ZNE
            tra_hy[0].data = N
            tra_hy[0].stats.channel = code+"HN"
            tra_hx[0].data = E
            tra_hx[0].stats.channel = code+"HE"

    except Exception as e:
        print(f"Error rotating translation data: {e}")
        return False
    
    # Combine data
    st = obs.Stream()
    st += rot_st
    st += tra_st
    
    # Set location code to empty
    for tr in st:
        tr.stats.location = ""

    # print(st.__str__(extended=True))

    # Write Data
    print(f"Writing output to {path_to_output}...")
    if not os.path.exists(path_to_output):
        os.makedirs(path_to_output)
    
    output_file = os.path.join(path_to_output, f"{date.replace('-', '')}.mseed")
    st.write(output_file)
    print(f"Successfully processed {date} -> {output_file}")
    
    return True


def main():
    """Main function to process all dates in a given folder."""
    parser = argparse.ArgumentParser(
        description="Prepare waveform data for analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python 00_brs_prepare_data_v2.py --path_to_data ./data/downloads/miniseed_LIGO/
  python 00_brs_prepare_data_v2.py --path_to_data ./data/downloads/miniseed_LIGO/ --stationxml ./metadata/LIGO_HAN_LIV_AB.xml
        """
    )
    
    parser.add_argument(
        '--path_to_data',
        type=str,
        required=True,
        help='Path to the folder containing date subfolders (e.g., ./data/downloads/miniseed_LIGO/)'
    )
    
    parser.add_argument(
        '--stationxml',
        type=str,
        default='./metadata/LIGO_HAN_LIV_AB.xml',
        help='Path to StationXML file (default: ./metadata/LIGO_HAN_LIV_AB.xml)'
    )
    
    parser.add_argument(
        '--output_base',
        type=str,
        default='./output/',
        help='Base path for output directories (default: ./output/)'
    )
    
    parser.add_argument(
        '--code',
        type=str,
        default='B',
        help='Channel code as in .xml metadata file (default: B)'
    )
    
    args = parser.parse_args()
    
    # Configuration
    path_to_data_base = args.path_to_data
    stationxml = args.stationxml
    output_base = args.output_base
    code = args.code

    # Specify seed names
    rot_seed = [f"LI.LIV.IX.{code}JY", f"LI.LIV.IY.{code}JX"]
    tra_seed = [f"LI.LIV.IX.{code}HZ", f"LI.LIV.IX.{code}HY", f"LI.LIV.IX.{code}HX"]
    
    # Check if base path exists
    if not os.path.exists(path_to_data_base):
        print(f"Error: Path does not exist: {path_to_data_base}")
        sys.exit(1)
    
    # Check if StationXML exists
    if not os.path.exists(stationxml):
        print(f"Error: StationXML file does not exist: {stationxml}")
        sys.exit(1)
    
    # Scan data path and extract dates from folder names
    date_folders = scan_data_path(path_to_data_base)
    
    if len(date_folders) == 0:
        print(f"\nNo date folders found in {path_to_data_base}")
        print("Expected folder names containing dates in formats like:")
        print("  - YYYY-MM-DD (e.g., 2020-07-22)")
        print("  - YYYY-MM-DD_L (e.g., 2020-07-22_L)")
        print("  - YYYYMMDD (e.g., 20200722)")
        print("  - YYYY_MM_DD (e.g., 2020_07_22)")
        print("  - Or any folder name containing a date pattern")
        sys.exit(1)
    
    # Assign _01, _02, ... when several folders share the same calendar date
    dated_events = assign_output_suffixes(date_folders)
    dated_events.sort(key=lambda x: (x[1], x[0]))
    print(f"\nFound {len(dated_events)} folder(s) to process")
    for date_folder, date, suffix in dated_events:
        print(f"  {date_folder} -> {os.path.basename(event_output_dir(output_base, date, suffix))}")
    
    # Process each date folder
    successful = 0
    failed = 0
    
    for date_folder, date, suffix in dated_events:
        # Construct paths
        path_to_data = os.path.join(path_to_data_base, date_folder)
        path_to_output = event_output_dir(output_base, date, suffix)
        
        try:
            if process_date(path_to_data, date, stationxml, rot_seed, tra_seed, path_to_output, code):
                successful += 1
            else:
                failed += 1
        except Exception as e:
            print(f"Error processing {date_folder} ({date}{suffix}): {e}")
            import traceback
            traceback.print_exc()
            failed += 1
    
    # Summary
    print(f"\n{'='*60}")
    print(f"Processing complete!")
    print(f"  Successful: {successful}")
    print(f"  Failed: {failed}")
    print(f"  Total: {len(dated_events)}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
