#!/usr/bin/env python3
"""
Combine BRS Analysis Results

This script reads all JSON result files from output/event* folders and combines them:
- Backazimuth results into a CSV file
- Dispersion curves into a YAML dictionary file
"""

import os
import sys
import argparse
import json
import glob
import yaml
import pandas as pd
from pathlib import Path


def find_result_files(output_base):
    """
    Find all JSON result files in event folders.
    
    Parameters:
    -----------
    output_base : str
        Base path to output directory
        
    Returns:
    --------
    list
        List of paths to JSON result files
    """
    result_files = []
    
    # Find all event folders
    event_folders = [
        os.path.join(output_base, d) 
        for d in os.listdir(output_base) 
        if os.path.isdir(os.path.join(output_base, d)) and d.startswith('event_')
    ]
    
    # Find JSON result files in each folder
    for event_folder in event_folders:
        json_files = glob.glob(os.path.join(event_folder, "*results.json"))
        result_files.extend(json_files)
    
    return sorted(result_files)


def extract_backazimuth_data(json_file):
    """
    Extract backazimuth data from JSON file.
    
    Parameters:
    -----------
    json_file : str
        Path to JSON result file
        
    Returns:
    --------
    dict
        Dictionary with backazimuth data
    """
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    # Extract basic event info
    result = {
        'event': data.get('event', 'unknown'),
        'date': data.get('date', 'unknown'),
        'magnitude': data.get('magnitude', None),
        'epicentral_distance_km': data.get('epicentral_distance_km', None),
        'depth_km': data.get('depth_km', None),
    }
    
    # Extract velocity backazimuth data
    baz_vel = data.get('backazimuth_velocity', {})
    if baz_vel and isinstance(baz_vel, dict):
        # Handle max method (may be None)
        max_vel = baz_vel.get('max')
        if max_vel and isinstance(max_vel, dict):
            result['baz_vel_max_peak'] = max_vel.get('kde_peak', None)
            result['baz_vel_max_uncertainty'] = max_vel.get('kde_uncertainty', None)
            result['baz_vel_max_count'] = max_vel.get('count', None)
        
        # Handle mid method (may be None)
        mid_vel = baz_vel.get('mid')
        if mid_vel and isinstance(mid_vel, dict):
            result['baz_vel_mid_peak'] = mid_vel.get('kde_peak', None)
            result['baz_vel_mid_uncertainty'] = mid_vel.get('kde_uncertainty', None)
            result['baz_vel_mid_count'] = mid_vel.get('count', None)

        result['baz_vel_theoretical'] = baz_vel.get('theoretical', None)
    
    # Extract acceleration backazimuth data
    baz_acc = data.get('backazimuth_acceleration', {})
    if baz_acc and isinstance(baz_acc, dict):
        # Handle max method (may be None)
        max_acc = baz_acc.get('max')
        if max_acc and isinstance(max_acc, dict):
            result['baz_acc_max_peak'] = max_acc.get('kde_peak', None)
            result['baz_acc_max_uncertainty'] = max_acc.get('kde_uncertainty', None)
            result['baz_acc_max_count'] = max_acc.get('count', None)

        # Handle mid method (may be None)
        mid_acc = baz_acc.get('mid')
        if mid_acc and isinstance(mid_acc, dict):
            result['baz_acc_mid_peak'] = mid_acc.get('kde_peak', None)
            result['baz_acc_mid_uncertainty'] = mid_acc.get('kde_uncertainty', None)
            result['baz_acc_mid_count'] = mid_acc.get('count', None)

        result['baz_acc_theoretical'] = baz_acc.get('theoretical', None)
    
    # Extract tangent backazimuth data
    tangent_baz = data.get('tangent_backazimuth', {})
    if tangent_baz and isinstance(tangent_baz, dict):
        result['tangent_vel_mid'] = tangent_baz.get('velocity_mid', None)
        result['tangent_vel_mid_uncertainty'] = tangent_baz.get('velocity_mid_uncertainty', None)
        result['tangent_vel_max'] = tangent_baz.get('velocity_max', None)
        result['tangent_vel_max_uncertainty'] = tangent_baz.get('velocity_max_uncertainty', None)
        result['tangent_acc_mid'] = tangent_baz.get('acceleration_mid', None)
        result['tangent_acc_mid_uncertainty'] = tangent_baz.get('acceleration_mid_uncertainty', None)
        result['tangent_acc_max'] = tangent_baz.get('acceleration_max', None)
        result['tangent_acc_max_uncertainty'] = tangent_baz.get('acceleration_max_uncertainty', None)


    # Extract cross-correlation data
    cc = data.get('cross_correlation', {})
    if cc:
        result['cc_lag_samples'] = cc.get('lag_samples', None)
        result['cc_max'] = cc.get('cc_max', None)
        result['cc_zero_lag_before'] = cc.get('cc_zero_lag_before', None)
        result['cc_zero_lag_after'] = cc.get('cc_zero_lag_after', None)
        result['polarity_flip'] = cc.get('polarity_flip', None)

    return result


def extract_dispersion_data(json_file):
    """
    Extract dispersion curve data from JSON file.
    
    Parameters:
    -----------
    json_file : str
        Path to JSON result file
        
    Returns:
    --------
    dict
        Dictionary with dispersion data organized by event
    """
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    # Extract basic event info
    event_name = data.get('event', 'unknown')
    
    result = {
        'event': event_name,
        'date': data.get('date', 'unknown'),
        'magnitude': data.get('magnitude', None),
        'epicentral_distance_km': data.get('epicentral_distance_km', None),
        'depth_km': data.get('depth_km', None),
        'velocity': None,
        'acceleration': None
    }
    
    # Extract velocity dispersion curve
    disp_vel = data.get('dispersion_curve_velocity', None)
    if disp_vel and isinstance(disp_vel, dict):
        frequencies = disp_vel.get('frequencies', [])
        velocities = disp_vel.get('velocities', [])
        backazimuths = disp_vel.get('backazimuths', [])
        
        if len(frequencies) > 0:
            result['velocity'] = {
                'frequencies': frequencies,
                'velocities': velocities,
                'backazimuths': backazimuths if len(backazimuths) > 0 else None
            }
    
    # Extract acceleration dispersion curve
    disp_acc = data.get('dispersion_curve_acceleration', None)
    if disp_acc and isinstance(disp_acc, dict):
        frequencies = disp_acc.get('frequencies', [])
        velocities = disp_acc.get('velocities', [])
        backazimuths = disp_acc.get('backazimuths', [])
        
        if len(frequencies) > 0:
            result['acceleration'] = {
                'frequencies': frequencies,
                'velocities': velocities,
                'backazimuths': backazimuths if len(backazimuths) > 0 else None
            }
    
    return result


def main():
    """Main function to combine results."""
    parser = argparse.ArgumentParser(
        description="Combine BRS Analysis Results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python 00_brs_combine_results.py --output_base ./output/
  python 00_brs_combine_results.py --output_base ./output/ --output_dir ./results/
        """
    )
    
    parser.add_argument(
        '--output_base',
        type=str,
        default='./output/',
        help='Base path to output directory containing event folders (default: ./output/)'
    )
    
    parser.add_argument(
        '--output_dir',
        type=str,
        default='./output/',
        help='Directory to save output files (default: ./output/)'
    )
    
    args = parser.parse_args()
    
    output_base = args.output_base
    output_dir = args.output_dir
    
    # Check if output base exists
    if not os.path.exists(output_base):
        print(f"Error: Output base path does not exist: {output_base}")
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")
    
    # Find all result files
    print(f"Scanning for result files in: {output_base}")
    result_files = find_result_files(output_base)
    
    if len(result_files) == 0:
        print(f"No result JSON files found in {output_base}")
        sys.exit(1)
    
    print(f"Found {len(result_files)} result file(s)")
    
    # Extract backazimuth data
    print("\nExtracting backazimuth data...")
    backazimuth_data = []
    for json_file in result_files:
        try:
            data = extract_backazimuth_data(json_file)
            backazimuth_data.append(data)
        except Exception as e:
            print(f"Error reading {json_file}: {e}")
            continue
    
    # Create DataFrame and save to CSV
    if len(backazimuth_data) > 0:
        df_baz = pd.DataFrame(backazimuth_data)
        baz_csv_file = os.path.join(output_dir, 'backazimuth_results.csv')
        df_baz.to_csv(baz_csv_file, index=False)
        print(f"Saved backazimuth results to: {baz_csv_file}")
        print(f"  Rows: {len(df_baz)}, Columns: {len(df_baz.columns)}")
    else:
        print("No backazimuth data found")
    
    # Extract dispersion curve data
    print("\nExtracting dispersion curve data...")
    dispersion_data = {}
    for json_file in result_files:
        try:
            event_data = extract_dispersion_data(json_file)
            event_name = event_data['event']
            if event_name != 'unknown':
                dispersion_data[event_name] = event_data
        except Exception as e:
            print(f"Error reading {json_file}: {e}")
            continue
    
    # Save dispersion curves as YAML
    if len(dispersion_data) > 0:
        # Clean up None values and empty entries
        cleaned_data = {}
        for event_name, event_data in dispersion_data.items():
            cleaned_event = {
                'date': event_data.get('date'),
                'magnitude': event_data.get('magnitude'),
                'epicentral_distance_km': event_data.get('epicentral_distance_km'),
                'depth_km': event_data.get('depth_km'),
            }
            
            # Add velocity data if available
            if event_data.get('velocity'):
                vel_data = event_data['velocity']
                cleaned_event['velocity'] = {
                    'frequencies': vel_data.get('frequencies', []),
                    'velocities': vel_data.get('velocities', [])
                }
                if vel_data.get('backazimuths'):
                    cleaned_event['velocity']['backazimuths'] = vel_data['backazimuths']
            
            # Add acceleration data if available
            if event_data.get('acceleration'):
                acc_data = event_data['acceleration']
                cleaned_event['acceleration'] = {
                    'frequencies': acc_data.get('frequencies', []),
                    'velocities': acc_data.get('velocities', [])
                }
                if acc_data.get('backazimuths'):
                    cleaned_event['acceleration']['backazimuths'] = acc_data['backazimuths']
            
            cleaned_data[event_name] = cleaned_event
        
        disp_yaml_file = os.path.join(output_dir, 'dispersion_curves.yaml')
        with open(disp_yaml_file, 'w') as f:
            yaml.dump(cleaned_data, f, default_flow_style=False, sort_keys=True, allow_unicode=True)
        print(f"Saved dispersion curves to: {disp_yaml_file}")
        print(f"  Events: {len(cleaned_data)}")
        
        # Count data points
        total_points = 0
        for event_data in cleaned_data.values():
            if event_data.get('velocity'):
                total_points += len(event_data['velocity'].get('frequencies', []))
            if event_data.get('acceleration'):
                total_points += len(event_data['acceleration'].get('frequencies', []))
        print(f"  Total data points: {total_points}")
    else:
        print("No dispersion curve data found")
    
    print(f"\n✓ Processing complete!")
    print(f"  Backazimuth results: {len(backazimuth_data)} events")
    print(f"  Dispersion curves: {len(dispersion_data)} events")


if __name__ == "__main__":
    main()
