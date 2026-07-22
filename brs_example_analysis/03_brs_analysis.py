#!/usr/bin/env python3
"""
Analyse BRS Event Data

This script processes all event folders in output/ and creates analysis figures.
Figures are stored in output/event*/figures/
"""

import os
import sys
import argparse
import re
import glob
import gc
import copy
import yaml
import json
import numpy as np
import matplotlib.pyplot as plt
from obspy import UTCDateTime
from obspy.signal.rotate import rotate2zne
from obspy.signal.cross_correlation import correlate, xcorr_max
from obspy.signal.rotate import rotate_ne_rt
import scipy.stats as sts

from sixdegrees.sixdegrees import sixdegrees
from sixdegrees.plots.plot_waveform_cc import plot_waveform_cc
from sixdegrees.plots.plot_cwt_all import plot_cwt_all
from sixdegrees.plots.plot_velocities import plot_velocities
from sixdegrees.plots.plot_backazimuth_results import plot_backazimuth_results
from sixdegrees.plots.plot_spectra import plot_spectra
from sixdegrees.plots.plot_event_overview import plot_event_overview
from sixdegrees.plots.plot_dispersion_traces import plot_dispersion_traces
from sixdegrees.plots.plot_dispersion_curves import plot_dispersion_curves


DEFAULT_ANALYSIS_PARAMS = {
    # Bandpass and backazimuth / cross-correlation
    'fmin': 0.02,
    'fmax': 0.2,
    'Twindow': 100,
    'Toverlap': 0.75,
    'cc_threshold': 0.8,
    # Rayleigh dispersion (sd.compute_dispersion_curve)
    'dispersion': {
        'wave_type': 'rayleigh',
        'fmin': 0.01,
        'fmax': 0.5,
        'octave_fraction': 3,
        'window_factor': 6.0,
        'cc_threshold': None,  # None -> use top-level cc_threshold
        'baz_step': 1,
        'time_window_overlap': 0.75,
        'velocity_method': 'odr',
        'zero_intercept': False,
        'verbose': False,
        'n_jobs': 5,
    },
}


def merge_analysis_params(params=None):
    """Deep-merge user params into DEFAULT_ANALYSIS_PARAMS (nested dispersion supported)."""
    merged = copy.deepcopy(DEFAULT_ANALYSIS_PARAMS)
    if not params:
        return merged
    for key, value in params.items():
        if key == 'dispersion' and isinstance(value, dict):
            merged['dispersion'].update(value)
        elif key.startswith('disp_'):
            legacy = {
                'disp_fmin': 'fmin',
                'disp_fmax': 'fmax',
                'disp_fraction_of_octave': 'octave_fraction',
            }
            if key in legacy:
                merged['dispersion'][legacy[key]] = value
        else:
            merged[key] = value
    return merged


def dispersion_compute_kwargs(params, cc_threshold):
    """Build keyword arguments for sd.compute_dispersion_curve from params."""
    disp = params['dispersion']
    return {
        'wave_type': disp.get('wave_type', 'rayleigh'),
        'fmin': disp['fmin'],
        'fmax': disp['fmax'],
        'octave_fraction': disp['octave_fraction'],
        'window_factor': disp['window_factor'],
        'cc_threshold': (
            disp['cc_threshold'] if disp.get('cc_threshold') is not None else cc_threshold
        ),
        'baz_step': disp['baz_step'],
        'time_window_overlap': disp['time_window_overlap'],
        'velocity_method': disp['velocity_method'],
        'zero_intercept': disp['zero_intercept'],
        'verbose': disp['verbose'],
        'n_jobs': disp['n_jobs'],
    }


def find_yaml_file(event_folder):
    """
    Find YAML config file in event folder.
    
    Parameters:
    -----------
    event_folder : str
        Path to event folder
        
    Returns:
    --------
    str or None
        Path to YAML file, or None if not found
    """
    # Look for YYYYMMDD.yml files
    yaml_files = glob.glob(os.path.join(event_folder, "[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9].yml"))
    if len(yaml_files) > 0:
        return yaml_files[0]
    return None


def extract_backazimuth_stats(baz_results, cc_threshold, baz_theo):
    """
    Extract backazimuth statistics from baz_results using KDE estimation.
    
    Parameters:
    -----------
    baz_results : dict
        Results from compute_backazimuth (contains 'baz_max', 'baz_mid', 'cc_max', 'cc_mid' arrays)
    cc_threshold : float
        Cross-correlation threshold
    baz_theo : float
        Theoretical backazimuth
        
    Returns:
    --------
    dict
        Dictionary with backazimuth statistics for max and mid methods using KDE estimation
    """
    stats = {
        'max': {},
        'mid': {},
        'theoretical': round(float(baz_theo), 1)
    }
    
    # Extract max method stats using KDE
    if 'baz_max' in baz_results and baz_results['baz_max'] is not None:
        baz_max_array = baz_results['baz_max']
        cc_max_array = baz_results.get('cc_max', None)
        
        if len(baz_max_array) > 0:
            # Filter by cc_threshold if cc_max is available
            if cc_max_array is not None and len(cc_max_array) == len(baz_max_array):
                mask = cc_max_array >= cc_threshold
                filtered_baz = baz_max_array[mask]
                filtered_cc = cc_max_array[mask]
            else:
                filtered_baz = baz_max_array
                filtered_cc = np.ones_like(filtered_baz)  # Equal weights if no CC
            
            if len(filtered_baz) > 0:
                # Use KDE estimation if enough data points
                if len(filtered_baz) >= 5:
                    try:
                        kde_stats = sixdegrees.get_kde_stats(
                            filtered_baz, 
                            filtered_cc, 
                            _baz_steps=0.5, 
                            Ndegree=60, 
                            plot=False
                        )
                        stats['max']['kde_peak'] = float(kde_stats['baz_estimate'])
                        stats['max']['kde_uncertainty'] = float(kde_stats['kde_dev'])
                    except Exception as e:
                        # Fallback to simple statistics if KDE fails
                        stats['max']['kde_peak'] = float(np.nanmean(filtered_baz))
                        stats['max']['kde_uncertainty'] = float(np.nanstd(filtered_baz))
                else:
                    # Not enough points for KDE, use simple statistics
                    stats['max']['kde_peak'] = float(np.nanmean(filtered_baz))
                    stats['max']['kde_uncertainty'] = float(np.nanstd(filtered_baz))
                
                # Also include basic statistics
                stats['max']['count'] = int(len(filtered_baz))
            else:
                stats['max'] = None
        else:
            stats['max'] = None
    else:
        stats['max'] = None
    
    # Extract mid method stats using KDE
    if 'baz_mid' in baz_results and baz_results['baz_mid'] is not None:
        baz_mid_array = baz_results['baz_mid']
        cc_mid_array = baz_results.get('cc_mid', None)
        
        if len(baz_mid_array) > 0:
            # Filter by cc_threshold if cc_mid is available
            if cc_mid_array is not None and len(cc_mid_array) == len(baz_mid_array):
                mask = cc_mid_array >= cc_threshold
                filtered_baz = baz_mid_array[mask]
                filtered_cc = cc_mid_array[mask]
            else:
                filtered_baz = baz_mid_array
                filtered_cc = np.ones_like(filtered_baz)  # Equal weights if no CC
            
            if len(filtered_baz) > 0:
                # Use KDE estimation if enough data points
                if len(filtered_baz) >= 5:
                    try:
                        kde_stats = sixdegrees.get_kde_stats(
                            filtered_baz, 
                            filtered_cc, 
                            _baz_steps=0.5, 
                            Ndegree=60, 
                            plot=False
                        )
                        stats['mid']['kde_peak'] = float(kde_stats['baz_estimate'])
                        stats['mid']['kde_uncertainty'] = float(kde_stats['kde_dev'])
                    except Exception as e:
                        # Fallback to simple statistics if KDE fails
                        stats['mid']['kde_peak'] = float(np.nanmean(filtered_baz))
                        stats['mid']['kde_uncertainty'] = float(np.nanstd(filtered_baz))
                else:
                    # Not enough points for KDE, use simple statistics
                    stats['mid']['kde_peak'] = float(np.nanmean(filtered_baz))
                    stats['mid']['kde_uncertainty'] = float(np.nanstd(filtered_baz))
                
                # Also include basic statistics
                stats['mid']['count'] = int(len(filtered_baz))
            else:
                stats['mid'] = None
        else:
            stats['mid'] = None
    else:
        stats['mid'] = None
    
    return stats


def save_analysis_results(event_folder, config, cc_results, baz_stats_vel, baz_stats_acc, 
                         rayleigh_dispersion_vel, rayleigh_dispersion_acc, params, ev=None,
                         tangent_results_vel_mid=None, tangent_results_vel_max=None,
                         tangent_results_acc_mid=None, tangent_results_acc_max=None):
    """
    Save analysis results to a JSON file.
    
    Parameters:
    -----------
    event_folder : str
        Path to event folder
    config : dict
        Configuration dictionary
    cc_results : dict
        Cross-correlation results
    baz_stats_vel : dict
        Backazimuth statistics for velocity analysis
    baz_stats_acc : dict
        Backazimuth statistics for acceleration analysis
    rayleigh_dispersion_vel : dict
        Dispersion curve data for velocity (dict with 'frequency' and 'velocity' keys)
    rayleigh_dispersion_acc : dict
        Dispersion curve data for acceleration (dict with 'frequency' and 'velocity' keys)
    params : dict
        Analysis parameters
    ev : dict, optional
        Event information dictionary containing distance_km, depth_km, etc.
    tangent_results_vel_mid : dict, optional
        Results from compare_backazimuth_methods for velocity, mid method
    tangent_results_vel_max : dict, optional
        Results from compare_backazimuth_methods for velocity, max method
    tangent_results_acc_mid : dict, optional
        Results from compare_backazimuth_methods for acceleration, mid method
    tangent_results_acc_max : dict, optional
        Results from compare_backazimuth_methods for acceleration, max method
    """
    # Ensure cc_threshold is in parameters
    params_output = params.copy()
    if 'cc_threshold' not in params_output:
        params_output['cc_threshold'] = 0.8  # Default value
    
    # Helper function to convert numpy arrays to lists in dict
    def convert_dict_arrays(d):
        """Convert numpy arrays in dict to lists for JSON serialization."""
        if d is None:
            return None
        result = {}
        for key, value in d.items():
            if isinstance(value, np.ndarray):
                result[key] = value.tolist()
            elif isinstance(value, dict):
                result[key] = convert_dict_arrays(value)
            else:
                result[key] = value
        return result
    
    # Extract tangent backazimuth estimates and uncertainties from compare_backazimuth_methods results
    def extract_tangent_data(results_dict):
        """Extract tangent backazimuth estimate and uncertainty from compare_backazimuth_methods results.
        
        Returns:
        --------
        tuple: (estimate, uncertainty) or (None, None) if not available
        """
        if results_dict is None:
            return None, None
        
        estimate = None
        uncertainty = None
        
        # Try different possible keys for tangent estimate
        if 'tangent' in results_dict:
            tangent_data = results_dict['tangent']
            if isinstance(tangent_data, dict):
                # Look for baz_estimate
                if 'baz_estimate' in tangent_data:
                    estimate = round(float(tangent_data['baz_estimate']), 1)
                
                # Look for uncertainty (try different possible keys; sixdegrees uses 'baz_std')
                if 'uncertainty' in tangent_data:
                    uncertainty = round(float(tangent_data['uncertainty']), 1)
                elif 'baz_std' in tangent_data:
                    uncertainty = round(float(tangent_data['baz_std']), 1)
                elif 'kde_dev' in tangent_data:
                    uncertainty = round(float(tangent_data['kde_dev']), 1)
                elif 'std' in tangent_data:
                    uncertainty = round(float(tangent_data['std']), 1)
                elif 'std_dev' in tangent_data:
                    uncertainty = round(float(tangent_data['std_dev']), 1)
        
        # If tangent data is directly a number (legacy format)
        elif isinstance(results_dict.get('tangent'), (int, float)):
            estimate = round(float(results_dict['tangent']), 1)
        
        return estimate, uncertainty
    
    tangent_vel_mid, tangent_vel_mid_unc = extract_tangent_data(tangent_results_vel_mid)
    tangent_vel_max, tangent_vel_max_unc = extract_tangent_data(tangent_results_vel_max)
    tangent_acc_mid, tangent_acc_mid_unc = extract_tangent_data(tangent_results_acc_mid)
    tangent_acc_max, tangent_acc_max_unc = extract_tangent_data(tangent_results_acc_max)
    
    # Prepare output data structure
    output_data = {
        'event': config.get('project', 'unknown'),
        'date': config.get('origin_time', 'unknown'),
        'magnitude': config.get('magnitude', None),
        'epicentral_distance_km': None,
        'depth_km': None,
        'analysis_parameters': params_output,
        'cross_correlation': {
            'lag_samples': int(cc_results['lag']),
            'cc_max': round(float(cc_results['cc_max']), 2),
            'cc_zero_lag_before': round(float(cc_results['cc_before']), 2),
            'cc_zero_lag_after': round(float(cc_results['cc_after']), 2),
            'polarity_flip': cc_results['polarity_flip']
        },
        'backazimuth_velocity': baz_stats_vel,
        'backazimuth_acceleration': baz_stats_acc,
        'tangent_backazimuth': {
            'velocity_mid': tangent_vel_mid,
            'velocity_mid_uncertainty': tangent_vel_mid_unc,
            'velocity_max': tangent_vel_max,
            'velocity_max_uncertainty': tangent_vel_max_unc,
            'acceleration_mid': tangent_acc_mid,
            'acceleration_mid_uncertainty': tangent_acc_mid_unc,
            'acceleration_max': tangent_acc_max,
            'acceleration_max_uncertainty': tangent_acc_max_unc
        },
        'dispersion_curve_velocity': convert_dict_arrays(rayleigh_dispersion_vel) if rayleigh_dispersion_vel is not None else None,
        'dispersion_curve_acceleration': convert_dict_arrays(rayleigh_dispersion_acc) if rayleigh_dispersion_acc is not None else None
    }
    
    # Add event information if available
    if ev is not None:
        if 'distance_km' in ev:
            output_data['epicentral_distance_km'] = round(float(ev['distance_km']), 2)
        if 'depth_km' in ev:
            output_data['depth_km'] = round(float(ev['depth_km']), 2)
    
    # Save to JSON file and flush immediately
    output_file = os.path.join(event_folder, f"{config['project'].replace('_', '')}_results.json")
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2, default=str)
        f.flush()
        os.fsync(f.fileno())  # Force write to disk
    
    # Clear output_data to free memory
    del output_data
    
    print(f"Analysis results saved to: {output_file}")
    return output_file


def compute_cross_correlation_lag(sd, backazimuth, apply_shift=False):
    """
    Compute cross-correlation between transverse rotation (JT) and vertical translation (HZ)
    and find optimal lag.
    
    Parameters:
    -----------
    sd : sixdegrees.sixdegrees
        Sixdegrees instance with loaded data
    backazimuth : float
        Backazimuth angle in degrees
    apply_shift : bool
        If True, apply the lag shift to the HZ data in the stream (default: False)
        
    Returns:
    --------
    dict
        Dictionary containing:
        - 'lag': int, optimal lag in samples
        - 'cc_max': float, maximum cross-correlation coefficient
        - 'cc_before': float, cross-correlation at zero lag before shift
        - 'cc_after': float, cross-correlation at zero lag after shift
        - 'HZ_shifted': numpy array, shifted HZ data (if apply_shift=False)
    """
    # Get rotation data
    JN = sd.get_stream("rotation").select(channel="*N")[0].data
    JE = sd.get_stream("rotation").select(channel="*E")[0].data
    
    # Rotate to radial/transverse components
    JR, JT = rotate_ne_rt(JN, JE, backazimuth)
    
    # Get vertical translation data
    HZ = sd.get_stream("translation").select(channel="*Z")[0].data
    
    # Compute cross-correlation at zero lag before shift
    cc0_before = correlate(JT, HZ, 0)
    
    # Compute full cross-correlation and find optimal lag
    ccX = correlate(JT, HZ, len(JT))
    lag, cc_max = xcorr_max(ccX, abs_max=False)
    
    # Apply lag shift
    HZ_shifted = np.roll(HZ, lag)
    
    # Compute cross-correlation at zero lag after shift
    cc0_after = correlate(JT, HZ_shifted, 0)
    
    # Apply shift to stream if requested
    if apply_shift:
        sd.st.select(channel="*HZ")[0].data = HZ_shifted
    
    return {
        'lag': lag,
        'cc_max': cc_max,
        'cc_before': cc0_before[0],
        'cc_after': cc0_after[0],
        'HZ_shifted': HZ_shifted
    }


def clear_memory():
    """
    Clear memory by closing all matplotlib figures and forcing garbage collection.
    """
    # Close all matplotlib figures
    plt.close('all')
    # Clear matplotlib backend cache
    try:
        plt.clf()
    except:
        pass
    # Force garbage collection
    gc.collect()


def process_event(event_folder, yaml_file, save_figs=True, params=None):
    """
    Process a single event.
    
    Parameters:
    -----------
    event_folder : str
        Path to event folder
    yaml_file : str
        Path to YAML config file
    save_figs : bool
        Whether to save figures
    params : dict, optional
        Analysis parameters (merged with DEFAULT_ANALYSIS_PARAMS). Top-level keys:
        'fmin', 'fmax', 'Twindow', 'Toverlap', 'cc_threshold'.
        Nested 'dispersion' dict is passed to compute_dispersion_curve
        (wave_type, fmin, fmax, octave_fraction, window_factor, cc_threshold,
        baz_step, time_window_overlap, velocity_method, zero_intercept,
        verbose, n_jobs). Legacy keys disp_fmin / disp_fmax / disp_fraction_of_octave
        are still accepted and mapped into 'dispersion'.
        
    Returns:
    --------
    bool
        True if successful, False otherwise
    """
    params = merge_analysis_params(params)

    # Extract parameters for easier access
    fmin = params['fmin']
    fmax = params['fmax']
    Twindow = params['Twindow']
    Toverlap = params['Toverlap']
    cc_threshold = params['cc_threshold']
    dispersion_kwargs = dispersion_compute_kwargs(params, cc_threshold)

    event_name = os.path.basename(event_folder)
    print(f"\n{'='*60}")
    print(f"Processing event: {event_name}")
    print(f"{'='*60}")

    # Initialize variables to None for cleanup
    sd = None
    baz_results = None
    baz_results_rayleigh = None
    results_vel1 = None
    rayleigh_dispersion_vel = None
    rayleigh_dispersion_acc = None
    results_vel_mid = None
    results_vel_max = None
    results_acc_mid = None
    results_acc_max = None

    try:
        # Load config from YAML file
        print(f"Loading config from: {yaml_file}")
        with open(yaml_file, 'r') as f:
            config = yaml.load(f, Loader=yaml.FullLoader)
        
        # Check if mseed file exists in event folder, if not try the path from config
        mseed_file_from_config = config.get('path_to_mseed_file', '')
        mseed_basename = os.path.basename(mseed_file_from_config)
        mseed_in_event_folder = os.path.join(event_folder, mseed_basename)
        
        if os.path.exists(mseed_file_from_config):
            print(f"Using mseed file from config: {mseed_file_from_config}")
        else:
            print(f"Warning: Mseed file not found: {mseed_file_from_config}")
            return False
        
        # Create directory for plotting outputs
        if save_figs:
            figures_dir = os.path.join(event_folder, "figures")
            if not os.path.exists(figures_dir):
                os.makedirs(figures_dir)
            # Update the path to the figures
            config['path_to_figures'] = figures_dir + "/"
            print(f"Figures will be saved to: {config['path_to_figures']}")
        
        # Setup a sixdegrees instance
        print("Setting up sixdegrees instance...")
        sd = sixdegrees(conf=config)
        
        # Get key information about the event
        print("Getting event information...")
        ev = sd.get_event_info(
            origin_time=config['origin_time'],
            base_catalog="USGS",
            magnitude=config['magnitude'],
            time_margin=1000
        )
        print(f"Found event: {ev.get('magnitude', 'N/A')} at {ev.get('distance_km', 'N/A')} km")
        
        # Load data
        print("Loading data...")
        sd.load_data(config['tbeg'], config['tend'], resample_rate=config['resample_rate'])
        
        # Ensure that all traces of stream have the same sample size
        sd.trim_stream(set_common=True, set_interpolate=True)
        
        # Frequency band to filter
        print("Filtering data...")
        sd.filter_data(fmin=fmin, fmax=fmax, output=False)
        
        # Compute cross-correlation between JT and HZ and find optimal lag
        cc_results = compute_cross_correlation_lag(sd, ev['backazimuth'], apply_shift=False)
        cc_results['polarity_flip'] = 0
        print(f"Cross-correlation: zero lag: {round(cc_results['cc_before'], 2)}, "
              f"with lag {cc_results['lag']} samples: {round(cc_results['cc_after'], 2)}")

        # Apply polarity correction if needed
        if cc_results['cc_before'] < 0:
            print("Applying polarity correction...")
            # apply polarity correction
            sd.polarity_stream(pol_dict={"HZ": -1}, raw=True)

            # update cc_results after polarity correction
            cc_results = compute_cross_correlation_lag(sd, ev['backazimuth'], apply_shift=False)
            cc_results['polarity_flip'] = 1
            print(f"Cross-correlation: zero lag: {round(cc_results['cc_before'], 2)}, "
                f"with lag {cc_results['lag']} samples: {round(cc_results['cc_after'], 2)}")
        else:
            print("No polarity correction needed.")
        

        # ===== VELOCITY ANALYSIS =====
        print("\n--- Velocity Analysis ---")
        
        # Plot waveforms
        print("Plotting velocity waveforms...")
        fig = plot_waveform_cc(
            sd_object=sd,
            fmin=fmin,
            fmax=fmax,
            twin_sec=Twindow,
            twin_overlap=Toverlap,
            wave_type="rayleigh",
            runit=config['runit'],
            tunit=config['tunit'],
            data_type="velocity"
        )
        if save_figs:
            fig.savefig(f"{config['path_to_figures']}{config['project']}_VEL_waveform_cc.png", bbox_inches='tight')
            plt.close(fig)
        else:
            plt.close(fig)
        del fig
        
        # Plot spectra
        try:
            print("Plotting velocity spectra...")
            fig = plot_spectra(
                sd_object=sd,
                ylog=False,
                xlog=False,
                fill=True,
                data_type="velocity"
            )
            if save_figs:
                fig.savefig(f"{config['path_to_figures']}{config['project']}_VEL_spectra.png", bbox_inches='tight')
                plt.close(fig)
            else:
                plt.close(fig)
            del fig
        except:
            print("> Error plotting spectra.")
            pass

        # Compute backazimuth
        print("Computing backazimuth...")
        baz_results = sd.compute_backazimuth(
            wave_type='rayleigh',
            baz_step=1,
            baz_win_sec=Twindow,
            baz_win_overlap=Toverlap,
            out=True
        )
        
        # Extract backazimuth statistics
        baz_stats_vel = extract_backazimuth_stats(baz_results, cc_threshold, ev['backazimuth'])
        
        # Plot backazimuth results (max method)
        print("Plotting backazimuth results (max method)...")
        fig = plot_backazimuth_results(
            sd,
            baz_results=baz_results,
            wave_type='rayleigh',
            baz_theo=ev['backazimuth'],
            cc_threshold=cc_threshold,
            cc_method='max',
            unitscale='nano',
            data_type="velocity"
        )
        if save_figs:
            fig.savefig(f"{config['path_to_figures']}{config['project']}_VEL_backazimuth_max.png", bbox_inches='tight')
            plt.close(fig)
        else:
            plt.close(fig)
        del fig
        
        # Plot backazimuth results (mid method)
        print("Plotting backazimuth results (mid method)...")
        fig = plot_backazimuth_results(
            sd,
            baz_results,
            wave_type='rayleigh',
            baz_theo=ev['backazimuth'],
            cc_threshold=cc_threshold,
            cc_method='mid',
            unitscale='nano',
            data_type="velocity"
        )
        if save_figs:
            fig.savefig(f"{config['path_to_figures']}{config['project']}_VEL_backazimuth_mid.png", bbox_inches='tight')
            plt.close(fig)
        else:
            plt.close(fig)
        del fig
        
        # Compare backazimuth methods (mid method)
        print("Comparing backazimuth methods...")
        fig, results_vel_mid = sd.compare_backazimuth_methods(
            Twin=Twindow,
            Toverlap=Toverlap,
            cc_threshold=cc_threshold,
            baz_theo=ev['backazimuth'],
            cc_method='mid',
            plot=True,
            output=True,
            wave_types=['rayleigh', 'tangent']
        )
        if save_figs:
            fig.savefig(f"{config['path_to_figures']}{config['project']}_VEL_backazimuth_methods_mid.png", bbox_inches='tight')
            plt.close(fig)
        else:
            plt.close(fig)
        del fig

        # Compare backazimuth methods (max method)
        fig, results_vel_max = sd.compare_backazimuth_methods(
            Twin=Twindow,
            Toverlap=Toverlap,
            cc_threshold=cc_threshold,
            baz_theo=ev['backazimuth'],
            cc_method='max',
            plot=True,
            output=True,
            wave_types=['rayleigh', 'tangent']
        )
        if save_figs:
            fig.savefig(f"{config['path_to_figures']}{config['project']}_VEL_backazimuth_methods_max.png", bbox_inches='tight')
            plt.close(fig)
        else:
            plt.close(fig)
        del fig

        # Estimate phase velocities
        print("Estimating phase velocities...")
        results_vel1 = sd.compute_velocities(
            wave_type='rayleigh',
            win_time_s=Twindow,
            overlap=Toverlap,
            baz=ev['backazimuth'],
            method="ransac"
        )
        
        # Plot velocities
        print("Plotting velocities...")
        fig = plot_velocities(
            sd,
            results_vel1,
            cc_threshold=cc_threshold,
            vmax=5000,
            data_type="velocity"
        )
        if save_figs:
            fig.savefig(f"{config['path_to_figures']}{config['project']}_VEL_velocities.png", bbox_inches='tight')
            plt.close(fig)
        else:
            plt.close(fig)
        del fig

        # Compute dispersion curve for Rayleigh waves (velocity)
        print("Computing dispersion curve for Rayleigh waves (velocity)...")
        rayleigh_dispersion_vel = sd.compute_dispersion_curve(**dispersion_kwargs)
        fig_traces_rayleigh = plot_dispersion_traces(
            dispersion_results=rayleigh_dispersion_vel,
            unitscale="nano",
            data_type="velocity"
        )
        fig_curves_rayleigh = plot_dispersion_curves(
            rayleigh_results=rayleigh_dispersion_vel,
            figsize=(10, 5),
            xlog=False,                 # Logarithmic x-axis
            ylog=False,                 # Logarithmic y-axis
            title="Rayleigh Dispersion Curve",  # Optional
            show_errors=True            # Show error bars from KDE deviations
        )

        if save_figs:
            fig_traces_rayleigh.savefig(f"{config['path_to_figures']}{config['project']}_VEL_dispersion_traces.png", bbox_inches='tight')
            plt.close(fig_traces_rayleigh)
            del fig_traces_rayleigh
            fig_curves_rayleigh.savefig(f"{config['path_to_figures']}{config['project']}_VEL_dispersion_curve.png", bbox_inches='tight')
            plt.close(fig_curves_rayleigh)
            del fig_curves_rayleigh
        else:
            plt.close(fig_traces_rayleigh)
            plt.close(fig_curves_rayleigh)
            del fig_traces_rayleigh, fig_curves_rayleigh

        # Plot event overview
        print("Plotting event overview...")
        fig = plot_event_overview(
            sd, 
            baz_results,
            results_vel1,
            event_info=ev,
            wave_type='rayleigh',
            baz_theo=ev['backazimuth'],
            cc_method='max',
            cc_threshold=cc_threshold,
            vmax=5000,
            twin_sec=Twindow,
            twin_overlap=Toverlap,
            zoom_t1=30,
            zoom_t2=300,
            rot_scale_factor=2.1,
            data_type="velocity"
        )
        if save_figs:
            fig.savefig(f"{config['path_to_figures']}{config['project']}_VEL_event_overview_max.png", bbox_inches='tight')
            plt.close(fig)
        else:
            plt.close(fig)
        del fig
        
        # Clear memory after velocity analysis
        clear_memory()
        
        # ===== ACCELERATION ANALYSIS =====
        print("\n--- Acceleration Analysis ---")
        
        # Convert to acceleration and rotation rate
        print("Converting to acceleration and rotation rate...")
        east_rotation_rate = sd.st.select(station="LIV", channel="*JE").differentiate()[0].data
        north_rotation_rate = sd.st.select(station="LIV", channel="*JN").differentiate()[0].data
        sd.st.select(station="LIV", channel="*JE")[0].data = east_rotation_rate
        sd.st.select(station="LIV", channel="*JN")[0].data = north_rotation_rate
        
        vertical_acceleration = sd.st.select(station="LIV", channel="*HZ").differentiate()[0].data
        north_acceleration = sd.st.select(station="LIV", channel="*HN").differentiate()[0].data
        east_acceleration = sd.st.select(station="LIV", channel="*HE").differentiate()[0].data
        sd.st.select(station="LIV", channel="*HN")[0].data = north_acceleration
        sd.st.select(station="LIV", channel="*HE")[0].data = east_acceleration
        sd.st.select(station="LIV", channel="*HZ")[0].data = vertical_acceleration
        
        # Plot waveforms
        print("Plotting acceleration waveforms...")
        fig = plot_waveform_cc(
            sd_object=sd,
            fmin=fmin,
            fmax=fmax,
            twin_sec=Twindow,
            twin_overlap=Toverlap,
            wave_type="rayleigh",
            runit=config['runit'],
            tunit=config['tunit']
        )
        if save_figs:
            fig.savefig(f"{config['path_to_figures']}{config['project']}_ACC_waveform_cc.png", bbox_inches='tight')
            plt.close(fig)
        else:
            plt.close(fig)
        del fig
        
        # Compute backazimuth
        print("Computing backazimuth (acceleration)...")
        baz_results_rayleigh = sd.compute_backazimuth(
            wave_type='rayleigh',
            baz_step=1,
            baz_win_sec=Twindow,
            baz_win_overlap=Toverlap,
            out=True,
        )
        
        # Extract backazimuth statistics for acceleration
        baz_stats_acc = extract_backazimuth_stats(baz_results_rayleigh, cc_threshold, ev['backazimuth'])
        
        # Plot backazimuth results (max method)
        print("Plotting backazimuth results (max method, acceleration)...")
        fig = plot_backazimuth_results(
            sd,
            baz_results=baz_results_rayleigh,
            wave_type='rayleigh',
            baz_theo=ev['backazimuth'],
            cc_threshold=cc_threshold,
            cc_method='max',
            unitscale='nano'
        )
        if save_figs:
            fig.savefig(f"{config['path_to_figures']}{config['project']}_ACC_backazimuth_max.png", bbox_inches='tight')
            plt.close(fig)
        else:
            plt.close(fig)
        del fig
        
        # Plot backazimuth results (mid method)
        print("Plotting backazimuth results (mid method, acceleration)...")
        fig = plot_backazimuth_results(
            sd,
            baz_results_rayleigh,
            wave_type='rayleigh',
            baz_theo=ev['backazimuth'],
            cc_threshold=cc_threshold,
            cc_method='mid',
            unitscale='nano'
        )
        if save_figs:
            fig.savefig(f"{config['path_to_figures']}{config['project']}_ACC_backazimuth_mid.png", bbox_inches='tight')
            plt.close(fig)
        else:
            plt.close(fig)
        del fig
        
        # Compare backazimuth methods (mid method)
        print("Comparing backazimuth methods (acceleration)...")
        fig, results_acc_mid = sd.compare_backazimuth_methods(
            Twin=Twindow,
            Toverlap=Toverlap,
            cc_threshold=cc_threshold,
            cc_method='mid',
            baz_theo=ev['backazimuth'],
            plot=True,
            output=True,
            wave_types=['rayleigh', 'tangent'],
        )
        if save_figs:
            fig.savefig(f"{config['path_to_figures']}{config['project']}_ACC_backazimuth_methods_mid.png", bbox_inches='tight')
            plt.close(fig)
        else:
            plt.close(fig)
        del fig
        
        # Compare backazimuth methods (max method)
        fig, results_acc_max = sd.compare_backazimuth_methods(
            Twin=Twindow,
            Toverlap=Toverlap,
            cc_threshold=cc_threshold,
            cc_method='max',
            baz_theo=ev['backazimuth'],
            plot=True,
            output=True,
            wave_types=['rayleigh', 'tangent']
        )
        if save_figs:
            fig.savefig(f"{config['path_to_figures']}{config['project']}_ACC_backazimuth_methods_max.png", bbox_inches='tight')
            plt.close(fig)
        else:
            plt.close(fig)
        del fig
        
        # Estimate phase velocities
        print("Estimating phase velocities (acceleration)...")
        results_vel1 = sd.compute_velocities(
            wave_type='rayleigh',
            win_time_s=Twindow,
            overlap=Toverlap,
            baz=ev['backazimuth'],
            method="odr"
        )
        
        # Plot velocities
        print("Plotting velocities (acceleration)...")
        fig = plot_velocities(
            sd,
            results_vel1,
            cc_threshold=cc_threshold,
            vmax=5000
        )
        if save_figs:
            fig.savefig(f"{config['path_to_figures']}{config['project']}_ACC_velocities.png", bbox_inches='tight')
            plt.close(fig)
        else:
            plt.close(fig)
        del fig
        
        # Plot event overview
        print("Plotting event overview (acceleration)...")
        # fig = plot_event_overview(
        #     sd, 
        #     baz_results_rayleigh,
        #     results_vel1,
        #     event_info=ev,
        #     wave_type='rayleigh',
        #     baz_theo=ev['backazimuth'],
        #     cc_method='max',
        #     cc_threshold=cc_threshold,
        #     vmax=5000,
        #     twin_sec=Twindow,
        #     twin_overlap=Toverlap,
        #     zoom_t1=30,
        #     zoom_t2=300,
        #     rot_scale_factor=2.1
        # )
        # if save_figs:
        #     fig.savefig(f"{config['path_to_figures']}{config['project']}_ACC_event_overview_max.png", bbox_inches='tight')
        #     plt.close(fig)
        # else:
        #     plt.close(fig)
        # del fig

        # Compute dispersion curve for Rayleigh waves (acceleration)
        print("Computing dispersion curve for Rayleigh waves (acceleration)...")
        # rayleigh_dispersion_acc = sd.compute_dispersion_curve(**dispersion_kwargs)
        # fig_traces_rayleigh = plot_dispersion_traces(
        #     dispersion_results=rayleigh_dispersion_acc,
        #     unitscale="nano",
        #     data_type="acceleration"
        # )
        # fig_curves_rayleigh = plot_dispersion_curves(
        #     rayleigh_results=rayleigh_dispersion_acc,
        #     figsize=(10, 5),
        #     xlog=False,                 # Logarithmic x-axis
        #     ylog=False,                 # Logarithmic y-axis
        #     title="Rayleigh Dispersion Curve",  # Optional
        #     show_errors=False,            # Show error bars from KDE deviations
        # )
        # if save_figs:
        #     fig_traces_rayleigh.savefig(f"{config['path_to_figures']}{config['project']}_ACC_dispersion_traces.png", bbox_inches='tight')
        #     plt.close(fig_traces_rayleigh)
        #     del fig_traces_rayleigh
        #     fig_curves_rayleigh.savefig(f"{config['path_to_figures']}{config['project']}_ACC_dispersion_curve.png", bbox_inches='tight')
        #     plt.close(fig_curves_rayleigh)
        #     del fig_curves_rayleigh
        # else:
        #     plt.close(fig_traces_rayleigh)
        #     plt.close(fig_curves_rayleigh)
        #     del fig_traces_rayleigh, fig_curves_rayleigh
        
        # Save analysis results to file immediately after processing
        print("Saving analysis results...")
        save_analysis_results(
            event_folder=event_folder,
            config=config,
            cc_results=cc_results,
            baz_stats_vel=baz_stats_vel,
            baz_stats_acc=baz_stats_acc,
            rayleigh_dispersion_vel=rayleigh_dispersion_vel,
            rayleigh_dispersion_acc=rayleigh_dispersion_acc,
            params=params,
            ev=ev,
            tangent_results_vel_mid=results_vel_mid,
            tangent_results_vel_max=results_vel_max,
            tangent_results_acc_mid=results_acc_mid,
            tangent_results_acc_max=results_acc_max
        )
        
        # Clear large data structures immediately after saving
        print("Clearing memory...")
        del baz_results, baz_results_rayleigh, results_vel1
        del rayleigh_dispersion_vel, rayleigh_dispersion_acc
        del results_vel_mid, results_vel_max, results_acc_mid, results_acc_max
        del baz_stats_vel, baz_stats_acc, cc_results, ev, config
        
        # Clear sixdegrees instance and its data
        if sd is not None:
            try:
                # Clear stream data if it exists
                if hasattr(sd, 'st') and sd.st is not None:
                    sd.st.clear()
                if hasattr(sd, 'rotation') and sd.rotation is not None:
                    sd.rotation.clear()
                if hasattr(sd, 'translation') and sd.translation is not None:
                    sd.translation.clear()
            except:
                pass
            del sd
        
        # Force memory cleanup
        clear_memory()
        
        print(f"\n✓ Successfully processed {event_name}")
        return True
        
    except Exception as e:
        print(f"\n✗ Error processing {event_name}: {e}")
        import traceback
        traceback.print_exc()
        
        # Clean up on error as well
        clear_memory()
        if sd is not None:
            try:
                if hasattr(sd, 'st') and sd.st is not None:
                    sd.st.clear()
            except:
                pass
            del sd
        gc.collect()
        
        return False


def main():
    """Main function to process all events."""
    parser = argparse.ArgumentParser(
        description="Analyse BRS Event Data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python 00_brs_analysis_v2.py --output_base ./output/
        """
    )
    
    parser.add_argument(
        '--output_base',
        type=str,
        default='./output/',
        help='Base path to output directory containing event folders (default: ./output/)'
    )
    
    parser.add_argument(
        '--event',
        type=str,
        default=None,
        help='Process only a specific event folder (e.g., event_20200722). If not specified, processes all events.'
    )
    
    parser.add_argument(
        '--no-figures',
        action='store_true',
        help='Do not save figures (for testing)'
    )
    
    args = parser.parse_args()
    
    output_base = args.output_base
    save_figs = not args.no_figures
    
    # Analysis parameters (see DEFAULT_ANALYSIS_PARAMS for full defaults)
    params = merge_analysis_params({
        'fmin': 0.02,
        'fmax': 0.2,
        'Twindow': 50,
        'Toverlap': 0.75,
        'cc_threshold': 0.8,
        'dispersion': {
            'wave_type': 'rayleigh',
            'fmin': 0.02,
            'fmax': 0.2,
            'octave_fraction': 2,
            'window_factor': 6.0,
            'cc_threshold': 0.8,
            'baz_step': 1,
            'time_window_overlap': 0.75,
            'velocity_method': 'odr',
            'zero_intercept': False,
            'verbose': False,
            'n_jobs': 4,
        },
    })

    # Check if output base exists
    if not os.path.exists(output_base):
        print(f"Error: Output base path does not exist: {output_base}")
        sys.exit(1)
    
    # Find all event folders
    if args.event:
        # Process specific event
        event_folders = [os.path.join(output_base, args.event)]
        if not os.path.exists(event_folders[0]):
            print(f"Error: Event folder does not exist: {event_folders[0]}")
            sys.exit(1)
    else:
        # Find all event_* folders
        event_folders = [
            os.path.join(output_base, d) 
            for d in os.listdir(output_base) 
            if os.path.isdir(os.path.join(output_base, d)) and d.startswith('event_')
        ]
        event_folders.sort()
    
    if len(event_folders) == 0:
        print(f"No event folders found in {output_base}")
        sys.exit(1)
    
    print(f"Found {len(event_folders)} event folder(s) to process")
    
    # Process each event
    successful = 0
    failed = 0
    
    for event_folder in event_folders:
        # Find YAML file
        yaml_file = find_yaml_file(event_folder)
        if yaml_file is None:
            print(f"\n⚠ Skipping {os.path.basename(event_folder)}: No YAML config file found")
            failed += 1
            continue
        
        if process_event(event_folder, yaml_file, save_figs=save_figs, params=params):
            successful += 1
        else:
            failed += 1
        
        # Aggressive memory cleanup after each event
        # Note: process_event already handles cleanup internally, but we ensure
        # everything is cleared here as well
        clear_memory()
    
    # Summary
    print(f"\n{'='*60}")
    print(f"Processing complete!")
    print(f"  Successful: {successful}")
    print(f"  Failed: {failed}")
    print(f"  Total: {len(event_folders)}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
