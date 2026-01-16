import matplotlib.pyplot as plt
import numpy as np
from obspy import read, read_inventory


def plot_traces_offset(traces, label, offset=None):
    """
    Plots multiple Obspy traces on the same subplot with vertical offsets so that all traces are visible.

    Parameters
    ----------
    traces : list
        A list of Obspy Trace objects.
    offset : float, optional
        The vertical offset between consecutive traces. If None, an offset is computed automatically.
    """
    # Automatically compute an offset if not provided.
    if offset is None:
        # Determine the maximum absolute amplitude among all traces.
        max_amp = max(np.max(np.abs(tr[0].data)) for tr in traces)
        # Multiply by a factor (e.g., 1.2) to ensure separation between traces.
        offset = max_amp * 1.2

    fig, ax = plt.subplots(figsize=(10, 6))

    for i, tr in enumerate(traces):
        # Create a time vector using the trace's sampling interval (delta).
        t = np.arange(len(tr[0].data)) * tr[0].stats.delta
        # Plot each trace with an offset. The label helps to identify the trace.
        ax.plot(tr[0].times('matplotlib'), tr[0].data + i * offset, label=label[i])

    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Amplitude (with offsets)")
    ax.legend()
    ax.xaxis_date()
    plt.show()

import numpy
from numpy import cos, sin
def rot_vec(phi, theta, psi, data):
    """
        Applies a rotation matrix to a given vector based on Euler angles.

        Parameters:
        phi (float): The roll angle (rotation about the x-axis) in radians.
        theta (float): The pitch angle (rotation about the y-axis) in radians.
        psi (float): The yaw angle (rotation about the z-axis) in radians.
        data (numpy.ndarray): The input vector to be rotated. Shape: (3,).

        Returns:
        numpy.ndarray: The rotated vector after applying the rotation matrix. Shape: (3,).

        Notes:
        - The rotation matrix R is computed based on the given Euler angles (phi, theta, psi).
        - The function returns the result of applying the rotation matrix R to the input vector data.
        """
    R = numpy.array([[cos(theta) * cos(psi), cos(theta) * sin(psi), -sin(theta)],
                     [sin(phi) * sin(theta) * cos(psi) - cos(phi) * sin(psi),
                      sin(phi) * sin(theta) * sin(psi) + cos(phi) * cos(psi), cos(theta) * sin(phi)],
                     [cos(phi) * sin(theta) * cos(psi) + sin(phi) * sin(psi),
                      cos(phi) * sin(theta) * sin(psi) - sin(phi) * cos(psi), cos(theta) * cos(phi)]])

    return R.dot(data)

def circular_positive_midpoint(corr_values):
    corr_values = numpy.asarray(corr_values)
    n = len(corr_values)

    # Step 1: Double the array to handle wrap-around
    doubled = numpy.concatenate([corr_values, corr_values])

    # Step 2: Identify where values are positive
    positive_mask = doubled > 0

    # Step 3: Find contiguous positive segments
    starts_ends = []
    in_segment = False
    for i, val in enumerate(positive_mask):
        if val and not in_segment:
            start = i
            in_segment = True
        elif not val and in_segment:
            end = i
            starts_ends.append((start, end))
            in_segment = False
    # Close final segment if it ends at the array end
    if in_segment:
        starts_ends.append((start, len(positive_mask)))

    # Step 4: Find longest segment that starts within the first 360 values
    max_len = 0
    best_mid = None
    for start, end in starts_ends:
        if start < n:  # only consider segments starting in the original range
            length = end - start
            if length > max_len:
                max_len = length
                best_mid = start + length // 2

    if best_mid is None:
        return None  # or numpy.nan if no positive segments

    # Step 5: Wrap back to original azimuth range
    return best_mid % n


def supylabel2(fig, s, **kwargs):
    defaults = {
        "x": 0.96,
        "y": 0.5,
        "horizontalalignment": "center",
        "verticalalignment": "center",
        "rotation": "vertical",
        "rotation_mode": "anchor",
        "size": plt.rcParams["figure.labelsize"],  # matplotlib >= 3.6
        "weight": plt.rcParams["figure.labelweight"],  # matplotlib >= 3.6
    }
    kwargs["s"] = s
    # kwargs = defaults | kwargs  # python >= 3.9
    kwargs = {**defaults, **kwargs}
    fig.text(**kwargs)

def get_ligo_data(folder, location):
    root_data = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/miniseed_LIGO' #zenodo_upload/miniseed_LIGO'
    root_xml = '/Users/yararossi/Documents/Work/Towards_Quantification/3_Projects/CRS_Collab/1_Data/metadata'

    if location == 'LIV':
        ## 1. READ DATA
        # BRS
        LIV_EY_DJX = read('%s/%s/LI_LIV_EY_DJX.ms' %(root_data,folder))
        LIV_EX_DJY = read('%s/%s/LI_LIV_EX_DJY.ms' %(root_data,folder))
        LIV_IY_DJX = read('%s/%s/LI_LIV_IY_DJX.ms' %(root_data,folder))
        LIV_IX_DJY = read('%s/%s/LI_LIV_IX_DJY.ms' %(root_data,folder))

        # seis
        LIV_EY_DH = read('%s/%s/LI_LIV_EY_DHX.ms' %(root_data,folder))
        LIV_EY_DH += read('%s/%s/LI_LIV_EY_DHY.ms' %(root_data,folder))
        LIV_EY_DH += read('%s/%s/LI_LIV_EY_DHZ.ms' %(root_data,folder))
        LIV_EX_DH = read('%s/%s/LI_LIV_EX_DHX.ms' %(root_data,folder))
        LIV_EX_DH += read('%s/%s/LI_LIV_EX_DHY.ms' %(root_data,folder))
        LIV_EX_DH += read('%s/%s/LI_LIV_EX_DHZ.ms' %(root_data,folder))

        LIV_IY_DH = read('%s/%s/LI_LIV_IY_DHX.ms' %(root_data,folder))
        LIV_IY_DH += read('%s/%s/LI_LIV_IY_DHY.ms' %(root_data,folder))
        LIV_IY_DH += read('%s/%s/LI_LIV_IY_DHZ.ms' %(root_data,folder))
        LIV_IX_DH = read('%s/%s/LI_LIV_IX_DHX.ms' %(root_data,folder))
        LIV_IX_DH += read('%s/%s/LI_LIV_IX_DHY.ms' %(root_data,folder))
        LIV_IX_DH += read('%s/%s/LI_LIV_IX_DHZ.ms' %(root_data,folder))

        if '_L' in folder:
            response = read_inventory('%s/LIGO_HAN_LIV_lowsps.xml' % root_xml)
        else:
            response = read_inventory('%s/LIGO_HAN_LIV.xml' % root_xml)

        ## 3. REMOVE RESPONSE
        # BRS
        LIV_EY_DJX.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
        LIV_EX_DJY.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
        LIV_IY_DJX.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
        LIV_IX_DJY.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)

        # SEIS
        LIV_EY_DH.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
        LIV_EX_DH.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
        LIV_IY_DH.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
        LIV_IX_DH.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)

        return LIV_EY_DJX, LIV_EX_DJY, LIV_IY_DJX, LIV_IX_DJY, LIV_IX_DH, LIV_IY_DH, LIV_EX_DH, LIV_EY_DH

    if location == 'HAN':
        # BRS
        HAN_EY_DJX = read('%s/%s/LI_HAN_EY_DJX.ms' % (root_data, folder))
        HAN_EX_DJY = read('%s/%s/LI_HAN_EX_DJY.ms' % (root_data, folder))

        # seis
        HAN_EY_DH = read('%s/%s/LI_HAN_EY_DHX.ms' % (root_data, folder))
        HAN_EY_DH += read('%s/%s/LI_HAN_EY_DHY.ms' % (root_data, folder))
        HAN_EY_DH += read('%s/%s/LI_HAN_EY_DHZ.ms' % (root_data, folder))
        HAN_EX_DH = read('%s/%s/LI_HAN_EX_DHX.ms' % (root_data, folder))
        HAN_EX_DH += read('%s/%s/LI_HAN_EX_DHY.ms' % (root_data, folder))
        HAN_EX_DH += read('%s/%s/LI_HAN_EX_DHZ.ms' % (root_data, folder))

        '''HAN_IY_DH = read('%s/%s/LI_HAN_IY_DHX.ms' % (root_data, folder))
        HAN_IY_DH += read('%s/%s/LI_HAN_IY_DHY.ms' % (root_data, folder))
        HAN_IY_DH += read('%s/%s/LI_HAN_IY_DHZ.ms' % (root_data, folder))'''
        try:
            HAN_VA_DH = read('%s/%s/LI_HAN_VA_DHX.ms' % (root_data, folder))
            HAN_VA_DH += read('%s/%s/LI_HAN_VA_DHY.ms' % (root_data, folder))
        except:
            HAN_VA_DH = None

        if '_L' in folder:
            response = read_inventory('%s/LIGO_HAN_LIV_lowsps.xml' % root_xml)
        else:
            response = read_inventory('%s/LIGO_HAN_LIV.xml' % root_xml)

        ## 3. REMOVE RESPONSE
        # BRS
        HAN_EY_DJX.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
        HAN_EX_DJY.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)

        # SEIS
        HAN_EY_DH.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
        HAN_EX_DH.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
        #HAN_IY_DH.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
        try:
            HAN_VA_DH.remove_response(inventory=response, pre_filt=None, output="DEF", water_level=60, plot=False)
        except:
            HAN_VA_DH = None

        return HAN_EY_DJX, HAN_EX_DJY, HAN_EY_DH, HAN_EX_DH, HAN_VA_DH # HAN_IY_DH,

