import matplotlib.pyplot as plt
import numpy as np


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
