"""Helpers to build dispersion curve tables from YAML/JSON (pandas-safe)."""
from __future__ import annotations

import numpy as np

try:
    import pandas as pd
except ImportError:  # pragma: no cover
    pd = None


def check_pandas_working() -> None:
    """Raise with fix instructions if pandas cannot build a simple DataFrame."""
    if pd is None:
        raise ImportError("pandas is required for dispersion plotting notebooks")
    try:
        pd.DataFrame({"velocity": [1.0]}, index=[0.01])
    except TypeError as exc:
        if "Cannot convert numpy.ndarray" in str(exc):
            raise RuntimeError(
                "pandas/numpy in this environment are incompatible (common when pandas>=2.2 "
                "was installed into the sixdegrees env). Fix with:\n"
                "  conda activate sixdegrees\n"
                "  conda install 'pandas>=1.5,<2.2' 'numpy>=1.24,<2' --force-reinstall\n"
                "Then restart the Jupyter kernel."
            ) from exc
        raise


def dispersion_to_dataframe(
    frequencies,
    velocities,
    value_col: str = "velocity",
):
    """
    One-event dispersion table: frequency index, single velocity column.

    Coerces inputs to 1D float arrays and uses Python scalars for pandas compatibility.
    """
    check_pandas_working()

    f = np.asarray(frequencies, dtype=float).ravel()
    v = np.asarray(velocities, dtype=float).ravel()
    if f.size != v.size:
        n = min(f.size, v.size)
        f, v = f[:n], v[:n]

    valid = np.isfinite(f) & (f > 0) & np.isfinite(v) & (v > 0)
    f, v = f[valid], v[valid]
    if f.size == 0:
        return None

    # Deduplicate by frequency (first wins), sorted
    seen = {}
    for freq, vel in sorted(zip(f.tolist(), v.tolist()), key=lambda x: x[0]):
        if freq not in seen:
            seen[float(freq)] = float(vel)

    index = list(seen.keys())
    values = [seen[k] for k in index]

    df = pd.DataFrame({value_col: values}, index=index)
    df.index.name = "frequency"
    return df


def merge_dispersion_tables(df_dict, target_frequencies):
    """Align per-event tables onto a common frequency index (NaN where missing)."""
    check_pandas_working()

    target = np.asarray(target_frequencies, dtype=float)
    target = np.unique(np.sort(target))

    merged = pd.DataFrame(index=target)
    merged.index.name = "frequency"

    for event_name, df_event in df_dict.items():
        col = df_event.iloc[:, 0]
        merged[event_name] = col.reindex(target)

    return merged
