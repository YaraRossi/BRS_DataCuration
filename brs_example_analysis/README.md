# BRS Analysis (LIGO LIV)

Analysis workflow for **broadband rotational seismology (BRS)** at the LIGO Hanford **LIV** station. The project combines co-located **translation** (IX) and **rotation** (tilt) MiniSEED data, then uses the [sixdegrees](https://github.com/andbrocode/SixDegreesOfFreedom) package to estimate event **backazimuth** and **Rayleigh-wave dispersion** from rotation–translation correlations.

Processing targets teleseismic events: raw downloads are prepared into per-event combined MiniSEED files, configured in YAML, analysed interactively or in batch, and summarized for multi-event comparison plots.

## Environment setup

The project uses a conda environment named `sixdegrees` (see notebook kernels). Dependencies are listed in `[environment.yml](environment.yml)`.

### Create the environment (new install)

From the repository root:

```bash
conda env create -f environment.yml
conda activate sixdegrees
```

This installs Python 3.10–3.12, **ObsPy**, **NumPy** (< 2), **SciPy**, **Matplotlib**, **pandas** (< 2.2), **PyYAML**, **Cartopy**, Jupyter, **sixdegrees**, and **acoustics** (via pip).

### Update an existing environment

```bash
conda env update -f environment.yml --prune
conda activate sixdegrees
```



### pandas / NumPy compatibility

The dispersion plotting notebooks build pandas DataFrames from YAML. If you see `TypeError: Cannot convert numpy.ndarray to numpy.ndarray`, the kernel likely has **pandas ≥ 2.2** with an incompatible NumPy build. Reinstall pinned versions and restart the kernel:

```bash
conda activate sixdegrees
conda install "pandas>=1.5,<2.2" "numpy>=1.24,<2" --force-reinstall
```

Helpers in `[brs_dispersion_utils.py](brs_dispersion_utils.py)` (`check_pandas_working`, `dispersion_to_dataframe`, `merge_dispersion_tables`) are also inlined in `05_brs_plot_dispersion_results.ipynb` imports.

### Jupyter kernel

Register the environment so notebooks can select the `sixdegrees` kernel:

```bash
conda activate sixdegrees
python -m ipykernel install --user --name sixdegrees --display-name "sixdegrees"
```

In Jupyter / VS Code / Cursor, choose the **sixdegrees** kernel when opening `.ipynb` files.

### Install sixdegrees from source (developers)

If you work on the [SixDegreesOfFreedom](https://github.com/andbrocode/SixDegreesOfFreedom) package locally, install it in editable mode after creating the env:

```bash
conda activate sixdegrees
pip install -e /path/to/SixDegreesOfFreedom
```

Remove any `prefix:` line before committing (conda adds your local path).

### Other requirements

- Station metadata: `metadata/LIGO_HAN_LIV_AB.xml` (or equivalent StationXML)
- Raw per-channel MiniSEED files under `data/downloads/` (from FDSN or local copies)
- Create `./figures/` before saving multi-event plots if it does not already exist



## Notebooks


| Notebook                                 | Description                                                                                                                                                                                                                                                                                                            |
| ---------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `00_brs_check_data_manually.ipynb`       | Optional manual QC before the automated pipeline. Loads translation and rotation MiniSEED from a download folder, applies downsampling and instrument correction, rotates to NEZ and radial/transverse, and plots overlays (e.g. vertical velocity vs. transverse rotation) to check timing, amplitude, and coherency. |
| `01_brs_make_config.ipynb`               | Builds the per-event YAML config for sixdegrees: time window, origin time, magnitude, station coordinates, channel SEED IDs, paths to MiniSEED and inventories, and processing options. Writes `./config/<project>.yml`.                                                                                               |
| `02_brs_prepare_data.ipynb`              | Prepares one event at a time: loads rotation and translation `.ms` files, corrects instruments, rotates to NEZ, merges into a single stream, and writes `./output/event_<YYYYMMDD>/<YYYYMMDD>.mseed`. Same logic as the batch script `02_brs_prepare_data.py`.                                                         |
| `03_brs_analysis.ipynb`                  | Interactive single-event analysis. Loads YAML from an event folder, runs sixdegrees (waveform CC, CWT, backazimuth, dispersion), saves figures and a JSON results file. Mirrors `03_brs_analysis.py` for batch runs.                                                                                                   |
| `05_brs_plot_backazimuth_results.ipynb`  | Multi-event plots from combined backazimuth results: cross-correlation scores, backazimuth vs. theoretical values, deviation histograms/KDE, and summary maps. Default input: `./output/backazimuth_results.csv`; figures: `./figures/`.                                                                               |
| `05_brs_plot_dispersion_results.ipynb`   | Multi-event dispersion plots (per-event curves, medians, percentiles). Default input: `./output/dispersion_curves.yaml`; figures: `./figures/`. Optional `EXCLUDE_EVENT_DATES` to drop events by origin date; optional renames for duplicate-date events (e.g. `event_20230206_01` → `event_20230206a`).               |
| `05_brs_plot_dispersion_waveforms.ipynb` | Single-event dispersion waveform panel (per-band + broadband overlay) from `output/event_*/` `*_results.json` and cached or recomputed dispersion output; figures under `./figures/`.                                                                                                                                  |




## Scripts


| Script                      | Description                                                                                                                                                                                                                |
| --------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `02_brs_prepare_data.py`    | Batch version of notebook 02: scans `data/downloads/` for date-named folders and writes combined MiniSEED under `output/event_<YYYYMMDD>/` (or `event_<YYYYMMDD>_01`, `_02`, … when multiple folders share the same date). |
| `03_brs_analysis.py`        | Batch sixdegrees analysis for all `output/event_*` folders that contain a `YYYYMMDD.yml` config and a valid MiniSEED path. Writes per-event `*_results.json` and optional figures under `output/event_*/figures/`.         |
| `04_brs_combine_results.py` | Aggregates per-event JSON results into `backazimuth_results.csv` and `dispersion_curves.yaml`.                                                                                                                             |
| `brs_dispersion_utils.py`   | Shared helpers for dispersion plotting notebooks (pandas-safe DataFrame build/merge).                                                                                                                                      |




## Workflow

Follow these steps for each new event (interactive path). Batch alternatives are noted where scripts replace notebooks.

### 1. Obtain data

Place per-channel MiniSEED files (e.g. `LI_LIV_IX_BHZ.ms`) in a folder under `data/downloads/`, such as `data/downloads/miniseed_LIGO/<date>/` or an event-named folder (e.g. `Jun28.2024/`).

### 2. Manual quality check (optional)

Open `00_brs_check_data_manually.ipynb`, set `path_to_data`, `tbeg`/`tend`, channel code, and `ev_back`, then run through the plots to confirm the waveforms look reasonable.

### 3. Prepare combined MiniSEED

**Notebook:** `02_brs_prepare_data.ipynb` — set `date`, `path_to_data`, `code` (`B` or `D`), and `path_to_output`, then run all cells.

**Batch:** prepare all dated folders in one run (run once per band code):

```bash
python 02_brs_prepare_data.py --path_to_data ./data/downloads/miniseed_LIGO/ --stationxml ./metadata/LIGO_HAN_LIV_AB.xml --output_base ./output/ --code D
python 02_brs_prepare_data.py --path_to_data ./data/downloads/miniseed_LIGO/ --stationxml ./metadata/LIGO_HAN_LIV_AB.xml --output_base ./output/ --code B
```

Output examples:

- `output/event_20240628/20240628.mseed`
- `output/event_20230206_01/20230206.mseed` and `output/event_20230206_02/20230206.mseed` when two download folders map to the same calendar date



### 4. Create event configuration

**Notebook:** `01_brs_make_config.ipynb` — edit the `config` dict (times, magnitude, `project`, SEED IDs, `path_to_mseed_file`, inventories), then run to write `./config/<project>.yml`.

Point `path_to_mseed_file` at the combined file from step 3, for example:

```text
./output/event_20240628/20240628.mseed
```

Copy (or symlink) the config into the matching event folder under `output/`, named as eight digits plus `.yml`:

```bash
cp config/2024_06_28.yml output/event_20240628/20240628.yml
```

Edit the copied YAML if paths differ. Analysis step 5 looks for `output/event_<name>/<YYYYMMDD>.yml` (use the same eight-digit date in the filename even when the folder has a `_01` / `_02` suffix).

### 5. Run analysis

**Notebook:** `03_brs_analysis.ipynb` — set `date` (and postfix if needed, e.g. `_01`) and load the YAML from the event folder.

**Batch:** process every configured event:

```bash
python 03_brs_analysis.py --output_base ./output/
```

Single event:

```bash
python 03_brs_analysis.py --output_base ./output/ --event event_20240628
```

Each successful run produces `<project>_results.json` (and figures under `output/event_*/figures/` when saving is enabled).

### 6. Combine results across events

By default, combined tables are written next to the per-event folders under `./output/`:

```bash
python 04_brs_combine_results.py --output_base ./output/
```

To write them elsewhere (e.g. a separate summary directory):

```bash
python 04_brs_combine_results.py --output_base ./output/ --output_dir ./results/
```

This creates `backazimuth_results.csv` and `dispersion_curves.yaml` in the chosen `--output_dir`.

### 7. Multi-event visualization

1. Run `04_brs_combine_results.py` (step 6).
2. Open `05_brs_plot_backazimuth_results.ipynb` and `05_brs_plot_dispersion_results.ipynb`.
3. Confirm input paths (defaults load from `./output/`):
  - backazimuth: `./output/backazimuth_results.csv`
  - dispersion: `./output/dispersion_curves.yaml`
   If you used `--output_dir ./results/`, point the notebooks at that directory instead.
4. In the dispersion notebook, edit `EXCLUDE_EVENT_DATES` (e.g. `"20190706"` or `"2019-07-06"`) to remove outliers from `df_vel` before statistics and plots; optionally rename duplicate-date columns (e.g. `event_20230206_01` → `event_20230206a`).
5. For per-event dispersion waveform stacks, use `05_brs_plot_dispersion_waveforms.ipynb`.

Summary figures from the `05_*` plotting notebooks are saved under `./figures/` (create the directory if needed). Per-event analysis figures remain under `output/event_*/figures/` when enabled in the event config.

---



## File structure

```text
BRS_analysis/
├── README.md
├── environment.yml                    # Conda env spec (name: sixdegrees)
├── brs_dispersion_utils.py            # Dispersion DataFrame helpers for plotting
│
├── 00_brs_check_data_manually.ipynb   # Optional manual waveform QC
├── 01_brs_make_config.ipynb           # Per-event YAML configuration
├── 02_brs_prepare_data.ipynb          # Prepare one combined MiniSEED file
├── 02_brs_prepare_data.py             # Batch prepare (all download folders)
├── 03_brs_analysis.ipynb              # Interactive sixdegrees analysis
├── 03_brs_analysis.py                 # Batch analysis
├── 04_brs_combine_results.py          # Merge JSON → CSV / YAML
├── 05_brs_plot_backazimuth_results.ipynb
├── 05_brs_plot_dispersion_results.ipynb
├── 05_brs_plot_dispersion_waveforms.ipynb
│
├── config/                            # Event YAML configs (from notebook 01)
│   └── <project>.yml
│
├── metadata/                          # StationXML inventories
│   └── LIGO_HAN_LIV_AB.xml
│
├── data/
│   ├── downloads/                     # Raw per-channel MiniSEED (.ms)
│   │   ├── miniseed_LIGO/           # Date-named subfolders (batch prepare)
│   │   └── <event folders>/         # e.g. Jun28.2024/, Dec08.2025/ (optional)
│
├── output/
│   ├── backazimuth_results.csv        # From 04_brs_combine_results (default dir)
│   ├── dispersion_curves.yaml
│   └── event_<YYYYMMDD>/            # Per-event workspace (optional _01, _02, …)
│       ├── <YYYYMMDD>.mseed         # Combined prepared data (step 3)
│       ├── <YYYYMMDD>.yml           # Config copy for analysis (step 4)
│       ├── <project>_results.json   # Analysis output (step 5)
│       ├── dispersion_*.pkl         # Cached dispersion (optional)
│       └── figures/                 # Per-event plots (optional)
│
├── results/                         # Optional; use --output_dir ./results/ in step 6
│   ├── backazimuth_results.csv
│   └── dispersion_curves.yaml
│
└── figures/                         # Multi-event summary plots (05_* notebooks)
```



## Disclaimer

**Disclaimer**: This documentation was created with the assistance of AI, yet verified by human.