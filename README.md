# Validated Metrics — scalp EEG computational metrics

This repository contains the computational metrics used at Beth Lopour's lab to
extract features from scalp EEG signals. Each metric has been rigorously
validated and cross-validated across datasets and simulations, and discussed
with clinical professionals to make sure the clinical interpretations drawn
from them are sound. They have been used in papers supporting the diagnosis of
infants with epilepsy (infantile epileptic spasms) and the assessment of their
response to hormonal therapies.

## The metrics

| Metric | What it measures | Main functions |
|---|---|---|
| Amplitude range | Overall EEG amplitude per epoch (median over 1-s sub-epochs) | `Calc_Amplitude_Range.m` |
| Power spectral density | Band power (delta, theta, alpha, beta, broadband), in dB, per epoch (median over 5-s sub-epochs) | `Calc_SEF_SpectralPower.m` |
| Spectral edge frequency (SEF) | Frequency below which 95% of the spectral power lies | `Calc_SEF_SpectralPower.m` |
| Shannon entropy | Signal complexity from the amplitude distribution, per band | `Calc_ShannonEntropy.m` |
| Permutation entropy | Signal complexity from the temporal ordering of samples, per band | `Calc_PermutationEntropy.m`, `PermutationEntropy.m` |
| Phase-amplitude coupling (PAC) | Coupling between delta phase and gamma amplitude: modulation index (Tort) and mean vector length (Canolty) | `Calc_PAC.m`, `Tort_MI.m`, `Canolty_MVL.m` |
| Long-range temporal correlations | Detrended fluctuation analysis (DFA): exponent and intercept, per band | `nbt_doDFA_Venus.m` |
| Functional connectivity | Cross-correlation networks; number of strong connections (> 0.1) and mean of the top 10% of connections | `calculate_crossCorrelation_v2.m`, `crossCorrFn_v2.m` |

## Repository layout

**Driver scripts for the NIMBIS clipped files** (one value set per clean 30-s
epoch / block; these are the scripts to run):

* `NIMBIS_Get_Metrics_clipped_files.m` — amplitude, PSD, SEF, Shannon and
  permutation entropy, and PAC, per clean 30-second epoch (channel Cz,
  linked-ears reference).
* `NIMBIS_DFA_clipped_files.m` — DFA exponent and intercept over independent,
  non-overlapping 5-minute blocks (10 consecutive clean epochs within one
  continuous piece of signal).
* `NIMBIS_DFA_pooled_short_clips.m` — DFA for the patients whose clips are too
  short for a single 5-minute block: the clean-epoch envelopes of all of a
  patient's clips (same state) are pooled first, then the 5-minute blocks are
  taken over the pooled epochs. Results are labeled `_pooled` so the table
  scripts can keep them apart.
* `NIMBIS_Get_XCorr_perBlock_clipped_files.m` — connectivity over 3-minute
  blocks (6 consecutive clean epochs); strength = number of connections above
  0.1, mean10 = mean of the strongest 10% of connections.
* `NIMBIS_Get_XCorr_perEpoch_clipped_files.m` — the same connectivity values
  per single clean 30-second epoch.
* `fable_NIMBIS_concat_mainMetrics.m` — merges two partial main-metrics runs
  of one site into `all_<state>_MainMetrics_<SITE>_...` files.

All five drivers share the same conventions: the epoch grid restarts at every
file seam (each source file keeps its own original 30-s epochs, so no epoch
ever contains a stitch), the artifact masks are sliced per piece using the
`clip_pieces` bookkeeping saved by the clipper, and `pre`/`post` is derived
from the clip name (`_DX_` = pre-treatment, `_FU_` = post-treatment/follow-up).

**Artifact detection and artifact-info tools:**

* `Get_All_Metrics_10min_artifact.m` — artifact detection over the source
  files (10-minute framework); writes the `<Pre|Post>_artifact_info_EAR.mat`
  files the drivers read. Signals shorter than the filter requirement are
  skipped and leave empty entries; the drivers drop those empties with both
  lists compacted together, so names and masks always stay aligned.
* `Detect_Artifacts.m` — the sample-level artifact detector itself.
* `get_long_short_art_info.m` — per-file summary of long/short artifacts.
* `binary_arts_to_times.m` — converts binary artifact masks to start/stop times.

**Core calculation functions** (called by the drivers, usually not run by
hand): `Rereference_EEG.m`, `Filter_EEG.m` + `Pick_Filter.m`,
`Find_Clean_Indices.m`, `Calc_Epoch_Indices.m`, `Check_EEG_Standardization.m`,
and the per-metric functions listed in the table above.

**Older drivers kept for reference:** `Get_All_Metrics.m`,
`Get_Metrics_clipped_files.m`, `Get_Cross_corr.m` (superseded by the
`NIMBIS_*` versions), plus the `others/` folder.

**`table_making_codes/`** — turns the result mats into long-format csv tables
(one row per epoch/block), ready for the python statistics and plotting
scripts:

* `fable_NIMBIS_make_mainMetrics_table.m` → `<site>_NIMBIS_clipped_mainMetrics_Cz.csv`
* `fable_NIMBIS_make_DFA_table.m` → `<site>_NIMBIS_clipped_DFA_table.csv`
* `fable_NIMBIS_add_pooled_DFA_to_table.m` → `..._DFA_table_withPooled.csv`
  (adds the pooled-clips DFA rows with a `pooled` true/false column)
* `fable_NIMBIS_make_connectivity_table_3min.m` → `..._connectivity_table_3min.csv`
  (one row per 3-minute block)
* `fable_NIMBIS_make_connectivity_table.m` → the older per-clip connectivity
  table (kept as is)

## Order of operations for one site

1. Convert the EDFs to mat files and fix channel order / sampling-frequency
   headers if needed (scripts in `EDF to mat codes` and
   `fable_segment_detection`).
2. Run the artifact detection (`Get_All_Metrics_10min_artifact.m`) on the
   Pre and Post folders → `Pre_artifact_info_EAR.mat`, `Post_artifact_info_EAR.mat`.
3. Detect sleep/wake segments and cut the clips
   (`fable_multiSite_flexible_segments.py`, `fable_clip_flexible_segments.m`
   in `fable_segment_detection`) → `clipped_<SITE>_flexible\sleep` and `\wake`.
4. Run the five `NIMBIS_*` drivers on the clipped folders. Set `site_tag` and
   `site_folder` at the top of each script. **Make sure every file loop reads
   `for f = 1:length(fileNames)`** — a leftover resume index (e.g. `f = 43:`)
   skips files silently, with no error anywhere downstream.
5. Run the `table_making_codes` scripts → the per-site csv tables.
6. Run the python layer (in `mat files\fable_segment_detection`):
   `fable_NIMBIS_patient_median_boxplots.py` (per-patient medians + pooled
   boxplots), then the per-site boxplots, Mann-Whitney workbooks, effect-size
   workbook, forest plots and caterpillar plots, which all read
   `NIMBIS_patient_medians.csv`.

## Requirements

MATLAB (scripts avoid post-R2017b-only functions such as `maxk`), Signal
Processing Toolbox for the filtering, and Python 3 with numpy, pandas, scipy,
matplotlib and openpyxl for the statistics/plotting layer.
