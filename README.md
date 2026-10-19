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


All five drivers share the same conventions: the epoch grid restarts at every
file seam (each source file keeps its own original 30-s epochs, so no epoch
ever contains a stitch), the artifact masks are sliced per piece using the
`clip_pieces` bookkeeping saved by the clipper, and `pre`/`post` is derived
from the clip name (`_DX_` = pre-treatment, `_FU_` = post-treatment/follow-up).

**Artifact detection and artifact-info tools:**

* `Detect_Artifacts.m` — the sample-level artifact detector itself.


**Core calculation functions** (called by the drivers, usually not run by
hand): `Rereference_EEG.m`, `Filter_EEG.m` + `Pick_Filter.m`,
`Find_Clean_Indices.m`, `Calc_Epoch_Indices.m`, `Check_EEG_Standardization.m`,
and the per-metric functions listed in the table above.


## Requirements

MATLAB (scripts avoid post-R2017b-only functions such as `maxk`), Signal
Processing Toolbox for the filtering, and Python 3 with numpy, pandas, scipy,
matplotlib and openpyxl for the statistics/plotting layer.
