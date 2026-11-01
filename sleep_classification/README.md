# Sleep/wake probability demo

This example uses the provided logistic regression model and saved
MinMax scalers to calculate sleep and wake probabilities from precomputed
EEG metrics. It writes one file: **`epoch_probabilities.csv`**.

## Run

Use Python 3.9. From this folder, run:

```bash
python -m pip install -r requirements.txt
python sleep_wake_demo.py
```

The two provided `.joblib` files must be beside the script. The included
`metrics/example_metrics.csv` lets you try it immediately; it contains
simulated data with artifacts, not participant data or a model-accuracy test.

## Use your own metrics

Replace the example CSV in `metrics/` with your metric CSVs, using one file
per recording and one row per epoch. Each CSV needs:

- `epoch_start_index` and `epoch_stop_index`.
- The 15 numeric, unscaled feature columns listed in `FEATURES` in the script.
  `amplitude` is also accepted instead of `amplitude_range`.

Calculate these features using the same methods and units as the training
pipeline. The script applies the saved scalers; do not scale the inputs first.
You can change the folder and file paths at the top of the script.

## Output

| Column | Meaning |
| --- | --- |
| `file_name` | Input CSV filename |
| `epoch_start_index` | Original epoch start sample index |
| `epoch_stop_index` | Original epoch stop sample index |
| `sleep_probability` | Model probability for class 1 (sleep) |
| `wake_probability` | Model probability for class 0 (wake), calculated as `1 - sleep_probability` |

A probability of 0.5 or higher was interpreted as sleep; lower values
as wake, as for the example showcase in this code. Note that the paper (Mostaghimi et al.,) has adaptive threshold which is not introduced in this example code. The CSV contains
probabilities only.

The provided model uses `liblinear` logistic regression with an L1 penalty.
It was saved with scikit-learn 1.0.2, which is pinned in `requirements.txt`.
Only load trusted joblib files.
