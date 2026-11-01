"""Use the saved sleep/wake model to create epoch_probabilities.csv."""

from pathlib import Path

import joblib
import pandas as pd


# Put your precomputed metric CSVs in the metrics folder (one recording per CSV).
FOLDER = Path(__file__).resolve().parent
METRICS_FOLDER = FOLDER / "metrics"
MODEL_FILE = FOLDER / "lr_liblinear_sleep_model.joblib"
SCALER_FILE = FOLDER / "sleep_model_feature_scaler.joblib"
OUTPUT_FILE = FOLDER / "epoch_probabilities.csv"

# Keep the same feature order used when training the model.
FEATURES = [
    "permutation_entropy_delta", "shannon_entropy_delta",
    "permutation_entropy_theta", "shannon_entropy_theta",
    "permutation_entropy_alpha", "shannon_entropy_alpha",
    "permutation_entropy_beta", "shannon_entropy_beta",
    "power_spectral_density_delta", "power_spectral_density_theta",
    "power_spectral_density_alpha", "power_spectral_density_beta",
    "power_spectral_density_broadband", "amplitude_range",
    "spectral_edge_frequency",
]


def main():
    # 1. Load the trained model and its saved feature scalers.
    model = joblib.load(MODEL_FILE)
    scalers = joblib.load(SCALER_FILE)
    sleep_column = list(model.classes_).index(1)  # Class 1 = sleep; class 0 = wake.
    results = []

    # 2. Read the precomputed metrics for each recording.
    for csv_file in sorted(METRICS_FOLDER.glob("*.csv")):
        data = pd.read_csv(csv_file)
        if data.empty:
            continue
        if "amplitude" in data.columns and "amplitude_range" not in data.columns:
            data = data.rename(columns={"amplitude": "amplitude_range"})

        # 3. Apply the training scalers without fitting them again.
        scaled = data[FEATURES].copy()
        for feature in FEATURES:
            scaled[feature] = scalers[feature].transform(data[[feature]]).ravel()

        # 4. Predict the sleep probability for every epoch.
        sleep_probability = model.predict_proba(scaled.to_numpy())[:, sleep_column]
        result = data[["epoch_start_index", "epoch_stop_index"]].copy()
        result.insert(0, "file_name", csv_file.name)
        result["sleep_probability"] = sleep_probability
        result["wake_probability"] = 1 - sleep_probability
        results.append(result)

    # 5. Combine all recordings into one probability CSV.
    if not results:
        raise ValueError(f"No non-empty metric CSVs found in {METRICS_FOLDER}")
    output = pd.concat(results, ignore_index=True)
    output.to_csv(OUTPUT_FILE, index=False)
    print(f"Saved {len(output)} epochs to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
