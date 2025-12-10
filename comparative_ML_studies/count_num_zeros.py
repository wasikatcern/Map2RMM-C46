import pandas as pd
import numpy as np

# -----------------------------
# Input file paths
# -----------------------------
FULL_RMM_TTBAR = "out/tev13.6pp_pythia8_ttbar_2lep_data10percent.csv.gz"
FULL_RMM_SH    = "out/pythia8_X1500GeV_SH2bbll_data100percent.csv.gz"
FULL_RMM_WZ    = "out/tev13.6pp_pythia8_wzjet_2lep_data10percent.csv.gz"

C46_TTBAR = "out_C46/tt_c46_frob_10.csv"
C46_SH    = "out_C46/sh_1500_c46_frob.csv"
C46_WZ    = "out_C46/wzjet_c46_frob.csv"


# -----------------------------
# Loaders
# -----------------------------
def load_full_rmm(path):
    df = pd.read_csv(path, compression="infer")
    return df.iloc[:, 4:].values  # remove 4 metadata columns


def load_c46(path):
    df = pd.read_csv(path, compression="infer")
    return df.iloc[:, 1:].values  # remove index column


# -----------------------------
# Zero statistics
# -----------------------------
def zero_stats(X):
    """
    Returns:
    - n_events      : number of rows (events)
    - n_features    : number of columns (values per event)
    - total_values  : n_events * n_features
    - total_zeros   : number of zeros in the full array
    - avg_zeros     : average zeros per event
    """
    n_events = X.shape[0]
    n_features = X.shape[1]
    total_values = X.size
    total_zeros  = np.sum(X == 0)
    avg_zeros    = np.mean(np.sum(X == 0, axis=1))

    return n_events, n_features, total_values, total_zeros, avg_zeros


# -----------------------------
# Load datasets
# -----------------------------
datasets = {
    "Full RMM ttbar sample": load_full_rmm(FULL_RMM_TTBAR),
    "Full RMM X→SH sample":  load_full_rmm(FULL_RMM_SH),
    "Full RMM WZJet sample": load_full_rmm(FULL_RMM_WZ),

    "C46 ttbar sample":      load_c46(C46_TTBAR),
    "C46 X→SH sample":       load_c46(C46_SH),
    "C46 WZJet sample":      load_c46(C46_WZ),
}

# -----------------------------
# Print results
# -----------------------------
print("===== Zero Statistics for Full RMM and C46 Representations =====\n")

for name, X in datasets.items():

    n_events, n_features, total_vals, total_zeros, avg_zeros = zero_stats(X)

    print(f"{name}")
    print(f"  Total events             : {n_events:,}")
    print(f"  Values per event         : {n_features:,}")
    print(f"  Total values             : {total_vals:,}")
    print(f"  Total zeros              : {total_zeros:,}")
    print(f"  Avg zeros per event      : {avg_zeros:.2f}")
    print()
