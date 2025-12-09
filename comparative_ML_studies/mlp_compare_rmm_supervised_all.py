#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import roc_auc_score, roc_curve

RANDOM_STATE = 42

# ----------------------------------------------------------------------
# File paths (adjust if needed)
# ----------------------------------------------------------------------

# Full RMM (original, high-dimensional)
FULL_RMM_BKG = "out/tev13.6pp_pythia8_ttbar_2lep_data10percent.csv.gz"
FULL_RMM_SIG = "out/pythia8_X1500GeV_SH2bbll_data100percent.csv.gz"

# RMM-C20 (20D compact)
C20_BKG = "out_C20/ttbar_10percent_C20.csv.gz"
C20_SIG = "out_C20/SH_1500GeV_C20.csv.gz"

# RMM-C45 (45D compact)
C45_BKG = "out_C45/tt_C45.csv"
C45_SIG = "out_C45/SH_1500_C45.csv"

# RMM-C46 (frob) (46D compact)
#C46_BKG = "tt_c46_frob.csv"
C46_BKG = "out_C46/tt_c46_frob_10.csv"
C46_SIG = "out_C46/sh_1500_c46_frob.csv"


# ----------------------------------------------------------------------
# Helper: load each representation
# ----------------------------------------------------------------------

def load_full_rmm(bkg_path, sig_path):
    """
    Load the full RMM representation.
    Assumes the first 4 columns are metadata (Run, Event, Weight, Label, ...)
    and the remaining columns are the RMM entries.
    """
    bkg = pd.read_csv(bkg_path, compression="infer")
    sig = pd.read_csv(sig_path, compression="infer")

    # Skip first 4 metadata columns, use remaining as features
    X_bkg = bkg.iloc[:, 4:].values
    X_sig = sig.iloc[:, 4:].values

    return X_bkg, X_sig


def load_compact_csv(bkg_path, sig_path, skip_cols=1):
    """
    Generic loader for compact representations (C20, C45, C46, etc).
    Assumes the first column is an event index; we skip it and use
    the remaining columns as features.
    """
    bkg = pd.read_csv(bkg_path)
    sig = pd.read_csv(sig_path)

    X_bkg = bkg.iloc[:, skip_cols:].values
    X_sig = sig.iloc[:, skip_cols:].values

    return X_bkg, X_sig


# ----------------------------------------------------------------------
# Helper: train MLP + compute ROC/AUC
# ----------------------------------------------------------------------

def train_mlp_and_roc(X_bkg, X_sig, label_name):
    """
    Build labels (0 = bkg, 1 = sig), split into train/test,
    standardize features, train an MLP, and compute ROC + AUC.
    """
    X = np.vstack([X_bkg, X_sig])
    y = np.concatenate([
        np.zeros(len(X_bkg), dtype=int),
        np.ones(len(X_sig), dtype=int),
    ])

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=RANDOM_STATE, stratify=y
    )

    # Standardize features based on training set only
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    clf = MLPClassifier(
        hidden_layer_sizes=(64, 32),
        activation="relu",
        solver="adam",
        max_iter=100,
        random_state=RANDOM_STATE,
        early_stopping=True,
        n_iter_no_change=5,
        validation_fraction=0.1,
    )

    clf.fit(X_train_scaled, y_train)
    y_score = clf.predict_proba(X_test_scaled)[:, 1]

    auc = roc_auc_score(y_test, y_score)
    fpr, tpr, _ = roc_curve(y_test, y_score)

    print(f"AUC ({label_name}) : {auc:.4f}")
    return fpr, tpr, auc


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    # 1) Load data for each representation
    print("=== Loading FULL RMM ===")
    X_full_bkg, X_full_sig = load_full_rmm(FULL_RMM_BKG, FULL_RMM_SIG)
    print("  Full RMM shapes:", X_full_bkg.shape, X_full_sig.shape)

    print("\n=== Loading RMM-C20 ===")
    X_c20_bkg, X_c20_sig = load_compact_csv(C20_BKG, C20_SIG, skip_cols=1)
    print("  C20 shapes:", X_c20_bkg.shape, X_c20_sig.shape)

    print("\n=== Loading RMM-C45 ===")
    X_c45_bkg, X_c45_sig = load_compact_csv(C45_BKG, C45_SIG, skip_cols=1)
    print("  C45 shapes:", X_c45_bkg.shape, X_c45_sig.shape)

    print("\n=== Loading RMM-C46 (frob) ===")
    X_c46_bkg, X_c46_sig = load_compact_csv(C46_BKG, C46_SIG, skip_cols=1)
    print("  C46 shapes:", X_c46_bkg.shape, X_c46_sig.shape)

    # 2) Train + evaluate each representation
    print("\n========== TRAINING MLPs AND COMPUTING AUC ==========\n")

    fpr_full, tpr_full, auc_full = train_mlp_and_roc(
        X_full_bkg, X_full_sig, "Full RMM"
    )
    fpr_c20, tpr_c20, auc_c20 = train_mlp_and_roc(
        X_c20_bkg, X_c20_sig, "RMM-C20"
    )
    fpr_c45, tpr_c45, auc_c45 = train_mlp_and_roc(
        X_c45_bkg, X_c45_sig, "RMM-C45"
    )
    fpr_c46, tpr_c46, auc_c46 = train_mlp_and_roc(
        X_c46_bkg, X_c46_sig, "RMM-C46 (frob)"
    )

    print("\n========== SUMMARY ==========")
    print(f"AUC (Full RMM)      : {auc_full:.4f}")
    print(f"AUC (RMM-C20)       : {auc_c20:.4f}")
    print(f"AUC (RMM-C45)       : {auc_c45:.4f}")
    print(f"AUC (RMM-C46, frob) : {auc_c46:.4f}")

    # 3) Plot ROC curves together
    plt.figure(figsize=(6.5, 5.5))
    plt.plot(fpr_full, tpr_full, label=f"Full RMM (AUC={auc_full:.4f})")
    plt.plot(fpr_c20,  tpr_c20,  label=f"RMM-C20 (AUC={auc_c20:.4f})")
    plt.plot(fpr_c45,  tpr_c45,  label=f"RMM-C45 (AUC={auc_c45:.4f})")
    plt.plot(fpr_c46,  tpr_c46,  label=f"RMM-C46 (frob) (AUC={auc_c46:.4f})")

    plt.plot([0, 1], [0, 1], "k--", alpha=0.5)
    plt.xlabel("False positive rate")
    plt.ylabel("True positive rate")
    plt.title("MLP ROC comparison: Full RMM vs RMM-C20/C45/C46")
    plt.legend(loc="lower right", fontsize=9)
    plt.tight_layout()
    plt.savefig("out_C46/roc_compare_rmm_all.pdf")
    plt.close()

    print("\nSaved ROC plot to: out_C46/roc_compare_rmm_all.pdf")


if __name__ == "__main__":
    main()
