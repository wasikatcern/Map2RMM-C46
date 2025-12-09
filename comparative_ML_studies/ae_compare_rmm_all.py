#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, roc_curve

import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset

RANDOM_STATE = 42
np.random.seed(RANDOM_STATE)
torch.manual_seed(RANDOM_STATE)

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
C46_BKG = "out_C46/tt_c46_frob_10.csv"
C46_SIG = "out_C46/sh_1500_c46_frob.csv"

# ----------------------------------------------------------------------
# Additional samples for multi-loss plot (second figure)
# ----------------------------------------------------------------------

# Full RMM
FULL_EXTRA_SAMPLES = {
    "ttbar (1%)": "out/tev13.6pp_pythia8_ttbar_2lep_data1percent.csv.gz",
    "WZJets": "out/tev13.6pp_pythia8_wzjet_2lep_data10percent.csv.gz",

    "HH 500 GeV":  "out/pythia8_X500GeV_HH2bbll_data100percent.csv.gz",
    "HH 700 GeV":  "out/pythia8_X700GeV_HH2bbll_data100percent.csv.gz",
    "HH 1000 GeV": "out/pythia8_X1000GeV_HH2bbll_data100percent.csv.gz",
    "HH 1500 GeV": "out/pythia8_X1500GeV_HH2bbll_data100percent.csv.gz",
    "HH 2000 GeV": "out/pythia8_X2000GeV_HH2bbll_data100percent.csv.gz",

    "SH 500 GeV":  "out/pythia8_X500GeV_SH2bbll_data100percent.csv.gz",
    "SH 700 GeV":  "out/pythia8_X700GeV_SH2bbll_data100percent.csv.gz",
    "SH 1000 GeV": "out/pythia8_X1000GeV_SH2bbll_data100percent.csv.gz",
    "SH 1500 GeV": "out/pythia8_X1500GeV_SH2bbll_data100percent.csv.gz",
    "SH 2000 GeV": "out/pythia8_X2000GeV_SH2bbll_data100percent.csv.gz",
}

# C20
C20_EXTRA_SAMPLES = {
    "ttbar (1%)": "out_C20/ttbar_1percent_C20.csv.gz",
    "WZJets": "out_C20/wjet_C20.csv.gz",

    "HH 500 GeV":  "out_C20/hh_500GeV_C20.csv.gz",
    "HH 700 GeV":  "out_C20/hh_700GeV_C20.csv.gz",
    "HH 1000 GeV": "out_C20/hh_1000GeV_C20.csv.gz",
    "HH 1500 GeV": "out_C20/hh_1500GeV_C20.csv.gz",
    "HH 2000 GeV": "out_C20/hh_2000GeV_C20.csv.gz",

    "SH 500 GeV":  "out_C20/SH_500GeV_C20.csv.gz",
    "SH 700 GeV":  "out_C20/SH_700GeV_C20.csv.gz",
    "SH 1000 GeV": "out_C20/SH_1000GeV_C20.csv.gz",
    "SH 1500 GeV": "out_C20/SH_1500GeV_C20.csv.gz",
    "SH 2000 GeV": "out_C20/SH_2000GeV_C20.csv.gz",
}

# C45
C45_EXTRA_SAMPLES = {
    "ttbar (1%)": "out_C45/tt_1p_C45.csv",
    "WZJets": "out_C45/wzjet_C45.csv",

    "HH 500 GeV":  "out_C45/HH_500_C45.csv",
    "HH 700 GeV":  "out_C45/HH_700_C45.csv",
    "HH 1000 GeV": "out_C45/HH_1000_C45.csv",
    "HH 1500 GeV": "out_C45/HH_1500_C45.csv",
    "HH 2000 GeV": "out_C45/HH_2000_C45.csv",

    "SH 500 GeV":  "out_C45/SH_500_C45.csv",
    "SH 700 GeV":  "out_C45/SH_700_C45.csv",
    "SH 1000 GeV": "out_C45/SH_1000_C45.csv",
    "SH 1500 GeV": "out_C45/SH_1500_C45.csv",
    "SH 2000 GeV": "out_C45/SH_2000_C45.csv",
}

# C46 (frob)
C46_EXTRA_SAMPLES = {
    "ttbar (1%)": "out_C46/tt_c46_frob.csv",
    "WZJets": "out_C46/wzjet_c46_frob.csv",

    "HH 500 GeV":  "out_C46/hh_500_c46_frob.csv",
    "HH 700 GeV":  "out_C46/hh_700_c46_frob.csv",
    "HH 1000 GeV": "out_C46/hh_1000_c46_frob.csv",
    "HH 1500 GeV": "out_C46/hh_1500_c46_frob.csv",
    "HH 2000 GeV": "out_C46/hh_2000_c46_frob.csv",

    "SH 500 GeV":  "out_C46/sh_500_c46_frob.csv",
    "SH 700 GeV":  "out_C46/sh_700_c46_frob.csv",
    "SH 1000 GeV": "out_C46/sh_1000_c46_frob.csv",
    "SH 1500 GeV": "out_C46/sh_1500_c46_frob.csv",
    "SH 2000 GeV": "out_C46/sh_2000_c46_frob.csv",
}

# ----------------------------------------------------------------------
# Per-sample event limits for the multi-loss plot
# ----------------------------------------------------------------------
# If value is None -> use all events
# If value is an int -> use only the first N events from that sample
EVENT_LIMITS = {
    "WZJets": 44000,
    # You can add more if you want:
    # "HH 500 GeV": 10000,
    # "SH 1500 GeV": 5000,
}

# ----------------------------------------------------------------------
# Helper: load each representation
# ----------------------------------------------------------------------

def load_full_rmm_single(path):
    """
    Load one full RMM CSV/CSV.GZ.
    Assumes the first 4 columns are metadata and the remaining are features.
    """
    df = pd.read_csv(path, compression="infer")
    X = df.iloc[:, 4:].values
    return X

def load_full_rmm_pair(bkg_path, sig_path):
    X_bkg = load_full_rmm_single(bkg_path)
    X_sig = load_full_rmm_single(sig_path)
    return X_bkg, X_sig

def load_compact_single(path, skip_cols=1):
    """
    Generic loader for compact representations (C20, C45, C46, etc).
    Assumes the first column is an event index; we skip it.
    """
    df = pd.read_csv(path, compression="infer")
    X = df.iloc[:, skip_cols:].values
    return X

def load_compact_pair(bkg_path, sig_path, skip_cols=1):
    X_bkg = load_compact_single(bkg_path, skip_cols=skip_cols)
    X_sig = load_compact_single(sig_path, skip_cols=skip_cols)
    return X_bkg, X_sig

# ----------------------------------------------------------------------
# Simple PyTorch Autoencoder
# ----------------------------------------------------------------------

class SimpleAE(nn.Module):
    def __init__(self, input_dim, bottleneck_dim=16):
        super(SimpleAE, self).__init__()
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 64),
            nn.ReLU(),
            nn.Linear(64, 32),
            nn.ReLU(),
            nn.Linear(32, bottleneck_dim),
            nn.ReLU(),
        )
        self.decoder = nn.Sequential(
            nn.Linear(bottleneck_dim, 32),
            nn.ReLU(),
            nn.Linear(32, 64),
            nn.ReLU(),
            nn.Linear(64, input_dim),
        )

    def forward(self, x):
        z = self.encoder(x)
        out = self.decoder(z)
        return out

def train_autoencoder(X_train, X_val, input_dim,
                      bottleneck_dim=16,
                      lr=1e-3,
                      batch_size=512,
                      num_epochs=50,
                      device="cpu"):
    model = SimpleAE(input_dim, bottleneck_dim=bottleneck_dim).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    criterion = nn.MSELoss(reduction="mean")

    train_dataset = TensorDataset(torch.from_numpy(X_train.astype(np.float32)))
    val_dataset   = TensorDataset(torch.from_numpy(X_val.astype(np.float32)))

    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader   = DataLoader(val_dataset,   batch_size=batch_size, shuffle=False)

    best_val_loss = np.inf
    best_state = None
    patience = 5
    wait = 0

    for epoch in range(num_epochs):
        model.train()
        train_loss = 0.0
        n_train    = 0

        for (xb,) in train_loader:
            xb = xb.to(device)
            optimizer.zero_grad()
            recon = model(xb)
            loss = criterion(recon, xb)
            loss.backward()
            optimizer.step()

            train_loss += loss.item() * xb.size(0)
            n_train    += xb.size(0)

        train_loss /= max(n_train, 1)

        model.eval()
        val_loss = 0.0
        n_val    = 0
        with torch.no_grad():
            for (xb,) in val_loader:
                xb = xb.to(device)
                recon = model(xb)
                loss = criterion(recon, xb)
                val_loss += loss.item() * xb.size(0)
                n_val    += xb.size(0)
        val_loss /= max(n_val, 1)

        print(f"  Epoch {epoch+1:03d}/{num_epochs} - train_loss={train_loss:.4e}, val_loss={val_loss:.4e}")

        if val_loss < best_val_loss - 1e-6:
            best_val_loss = val_loss
            best_state = model.state_dict()
            wait = 0
        else:
            wait += 1
            if wait >= patience:
                print("  Early stopping.")
                break

    if best_state is not None:
        model.load_state_dict(best_state)
    return model

def compute_reco_losses(model, X, device="cpu"):
    model.eval()
    dataset = TensorDataset(torch.from_numpy(X.astype(np.float32)))
    loader  = DataLoader(dataset, batch_size=1024, shuffle=False)

    losses = []
    criterion = nn.MSELoss(reduction="none")

    with torch.no_grad():
        for (xb,) in loader:
            xb = xb.to(device)
            recon = model(xb)
            per_elem = criterion(recon, xb)
            per_event = per_elem.view(per_elem.size(0), -1).mean(dim=1)
            losses.append(per_event.cpu().numpy())

    return np.concatenate(losses, axis=0)

# ----------------------------------------------------------------------
# Train AE and compute ROC for one representation
# ----------------------------------------------------------------------

def ae_train_and_roc(X_bkg, X_sig, label_name, device="cpu"):
    """
    Unsupervised AE: train only on background (ttbar), then evaluate on
    (held-out bkg + signal) to build ROC.
    """
    # Split background into train/val/test
    X_bkg_train, X_bkg_tmp = train_test_split(
        X_bkg, test_size=0.4, random_state=RANDOM_STATE
    )
    X_bkg_val, X_bkg_test = train_test_split(
        X_bkg_tmp, test_size=0.5, random_state=RANDOM_STATE
    )

    # Scale using train background only
    scaler = StandardScaler()
    X_bkg_train_scaled = scaler.fit_transform(X_bkg_train)
    X_bkg_val_scaled   = scaler.transform(X_bkg_val)
    X_bkg_test_scaled  = scaler.transform(X_bkg_test)
    X_sig_scaled       = scaler.transform(X_sig)

    input_dim = X_bkg_train_scaled.shape[1]
    print(f"\n[{label_name}] Training AE with input_dim={input_dim}")

    model = train_autoencoder(
        X_bkg_train_scaled,
        X_bkg_val_scaled,
        input_dim=input_dim,
        bottleneck_dim=16,
        lr=1e-3,
        batch_size=512,
        num_epochs=50,
        device=device,
    )

    # Reconstruction losses
    loss_bkg = compute_reco_losses(model, X_bkg_test_scaled, device=device)
    loss_sig = compute_reco_losses(model, X_sig_scaled, device=device)

    # Higher loss = more anomalous (signal-like)
    y_true  = np.concatenate([np.zeros_like(loss_bkg), np.ones_like(loss_sig)])
    scores  = np.concatenate([loss_bkg, loss_sig])

    auc = roc_auc_score(y_true, scores)
    fpr, tpr, _ = roc_curve(y_true, scores)

    print(f"[{label_name}] ROC AUC: {auc:.4f}")

    return {
        "label": label_name,
        "model": model,
        "scaler": scaler,
        "fpr": fpr,
        "tpr": tpr,
        "auc": auc,
        "loss_bkg_test": loss_bkg,
        "loss_sig_test": loss_sig,
    }

# ----------------------------------------------------------------------
# Multi-loss distributions for many samples (second figure)
# ----------------------------------------------------------------------

def compute_multi_losses(name, model, scaler, loader_single, samples_paths, device="cpu"):
    """
    For a given representation and trained AE+scaler, compute reconstruction
    losses for many different samples. Uses EVENT_LIMITS to optionally
    restrict the number of events per sample.
    """
    samples_losses = {}
    print(f"\n[{name}] Computing multi-sample losses...")

    for lab, path in samples_paths.items():
        if not os.path.exists(path):
            print(f"  [WARN] {name}: file not found: {path} (skipping {lab})")
            continue

        X_raw = loader_single(path)

        # Apply per-sample event limits, if defined
        limit = EVENT_LIMITS.get(lab, None)
        if limit is not None and X_raw.shape[0] > limit:
            X_raw = X_raw[:limit]

        X_scaled = scaler.transform(X_raw)
        losses = compute_reco_losses(model, X_scaled, device=device)
        samples_losses[lab] = losses

        print(f"  {lab}: loaded {X_raw.shape[0]} events, computed losses")

    return samples_losses


def plot_multi_loss_distributions(all_losses_dict, out_pdf):
    """
    all_losses_dict: {
        "Full RMM": {"ttbar (1%)": losses, "WZJets": losses, ...},
        "RMM-C20": {...},
        ...
    }
    """
    plt.figure(figsize=(10, 8))

    # 2x2 grid assuming 4 representations
    grid_positions = {
        "Full RMM": 221,
        "RMM-C20": 222,
        "RMM-C45": 223,
        "RMM-C46": 224,
    }

    colors = [
        "C0", "C1", "C2", "C3", "C4", "C5", "C6",
        "C7", "C8", "C9", "tab:gray", "tab:brown"
    ]

    for rep_name, losses_dict in all_losses_dict.items():
        if rep_name not in grid_positions:
            continue
        pos = grid_positions[rep_name]
        ax = plt.subplot(pos)

        for i, (lab, loss_arr) in enumerate(losses_dict.items()):
            if len(loss_arr) == 0:
                continue
            c = colors[i % len(colors)]
            log_loss = np.log(loss_arr + 1e-12)
            ax.hist(
                log_loss,
                bins=60,
                histtype="step",
                linewidth=1.5,
                density=True,
                label=lab,
                color=c,
                alpha=0.9,
            )

        ax.set_xlabel("log(Reconstruction loss)")
        ax.set_ylabel("Density")
        ax.set_title(rep_name)
        ax.legend(fontsize=7, loc="upper right")

    plt.tight_layout()
    plt.savefig(out_pdf)
    plt.close()
    print(f"Saved multi-loss distributions to: {out_pdf}")



# ----------------------------------------------------------------------
# NEW: AUC vs mass for HH/SH for each representation
# ----------------------------------------------------------------------

def compute_and_plot_auc_vs_mass(all_losses_dict, out_pdf):
    """
    Using the already computed reconstruction losses in all_losses_dict,
    compute AUC(ttbar (1%) vs HH mX) and AUC(ttbar (1%) vs SH mX)
    for mX = 500, 700, 1000, 1500, 2000 GeV,
    for each representation (Full RMM, C20, C45, C46),
    and plot AUC vs mass for HH and SH in a single figure with two subplots.
    """
    mass_points = [500, 700, 1000, 1500, 2000]
    #reps = ["Full RMM", "RMM-C20", "RMM-C45", "RMM-C46 (frob)"]
    reps = ["Full RMM", "RMM-C45", "RMM-C46"]

    # Store AUCs
    hh_aucs = {rep: [] for rep in reps}
    hh_ms   = {rep: [] for rep in reps}
    sh_aucs = {rep: [] for rep in reps}
    sh_ms   = {rep: [] for rep in reps}

    print("\n========== AUC vs mass (HH, SH) using ttbar (1%) as background ==========\n")

    for rep in reps:
        if rep not in all_losses_dict:
            print(f"[WARN] Representation {rep} not found in all_losses_dict, skipping.")
            continue

        rep_losses = all_losses_dict[rep]
        if "ttbar (1%)" not in rep_losses:
            print(f"[WARN] 'ttbar (1%)' not present for {rep}, skipping AUC for this rep.")
            continue

        bkg_losses = rep_losses["ttbar (1%)"]

        print(f"\n=== {rep} ===")
        for m in mass_points:
            # HH
            hh_key = f"HH {m} GeV"
            sh_key = f"SH {m} GeV"

            if hh_key in rep_losses:
                sig_hh = rep_losses[hh_key]
                y_true = np.concatenate([
                    np.zeros_like(bkg_losses),
                    np.ones_like(sig_hh),
                ])
                scores = np.concatenate([bkg_losses, sig_hh])
                auc_hh = roc_auc_score(y_true, scores)
                hh_aucs[rep].append(auc_hh)
                hh_ms[rep].append(m)
                print(f"  AUC(ttbar vs {hh_key}) = {auc_hh:.4f}")
            else:
                print(f"  [WARN] Missing sample {hh_key} for {rep}")

            # SH
            if sh_key in rep_losses:
                sig_sh = rep_losses[sh_key]
                y_true = np.concatenate([
                    np.zeros_like(bkg_losses),
                    np.ones_like(sig_sh),
                ])
                scores = np.concatenate([bkg_losses, sig_sh])
                auc_sh = roc_auc_score(y_true, scores)
                sh_aucs[rep].append(auc_sh)
                sh_ms[rep].append(m)
                print(f"  AUC(ttbar vs {sh_key}) = {auc_sh:.4f}")
            else:
                print(f"  [WARN] Missing sample {sh_key} for {rep}")

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    colors = {
        "Full RMM": "C0",
        #"RMM-C20": "C1",
        "RMM-C45": "C2",
        "RMM-C46": "C3",
    }

    # HH subplot
    ax_hh = axes[0]
    for rep in reps:
        if len(hh_ms[rep]) == 0:
            continue
        ax_hh.plot(
            hh_ms[rep],
            hh_aucs[rep],
            marker="o",
            label=rep,
            color=colors.get(rep, None),
        )
    ax_hh.set_xlabel(r"$m_X$ [GeV]")
    ax_hh.set_ylabel("AUC")
    ax_hh.set_title(r"X $\to$ HH vs ttbar")
    ax_hh.set_xticks(mass_points)
    ax_hh.grid(True, alpha=0.3)

    # SH subplot
    ax_sh = axes[1]
    for rep in reps:
        if len(sh_ms[rep]) == 0:
            continue
        ax_sh.plot(
            sh_ms[rep],
            sh_aucs[rep],
            marker="o",
            label=rep,
            color=colors.get(rep, None),
        )
    ax_sh.set_xlabel(r"$m_X$ [GeV]")
    ax_sh.set_title(r"X $\to$ SH vs ttbar")
    ax_sh.set_xticks(mass_points)
    ax_sh.grid(True, alpha=0.3)

    # Put legend only once
    handles, labels = ax_hh.get_legend_handles_labels()
    if not handles:
        handles, labels = ax_sh.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=8)

    plt.tight_layout(rect=[0, 0.08, 1, 1])
    os.makedirs(os.path.dirname(out_pdf), exist_ok=True)
    plt.savefig(out_pdf)
    plt.close()
    print(f"\nSaved AUC vs mass figure to: {out_pdf}")

# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Using device: {device}")

    # --------------------------------------------------
    # 1) Load base train/test sets for ROC comparison
    # --------------------------------------------------
    print("=== Loading FULL RMM (train/test) ===")
    X_full_bkg, X_full_sig = load_full_rmm_pair(FULL_RMM_BKG, FULL_RMM_SIG)
    print("  Full RMM shapes:", X_full_bkg.shape, X_full_sig.shape)

    print("\n=== Loading RMM-C20 (train/test) ===")
    X_c20_bkg, X_c20_sig = load_compact_pair(C20_BKG, C20_SIG, skip_cols=1)
    print("  C20 shapes:", X_c20_bkg.shape, X_c20_sig.shape)

    print("\n=== Loading RMM-C45 (train/test) ===")
    X_c45_bkg, X_c45_sig = load_compact_pair(C45_BKG, C45_SIG, skip_cols=1)
    print("  C45 shapes:", X_c45_bkg.shape, X_c45_sig.shape)

    print("\n=== Loading RMM-C46 (frob, train/test) ===")
    X_c46_bkg, X_c46_sig = load_compact_pair(C46_BKG, C46_SIG, skip_cols=1)
    print("  C46 shapes:", X_c46_bkg.shape, X_c46_sig.shape)

    # --------------------------------------------------
    # 2) Train AEs and compute ROC for each representation
    # --------------------------------------------------
    print("\n========== TRAINING AEs AND COMPUTING AUC ==========\n")

    results = []

    res_full = ae_train_and_roc(X_full_bkg, X_full_sig, "Full RMM", device=device)
    results.append(res_full)

    res_c20 = ae_train_and_roc(X_c20_bkg, X_c20_sig, "RMM-C20", device=device)
    results.append(res_c20)

    res_c45 = ae_train_and_roc(X_c45_bkg, X_c45_sig, "RMM-C45", device=device)
    results.append(res_c45)

    res_c46 = ae_train_and_roc(X_c46_bkg, X_c46_sig, "RMM-C46", device=device)
    results.append(res_c46)

    # --------------------------------------------------
    # 3) ROC plot comparison
    # --------------------------------------------------
    print("\n========== SUMMARY (ROC AUC) ==========")
    for res in results:
        print(f"AUC ({res['label']}) : {res['auc']:.4f}")

    os.makedirs("out_C46", exist_ok=True)

    plt.figure(figsize=(6.5, 5.5))
    for res in results:
        plt.plot(
            res["fpr"],
            res["tpr"],
            label=f"{res['label']} (AUC={res['auc']:.4f})"
        )

    plt.plot([0, 1], [0, 1], "k--", alpha=0.5)
    plt.xlabel("False positive rate")
    plt.ylabel("True positive rate")
    #plt.title("Autoencoder ROC: Full RMM vs RMM-C20/C45/C46")
    plt.title("Autoencoder ROC: Full RMM vs RMM-C45/C46")
    plt.legend(loc="lower right", fontsize=9)
    plt.tight_layout()
    plt.savefig("out_C46/roc_compare_rmm_all_AE.pdf")
    plt.close()
    print("Saved ROC plot to: out_C46/roc_compare_rmm_all_AE.pdf")

    # --------------------------------------------------
    # 4) Multi-loss distribution plots for many samples
    # --------------------------------------------------
    print("\n========== MULTI-SAMPLE LOSS DISTRIBUTIONS ==========")

    all_losses_dict = {}

    # Full RMM
    full_losses = compute_multi_losses(
        "Full RMM",
        res_full["model"],
        res_full["scaler"],
        loader_single=load_full_rmm_single,
        samples_paths=FULL_EXTRA_SAMPLES,
        device=device,
    )
    all_losses_dict["Full RMM"] = full_losses
    
    
    # C20
    c20_losses = compute_multi_losses(
        "RMM-C20",
        res_c20["model"],
        res_c20["scaler"],
        loader_single=lambda p: load_compact_single(p, skip_cols=1),
        samples_paths=C20_EXTRA_SAMPLES,
        device=device,
    )
    all_losses_dict["RMM-C20"] = c20_losses
    

    # C45
    c45_losses = compute_multi_losses(
        "RMM-C45",
        res_c45["model"],
        res_c45["scaler"],
        loader_single=lambda p: load_compact_single(p, skip_cols=1),
        samples_paths=C45_EXTRA_SAMPLES,
        device=device,
    )
    all_losses_dict["RMM-C45"] = c45_losses

    # C46
    c46_losses = compute_multi_losses(
        "RMM-C46",
        res_c46["model"],
        res_c46["scaler"],
        loader_single=lambda p: load_compact_single(p, skip_cols=1),
        samples_paths=C46_EXTRA_SAMPLES,
        device=device,
    )
    all_losses_dict["RMM-C46"] = c46_losses

    plot_multi_loss_distributions(
        all_losses_dict,
        out_pdf="out_C46/reco_lossx_distributions_rmm_all_AE_multi.pdf"
    )

    # --------------------------------------------------
    # 5) NEW: AUC vs mass for HH/SH for each representation
    # --------------------------------------------------
    compute_and_plot_auc_vs_mass(
        all_losses_dict,
        out_pdf="out_C46/auc_vs_mass_all_AE.pdf"
    )

if __name__ == "__main__":
    main()
