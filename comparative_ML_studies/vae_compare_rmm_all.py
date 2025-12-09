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
EVENT_LIMITS = {
    "WZJets": 44000,
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
    Generic loader for compact representations.
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
# Variational Autoencoder
# ----------------------------------------------------------------------

class VariationalAE(nn.Module):
    def __init__(self, input_dim, latent_dim=16):
        super(VariationalAE, self).__init__()
        # Encoder
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 64),
            nn.ReLU(),
            nn.Linear(64, 32),
            nn.ReLU(),
        )
        self.fc_mu     = nn.Linear(32, latent_dim)
        self.fc_logvar = nn.Linear(32, latent_dim)

        # Decoder
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, 32),
            nn.ReLU(),
            nn.Linear(32, 64),
            nn.ReLU(),
            nn.Linear(64, input_dim),
        )

    def encode(self, x):
        h = self.encoder(x)
        mu     = self.fc_mu(h)
        logvar = self.fc_logvar(h)
        # numeric safety: prevent exp(logvar) overflow
        logvar = torch.clamp(logvar, min=-10.0, max=10.0)
        return mu, logvar

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def decode(self, z):
        return self.decoder(z)

    def forward(self, x):
        """
        Returns:
          recon_x: reconstructed input
          mu, logvar: latent parameters for KL term
        """
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        recon_x = self.decode(z)
        return recon_x, mu, logvar

def vae_loss_function(recon_x, x, mu, logvar, kl_weight=1.0):
    """
    Standard VAE loss: reconstruction + KL.
    Inputs are already standardized (Gaussian-ish), so we use MSE for recon.
    """
    recon_loss = nn.functional.mse_loss(recon_x, x, reduction="mean")
    kl_loss = -0.5 * torch.mean(1 + logvar - mu.pow(2) - logvar.exp())
    loss = recon_loss + kl_weight * kl_loss
    return loss, recon_loss, kl_loss

def train_vae(
    X_train,
    X_val,
    input_dim,
    latent_dim=16,
    lr=5e-4,          # smaller LR for stability
    batch_size=512,
    num_epochs=50,
    kl_weight=1.0,
    device="cpu",
):
    model = VariationalAE(input_dim, latent_dim=latent_dim).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)

    train_dataset = TensorDataset(torch.from_numpy(X_train.astype(np.float32)))
    val_dataset   = TensorDataset(torch.from_numpy(X_val.astype(np.float32)))

    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader   = DataLoader(val_dataset,   batch_size=batch_size, shuffle=False)

    best_val_loss = np.inf
    best_state = None
    patience = 5
    wait = 0

    for epoch in range(num_epochs):
        # KL warm-up: ramp from 0 → kl_weight over first 10 epochs
        warmup_factor = min(1.0, (epoch + 1) / 10.0)
        current_kl_weight = kl_weight * warmup_factor

        # -------------------------
        # Training
        # -------------------------
        model.train()
        train_loss = 0.0
        n_train = 0

        for (xb,) in train_loader:
            xb = xb.to(device)
            optimizer.zero_grad()

            recon, mu, logvar = model(xb)
            loss, recon_l, kl_l = vae_loss_function(
                recon, xb, mu, logvar, kl_weight=current_kl_weight
            )
            loss.backward()

            # gradient clipping to avoid explosion
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)

            optimizer.step()

            train_loss += loss.item() * xb.size(0)
            n_train    += xb.size(0)

        train_loss /= max(n_train, 1)

        # -------------------------
        # Validation
        # -------------------------
        model.eval()
        val_loss = 0.0
        n_val    = 0
        with torch.no_grad():
            for (xb,) in val_loader:
                xb = xb.to(device)
                recon, mu, logvar = model(xb)
                loss, recon_l, kl_l = vae_loss_function(
                    recon, xb, mu, logvar, kl_weight=current_kl_weight
                )
                val_loss += loss.item() * xb.size(0)
                n_val    += xb.size(0)

        val_loss /= max(n_val, 1)

        print(
            f"  Epoch {epoch+1:03d}/{num_epochs} "
            f"- train_loss={train_loss:.4e}, val_loss={val_loss:.4e}, "
            f"kl_w={current_kl_weight:.3f}"
        )

        # Early stopping on total loss (recon + KL)
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
    """
    Compute per-event reconstruction loss for a trained VAE.
    We use the MSE part only (no KL) as anomaly score.
    NaN/Inf are replaced by a large finite number.
    """
    model.eval()
    dataset = TensorDataset(torch.from_numpy(X.astype(np.float32)))
    loader  = DataLoader(dataset, batch_size=1024, shuffle=False)

    losses = []
    criterion = nn.MSELoss(reduction="none")

    with torch.no_grad():
        for (xb,) in loader:
            xb = xb.to(device)
            recon, _, _ = model(xb)
            per_elem = criterion(recon, xb)  # shape (batch, input_dim)
            per_event = per_elem.view(per_elem.size(0), -1).mean(dim=1)

            # Safety: replace any NaN/inf with a large finite value
            per_event = torch.nan_to_num(
                per_event,
                nan=1e6,
                posinf=1e6,
                neginf=1e6,
            )

            losses.append(per_event.cpu().numpy())

    return np.concatenate(losses, axis=0)

# ----------------------------------------------------------------------
# Train VAE and compute ROC for one representation
# ----------------------------------------------------------------------

def vae_train_and_roc(X_bkg, X_sig, label_name, device="cpu"):
    """
    Unsupervised VAE: train only on background (ttbar), then evaluate on
    (held-out bkg + signal) to build ROC using reconstruction loss.
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
    print(f"\n[{label_name}] Training VAE with input_dim={input_dim}")

    # Slightly smaller KL weight for very high-dimensional input
    kl_w = 0.2 if input_dim > 100 else 1.0

    model = train_vae(
        X_bkg_train_scaled,
        X_bkg_val_scaled,
        input_dim=input_dim,
        latent_dim=16,
        lr=5e-4,
        batch_size=512,
        num_epochs=50,
        kl_weight=kl_w,
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

    print(f"[{label_name}] ROC AUC (VAE): {auc:.3f}")

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
    For a given representation and trained VAE+scaler, compute reconstruction
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
        "RMM-C46": {...},
    }
    Only "Full RMM" and "RMM-C46" are plotted here.
    """
    plt.figure(figsize=(10, 4))

    # 1x2 grid: Full RMM (left), RMM-C46 (right)
    grid_positions = {
        "Full RMM": 121,
        "RMM-C46": 122,
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
# AUC vs mass for HH/SH for each representation
# ----------------------------------------------------------------------

def compute_and_plot_auc_vs_mass(all_losses_dict, out_pdf):
    """
    Using the already computed reconstruction losses in all_losses_dict,
    compute AUC(ttbar (1%) vs HH mX) and AUC(ttbar (1%) vs SH mX)
    for mX = 500, 700, 1000, 1500, 2000 GeV,
    for each representation (Full RMM, C46),
    and plot AUC vs mass for HH and SH in a single figure with two subplots.
    """
    mass_points = [500, 700, 1000, 1500, 2000]
    reps = ["Full RMM", "RMM-C46"]

    # Store AUCs
    hh_aucs = {rep: [] for rep in reps}
    hh_ms   = {rep: [] for rep in reps}
    sh_aucs = {rep: [] for rep in reps}
    sh_ms   = {rep: [] for rep in reps}

    print("\n========== AUC vs mass (HH, SH) using ttbar (1%) as background (VAE) ==========\n")

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
            hh_key = f"HH {m} GeV"
            sh_key = f"SH {m} GeV"

            # HH
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
                print(f"  AUC(ttbar vs {hh_key}) = {auc_hh:.3f}")
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
                print(f"  AUC(ttbar vs {sh_key}) = {auc_sh:.3f}")
            else:
                print(f"  [WARN] Missing sample {sh_key} for {rep}")

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    colors = {
        "Full RMM": "C0",
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
    ax_hh.set_title(r"X $\to$ HH vs ttbar (VAE)")
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
    ax_sh.set_title(r"X $\to$ SH vs ttbar (VAE)")
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

    print("\n=== Loading RMM-C46 (frob, train/test) ===")
    X_c46_bkg, X_c46_sig = load_compact_pair(C46_BKG, C46_SIG, skip_cols=1)
    print("  C46 shapes:", X_c46_bkg.shape, X_c46_sig.shape)

    # --------------------------------------------------
    # 2) Train VAEs and compute ROC for each representation
    # --------------------------------------------------
    print("\n========== TRAINING VAEs AND COMPUTING AUC ==========\n")

    results = []

    res_full = vae_train_and_roc(X_full_bkg, X_full_sig, "Full RMM", device=device)
    results.append(res_full)

    res_c46 = vae_train_and_roc(X_c46_bkg, X_c46_sig, "RMM-C46", device=device)
    results.append(res_c46)

    # --------------------------------------------------
    # 3) ROC plot comparison
    # --------------------------------------------------
    print("\n========== SUMMARY (ROC AUC, VAE) ==========")
    for res in results:
        print(f"AUC ({res['label']}, VAE) : {res['auc']:.3f}")

    os.makedirs("out_C46", exist_ok=True)

    plt.figure(figsize=(6.5, 5.5))
    for res in results:
        plt.plot(
            res["fpr"],
            res["tpr"],
            label=f"{res['label']} (AUC={res['auc']:.3f})"
        )

    plt.plot([0, 1], [0, 1], "k--", alpha=0.5)
    plt.xlabel("False positive rate")
    plt.ylabel("True positive rate")
    plt.title("VAE ROC: Full RMM vs RMM-C46")
    plt.legend(loc="lower right", fontsize=9)
    plt.tight_layout()
    plt.savefig("out_C46/roc_compare_rmm_all_VAE.pdf")
    plt.close()
    print("Saved ROC plot to: out_C46/roc_compare_rmm_all_VAE.pdf")

    # --------------------------------------------------
    # 4) Multi-loss distribution plots for many samples
    # --------------------------------------------------
    print("\n========== MULTI-SAMPLE LOSS DISTRIBUTIONS (VAE) ==========")

    all_losses_dict = {}

    # Full RMM
    full_losses = compute_multi_losses(
        "Full RMM (VAE)",
        res_full["model"],
        res_full["scaler"],
        loader_single=load_full_rmm_single,
        samples_paths=FULL_EXTRA_SAMPLES,
        device=device,
    )
    all_losses_dict["Full RMM"] = full_losses
    
    # C46
    c46_losses = compute_multi_losses(
        "RMM-C46 (VAE)",
        res_c46["model"],
        res_c46["scaler"],
        loader_single=lambda p: load_compact_single(p, skip_cols=1),
        samples_paths=C46_EXTRA_SAMPLES,
        device=device,
    )
    all_losses_dict["RMM-C46"] = c46_losses

    plot_multi_loss_distributions(
        all_losses_dict,
        out_pdf="out_C46/reco_loss_distributions_rmm_all_VAE_multi.pdf"
    )

    # --------------------------------------------------
    # 5) AUC vs mass for HH/SH for each representation
    # --------------------------------------------------
    compute_and_plot_auc_vs_mass(
        all_losses_dict,
        out_pdf="out_C46/auc_vs_mass_all_VAE.pdf"
    )

if __name__ == "__main__":
    main()


