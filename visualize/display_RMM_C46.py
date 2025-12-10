#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE

# --------------------------------------------------
# Config: input files and labels
# --------------------------------------------------

'''
FILE_MAP_ADD = {
    "tt":    "out_C46/tt_c46_add.csv",
    "WZJet": "out_C46/wzjet_c46_add.csv",
    #"HH":    "out_C46/hh_1500_c46_add.csv",
    #"SH":    "out_C46/sh_1500_c46_add.csv",
}
'''
FILE_MAP_FROB = {
    "tt":    "out_C46/tt_c46_frob.csv",
    "WZJet": "out_C46/wzjet_c46_frob.csv",
    "HH":    "out_C46/hh_1500_c46_frob.csv",
    "SH":    "out_C46/sh_1500_c46_frob.csv",
}

OUTFIG_ADD_HEAT  = "out_C46/C46_mean_heatmap_add.pdf"
OUTFIG_FROB_HEAT = "out_C46/C46_mean_heatmap_frob.png"

OUTFIG_ADD_TSNE  = "out_C46/C46_tsne_add.pdf"
OUTFIG_FROB_TSNE = "out_C46/C46_tsne_frob.pdf"

# Correlation matrix output
OUTDIR_CORR = "out_C46"
# e.g. out_C46/C46_corr_tt_add.pdf, C46_corr_tt_frob.pdf

# Preferred order of samples (columns in heatmap, labels in t-SNE)
SAMPLE_ORDER = ["tt", "WZJet", "HH", "SH"]

# --------------------------------------------------
# Per-sample maximum events for t-SNE (None = use all)
# --------------------------------------------------
MAX_EVENTS = {
    "tt":    15000,
    "WZJet": 15000,
    "HH":    15000,
    "SH":    15000,
}

# ------------------------- mean heatmap tools -----------------------------

def build_mean_df(file_map):
    """
    Build a DataFrame of mean C46 values:
      rows = features (46),
      cols = samples in SAMPLE_ORDER (if present in file_map).
    """
    mean_vectors = {}
    feature_names = None

    for label, path in file_map.items():
        if not os.path.isfile(path):
            raise FileNotFoundError(path)

        df = pd.read_csv(path)

        # Assume first column is event index, remaining are the 46 C46 variables
        if feature_names is None:
            feature_names = df.columns[1:]
        else:
            if not np.array_equal(feature_names, df.columns[1:]):
                raise ValueError(f"Column mismatch in {path}")

        mean_vectors[label] = df.iloc[:, 1:].mean(axis=0)

    cols = [x for x in SAMPLE_ORDER if x in mean_vectors]

    mean_df = pd.DataFrame(
        {lab: mean_vectors[lab].values for lab in cols},
        index=feature_names
    )
    return mean_df, cols


def plot_heatmap(mean_df, sample_order, outfig, title):
    """
    Plot a 46×N heatmap:
      rows = C46 features,
      columns = samples.
    """
    feature_names = mean_df.index.values
    fig, ax = plt.subplots(figsize=(6.0, 10.0))

    im = ax.imshow(mean_df.values, aspect="auto")

    # Y-axis: feature names
    ax.set_yticks(np.arange(len(feature_names)))
    ax.set_yticklabels(feature_names, fontsize=7)

    # X-axis: sample labels on top
    ax.set_xticks(np.arange(len(sample_order)))
    ax.set_xticklabels(sample_order, fontsize=10)
    ax.xaxis.tick_top()

    ax.set_title(title, pad=30)

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Mean value", rotation=90)

    plt.tight_layout()
    plt.savefig(outfig)
    plt.close(fig)
    print(f"Saved heatmap: {outfig}")
    print("  Matrix shape:", mean_df.shape)


# ---------------------- correlation / similarity matrix tools --------------------------

def compute_corr_matrix_from_file(path, n_events=None):
    """
    Load a single C46 CSV file and compute the 46×46 correlation matrix.

    Parameters
    ----------
    path : str
        Path to the CSV file.
    n_events : int or None
        If not None, only the first `n_events` rows (events) are used
        to compute the correlation. If None, all events are used.

    Returns
    -------
    corr_df : pandas.DataFrame
        46×46 correlation matrix (features × features).

    Note
    ----
    If a feature has zero variance (e.g. always 0) in the subset used,
    its correlation with all others is undefined and will appear as NaN.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(path)

    df = pd.read_csv(path)

    # Restrict to the first n_events if requested
    if n_events is not None:
        df = df.iloc[:n_events, :]

    # Drop the first column (event index) and keep only C46 features
    features_df = df.iloc[:, 1:]

    # Pearson correlation between features; constant columns -> NaN correlations
    corr_df = features_df.corr()
    return corr_df


def compute_single_event_similarity(row):
    """
    Given a single event (row of 46 features), build a 46×46
    'correlation-like' similarity matrix from its normalized outer product.

        u = v / ||v||
        M_ij = u_i * u_j

    This is NOT a statistical correlation (needs many events),
    but a per-event pattern of co-occurring feature strengths.
    """
    # row is a pandas Series of length 46
    v = row.to_numpy(dtype=float)
    norm = np.linalg.norm(v)

    if norm == 0.0:
        # All features zero => similarity undefined; fill with NaNs
        mat = np.full((len(v), len(v)), np.nan)
    else:
        u = v / norm
        mat = np.outer(u, u)  # values in [-1, 1]

    sim_df = pd.DataFrame(mat, index=row.index, columns=row.index)
    return sim_df

def plot_corr_matrix(corr_df, outfig, title):
    """
    Plot a 46×46 matrix as a heatmap.

    Any cell equal to 0.0 is shown as white (empty/no structure).
    NaN entries (undefined correlation) are also shown as white.
    """
    feature_names = corr_df.columns.values
    corr_values = corr_df.values

    # Mask zeros AND NaNs so they appear as white.
    corr_masked = np.ma.masked_where(corr_values == 0.0, corr_values)
    corr_masked = np.ma.masked_invalid(corr_masked)

    # Use coolwarm colormap, but set masked entries ("bad") to white.
    cmap = plt.cm.coolwarm.copy()
    cmap.set_bad(color='white')

    fig, ax = plt.subplots(figsize=(8.0, 7.0))
    im = ax.imshow(
        corr_masked,
        vmin=-1.0,
        vmax=1.0,
        cmap=cmap,
        aspect="equal"
    )

    ax.set_xticks(np.arange(len(feature_names)))
    ax.set_yticks(np.arange(len(feature_names)))

    ax.set_xticklabels(feature_names, fontsize=5, rotation=90)
    ax.set_yticklabels(feature_names, fontsize=5)

    ax.set_title(title, pad=20)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Correlation / similarity", rotation=90)

    plt.tight_layout()
    plt.savefig(outfig)
    plt.close(fig)
    print(f"Saved matrix plot: {outfig}")
    print("  Matrix shape:", corr_df.shape)


def build_and_plot_corr_matrices(file_map, tag, n_single_events=5):
    """
    For a given file_map (ADD or FROB), build and plot matrices for each sample.

    For each sample (tt, WZJet, HH, SH):

      (A) Per-event similarity (outer product):
          - For the first `n_single_events` events (rows 0..n_single_events-1),
            compute a 46×46 similarity matrix from that single event.
          - Save each as:
                C46_corr_<sample>_<tag>_event<i>.pdf

      (B) Full-sample correlation:
          - Compute the usual 46×46 Pearson correlation matrix over ALL events.
          - Save as:
                C46_corr_<sample>_<tag>.pdf

    'tag' is a short label for the variant ('add' or 'frob').
    """
    os.makedirs(OUTDIR_CORR, exist_ok=True)

    for sample in SAMPLE_ORDER:
        if sample not in file_map:
            continue

        path = file_map[sample]
        print(f"\n=== Sample: {sample} ({tag}) ===")
        print(f"File: {path}")

        df = pd.read_csv(path)
        features_df = df.iloc[:, 1:]  # drop event index
        feature_names = features_df.columns

        # --- (A) Per-event similarity matrices for first n_single_events ---
        n_events_available = len(features_df)
        n_events_to_use = min(n_single_events, n_events_available)
        print(f"-> Building per-event similarity matrices for first {n_events_to_use} events...")

        for i_evt in range(n_events_to_use):
            row = features_df.iloc[i_evt]
            sim_df = compute_single_event_similarity(row)

            outfig_evt = os.path.join(
                OUTDIR_CORR,
                f"C46_corr_{sample}_{tag}_event{i_evt+1}.png"
            )
            #title_evt = f"RMM-C46({tag}) per-event similarity: {sample}, event {i_evt+1}"
            title_evt = f"Compressed RMM refgions - per-event correlation: {sample}, event {i_evt+1}"
            plot_corr_matrix(sim_df, outfig_evt, title_evt)

        # --- (B) Full-sample correlation matrix (over all events) ---
        print("-> Computing full-sample correlation matrix over ALL events...")
        corr_full = features_df.corr()

        outfig_full = os.path.join(
            OUTDIR_CORR,
            f"C46_corr_{sample}_{tag}.png"
        )
        title_full = f"RMM-C46({tag}) feature correlations: {sample} (all events)"
        plot_corr_matrix(corr_full, outfig_full, title_full)


# ------------------------------- main ------------------------------------

def main():
    '''
    # 1) ADD version heatmap
    print("=== Building ADD-based C46 heatmap ===")
    mean_df_add, cols_add = build_mean_df(FILE_MAP_ADD)
    plot_heatmap(
        mean_df_add,
        cols_add,
        OUTFIG_ADD_HEAT,
        "Mean RMM-C46(add) values per sample"
    )
    '''

    # 2) FROB version heatmap
    print("\n=== Building FROB-based C46 heatmap ===")
    mean_df_frob, cols_frob = build_mean_df(FILE_MAP_FROB)
    plot_heatmap(
        mean_df_frob,
        cols_frob,
        OUTFIG_FROB_HEAT,
        #"Mean RMM-C46(frob) values per sample"
        "Mean values on compressed RMM regions per sample"

    )

    # 3) Correlation / similarity matrices (FROB)
    #    - per-event outer-product similarities for first 5 events
    #    - full-sample correlations
    print("\n=== Building correlation / similarity matrices for RMM-C46(frob) ===")
    build_and_plot_corr_matrices(FILE_MAP_FROB, tag="frob", n_single_events=5)


if __name__ == "__main__":
    main()

