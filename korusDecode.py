#!/usr/bin/env python3
# RLW 04/2026

import argparse
import pandas as pd
import numpy as np

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

import matplotlib.pyplot as plt
import sys
import gzip
import io
import os
import plotly.express as px
import plotly.graph_objects as go

import matplotlib.cm as cm
import matplotlib.colors as mcolors

norm = mcolors.Normalize(vmin=0, vmax=1)

KORUSDECODE_VERSION = "v0.2.5"

# Define cluster colors
UNSUPERVISED_PALETTE = [
    "#4E79A7",  # blue
    "#F28E2B",  # orange
    "#E15759",  # red
    "#76B7B2",  # teal
    "#59A14F",  # green
    "#EDC948",  # yellow
    "#B07AA1",  # purple
    "#FF9DA7",  # pink
    "#9C755F",  # brown
    "#BAB0AC",  # gray
    "#8dd3c7",  # pastels V
    "#ffffb3",
    "#bebada",
    "#fb8072",
    "#80b1d3",
    "#fdb462",
    "#b3de69",
    "#fccde5",
    "#d9d9d9",
    "#bc80bd"
]

SPECIAL_CLASSES = {"NA", "REF"}

# -------------------------
# Define unsupervised colour
# -------------------------
def get_cluster_color(cluster_id, unique_clusters):
    palette = UNSUPERVISED_PALETTE
    #idx = list(unique_clusters).index(cluster_id)
    cluster_to_idx = {c: i for i, c in enumerate(unique_clusters)}
    idx = cluster_to_idx[cluster_id]
    return palette[idx % len(palette)]

# -------------------------
# Define dynamic colour
# -------------------------
def build_dynamic_ancestry_colors(ancestries, color_map, palette):
    dynamic = {}
    used = set(color_map.keys())

    # stable ordering prevents nondeterministic coloring
    unknowns = sorted([a for a in ancestries if a not in used])

    for i, anc in enumerate(unknowns):
        dynamic[anc] = palette[i % len(palette)]

    return dynamic

# -------------------------
# Tee logger: duplicate stdout to file
# -------------------------
class Tee:
    def __init__(self, fname):
        # Open log file and keep reference to original stdout
        self.file = open(fname, "w")
        self.stdout = sys.stdout

    def write(self, data):
        # Write to both file and terminal
        self.file.write(data)
        self.stdout.write(data)

    def flush(self):
        # Ensure both outputs are flushed
        self.file.flush()
        self.stdout.flush()

# -------------------------
# Support gzip
# -------------------------
def open_maybe_gz(path):
    """
    Open plain text or gz/bgz files as a text stream.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"Input file does not exist: {path}")

    if path.endswith((".gz", ".bgz", ".gzip")):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")

    return open(path, "r")

# -------------------------
# Entropy calcs
# -------------------------
def compute_entropy(row):
    p = row.values.astype(float)
    p = p[p > 0]
    return -np.sum(p * np.log(p))

# -------------------------
# PCA plot titles
# -------------------------
def make_pca_title(mode, dim, k=None, hmm=False, interactive=False):
    base = "PCA"

    if dim == 3 and interactive:
        base = "PCA"
    elif dim == 3:
        base = "3D PCA"
    else:
        base = "PCA"

    if mode == "fractions":
        title = f"Soft-supervised ancestry {base} projection"
    elif mode == "labels":
        title = f"Hard-supervised ancestry {base} projection"
    else:
        title = f"Unsupervised ancestry clustering (KMeans k={k})"

    if hmm:
        title += " [HMM-smoothed]"

    return title

# -------------------------
# Load reference matrix
# -------------------------
def load_reference_matrix(path):
    # Load TSV with chromosome as string (prevents numeric coercion)
    df = pd.read_csv(open_maybe_gz(path), sep="\t", dtype={"chr": str})

    # Extract window size from header (format: window_start:<size>)
    window_col = [c for c in df.columns if c.startswith("window_start")][0]
    window_size = int(window_col.split(":")[1])

    # Normalize column name
    df = df.rename(columns={window_col: "window_start"})

    # Separate genomic coordinates
    coords = df[["chr", "window_start"]]

    # Transpose: rows = samples, columns = genomic windows
    X = df.drop(columns=["chr", "window_start"]).T

    return X, coords, window_size

# -------------------------
# Load sample matrix
# -------------------------
def load_sample_matrix(path):
    # Same logic as reference loader
    df = pd.read_csv(open_maybe_gz(path), sep="\t", dtype={"chr": str})

    window_col = [c for c in df.columns if c.startswith("window_start")][0]
    window_size = int(window_col.split(":")[1])

    df = df.rename(columns={window_col: "window_start"})

    coords = df[["chr", "window_start"]]
    X = df.drop(columns=["chr", "window_start"]).T

    return X, coords, window_size

# -------------------------
# Load ancestry labels
# -------------------------
def load_labels(path):
    # Load label file (expects sample ↔ ancestry mapping)
    labels = pd.read_csv(open_maybe_gz(path), sep="\t")

    # Normalize column names for robustness
    labels.columns = [c.strip().lower() for c in labels.columns]

    # Validate required columns
    if "sample" not in labels.columns or "ancestry" not in labels.columns:
        raise ValueError("Labels file must contain 'sample' and 'ancestry' columns")

    # Use sample IDs as index
    labels = labels.set_index("sample")

    return labels

# -------------------------
# Load ancestry fractions
# -------------------------
def load_fraction_labels(path):
    """
    Load ancestry fraction table.

    Expected format:
    sample  AFR  EUR  EAS ...
    HG001   0.7  0.3  0.0

    Returns:
        frac_df : DataFrame indexed by sample
    """

    frac_df = pd.read_csv(open_maybe_gz(path), sep="\t")

    # normalize column names
    frac_df.columns = [c.strip() for c in frac_df.columns]

    # require sample column
    if "sample" not in frac_df.columns:
        raise ValueError(
            "Fractions file must contain a 'sample' column"
        )

    frac_df = frac_df.set_index("sample")

    # ensure numeric ancestry columns
    for col in frac_df.columns:
        frac_df[col] = pd.to_numeric(frac_df[col], errors="raise")

    # optional normalization safeguard
    row_sums = frac_df.sum(axis=1)

    if np.any(np.abs(row_sums - 1.0) > 1e-3):
        print("WARNING: Fraction rows do not sum to 1. Renormalizing.")
        frac_df = frac_df.div(row_sums, axis=0)

    return frac_df

# -------------------------
# Weighted Gaussian estimation
# -------------------------
def estimate_weighted_gaussians(Z_train, frac_df, classes, reg=1e-3):
    """
    Estimate weighted ancestry-specific Gaussian parameters.

    Prints diagnostic information about the fractional training
    populations but otherwise preserves the original behavior.
    """

    means = {}
    covs = {}

    print("\nWeighted Gaussian diagnostics:")

    for c in classes:

        weights = frac_df[c].values.astype(float)

        weight_sum = np.sum(weights)

        if weight_sum <= 0:
            continue

        # Effective sample size
        neff = (
            weight_sum ** 2 /
            np.sum(weights ** 2)
        )

        # Weighted mean
        mu = np.average(
            Z_train,
            axis=0,
            weights=weights
        )

        # Weighted covariance
        diff = Z_train - mu

        cov = np.einsum(
            "i,ij,ik->jk",
            weights,
            diff,
            diff
        ) / weight_sum

        # Regularization
        cov += np.eye(cov.shape[0]) * reg

        # Diagnostics
        mean_variance = np.mean(np.diag(cov))

        sign, logdet = np.linalg.slogdet(cov)

        if sign <= 0:
            raise ValueError(
                f"Covariance matrix for {c} "
                f"is not positive definite"
            )

        print(
            f"  {c:4s} "
            f"sum={weight_sum:10.2f} "
            f"Neff={neff:10.2f} "
            f"meanVar={mean_variance:10.5f} "
            f"logdet={logdet:10.3f}"
        )

        means[c] = mu
        covs[c] = np.linalg.inv(cov)

    return means, covs

# -------------------------
# Fraction-supervised assignment
# -------------------------
def window_fractional_assignment(
    X_ref,
    X_sample,
    frac_labels,
    local_pcs=10,
    knn_mode="window",
    window_radius=3,
    require_full_window=False
):

    preds_all = [[] for _ in range(X_sample.shape[0])]
    probs_all = [[] for _ in range(X_sample.shape[0])]

    cols_list = list(X_ref.columns)
    classes = np.array(list(frac_labels.columns))

    for i in range(len(cols_list)):

        # -------------------------
        # Window selection
        # -------------------------

        if knn_mode == "raw":
            cols = [cols_list[i]]

        elif knn_mode == "window":
            start = max(0, i - window_radius)
            end = min(len(cols_list), i + window_radius + 1)
            cols = cols_list[start:end]

            if require_full_window and len(cols) != (2 * window_radius + 1):
                continue

        else:
            raise ValueError(f"Unknown knn_mode: {knn_mode}")

        if len(cols) == 0:
            continue

        # -------------------------
        # Data extraction
        # -------------------------

        X_train = X_ref[cols].values
        X_test = X_sample[cols].values

        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)

        # -------------------------
        # PCA stabilization
        # -------------------------

        n_comp = min(
            local_pcs,
            X_train.shape[1],
            X_train.shape[0] - 1
        )

        pca = PCA(n_components=n_comp)

        Z_train = pca.fit_transform(X_train)
        Z_test = pca.transform(X_test)

        # -------------------------
        # Weighted Gaussian estimation
        # -------------------------

        means, covs = estimate_weighted_gaussians(
            Z_train,
            frac_labels,
            classes
        )

        # -------------------------
        # Calculate log determinants
        # -------------------------

        logdets = {}

        for c in classes:
            # estimate_weighted_gaussians() returns the inverse
            # covariance matrix, so recover the covariance matrix
            # for calculation of log|Sigma|.
            cov_inv = covs[c]

            sign, logdet_inv = np.linalg.slogdet(cov_inv)

            if sign <= 0:
                raise ValueError(
                    f"Inverse covariance matrix for {c} "
                    f"is not positive definite"
                )

            # log|Sigma| = -log|Sigma^-1|
            logdets[c] = -logdet_inv

        # -------------------------
        # Classify samples
        # -------------------------

        for s in range(Z_test.shape[0]):

            z = Z_test[s]
            scores = []

            for c in classes:

                diff = z - means[c]

                # Mahalanobis distance
                dist = diff @ covs[c] @ diff

                # Full Gaussian log-likelihood
                # (common -d/2*log(2*pi) term omitted)
                score = -0.5 * (
                    dist + logdets[c]
                )

                scores.append(score)

            scores = np.array(scores)

            # Softmax normalization
            exp_scores = np.exp(
                scores - np.max(scores)
            )

            probs = exp_scores / exp_scores.sum()

            pred = classes[np.argmax(probs)]

            preds_all[s].append(pred)
            probs_all[s].append(probs)

    return preds_all, probs_all, classes

# -------------------------
# Align windows between ref and sample
# -------------------------
def align_matrices(X_ref, coords_ref, X_sample, coords_sample):
    # Construct unique window IDs (chr:start)
    ref_ids = coords_ref["chr"].astype(str) + ":" + coords_ref["window_start"].astype(str)
    sample_ids = coords_sample["chr"].astype(str) + ":" + coords_sample["window_start"].astype(str)

    # Assign these as column names
    X_ref.columns = ref_ids
    X_sample.columns = sample_ids

    # Find shared genomic windows
    common = X_ref.columns.intersection(X_sample.columns)

    if len(common) == 0:
        raise ValueError("No overlapping windows between reference and sample")

    # Subset both matrices to shared windows
    X_ref = X_ref[common]
    X_sample = X_sample[common]

    print(f"Using {len(common)} shared windows")
    return X_ref, X_sample

# -------------------------
# Global PCA (for visualization only)
# -------------------------
def run_pca(X_ref, X_sample, n_components=2):
    # Standardize features (per window)
    scaler = StandardScaler()

    X_ref_scaled = scaler.fit_transform(X_ref)

    # Fit PCA on reference only
    pca = PCA(n_components=n_components, svd_solver="full")
    X_ref_pca = pca.fit_transform(X_ref_scaled)

    # Project samples into same PCA space
    X_sample_scaled = scaler.transform(X_sample)
    X_sample_pca = pca.transform(X_sample_scaled)

    return X_ref_pca, X_sample_pca

# -------------------------
# Inference PCA
# -------------------------
def run_inference_pca(X_ref, X_sample, n_components):

    scaler = StandardScaler()

    X_ref_scaled = scaler.fit_transform(X_ref)

    n_pc = min(
        n_components,
        X_ref.shape[0] - 1,
        X_ref.shape[1]
    )

    pca = PCA(
        n_components=n_pc,
        svd_solver="full"
    )

    X_ref_pca = pca.fit_transform(X_ref_scaled)

    X_sample_scaled = scaler.transform(X_sample)

    X_sample_pca = pca.transform(X_sample_scaled)

    return X_ref_pca, X_sample_pca

# -------------------------
# KMeans clustering
# -------------------------
def run_kmeans(X, k):
    # Standard KMeans clustering with fixed seed for reproducibility
    kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
    clusters = kmeans.fit_predict(X)
    return clusters, kmeans

# -------------------------
# PCA Plot (2D)
# -------------------------
def plot_pca(X_ref_pca, labels, clusters, X_sample_pca, X_sample_names, outprefix, title):
    plt.figure(figsize=(8, 6))

    # Predefined ancestry color scheme
    color_map = {
        "AFR": "orange",
        "EUR": "green",
        "EAS": "red",
        "AMR": "black",
        "SAS": "blue",
        "CAS": "violet",
        "OCE": "lightgreen",
        "NA": "#c000c0",
        "REF": "#d3d3d3" 
    }

    all_anc = labels["ancestry"].unique() if labels is not None else []
    dynamic_color_map = build_dynamic_ancestry_colors(
        all_anc,
        color_map,
        UNSUPERVISED_PALETTE
    )

    # Determine ancestry categories
    if labels is not None:
        unique_anc = sorted(labels["ancestry"].unique())
    else:
        unique_anc = []

    # -------------------------
    # Plot reference samples
    # -------------------------
    if labels is not None:
        for i, sample in enumerate(labels.index):
            anc = labels.loc[sample, "ancestry"]

            if anc in color_map:
                color = color_map[anc]
            elif anc in dynamic_color_map:
                color = dynamic_color_map[anc]
            elif anc in SPECIAL_CLASSES:
                color = color_map.get(anc, "grey")
            else:
                color = dynamic_color_map.get(anc, "#999999")

            if clusters is not None and clusters[i] >= 0:
                label = str(clusters[i] + 1)
            else:
                label = ""

            plt.scatter(X_ref_pca[i, 0], X_ref_pca[i, 1], c=color, s=20, alpha=0.6)
            plt.text(X_ref_pca[i, 0], X_ref_pca[i, 1], label, fontsize=6, alpha=0.5)

    # -------------------------
    # UNSUPERVISED MODE
    # -------------------------
    else:
        unique_clusters = np.unique(clusters)

        # stable cluster → color mapping (palette-based)
        cluster_to_color = {
            c: get_cluster_color(c, unique_clusters)
            for c in unique_clusters
        }

        for i in range(X_ref_pca.shape[0]):

            cluster_id = clusters[i]
            color = cluster_to_color[cluster_id]

            label = f"C{cluster_id + 1}"

            plt.scatter(
                X_ref_pca[i, 0],
                X_ref_pca[i, 1],
                c=color,
                s=12,
                alpha=0.7
            )

            plt.text(
                X_ref_pca[i, 0],
                X_ref_pca[i, 1],
                label,
                fontsize=5,
                alpha=0.6
            )

        # cluster legend (consistent with points)
        for c in unique_clusters:
            plt.scatter([], [], c=cluster_to_color[c], label=f"C{c + 1}")

    # -------------------------
    # Plot query samples
    # -------------------------
    for i, name in enumerate(X_sample_names):
        plt.scatter(X_sample_pca[i, 0], X_sample_pca[i, 1], c="black", marker="x", s=100)
        plt.text(X_sample_pca[i, 0], X_sample_pca[i, 1], name, fontsize=7, fontweight="bold")

    # -------------------------
    # Build legend (ancestries)
    # -------------------------
    for anc in unique_anc:
        if anc in color_map:
            color = color_map[anc]
        elif anc in dynamic_color_map:
            color = dynamic_color_map[anc]
        elif anc in SPECIAL_CLASSES:
            color = color_map.get(anc, "grey")
        else:
            color = dynamic_color_map.get(anc, "#999999")
        plt.scatter([], [], c=color, label=anc)

    plt.scatter([], [], c="black", marker="x", s=100, label="query")

    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title(title)
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(f"{outprefix}.pca.png", dpi=300)
    print(f"\tWriting PCA plot to: {outprefix}.pca.png")

# -------------------------
# PCA Plot (3D)
# -------------------------
def plot_pca_3d(X_ref, X_sample, labels, clusters, sample_names, outprefix, title):
    from mpl_toolkits.mplot3d import Axes3D

    # Standardize before PCA
    scaler = StandardScaler()
    X_ref_scaled = scaler.fit_transform(X_ref)
    X_sample_scaled = scaler.transform(X_sample)

    # 3D PCA
    pca = PCA(n_components=3)
    X_ref_pca = pca.fit_transform(X_ref_scaled)
    X_sample_pca = pca.transform(X_sample_scaled)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    # -------------------------
    # Color map (supervised)
    # -------------------------
    color_map = {
        "AFR": "orange",
        "EUR": "green",
        "EAS": "red",
        "AMR": "black",
        "SAS": "blue",
        "CAS": "violet",
        "OCE": "lightgreen",
        "NA": "#c000c0",
        "REF": "#d3d3d3"
    }

    all_anc = labels["ancestry"].unique() if labels is not None else []
    dynamic_color_map = build_dynamic_ancestry_colors(
        all_anc,
        color_map,
        UNSUPERVISED_PALETTE
    )

    if labels is not None:
        unique_anc = sorted(labels["ancestry"].unique())
    else:
        unique_anc = []

    # =========================================================
    # UNSUPERVISED MODE
    # =========================================================
    if labels is None:

        unique_clusters = np.unique(clusters)

        cluster_to_color = {
            c: get_cluster_color(c, unique_clusters)
            for c in unique_clusters
        }

        for i in range(X_ref_pca.shape[0]):

            cluster_id = clusters[i]
            color = cluster_to_color[cluster_id]
            label = f"C{cluster_id + 1}"

            ax.scatter(
                X_ref_pca[i, 0],
                X_ref_pca[i, 1],
                X_ref_pca[i, 2],
                c=color,
                s=12,
                alpha=0.7,
                depthshade=False,
                zorder=1
            )

            ax.text(
                X_ref_pca[i, 0],
                X_ref_pca[i, 1],
                X_ref_pca[i, 2],
                label,
                fontsize=5,
                alpha=0.6
            )

        # cluster legend (consistent)
        for c in unique_clusters:
            ax.scatter([], [], [], c=cluster_to_color[c], label=f"C{c + 1}")

    # =========================================================
    # SUPERVISED MODE
    # =========================================================
    else:

        for i, sample in enumerate(labels.index):
            anc = labels.loc[sample, "ancestry"]

            if anc in color_map:
                color = color_map[anc]
            elif anc in dynamic_color_map:
                color = dynamic_color_map[anc]
            elif anc in SPECIAL_CLASSES:
                color = color_map.get(anc, "grey")
            else:
                color = dynamic_color_map.get(anc, "#999999")

            if clusters is not None and clusters[i] >= 0:
                label = str(clusters[i] + 1)
            else:
                label = ""

            ax.scatter(
                X_ref_pca[i, 0],
                X_ref_pca[i, 1],
                X_ref_pca[i, 2],
                c=color,
                s=20,
                alpha=0.6,
                depthshade=False,
                zorder=1
            )

            ax.text(
                X_ref_pca[i, 0],
                X_ref_pca[i, 1],
                X_ref_pca[i, 2],
                label,
                fontsize=6,
                alpha=0.5
            )

        # ancestry legend (SINGLE PASS — FIXED)
        for anc in unique_anc:
            if anc in color_map:
                color = color_map[anc]
            elif anc in dynamic_color_map:
                color = dynamic_color_map[anc]
            elif anc in SPECIAL_CLASSES:
                color = color_map.get(anc, "grey")
            else:
                color = dynamic_color_map.get(anc, "#999999")
            ax.scatter([], [], [], c=color, label=anc)

    # -------------------------
    # Query samples
    # -------------------------
    for i, name in enumerate(sample_names):
        ax.scatter(
            X_sample_pca[i, 0],
            X_sample_pca[i, 1],
            X_sample_pca[i, 2],
            c="black",
            marker="x",
            s=100,
            depthshade=False,
            zorder=1000
        )
        ax.text(
            X_sample_pca[i, 0],
            X_sample_pca[i, 1],
            X_sample_pca[i, 2],
            name,
            fontsize=7,
            fontweight="bold",
            zorder=1001,
            clip_on=False
        )

    # -------------------------
    # FINAL LEGEND (ONLY ONCE)
    # -------------------------
    ax.scatter([], [], [], c="black", marker="x", s=100, label="query")

    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_zlabel("PC3")
    ax.set_title(title)

    ax.legend(loc='center left', bbox_to_anchor=(1.07, 0.5))
    plt.tight_layout(rect=[0, 0, 0.8, 1])
    plt.savefig(f"{outprefix}.pca3d.png", dpi=300)

    print(f"\tWriting 3D PCA plot to: {outprefix}.pca3d.png")

# -------------------------
# Interactive PCA Plot (3D)
# -------------------------
def plot_pca_3dint(X_ref, X_sample, labels, clusters, sample_names, outprefix, frac_df, title):

    # -------------------------
    # Helper: bar-style fractions
    # -------------------------
    def format_bar(fractions, width=14, min_frac=0.01):

        items = [
            (k, v)
            for k, v in fractions.items()
            if v > min_frac
        ]

        items.sort(key=lambda x: -x[1])

        lines = []

        for pop, val in items:

            try:
                pop_label = f"C{int(pop) + 1}"
            except:
                pop_label = str(pop)

            n = int(val * width)
            bar = "█" * n + " " * (width - n)

            lines.append(f"{pop_label:<4} {val:.2f} {bar}")

        return "<br>".join(lines) if lines else "NA"

    # -------------------------
    # Color map (base)
    # -------------------------
    color_map = {
        "AFR": "orange",
        "EUR": "green",
        "EAS": "red",
        "AMR": "black",
        "SAS": "blue",
        "CAS": "violet",
        "query": "cyan",
        "OCE": "lightgreen",
        "NA": "#c000c0",
        "REF": "#d3d3d3"
    }

    all_anc = labels["ancestry"].unique() if labels is not None else []
    dynamic_color_map = build_dynamic_ancestry_colors(
        all_anc,
        color_map,
        UNSUPERVISED_PALETTE
    )

    # -------------------------
    # Ensure numeric only
    # -------------------------
    X_ref = X_ref.select_dtypes(include=[np.number])
    X_sample = X_sample.select_dtypes(include=[np.number])

    scaler = StandardScaler()
    X_ref_scaled = scaler.fit_transform(X_ref)
    X_sample_scaled = scaler.transform(X_sample)

    pca = PCA(n_components=3)
    X_ref_pca = pca.fit_transform(X_ref_scaled)
    X_sample_pca = pca.transform(X_sample_scaled)

    # -------------------------
    # Build dataframe
    # -------------------------
    df_list = []

    # =========================================================
    # REFERENCE POINTS
    # =========================================================
    if labels is not None:

        unique_anc = sorted(labels["ancestry"].unique())

        for i in range(X_ref_pca.shape[0]):

            plot_label = labels.iloc[i]["ancestry"]
            ref_name = labels.index[i]

            df_list.append({
                "PC1": X_ref_pca[i, 0],
                "PC2": X_ref_pca[i, 1],
                "PC3": X_ref_pca[i, 2],
                "label": plot_label,
                "type": "reference",
                "name": ref_name,
                "cluster": ""
            })

    else:
        unique_clusters = np.unique(clusters)
        unique_anc = [f"C{i+1}" for i in unique_clusters]

        cluster_to_color_label = {
            c: f"C{c + 1}" for c in unique_clusters
        }

        for i in range(X_ref_pca.shape[0]):

            cluster_id = clusters[i]
            plot_label = cluster_to_color_label[cluster_id]

            df_list.append({
                "PC1": X_ref_pca[i, 0],
                "PC2": X_ref_pca[i, 1],
                "PC3": X_ref_pca[i, 2],
                "label": plot_label,
                "type": "reference",
                #"name": f"ref_{i}",
                "name": X_ref.index[i],
                "cluster": plot_label
            })

    # =========================================================
    # QUERY POINTS
    # =========================================================
    query_labels = []

    for i in range(X_sample_pca.shape[0]):

        sample = sample_names[i]

        if sample in frac_df.index:
            fractions = frac_df.loc[sample].to_dict()
            frac_str = format_bar(fractions)
        else:
            frac_str = "NA"

        df_list.append({
            "PC1": X_sample_pca[i, 0],
            "PC2": X_sample_pca[i, 1],
            "PC3": X_sample_pca[i, 2],
            "label": "query",
            "type": "query",
            "name": sample,
            "fractions": frac_str,
            "cluster": ""
        })

        # IMPORTANT: label text only for query points
        query_labels.append(sample)

    df = pd.DataFrame(df_list)

    # -------------------------
    # COMPLETE color map
    # -------------------------
    color_discrete_map = dict(color_map)

    if labels is not None:
        for anc in unique_anc:
            if anc not in color_discrete_map:
                color_discrete_map[anc] = dynamic_color_map.get(anc, "#999999")
    else:
        unique_clusters = np.unique(clusters)

        cluster_to_color = {
            f"C{c+1}": get_cluster_color(c, unique_clusters)
            for c in unique_clusters
        }

        color_discrete_map.update(cluster_to_color)

    color_discrete_map["query"] = "cyan"

    # -------------------------
    # BUILD FIG (GRAPH OBJECTS)
    # -------------------------
    fig = go.Figure()

    # -------------------------
    # References
    # -------------------------
    for lab in sorted(df[df["type"] == "reference"]["label"].unique()):

        sub = df[(df["type"] == "reference") & (df["label"] == lab)]

        fig.add_trace(go.Scatter3d(
            x=sub["PC1"],
            y=sub["PC2"],
            z=sub["PC3"],
            mode="markers",
            name=lab,
            marker=dict(size=5, color=color_discrete_map.get(lab, "#999999"), opacity=0.6),
            text=sub["name"],
            customdata=sub["cluster"],
            hovertemplate=
                "<b>%{text}</b><br>"
                "PC1: %{x:.3f}<br>"
                "PC2: %{y:.3f}<br>"
                "PC3: %{z:.3f}<br><br>"
                "%{customdata}<br>"
                "<extra></extra>"
        ))

    # -------------------------
    # Queries (WITH LABELS ADDED HERE)
    # -------------------------
    qdf = df[df["type"] == "query"]

    fig.add_trace(go.Scatter3d(
        x=qdf["PC1"],
        y=qdf["PC2"],
        z=qdf["PC3"],
        mode="markers+text",
        name="query",
        text=query_labels,
        textposition="top center",
        textfont=dict(size=9),
        marker=dict(size=5, color="cyan", symbol="x"),
        customdata=qdf["fractions"],
        hovertemplate=
            "<b>%{text}</b><br>"
            "PC1: %{x:.3f}<br>"
            "PC2: %{y:.3f}<br>"
            "PC3: %{z:.3f}<br><br>"
            "%{customdata}"
            "<extra></extra>"
    ))

    # -------------------------
    # layout
    # -------------------------
    fig.update_layout(
        autosize=True,
        scene=dict(
            xaxis_title="PC1",
            yaxis_title="PC2",
            zaxis_title="PC3"
        ),
        title=dict(
            text=title,
            x=0.5,
            xanchor="center",
            y=0.95,
            yanchor="top",
            pad=dict(t=0, b=0),
            font=dict(size=18)
        ),
        legend=dict(
        #    title=dict(text="Samples"),
            font=dict(size=12),
            itemsizing="constant"
        ),
        margin=dict(l=40, r=140, t=60, b=40)
    )

    outfile = f"{outprefix}.pca3d.html"
    fig.write_html(outfile, config={"responsive": True})

    print(f"\tWriting interactive PCA plot to: {outfile}")

# -------------------------
# PCA Plot (2D) -- fractions mode
# -------------------------
def plot_pca_fraction_2d(X_ref_pca, entropy_vals, X_sample_pca, X_sample_names, outprefix, title):

    fig, ax = plt.subplots(figsize=(8, 6))

    cmap = plt.cm.viridis

    # -------------------------
    # Reference samples
    # -------------------------
    sc = ax.scatter(
        X_ref_pca[:, 0],
        X_ref_pca[:, 1],
        c=entropy_vals.values,
        cmap=cmap,
        norm=norm,
        s=15,
        alpha=0.7
    )

    # -------------------------
    # Query samples
    # -------------------------
    ax.scatter(
        X_sample_pca[:, 0],
        X_sample_pca[:, 1],
        c="black",
        marker="x",
        s=100,
        label="query"
    )

    for i, name in enumerate(X_sample_names):
        ax.text(
            X_sample_pca[i, 0],
            X_sample_pca[i, 1],
            name,
            fontsize=7,
            fontweight="bold"
        )

    # -------------------------
    # Colorbar legend
    # -------------------------
    cbar = plt.colorbar(
        sc,
        ax=ax,
        orientation="horizontal",
        pad=0.12,
        fraction=0.05
    )

    cbar.set_label("Normalized ancestry entropy")

    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")

    ax.set_title(title)

    ax.legend(loc="lower right")

    plt.tight_layout()

    outfile = f"{outprefix}.fractions.pca.png"

    plt.savefig(outfile, dpi=300)

    print(f"\tWriting entropy PCA plot to: {outfile}")

# -------------------------
# PCA Plot (3D) -- fractions mode
# -------------------------
def plot_pca_fraction_3d(X_ref, X_sample, entropy_vals, sample_names, outprefix, title):

    scaler = StandardScaler()

    X_ref_scaled = scaler.fit_transform(X_ref)
    X_sample_scaled = scaler.transform(X_sample)

    pca = PCA(n_components=3)

    X_ref_pca = pca.fit_transform(X_ref_scaled)
    X_sample_pca = pca.transform(X_sample_scaled)

    fig = plt.figure(figsize=(8, 6))

    ax = fig.add_subplot(111, projection='3d')

    cmap = plt.cm.viridis

    # -------------------------
    # Reference samples
    # -------------------------
    sc = ax.scatter(
        X_ref_pca[:, 0],
        X_ref_pca[:, 1],
        X_ref_pca[:, 2],
        c=entropy_vals.values,
        cmap=cmap,
        norm=norm,
        s=15,
        alpha=0.8,
        depthshade=False,
        zorder=1
    )

    # -------------------------
    # Query samples
    # -------------------------
    ax.scatter(
        X_sample_pca[:, 0],
        X_sample_pca[:, 1],
        X_sample_pca[:, 2],
        c="black",
        marker="x",
        s=100,
        label="query",
        depthshade=False,
        zorder=1000
    )

    for i, name in enumerate(sample_names):

        ax.text(
            X_sample_pca[i, 0],
            X_sample_pca[i, 1],
            X_sample_pca[i, 2],
            name,
            fontsize=7,
            fontweight="bold",
            zorder=1001,
            clip_on=False
        )

    # -------------------------
    # Colorbar
    # -------------------------
    cbar = fig.colorbar(
        sc,
        ax=ax,
        orientation="horizontal",
        shrink=0.7,
        pad=0.08,
        fraction=0.05
    )

    cbar.set_label("Normalized ancestry entropy")

    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_zlabel("PC3")

    ax.set_title(title)

    ax.legend(loc="center left",bbox_to_anchor=(1.07, 0.5))
    plt.tight_layout(rect=[0, 0.05, 1, 1])

    outfile = f"{outprefix}.fractions.pca3d.png"

    plt.savefig(outfile, dpi=300)

    print(f"\tWriting entropy 3D PCA plot to: {outfile}")

# -------------------------
# Interactive PCA Plot (3D) -- fractions mode
# -------------------------
def plot_pca_fraction_3dint(X_ref, X_sample, entropy_vals, sample_names, outprefix, frac_df, title):

    # -------------------------
    # Helper: bar-style fractions
    # -------------------------
    def format_bar(fractions, width=14, min_frac=0.01):

        items = [
            (k, v)
            for k, v in fractions.items()
            if v > min_frac
        ]

        items.sort(key=lambda x: -x[1])

        lines = []

        for pop, val in items:

            n = int(val * width)

            bar = "█" * n + " " * (width - n)

            lines.append(
                f"{pop:<4} {val:.2f} {bar}"
            )

        return "<br>".join(lines) if lines else "NA"

    # -------------------------
    # PCA
    # -------------------------
    scaler = StandardScaler()

    X_ref_scaled = scaler.fit_transform(X_ref)
    X_sample_scaled = scaler.transform(X_sample)

    pca = PCA(n_components=3)

    X_ref_pca = pca.fit_transform(X_ref_scaled)
    X_sample_pca = pca.transform(X_sample_scaled)

    # -------------------------
    # Build reference dataframe
    # -------------------------
    ref_df = pd.DataFrame({
        "PC1": X_ref_pca[:, 0],
        "PC2": X_ref_pca[:, 1],
        "PC3": X_ref_pca[:, 2],
        "entropy": entropy_vals.values,
        "name": entropy_vals.index
    })

    # -------------------------
    # Build query dataframe
    # -------------------------
    query_rows = []

    for i, sample in enumerate(sample_names):

        if sample in frac_df.index:

            fractions = frac_df.loc[sample].to_dict()

            frac_str = format_bar(fractions)

        else:

            frac_str = "NA"

        query_rows.append({
            "PC1": X_sample_pca[i, 0],
            "PC2": X_sample_pca[i, 1],
            "PC3": X_sample_pca[i, 2],
            "name": sample,
            "fractions": frac_str
        })

    query_df = pd.DataFrame(query_rows)

    # -------------------------
    # Plot reference samples
    # -------------------------
    fig = px.scatter_3d(
        ref_df,
        x="PC1",
        y="PC2",
        z="PC3",
        color="entropy",
        color_continuous_scale="Viridis",
        range_color=[0, 1],
        hover_name="name",
        opacity=0.6
    )

    # Match original PCA_3Dint point sizing
    fig.update_traces(
        selector=dict(type="scatter3d", mode="markers"),
        marker=dict(size=5)
    )

    # -------------------------
    # Add query samples
    # -------------------------
    fig.add_scatter3d(
        x=query_df["PC1"],
        y=query_df["PC2"],
        z=query_df["PC3"],
        mode="markers+text",
        text=query_df["name"],
        textposition="top center",
        textfont=dict(size=9),

        marker=dict(
            color="cyan",
            size=5,
            symbol="x"
        ),

        customdata=query_df[["fractions"]].values,

        name="query",

        showlegend=True,

        hovertemplate=(
            "<b>%{text}</b><br>"
            "PC1: %{x:.3f}<br>"
            "PC2: %{y:.3f}<br>"
            "PC3: %{z:.3f}<br><br>"
            "%{customdata[0]}"
            "<extra></extra>"
        )
    )

    # -------------------------
    # Legend cleanup
    # -------------------------
    for tr in fig.data:

        if tr.name == "query":
            tr.showlegend = True
            tr.marker.size = 6
        else:
            tr.showlegend = False

    # -------------------------
    # Reference hover
    # -------------------------
    for tr in fig.data:

        if tr.name != "query":

            tr.hovertemplate = (
                "<b>%{hovertext}</b><br>"
                "PC1: %{x:.3f}<br>"
                "PC2: %{y:.3f}<br>"
                "PC3: %{z:.3f}<br><br>"
                "Entropy: %{marker.color:.3f}"
                "<extra></extra>"
            )

    # -------------------------
    # Layout
    # -------------------------
    fig.update_layout(
        autosize=True,
        title=dict(
            text=title,
            x=0.5,
            xanchor="center",
            y=0.95,
            yanchor="top",
            pad=dict(t=0, b=0),
            font=dict(size=18)
        ),
        #font=dict(size=14),
        coloraxis_colorbar=dict(
            title=dict(
                text="Normalized entropy",
                font=dict(size=18)
            ),
            tickfont=dict(size=15),
            orientation="h",
            x=0.5,
            xanchor="center",
            y=-0.18,
            len=0.6,
            thickness=18
        ),

        legend=dict(
        #    title=dict(text="Samples"),
            font=dict(size=12),
            itemsizing="constant"
        ),

        margin=dict(l=40, r=140, t=60, b=40),
    )

    outfile = f"{outprefix}.fractions.pca3d.html"

    fig.write_html(outfile, config={"responsive": True})

    print(f"\tWriting entropy interactive PCA plot to: {outfile}")

# -------------------------
# Mahalanobis-based LAI assignment
# -------------------------
def window_cluster_assignment(
    X_ref,
    X_sample,
    y,
    k,
    local_pcs=10,
    knn_mode="window",
    window_radius=3,
    require_full_window=False,
    kmeans_labels=None
):

    # Per-sample storage of predictions and probabilities
    preds_all = [[] for _ in range(X_sample.shape[0])]
    probs_all = [[] for _ in range(X_sample.shape[0])]

    cols_list = list(X_ref.columns)

    # Determine class labels (supervised vs unsupervised)
    if y is None:
        if kmeans_labels is None:
            raise ValueError("kmeans_labels must be provided when y is None")

        y = kmeans_labels
        classes = np.unique(y)
    else:
        classes = np.unique(y)

    # Iterate over genomic windows
    for i in range(len(cols_list)):

        # Window selection strategy
        if knn_mode == "raw":
            cols = [cols_list[i]]

        elif knn_mode == "window":
            start = max(0, i - window_radius)
            end = min(len(cols_list), i + window_radius + 1)
            cols = cols_list[start:end]

            # Optionally enforce full window size
            if require_full_window and len(cols) != (2 * window_radius + 1):
                continue
        else:
            raise ValueError(f"Unknown knn_mode: {knn_mode}")

        if len(cols) == 0:
            continue

        # Extract training and test data
        X_train = X_ref[cols].values
        X_test = X_sample[cols].values

        # Standardize
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)

        # PCA dimensionality reduction for numerical stability
        n_comp = min(local_pcs, X_train.shape[1], X_train.shape[0] - 1)
        pca = PCA(n_components=n_comp)
        Z_train = pca.fit_transform(X_train)
        Z_test = pca.transform(X_test)

        # Compute class-specific Gaussian parameters
        means = {}
        covs = {}

        for c in classes:
            Zc = Z_train[y == c]
            if len(Zc) < 2:
                continue

            mu = np.mean(Zc, axis=0)
            cov = np.cov(Zc, rowvar=False)

            # Ensure covariance is 2D
            if np.ndim(cov) == 0:
                cov = np.array([[cov]])
            elif cov.ndim == 1:
                cov = np.diag(cov)

            # Diagonal regularization (covariance shrinkage)
            cov += np.eye(cov.shape[0]) * 1e-3

            means[c] = mu
            covs[c] = np.linalg.inv(cov)

        # Classify each sample window
        for s in range(Z_test.shape[0]):
            z = Z_test[s]

            scores = []
            valid_classes = []

            # Compute Mahalanobis-based log-likelihood
            for c in classes:
                if c not in means:
                    continue

                diff = z - means[c]
                dist = diff @ covs[c] @ diff

                score = -0.5 * dist
                scores.append(score)
                valid_classes.append(c)

            scores = np.array(scores)
            if scores.size == 0:
               continue

            # Convert scores to probabilities (softmax)
            exp_scores = np.exp(scores - np.max(scores))
            probs = exp_scores / exp_scores.sum()

            # Expand to full class vector
            full_probs = np.zeros(len(classes))
            for idx, c in enumerate(valid_classes):
                full_probs[np.where(classes == c)[0][0]] = probs[idx]

            # Hard assignment
            pred = classes[np.argmax(full_probs)]

            preds_all[s].append(pred)
            probs_all[s].append(full_probs)

    # Convert classes to string labels
    # Convert predictions to string labels
    # Unsupervised mode
    if kmeans_labels is not None:
        classes = np.array([f"C{i+1}" for i in classes])

        for s in range(len(preds_all)):
            preds_all[s] = [f"C{p+1}" for p in preds_all[s]]

    # Supervised mode
    else:
        classes = np.array([str(i) for i in classes])

        for s in range(len(preds_all)):
            preds_all[s] = [str(p) for p in preds_all[s]]

    return preds_all, probs_all, classes

# -------------------------
# Aggregate probabilities into ancestry fractions
# -------------------------
def compute_fractions(preds_all, probs_all, classes):
    # Number of samples and classes
    n_samples = len(preds_all)
    n_classes = len(classes)

    # Initialize output matrices
    frac_matrix = np.zeros((n_samples, n_classes))  # hard assignment fractions
    prob_matrix = np.zeros((n_samples, n_classes))  # mean posterior probabilities

    for s in range(n_samples):
        # Skip samples with no predictions
        if len(preds_all[s]) == 0:
            continue

        # --- HARD COUNTS ---
        # Count occurrences of each class across windows
        for p in preds_all[s]:
            idx = np.where(classes == p)[0][0]
            frac_matrix[s, idx] += 1

        # Normalize to fractions
        frac_matrix[s] /= len(preds_all[s])

        # --- SOFT PROBABILITIES ---
        # Average posterior probabilities across windows
        probs_stack = np.array(probs_all[s])
        prob_matrix[s] = probs_stack.mean(axis=0)

    return frac_matrix, prob_matrix

# -------------------------
# Write per-window LAI outputs
# -------------------------
def write_lai_outputs(preds_all, probs_all, classes, coords, sample_names, outprefix, window_size):

    # Map numeric class labels to cluster-style labels (C1, C2, ...)
    class_to_cluster = {c: f"C{i+1}" for i, c in enumerate(classes)}

    # --- per-sample LAI ---
    for s_idx, sample in enumerate(sample_names):

        safe_sample = str(sample).replace("/", "-").replace(" ", "-")
        outfile = f"{safe_sample}_{outprefix}.lai.tsv"
        print(f"\tWriting LAI file: {outfile}")

        with open(outfile, "w") as f:
            # Header includes genomic coordinates, prediction, and per-class probabilities
            header = ["chrom", "start", "end", "ancestry_prediction"]
            header += [f"{c}_prob" for c in classes]
            f.write("\t".join(header) + "\n")

            # Iterate over genomic windows
            for i, (chr_, start_idx) in enumerate(zip(coords["chr"], coords["window_start"])):

                # Skip if predictions do not cover this window
                if i >= len(preds_all[s_idx]):
                    continue

                # Convert window index to genomic coordinates
                start = start_idx * window_size + 1
                end = start + window_size - 1

                pred = preds_all[s_idx][i]
                probs = probs_all[s_idx][i]

                # Combine cluster ID and class label
                cluster_label = class_to_cluster[pred]
                pred_label = pred 

                # Build output row
                row = [chr_, str(start), str(end), pred_label]
                row += [f"{p:.4f}" for p in probs]

                f.write("\t".join(row) + "\n")

# -------------------------
# HMM smoothing via Viterbi
# -------------------------
def viterbi_smooth(pred_probs, classes, stay_prob=0.9):
    """
    Minimal HMM Viterbi smoothing over window ancestry calls.

    pred_probs: list of arrays [window][class_probs]
    classes: np.array of class labels
    """

    n_states = len(classes)
    n_steps = len(pred_probs)

    # Transition matrix (log-space)
    # High probability of staying in same state, small probability of switching
    log_trans = np.full((n_states, n_states), np.log((1 - stay_prob) / (n_states - 1)))
    np.fill_diagonal(log_trans, np.log(stay_prob))

    # Dynamic programming table
    log_dp = np.zeros((n_steps, n_states))
    backptr = np.zeros((n_steps, n_states), dtype=int)

    # Initialization (first observation)
    log_dp[0] = np.log(pred_probs[0] + 1e-12)

    # Forward pass
    for t in range(1, n_steps):
        obs = np.log(pred_probs[t] + 1e-12)

        for j in range(n_states):
            # Evaluate all possible transitions into state j
            scores = log_dp[t - 1] + log_trans[:, j] + obs[j]
            backptr[t, j] = np.argmax(scores)
            log_dp[t, j] = np.max(scores)

    # Backtracking to recover optimal path
    path = np.zeros(n_steps, dtype=int)
    path[-1] = np.argmax(log_dp[-1])

    for t in range(n_steps - 2, -1, -1):
        path[t] = backptr[t + 1, path[t + 1]]

    return classes[path]

# -------------------------
# Write LAI summary (per sample)
# -------------------------
def write_lai_summary(preds_all, probs_all, classes, sample_names, outprefix):

    # Map class labels to cluster IDs
    class_to_cluster = {c: f"C{i+1}" for i, c in enumerate(classes)}

    outfile = f"{outprefix}.lai_summary.tsv"
    print(f"\tWriting summary: {outfile}")

    with open(outfile, "w") as f:

        # Write cluster-to-class mapping for interpretability
        f.write("--- Cluster Label Mapping ---\n")
        for c in classes:
            f.write(f"{class_to_cluster[c]}: {c}\n")
        f.write("-----------------------------\n\n")

        # Header: sample ID + cluster columns
        header = ["ID"] + [class_to_cluster[c] for c in classes]
        f.write("\t".join(header) + "\n")

        for s_idx, sample in enumerate(sample_names):

            preds = preds_all[s_idx]
            probs = np.array(probs_all[s_idx])

            if len(preds) == 0:
                continue

            row = [sample]

            for c in classes:
                idx = np.where(classes == c)[0][0]

                # Hard fraction
                frac = np.mean([p == c for p in preds])

                # Mean posterior probability
                prob = probs[:, idx].mean()

                row.append(f"{frac:.6f};p={prob:.6f}")

            f.write("\t".join(row) + "\n")

# -------------------------
# MAIN
# -------------------------
def main():

    parser = argparse.ArgumentParser(
        description="Ancestry inference with PCA + Mahalanobis window assignment"
    )

    parser.add_argument("--ref", required=True,
                        help="Reference matrix TSV")

    parser.add_argument("--sample", required=True,
                        help="Sample matrix TSV")

    parser.add_argument("--labels", required=False, default=None,
                        help="Discrete ancestry labels TSV")

    parser.add_argument("--fractions", required=False, default=None,
                        help="Reference ancestry fractions TSV")

    parser.add_argument("--c", type=int, default=5,
                        help="KMeans n_clusters")

    parser.add_argument("--outprefix", required=True,
                        help="Output prefix")

    parser.add_argument("--window_radius",
                        type=int,
                        default=3)

    parser.add_argument("--require_full_window",
                        action="store_true")

    parser.add_argument("--stay_prob",
                        type=float,
                        default=0.98,
                        help="HMM self-transition probability")

    parser.add_argument("--no_hmm",
                        action="store_true",
                        help="Disable HMM smoothing")

    parser.add_argument("--global_pcs",
                        type=int,
                        default=10,
                        help="Number of PCs used for unsupervised K-means clustering")

    parser.add_argument("--local_pcs",
                        type=int,
                        default=10,
                        help="Maximum number of PCs used for local Mahalanobis ancestry inference")

    parser.add_argument("-V", "--version",
                        action="version",
                        version=KORUSDECODE_VERSION)

    args = parser.parse_args()

    if args.fractions is not None and args.labels is not None:
        print("WARNING: --fractions provided; ignoring --labels")

    # -------------------------
    # smoothing logic
    # -------------------------
    if not args.no_hmm:
        args.window_radius = 1

    # -------------------------
    # logging
    # -------------------------
    sys.stdout = Tee(args.outprefix + ".log")

    # -------------------------
    # load matrices
    # -------------------------
    print("Loading data...")

    try:

        X_ref, coords_ref, window_size_ref = \
            load_reference_matrix(args.ref)

        X_sample, coords_sample, window_size_sample = \
            load_sample_matrix(args.sample)

    except FileNotFoundError as e:

        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    window_size = window_size_ref

    # -------------------------
    # load supervision
    # -------------------------
    labels = None
    frac_labels = None
    y_in = None

    if args.fractions is not None:
        # soft-supervised mode
        print("Using fractional ancestry supervision")

        frac_labels = load_fraction_labels(args.fractions)

        # align samples
        X_ref = X_ref.loc[frac_labels.index]

        frac_labels = frac_labels.loc[X_ref.index]

    elif args.labels is not None:
        # hard-supervised mode
        print("Using discrete ancestry supervision")

        labels = load_labels(args.labels)

        X_ref = X_ref.loc[labels.index]

        y_in = labels["ancestry"].values

    else:
        # unsupervised mode
        print("Using unsupervised mode")

    # -------------------------
    # align windows
    # -------------------------
    X_ref, X_sample = align_matrices(
        X_ref,
        coords_ref,
        X_sample,
        coords_sample
    )

    # -------------------------
    # assignment mode
    # -------------------------
    if frac_labels is not None:
        # soft-supervised mode
        print("Running fractional Mahalanobis assignment...")

        preds_all, probs_all, classes = \
            window_fractional_assignment(
                X_ref,
                X_sample,
                frac_labels,
                local_pcs=args.local_pcs,
                knn_mode="window",
                window_radius=args.window_radius,
                require_full_window=args.require_full_window
            )

        clusters = None

    elif labels is not None:
        # hard-supervised mode
        print("Running supervised Mahalanobis assignment...")

        preds_all, probs_all, classes = \
            window_cluster_assignment(
                X_ref,
                X_sample,
                y_in,
                args.c,
                local_pcs=args.local_pcs,
                knn_mode="window",
                window_radius=args.window_radius,
                require_full_window=args.require_full_window,
                kmeans_labels=None
            )

        clusters = None

    else:
        # unsupervised mode
        print(
            f"Running unsupervised KMeans mode "
            f"(k={args.c}, PCA dimensions={args.global_pcs})..."
        )

        X_ref_pca, _ = run_inference_pca(X_ref, X_sample, args.global_pcs)

        clusters, kmeans_model = run_kmeans(
            X_ref_pca,
            k=args.c
        )

        cluster_df = pd.DataFrame({
            "reference_id": X_ref.index,
            "cluster": [f"C{i+1}" for i in clusters]
        })

        cluster_out = args.outprefix + ".reference_clusters.tsv"

        cluster_df.to_csv(
            cluster_out,
            sep="\t",
            index=False
        )

        print(f"\tWriting reference clusters to: {cluster_out}")

        preds_all, probs_all, classes = \
            window_cluster_assignment(
                X_ref,
                X_sample,
                None,
                args.c,
                local_pcs=args.local_pcs,
                knn_mode="window",
                window_radius=args.window_radius,
                require_full_window=args.require_full_window,
                kmeans_labels=clusters
            )

    # -------------------------
    # optional HMM smoothing
    # -------------------------
    if not args.no_hmm:

        print("Running HMM smoothing...")

        for s in range(len(preds_all)):

            if len(probs_all[s]) == 0:
                continue

            smoothed = viterbi_smooth(
                probs_all[s],
                classes,
                stay_prob=args.stay_prob
            )

            preds_all[s] = list(smoothed)

    # -------------------------
    # aggregate fractions
    # -------------------------
    print("Aggregating fractions...")

    frac_matrix, prob_matrix = compute_fractions(
        preds_all,
        probs_all,
        classes
    )

    sample_names = list(X_sample.index)

    # -------------------------
    # write LAI outputs
    # -------------------------
    print("Writing LAI outputs...")

    write_lai_outputs(
        preds_all,
        probs_all,
        classes,
        coords_sample,
        sample_names,
        args.outprefix,
        window_size
    )

    write_lai_summary(
        preds_all,
        probs_all,
        classes,
        sample_names,
        args.outprefix
    )

    # -------------------------
    # output matrices
    # -------------------------
    frac_df = pd.DataFrame(
        frac_matrix,
        index=sample_names,
        columns=classes
    )

    # Compute entropy for coloring in --fractions mode
    if frac_labels is not None:
        entropy_vals = frac_labels.apply(compute_entropy, axis=1)
        entropy_vals = entropy_vals / np.log(len(classes))
    else:
        entropy_vals = None

    prob_df = pd.DataFrame(
        prob_matrix,
        index=sample_names,
        columns=classes
    )

    frac_out = args.outprefix + ".fractions.tsv"
    prob_out = args.outprefix + ".probabilities.tsv"

    frac_df.to_csv(frac_out, sep="\t")
    prob_df.to_csv(prob_out, sep="\t")

    print(f"\tWriting fractions to: {frac_out}")
    print(f"\tWriting probabilities to: {prob_out}")

    # -------------------------
    # PCA visualization
    # -------------------------
    print("Running PCA for visualization...")

    # Preparing the plot titles
    hmm_on = not args.no_hmm
    k = args.c if labels is None and frac_labels is None else None
    if frac_labels is not None:
        mode = "fractions"
    elif labels is not None:
        mode = "labels"
    else:
        mode = "unsupervised"
    title_2d = make_pca_title(mode, dim=2, k=k, hmm=hmm_on)
    title_3d = make_pca_title(mode, dim=3, k=k, hmm=hmm_on)
    title_3dint = make_pca_title(mode, dim=3, k=k, hmm=hmm_on, interactive=True)

    X_ref_pca, X_sample_pca = run_pca(
        X_ref,
        X_sample
    )

    if clusters is None:
        plot_clusters = np.full(X_ref.shape[0], -1, dtype=int)
    else:
        plot_clusters = clusters

    # -------------------------
    # Fraction-supervised plots
    # -------------------------
    if frac_labels is not None:

        plot_pca_fraction_2d(
            X_ref_pca,
            entropy_vals,
            X_sample_pca,
            sample_names,
            args.outprefix,
            title_2d
        )

        plot_pca_fraction_3d(
            X_ref,
            X_sample,
            entropy_vals,
            sample_names,
            args.outprefix,
            title_3d
        )

        plot_pca_fraction_3dint(
            X_ref,
            X_sample,
            entropy_vals,
            sample_names,
            args.outprefix,
            frac_df,
            title_3dint
        )

    # -------------------------
    # Standard plots
    # -------------------------
    else:

        plot_pca(
            X_ref_pca,
            labels,
            plot_clusters,
            X_sample_pca,
            sample_names,
            args.outprefix,
            title_2d
        )

        plot_pca_3d(
            X_ref,
            X_sample,
            labels,
            plot_clusters,
            sample_names,
            args.outprefix,
            title_3d
        )

        plot_pca_3dint(
            X_ref,
            X_sample,
            labels,
            plot_clusters,
            sample_names,
            args.outprefix,
            frac_df,
            title_3dint
        )

# -------------------------
# ENTRY POINT
# -------------------------
if __name__ == "__main__":
    main()
