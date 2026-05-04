#!/usr/bin/env python

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
    if path.endswith((".gz", ".bgz", ".gzip")):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return open(path, "r")

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
def plot_pca(X_ref_pca, labels, clusters, X_sample_pca, X_sample_names, outprefix):
    plt.figure(figsize=(8, 6))

    # Predefined ancestry color scheme
    color_map = {
        "AFR": "orange",
        "EUR": "green",
        "EAS": "red",
        "AMR": "black",
        "SAS": "blue",
        "NA": "grey"
    }

    # Determine ancestry categories
    if labels is not None:
        unique_anc = sorted(labels["ancestry"].unique())
    else:
        unique_anc = ["NA"]

    # Dynamically assign colors to unseen ancestries
    import matplotlib.colors as mcolors
    cmap = plt.get_cmap("tab20", len(unique_anc))

    dynamic_color_map = {}
    for i, anc in enumerate(unique_anc):
        if anc not in color_map:
            dynamic_color_map[anc] = mcolors.to_hex(cmap(i))

    # Plot reference samples
    if labels is not None:
        for i, sample in enumerate(labels.index):
            anc = labels.loc[sample, "ancestry"]

            if anc in color_map:
                color = color_map[anc]
            elif anc == "NA":
                color = color_map["NA"]
            else:
                color = dynamic_color_map.get(anc, "grey")

            plt.scatter(X_ref_pca[i, 0], X_ref_pca[i, 1], c=color, s=20, alpha=0.6)
            # Annotate with cluster ID
            plt.text(X_ref_pca[i, 0], X_ref_pca[i, 1], str(clusters[i] + 1), fontsize=6, alpha=0.5)
    else:
        # Unlabeled mode
        for i in range(X_ref_pca.shape[0]):
            plt.scatter(X_ref_pca[i, 0], X_ref_pca[i, 1], c="lightgrey", s=10, alpha=0.4)
            plt.text(X_ref_pca[i, 0], X_ref_pca[i, 1], str(clusters[i] + 1), fontsize=5, alpha=0.5)

    # Plot query samples
    for i, name in enumerate(X_sample_names):
        plt.scatter(X_sample_pca[i, 0], X_sample_pca[i, 1], c="black", marker="x", s=120)
        plt.text(X_sample_pca[i, 0], X_sample_pca[i, 1], name, fontsize=9, fontweight="bold")

    # Build legend
    for anc in unique_anc:
        if anc in color_map:
            c = color_map[anc]
        elif anc == "NA":
            c = color_map["NA"]
        else:
            c = dynamic_color_map.get(anc, "grey") 

        plt.scatter([], [], c=c, label=anc)

    plt.scatter([], [], c="black", marker="x", s=120, label="query")

    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("PCA projection with clusters")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{outprefix}.pca.png", dpi=300)
    print(f"\tWriting PCA plot to: {outprefix}.pca.png")

# -------------------------
# NEW: PCA Plot (3D)
# -------------------------
def plot_pca_3d(X_ref, X_sample, labels, clusters, sample_names, outprefix):
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

    # Same color logic as 2D plot
    color_map = {
        "AFR": "orange",
        "EUR": "green",
        "EAS": "red",
        "AMR": "black",
        "SAS": "blue",
        "NA": "grey"
    }

    if labels is not None:
        unique_anc = sorted(labels["ancestry"].unique())
    else:
        unique_anc = ["NA"]
    

    import matplotlib.colors as mcolors
    cmap = plt.get_cmap("tab20", len(unique_anc))

    dynamic_color_map = {}
    for i, anc in enumerate(unique_anc):
        if anc not in color_map:
            dynamic_color_map[anc] = mcolors.to_hex(cmap(i))

    # Plot reference points
    if labels is not None:
        for i, sample in enumerate(labels.index):
            anc = labels.loc[sample, "ancestry"]

            if anc in color_map:
                color = color_map[anc]
            elif anc == "NA":
                color = color_map["NA"]
            else:
                color = dynamic_color_map.get(anc, "grey")

            ax.scatter(X_ref_pca[i, 0], X_ref_pca[i, 1], X_ref_pca[i, 2],c=color, s=20, alpha=0.6)
            ax.text(X_ref_pca[i, 0], X_ref_pca[i, 1], X_ref_pca[i, 2], str(clusters[i] + 1), fontsize=6, alpha=0.5)
    else:
        for i in range(X_ref_pca.shape[0]):
            ax.scatter(X_ref_pca[i, 0], X_ref_pca[i, 1], X_ref_pca[i, 2],
                       c="lightgrey", s=10, alpha=0.4)
            ax.text(X_ref_pca[i, 0], X_ref_pca[i, 1], X_ref_pca[i, 2],
                    str(clusters[i] + 1), fontsize=5, alpha=0.5)

    # Plot query samples
    for i, name in enumerate(sample_names):
        ax.scatter(X_sample_pca[i, 0], X_sample_pca[i, 1], X_sample_pca[i, 2], c="black", marker="x", s=120)
        ax.text(X_sample_pca[i, 0], X_sample_pca[i, 1], X_sample_pca[i, 2], name, fontsize=9, fontweight="bold")

    # Legend
    for anc in unique_anc:
        if anc in color_map:
            c = color_map[anc]
        elif anc == "NA":
            c = color_map["NA"]
        else:
            c = dynamic_color_map.get(anc, "grey") 

        ax.scatter([], [], [], c=c, label=anc)

    ax.scatter([], [], [], c="black", marker="x", s=120, label="query")

    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_zlabel("PC3")
    ax.set_title("3D PCA projection with clusters")
    ax.legend(loc='center left', bbox_to_anchor=(1.07, 0.5))

    plt.tight_layout(rect=[0, 0, 0.8, 1])
    plt.savefig(f"{outprefix}.pca3d.png", dpi=300)
    print(f"\tWriting 3D PCA plot to: {outprefix}.pca3d.png")

# -------------------------
# Mahalanobis-based LAI assignment (REPLACES KNN)
# -------------------------
def window_cluster_assignment(X_ref, X_sample, y, k, knn_mode="window", window_radius=3, require_full_window=False, kmeans_labels=None):

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
        n_comp = min(10, X_train.shape[1], X_train.shape[0] - 1)
        #n_comp = max(1, min(10, X_train.shape[1], X_train.shape[0] - 1))
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
            #eps = 1e-3 * np.trace(cov) / cov.shape[0]
            #cov += np.eye(cov.shape[0]) * eps

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
    classes = np.array([f"{i}" for i in classes])

    # Convert predictions to string labels
    for s in range(len(preds_all)):
        preds_all[s] = [f"{p}" for p in preds_all[s]]

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

        outfile = f"{outprefix}.{sample}.lai.tsv"
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
                pred_label = f"{cluster_label}/{pred}"

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
    # Argument parsing
    parser = argparse.ArgumentParser(description="Ancestry inference with PCA + Mahalanobis window assignment")

    parser.add_argument("--ref", required=True, help="Reference matrix")
    parser.add_argument("--sample", required=True, help="Sample matrix")
    parser.add_argument("--labels", required=False, help="Reference labels", default=None)
    parser.add_argument("--k", type=int, help="KMeans n_clusters", default=5)
    parser.add_argument("--outprefix", required=True, help="Output prefix")

    parser.add_argument("--window_radius", type=int, default=3)
    parser.add_argument("--require_full_window", action="store_true")

    parser.add_argument("--stay_prob", type=float, default=0.99,help="HMM self-transition probability (smoothing strength)")
    parser.add_argument("--no_hmm", action="store_true", help="Disable HMM smoothing")

    args = parser.parse_args()

    # dual smoothing approach, but tune window_radius down if hmm enabled
    if not args.no_hmm:
        args.window_radius = 1 # dials down non-parametric smoother

    # Redirect stdout to both terminal and log file
    sys.stdout = Tee(args.outprefix + ".log")

    print("Loading data...")
    X_ref, coords_ref, window_size_ref = load_reference_matrix(args.ref)
    X_sample, coords_sample, window_size_sample = load_sample_matrix(args.sample)

    # Assume same window size for ref and sample
    window_size = window_size_ref

    # Load labels if provided (supervised mode)
    if args.labels is not None:
        labels = load_labels(args.labels)
        # Align reference samples to label order
        X_ref = X_ref.loc[labels.index]
        y = labels["ancestry"].values
    else:
        labels = None
        y = None

    # Align genomic windows between reference and sample
    X_ref, X_sample = align_matrices(X_ref, coords_ref, X_sample, coords_sample)

    print("Preparing labels / clusters...")
    if labels is not None:
        # Supervised mode
        y_in = labels["ancestry"].values
        clusters = None
    else:
        # Unsupervised mode: generate clusters from PCA
        print("Running KMeans (unsupervised mode)...")
        X_ref_pca, _ = run_pca(X_ref, X_sample)   # PCA for clustering stability
        clusters, kmeans_model = run_kmeans(X_ref_pca, k=args.k)
        y_in = None

    print("Running Mahalanobis window assignment...")

    preds_all, probs_all, classes = window_cluster_assignment(
        X_ref,
        X_sample,
        y_in,
        args.k,
        knn_mode="window",
        window_radius=args.window_radius,
        require_full_window=args.require_full_window,
        kmeans_labels=clusters
    )

    # -------------------------
    # HMM smoothing (optional)
    # -------------------------
    if not args.no_hmm:
        print("Running HMM smoothing...")
        for s in range(len(preds_all)):
            if len(probs_all[s]) == 0:
                continue

            # Apply Viterbi smoothing to per-window predictions
            smoothed = viterbi_smooth(
                probs_all[s],
                classes,
                stay_prob=args.stay_prob
            )

            preds_all[s] = list(smoothed)

    print("Aggregating fractions...")
    frac_matrix, prob_matrix = compute_fractions(preds_all, probs_all, classes)

    sample_names = list(X_sample.index)

    print("Running window-level assignment (LAI-style)...")

    # Write detailed per-window outputs
    write_lai_outputs(
        preds_all,
        probs_all,
        classes,
        coords_sample,
        sample_names,
        args.outprefix,
        window_size
    )

    # Write summary table
    write_lai_summary(
        preds_all,
        probs_all,
        classes,
        sample_names,
        args.outprefix
    )

    # Convert results to DataFrames
    frac_df = pd.DataFrame(frac_matrix, index=sample_names, columns=classes)
    prob_df = pd.DataFrame(prob_matrix, index=sample_names, columns=classes)

    frac_out = args.outprefix + ".fractions.tsv"
    prob_out = args.outprefix + ".probabilities.tsv"

    # Save outputs
    frac_df.to_csv(frac_out, sep="\t")
    prob_df.to_csv(prob_out, sep="\t")

    print(f"\tWriting fractions to: {frac_out}")
    print(f"\tWriting probabilities to: {prob_out}")

    print("Running PCA for visualization...")
    X_ref_pca, X_sample_pca = run_pca(X_ref, X_sample)

    # Recompute clusters for plotting (independent of assignment phase)
    clusters, kmeans_model = run_kmeans(X_ref_pca, k=args.k)

    # Generate plots
    plot_pca(X_ref_pca, labels, clusters, X_sample_pca, sample_names, args.outprefix)
    plot_pca_3d(X_ref, X_sample, labels, clusters, sample_names, args.outprefix)

# -------------------------
# ENTRY POINT
# -------------------------
if __name__ == "__main__":
    main()
