"""
Compare strain perturbations across probe sets using series-level real-vs-random metrics.

Expected parent folder layout:

Parent folder/
├── yWL333/
│   ├── MS2(ATP6)/
│   │   └── pooled_real_vs_random_comparisons/
│   ├── MS2(ATP8)/
│   │   └── pooled_real_vs_random_comparisons/
│   ├── MS2(ATP2)/
│   │   └── pooled_real_vs_random_comparisons/
│   ├── MS2(ATP3)/
│   │   └── pooled_real_vs_random_comparisons/
│   └── MS2(TIM50)/
│       └── pooled_real_vs_random_comparisons/
├── yMM002(ATP11)/
│   ├── MS2(ATP6_8)/
│   │   └── pooled_real_vs_random_comparisons/
│   ├── ATP2/
│   │   └── pooled_real_vs_random_comparisons/
│   ├── ATP3/
│   │   └── pooled_real_vs_random_comparisons/
│   └── TIM50/
│       └── pooled_real_vs_random_comparisons/
└── ...

The script also supports a repeated inner strain folder, e.g.
Parent/strain/probe/strain/pooled_real_vs_random_comparisons/.

Each pooled_real_vs_random_comparisons folder must contain:
    pooled_input_file_pairing_summary.csv

Main output:
    strain_probe_delta_median_grouped_boxplot.png

The plot uses one dot per series and a box/whisker summary per strain within each probe set.
ATP6, ATP8, and row labels such as MS2 are automatically combined into one ATP6/8 probe set.

Important behavior for mixed probe folders:
If a folder is named something like MS2 (ATP2), the script DOES NOT force all rows
from that folder into ATP2. Instead, it checks each row in the pairing CSV first.
For example, rows labeled mRNA = MS2 are assigned to ATP6/8, while rows labeled
mRNA = ATP2 are assigned to ATP2. Folder names are used only as a fallback.
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import ks_2samp, mannwhitneyu, kruskal
from tkinter import filedialog
import tkinter as tk


# ==========================================================
# Parent folder picker.
# Select the parent folder that contains the strain folders.
# ==========================================================
root = tk.Tk()
root.withdraw()
parent_folder = Path(filedialog.askdirectory(title="Select parent folder containing strain folders"))
print("Selected parent folder:", parent_folder)


# ==========================================================
# User-adjustable settings
# ==========================================================
output_folder_name = "strain_probe_comparison"

# Main metric for box/whisker plots and strain/probe comparisons.
# Options include:
#   "delta_median" = median(real) - median(random)
#   "delta_mean"   = mean(real) - mean(random)
#   "ratio_median" = median(real) / median(random)
#   "ks_statistic" = KS distance between real and random distributions
primary_metric = "delta_median"

# Distances below these thresholds will also be summarized.
close_distance_thresholds_um = [0.25, 0.5, 1.0]

# Display order for the probe sets. Any unrecognized probe set is appended after these.
probe_set_order = ["ATP6/8", "ATP2", "ATP3", "TIM50"]

# Optional manual renaming. Useful if a folder name should be displayed differently.
# Example: strain_name_overrides = {"yMM002(ATP11)": "yMM002 ATP11"}
strain_name_overrides = {}
probe_name_overrides = {}

# Display options for strain labels/order.
wt_strain_names = {"yWL333"}

# Plot options.
make_individual_probe_plots = True
# Cleaner default output set: keep the main plot and individual probe plots.
make_alternate_strain_grouped_plot = False
make_summary_heatmap = False
save_svg_copies = True
random_jitter_seed = 123

# Aesthetic options
base_font_size = 13
axis_label_font_size = 14
title_font_size = 16
tick_font_size = 12
legend_font_size = 11
main_xtick_rotation = 0
strain_xtick_rotation = 25
dot_size = 46

# Statistical annotation options.
annotate_wt_vs_ko = True
show_nonsignificant_annotations = False
minimum_replicates_for_stats = 2
# Options: "none", "fdr_bh_by_probe", or "fdr_bh_global".
wt_vs_ko_pvalue_correction = "fdr_bh_by_probe"
# Stars are based on adjusted p-values when correction is enabled.


# ==========================================================
# Helpers
# ==========================================================
def format_strain_label(name):
    name = str(name).strip()

    if name in strain_name_overrides:
        return str(strain_name_overrides[name])

    if name in wt_strain_names:
        return "WT"

    m = re.search(r"\(([^)]+)\)", name)
    if m:
        gene = m.group(1).strip().lower()
        if gene:
            return f"Δ{gene}"

    return name


def strain_sort_key(name):
    name = str(name).strip()
    if name in wt_strain_names:
        return (0, natural_key(name))
    return (1, natural_key(name))


def metric_axis_label(metric):
    labels = {
        "delta_median": "Δ median distance",
        "delta_mean": "Δ mean distance",
        "ratio_median": "Median ratio",
        "ratio_mean": "Mean ratio",
        "ks_statistic": "KS statistic",
    }
    return labels.get(metric, metric.replace("_", " "))


def sanitize_name(name):
    name = str(name).strip()
    if name == "":
        name = "unknown"
    for bad in ['\\', '/', ':', '*', '?', '"', '<', '>', '|']:
        name = name.replace(bad, "_")
    name = re.sub(r"\s+", "_", name)
    return name


def natural_key(text):
    """Sort strings naturally, e.g. strain2 before strain10."""
    parts = re.split(r"(\d+)", str(text))
    return [int(p) if p.isdigit() else p.casefold() for p in parts]


def clean_label(text, remove_ms2=True):
    """Normalize labels for robust matching while preserving biology terms."""
    text = str(text).strip()
    text = text.replace("＋", "+")
    if remove_ms2:
        text = re.sub(r"(?i)\bMS2\b", "", text)
    text = text.replace("(", " ").replace(")", " ")
    text = text.replace("[", " ").replace("]", " ")
    text = text.replace("-", " ").replace("_", " ").replace("/", " ").replace("+", " ")
    text = re.sub(r"\s+", " ", text)
    return text.strip().upper()


def normalize_probe_set(*candidate_labels, ms2_means_atp6_8=False, allow_fallback_label=True):
    """
    Convert labels into one of the intended probe sets.

    Rules:
      - ATP6 and ATP8 are deliberately collapsed into ATP6/8.
      - If ms2_means_atp6_8=True, a label that contains MS2 but no explicit
        ATP2/ATP3/TIM50 is assigned to ATP6/8. This is used for row-level
        pairing-summary labels where MS2 marks the ATP6/8 probe.
      - Folder names are used only as fallbacks, so mixed folders such as
        MS2 (ATP2) can contribute both MS2/ATP6/8 rows and ATP2 rows.
    """
    labels = [str(x) for x in candidate_labels if pd.notna(x) and str(x).strip() != ""]
    joined_raw = " ".join(labels)

    # Manual overrides are checked on raw labels first.
    for label in labels:
        if str(label) in probe_name_overrides:
            return probe_name_overrides[str(label)]

    cleaned_keep_ms2 = clean_label(joined_raw, remove_ms2=False)
    cleaned_no_ms2 = clean_label(joined_raw, remove_ms2=True)

    has_atp6_or_8 = bool(re.search(r"\bATP6\b", cleaned_keep_ms2) or re.search(r"\bATP8\b", cleaned_keep_ms2))
    has_atp2 = bool(re.search(r"\bATP2\b", cleaned_keep_ms2))
    has_atp3 = bool(re.search(r"\bATP3\b", cleaned_keep_ms2))
    has_tim50 = bool(re.search(r"\bTIM50\b", cleaned_keep_ms2))
    has_ms2 = bool(re.search(r"\bMS2\b", cleaned_keep_ms2))

    if has_atp6_or_8:
        return "ATP6/8"

    # In this dataset, row-level MS2 is the ATP6/8 marker.
    # Do this before ATP2/ATP3/TIM50 only when the caller is intentionally
    # passing a specific row label, not a mixed folder label.
    if ms2_means_atp6_8 and has_ms2 and not (has_atp2 or has_atp3 or has_tim50):
        return "ATP6/8"

    if has_atp2:
        return "ATP2"
    if has_atp3:
        return "ATP3"
    if has_tim50:
        return "TIM50"

    if allow_fallback_label and labels:
        fallback = cleaned_no_ms2 if cleaned_no_ms2 else cleaned_keep_ms2
        return fallback if fallback else "UNKNOWN"

    return "UNKNOWN"


def path_basename_any_platform(raw_path):
    parts = split_path_parts_any_platform(raw_path)
    return parts[-1] if parts else str(raw_path)


def infer_probe_set_for_row(row, folder_probe_set):
    """
    Assign ProbeSet for a single pairing-summary row.

    Priority:
      1. mRNA column, because this identifies the actual molecule/probe for the row.
      2. File basenames, in case mRNA is generic but filenames contain MS2/ATP2/etc.
      3. SeriesFolder basename.
      4. Folder-derived probe set as a fallback.

    This is the key fix for mixed folders such as MS2 (ATP2): an MS2 row becomes
    ATP6/8, while an ATP2 row from the same folder remains ATP2.
    """
    checks = [
        ("mRNA", row.get("mRNA", "")),
        ("RealFile basename", path_basename_any_platform(row.get("RealFile", ""))),
        ("RandomFile basename", path_basename_any_platform(row.get("RandomFile", ""))),
        ("SeriesFolder basename", path_basename_any_platform(row.get("SeriesFolder", ""))),
    ]

    for source, label in checks:
        probe = normalize_probe_set(label, ms2_means_atp6_8=True, allow_fallback_label=False)
        if probe != "UNKNOWN":
            return probe, source

    return folder_probe_set, "folder fallback"


def infer_strain_and_probe(condition_folder, experiment_root):
    """
    Infer strain and probe folder from a folder that contains pooled outputs.

    This script supports both of these layouts:

        experiment_root / strain / probe / pooled_real_vs_random_comparisons
        experiment_root / strain / probe / repeated_strain / pooled_real_vs_random_comparisons

    In both cases, the correct strain and probe are the FIRST TWO folders
    under the selected parent folder. This is important for datasets where
    each probe folder contains another strain-named folder before the pooled
    output folder, for example:

        Aim 2 / yMM002 (atp11) / MS2 (ATP2) / yMM002 (atp11) / pooled...
    """
    try:
        rel_parts = condition_folder.relative_to(experiment_root).parts
    except ValueError:
        rel_parts = condition_folder.parts

    if len(rel_parts) >= 2:
        # Correct for nested layouts with a repeated strain folder.
        strain = rel_parts[0]
        probe_folder = rel_parts[1]
    elif len(rel_parts) == 1:
        # Directly under root. Usually this is a strain-level folder and will
        # be skipped later if it lacks a valid pairing CSV.
        strain = rel_parts[0]
        probe_folder = rel_parts[0]
    else:
        strain = condition_folder.parent.name
        probe_folder = condition_folder.name

    strain = strain_name_overrides.get(strain, strain)
    probe_set = normalize_probe_set(probe_folder)
    return strain, probe_folder, probe_set


def find_condition_folders(root_folder):
    """Find folders containing pooled real-vs-random comparison outputs."""
    pooled_dirs = [
        p for p in root_folder.rglob("pooled_real_vs_random_comparisons")
        if p.is_dir()
    ]
    return sorted({p.parent for p in pooled_dirs}, key=lambda p: natural_key(str(p)))


def load_clean_npy(path):
    arr = np.load(path)
    arr = np.asarray(arr, dtype=float).ravel()
    arr = arr[np.isfinite(arr)]
    return arr


def read_pairing_summary(condition_folder):
    pairing_file = (
        condition_folder
        / "pooled_real_vs_random_comparisons"
        / "pooled_input_file_pairing_summary.csv"
    )

    if not pairing_file.exists():
        return None

    df = pd.read_csv(pairing_file)
    df.columns = [str(c).strip() for c in df.columns]

    required = {"mRNA", "SeriesFolder", "RealFile", "RandomFile", "Status"}
    if not required.issubset(set(df.columns)):
        print(f"  Pairing file is missing required columns: {pairing_file}")
        print(f"  Required columns: {sorted(required)}")
        print(f"  Found columns:    {list(df.columns)}")
        return None

    return df


def split_path_parts_any_platform(raw_path):
    """Split Windows or POSIX paths into comparable path parts."""
    raw = str(raw_path).strip().strip('"').strip("'")
    return [part for part in re.split(r"[\\/]+", raw) if part not in {"", "."}]


def build_npy_file_index(experiment_root):
    """Index .npy files by filename for robust recovery from stale absolute paths."""
    index = {}
    for path in experiment_root.rglob("*.npy"):
        index.setdefault(path.name, []).append(path)
    return index


def score_candidate_path(candidate, raw_parts, condition_folder):
    """Prefer candidates sharing more trailing path parts with the saved path."""
    cand_parts = [p.casefold() for p in candidate.parts]
    raw_parts = [p.casefold() for p in raw_parts]

    score = 0
    for c, r in zip(reversed(cand_parts), reversed(raw_parts)):
        if c == r:
            score += 1
        else:
            break

    # Prefer files inside the current condition subtree if there is a tie.
    try:
        candidate.relative_to(condition_folder)
        score += 100
    except ValueError:
        pass

    return score


def resolve_file_path(raw_path, condition_folder, experiment_root, file_index=None):
    """
    Resolve NPY paths stored as absolute paths or paths relative to likely roots.

    The pairing summary sometimes contains absolute paths from an earlier run.
    If those paths no longer exist, this function falls back to a filename/suffix
    search under the selected parent folder.
    """
    raw = str(raw_path).strip().strip('"').strip("'")
    path = Path(raw)

    if path.exists():
        return path

    candidates = []
    if not path.is_absolute():
        candidates.extend([
            condition_folder / path,
            condition_folder / "pooled_real_vs_random_comparisons" / path,
            experiment_root / path,
        ])

    raw_parts = split_path_parts_any_platform(raw)

    # If the saved path contains the selected parent folder name, rebuild the
    # suffix below that folder under the currently selected parent.
    parent_name = experiment_root.name.casefold()
    raw_parts_casefold = [part.casefold() for part in raw_parts]
    if parent_name in raw_parts_casefold:
        idx = len(raw_parts_casefold) - 1 - list(reversed(raw_parts_casefold)).index(parent_name)
        suffix_parts = raw_parts[idx + 1:]
        if suffix_parts:
            candidates.append(experiment_root.joinpath(*suffix_parts))

    for candidate in candidates:
        if candidate.exists():
            return candidate

    # Last-resort recovery: same filename somewhere under the selected parent.
    if file_index is not None and raw_parts:
        filename = raw_parts[-1]
        same_name = file_index.get(filename, [])
        if same_name:
            best = max(same_name, key=lambda p: score_candidate_path(p, raw_parts, condition_folder))
            if best.exists():
                return best

    # Return the original path so the missing-file status is informative.
    return path


def compute_metrics_for_pair(strain, probe_folder, folder_probe_set, condition_folder, row, experiment_root, file_index=None):
    mrna = str(row["mRNA"])
    series_folder = str(row["SeriesFolder"])

    probe_set, probe_assignment_source = infer_probe_set_for_row(row, folder_probe_set)

    real_file = resolve_file_path(row["RealFile"], condition_folder, experiment_root, file_index=file_index)
    random_file = resolve_file_path(row["RandomFile"], condition_folder, experiment_root, file_index=file_index)

    base = {
        "Strain": strain,
        "ProbeFolder": probe_folder,
        "FolderProbeSet": folder_probe_set,
        "ProbeSet": probe_set,
        "ProbeAssignmentSource": probe_assignment_source,
        "mRNA": mrna,
        "SeriesFolder": series_folder,
        "ConditionFolder": str(condition_folder),
        "RealFile": str(real_file),
        "RandomFile": str(random_file),
        "Condition": f"{strain} | {probe_set}",
    }

    if not real_file.exists():
        return {**base, "Status": "missing_real_file"}

    if not random_file.exists():
        return {**base, "Status": "missing_random_file"}

    real = load_clean_npy(real_file)
    random = load_clean_npy(random_file)

    if len(real) == 0 or len(random) == 0:
        return {
            **base,
            "Status": "empty_real_or_random",
            "real_n": len(real),
            "random_n": len(random),
        }

    real_mean = float(np.mean(real))
    random_mean = float(np.mean(random))
    real_median = float(np.median(real))
    random_median = float(np.median(random))

    out = {
        **base,
        "Status": "processed",
        "real_n": int(len(real)),
        "random_n": int(len(random)),
        "real_mean": real_mean,
        "random_mean": random_mean,
        "delta_mean": real_mean - random_mean,
        "real_median": real_median,
        "random_median": random_median,
        "delta_median": real_median - random_median,
        "ratio_mean": real_mean / random_mean if random_mean != 0 else np.nan,
        "ratio_median": real_median / random_median if random_median != 0 else np.nan,
        "real_std": float(np.std(real, ddof=1)) if len(real) > 1 else np.nan,
        "random_std": float(np.std(random, ddof=1)) if len(random) > 1 else np.nan,
    }

    try:
        ks = ks_2samp(real, random)
        out["ks_statistic"] = float(ks.statistic)
        out["ks_pvalue_distance_level"] = float(ks.pvalue)
    except Exception:
        out["ks_statistic"] = np.nan
        out["ks_pvalue_distance_level"] = np.nan

    try:
        mw = mannwhitneyu(real, random, alternative="two-sided")
        out["mannwhitney_pvalue_distance_level"] = float(mw.pvalue)
    except Exception:
        out["mannwhitney_pvalue_distance_level"] = np.nan

    for threshold in close_distance_thresholds_um:
        key = str(threshold).replace(".", "p")
        real_frac = float(np.mean(real <= threshold))
        random_frac = float(np.mean(random <= threshold))

        out[f"real_fraction_le_{key}_um"] = real_frac
        out[f"random_fraction_le_{key}_um"] = random_frac
        out[f"delta_fraction_le_{key}_um"] = real_frac - random_frac

    return out


def ordered_unique(values, preferred_order=None):
    values = [v for v in pd.Series(values).dropna().unique()]
    preferred_order = preferred_order or []
    ordered = [v for v in preferred_order if v in values]
    ordered.extend(sorted([v for v in values if v not in ordered], key=natural_key))
    return ordered


def summarize_by_strain_probe(series_metrics, metric):
    valid = series_metrics[series_metrics["Status"] == "processed"].copy()
    if len(valid) == 0:
        return pd.DataFrame()

    metric_cols = [
        "real_n",
        "random_n",
        "real_mean",
        "random_mean",
        "delta_mean",
        "real_median",
        "random_median",
        "delta_median",
        "ratio_mean",
        "ratio_median",
        "ks_statistic",
    ]

    for threshold in close_distance_thresholds_um:
        key = str(threshold).replace(".", "p")
        metric_cols.extend([
            f"real_fraction_le_{key}_um",
            f"random_fraction_le_{key}_um",
            f"delta_fraction_le_{key}_um",
        ])

    rows = []
    for (strain, probe_set), sub in valid.groupby(["Strain", "ProbeSet"], sort=False):
        row = {
            "Strain": strain,
            "ProbeSet": probe_set,
            "ProbeFoldersCombined": "; ".join(sorted(sub["ProbeFolder"].astype(str).unique(), key=natural_key)),
            "mRNAsIncluded": "; ".join(sorted(sub["mRNA"].astype(str).unique(), key=natural_key)),
            "n_series": int(len(sub)),
            "total_real_distances": int(sub["real_n"].sum()),
            "total_random_distances": int(sub["random_n"].sum()),
        }

        vals_for_metric = pd.to_numeric(sub[metric], errors="coerce").dropna()
        row[f"{metric}_mean_across_series"] = float(vals_for_metric.mean()) if len(vals_for_metric) else np.nan
        row[f"{metric}_median_across_series"] = float(vals_for_metric.median()) if len(vals_for_metric) else np.nan
        row[f"{metric}_sem_across_series"] = (
            float(vals_for_metric.std(ddof=1) / np.sqrt(len(vals_for_metric)))
            if len(vals_for_metric) > 1 else np.nan
        )

        for col in metric_cols:
            if col not in sub.columns:
                continue
            vals = pd.to_numeric(sub[col], errors="coerce").dropna()
            if len(vals) == 0:
                row[f"{col}_mean_across_series"] = np.nan
                row[f"{col}_median_across_series"] = np.nan
                row[f"{col}_sem_across_series"] = np.nan
            else:
                row[f"{col}_mean_across_series"] = float(vals.mean())
                row[f"{col}_median_across_series"] = float(vals.median())
                row[f"{col}_sem_across_series"] = (
                    float(vals.std(ddof=1) / np.sqrt(len(vals))) if len(vals) > 1 else np.nan
                )

        rows.append(row)

    out = pd.DataFrame(rows)
    if len(out):
        ordered_probes = ordered_unique(out["ProbeSet"], probe_set_order)
        ordered_strains = sorted(out["Strain"].unique(), key=strain_sort_key)
        probe_rank = {probe: i for i, probe in enumerate(ordered_probes)}
        strain_rank = {strain: i for i, strain in enumerate(ordered_strains)}
        out["_probe_rank"] = out["ProbeSet"].map(probe_rank)
        out["_strain_rank"] = out["Strain"].map(strain_rank)
        out = out.sort_values(["_probe_rank", "_strain_rank"]).drop(columns=["_probe_rank", "_strain_rank"])
    return out


def statistics_by_probe_across_strains(series_metrics, metric):
    """For each probe set, compare strains using series-level metric values as replicates."""
    valid = series_metrics[series_metrics["Status"] == "processed"].copy()
    valid[metric] = pd.to_numeric(valid[metric], errors="coerce")
    valid = valid.dropna(subset=[metric])

    rows = []
    for probe_set in ordered_unique(valid["ProbeSet"], probe_set_order):
        sub_probe = valid[valid["ProbeSet"] == probe_set]
        groups = []
        labels = []

        for strain in sorted(sub_probe["Strain"].unique(), key=strain_sort_key):
            vals = sub_probe.loc[sub_probe["Strain"] == strain, metric].dropna().values
            if len(vals) > 0:
                labels.append(strain)
                groups.append(vals)

        if len(groups) < 2:
            continue

        try:
            kw = kruskal(*groups)
            rows.append({
                "Scope": "within_probe_set_across_strains",
                "ProbeSet": probe_set,
                "Strain": "",
                "comparison": "all_strains",
                "metric": metric,
                "test": "Kruskal-Wallis on series-level metric",
                "group_a": "",
                "group_b": "",
                "n_a": "",
                "n_b": "",
                "statistic": float(kw.statistic),
                "pvalue": float(kw.pvalue),
            })
        except Exception as e:
            rows.append({
                "Scope": "within_probe_set_across_strains",
                "ProbeSet": probe_set,
                "Strain": "",
                "comparison": "all_strains",
                "metric": metric,
                "test": f"Kruskal-Wallis failed: {e}",
            })

        for i in range(len(groups)):
            for j in range(i + 1, len(groups)):
                a = groups[i]
                b = groups[j]
                label_a = labels[i]
                label_b = labels[j]
                try:
                    test = mannwhitneyu(a, b, alternative="two-sided")
                    rows.append({
                        "Scope": "within_probe_set_across_strains",
                        "ProbeSet": probe_set,
                        "Strain": "",
                        "comparison": f"{label_a} vs {label_b}",
                        "metric": metric,
                        "test": "Mann-Whitney U on series-level metric",
                        "group_a": label_a,
                        "group_b": label_b,
                        "n_a": int(len(a)),
                        "n_b": int(len(b)),
                        "median_a": float(np.median(a)),
                        "median_b": float(np.median(b)),
                        "difference_median_a_minus_b": float(np.median(a) - np.median(b)),
                        "statistic": float(test.statistic),
                        "pvalue": float(test.pvalue),
                    })
                except Exception as e:
                    rows.append({
                        "Scope": "within_probe_set_across_strains",
                        "ProbeSet": probe_set,
                        "Strain": "",
                        "comparison": f"{label_a} vs {label_b}",
                        "metric": metric,
                        "test": f"Mann-Whitney failed: {e}",
                        "group_a": label_a,
                        "group_b": label_b,
                        "n_a": int(len(a)),
                        "n_b": int(len(b)),
                    })

    return pd.DataFrame(rows)


def statistics_by_strain_across_probes(series_metrics, metric):
    """For each strain, compare probe sets using series-level metric values as replicates."""
    valid = series_metrics[series_metrics["Status"] == "processed"].copy()
    valid[metric] = pd.to_numeric(valid[metric], errors="coerce")
    valid = valid.dropna(subset=[metric])

    rows = []
    for strain in sorted(valid["Strain"].unique(), key=strain_sort_key):
        sub_strain = valid[valid["Strain"] == strain]
        groups = []
        labels = []

        for probe_set in ordered_unique(sub_strain["ProbeSet"], probe_set_order):
            vals = sub_strain.loc[sub_strain["ProbeSet"] == probe_set, metric].dropna().values
            if len(vals) > 0:
                labels.append(probe_set)
                groups.append(vals)

        if len(groups) < 2:
            continue

        try:
            kw = kruskal(*groups)
            rows.append({
                "Scope": "within_strain_across_probe_sets",
                "ProbeSet": "",
                "Strain": strain,
                "comparison": "all_probe_sets",
                "metric": metric,
                "test": "Kruskal-Wallis on series-level metric",
                "group_a": "",
                "group_b": "",
                "n_a": "",
                "n_b": "",
                "statistic": float(kw.statistic),
                "pvalue": float(kw.pvalue),
            })
        except Exception as e:
            rows.append({
                "Scope": "within_strain_across_probe_sets",
                "ProbeSet": "",
                "Strain": strain,
                "comparison": "all_probe_sets",
                "metric": metric,
                "test": f"Kruskal-Wallis failed: {e}",
            })

        for i in range(len(groups)):
            for j in range(i + 1, len(groups)):
                a = groups[i]
                b = groups[j]
                label_a = labels[i]
                label_b = labels[j]
                try:
                    test = mannwhitneyu(a, b, alternative="two-sided")
                    rows.append({
                        "Scope": "within_strain_across_probe_sets",
                        "ProbeSet": "",
                        "Strain": strain,
                        "comparison": f"{label_a} vs {label_b}",
                        "metric": metric,
                        "test": "Mann-Whitney U on series-level metric",
                        "group_a": label_a,
                        "group_b": label_b,
                        "n_a": int(len(a)),
                        "n_b": int(len(b)),
                        "median_a": float(np.median(a)),
                        "median_b": float(np.median(b)),
                        "difference_median_a_minus_b": float(np.median(a) - np.median(b)),
                        "statistic": float(test.statistic),
                        "pvalue": float(test.pvalue),
                    })
                except Exception as e:
                    rows.append({
                        "Scope": "within_strain_across_probe_sets",
                        "ProbeSet": "",
                        "Strain": strain,
                        "comparison": f"{label_a} vs {label_b}",
                        "metric": metric,
                        "test": f"Mann-Whitney failed: {e}",
                        "group_a": label_a,
                        "group_b": label_b,
                        "n_a": int(len(a)),
                        "n_b": int(len(b)),
                    })

    return pd.DataFrame(rows)


def benjamini_hochberg_adjust(pvalues):
    """Return Benjamini-Hochberg adjusted p-values, preserving NaNs."""
    pvalues = np.asarray(pvalues, dtype=float)
    adjusted = np.full_like(pvalues, np.nan, dtype=float)
    finite_mask = np.isfinite(pvalues)
    finite_p = pvalues[finite_mask]

    if len(finite_p) == 0:
        return adjusted

    order = np.argsort(finite_p)
    ranked = finite_p[order]
    n = len(ranked)
    ranked_adj = ranked * n / np.arange(1, n + 1)
    ranked_adj = np.minimum.accumulate(ranked_adj[::-1])[::-1]
    ranked_adj = np.clip(ranked_adj, 0, 1)

    tmp = np.empty_like(ranked_adj)
    tmp[order] = ranked_adj
    adjusted[finite_mask] = tmp
    return adjusted


def pvalue_to_stars(pvalue):
    if not np.isfinite(pvalue):
        return ""
    if pvalue < 0.0001:
        return "****"
    if pvalue < 0.001:
        return "***"
    if pvalue < 0.01:
        return "**"
    if pvalue < 0.05:
        return "*"
    return "ns"


def wt_vs_ko_statistics(series_metrics, metric):
    """Compare WT against each KO within each probe set using series-level values."""
    valid = series_metrics[series_metrics["Status"] == "processed"].copy()
    if metric not in valid.columns:
        return pd.DataFrame()

    valid[metric] = pd.to_numeric(valid[metric], errors="coerce")
    valid = valid.dropna(subset=[metric])

    rows = []
    for probe_set in ordered_unique(valid["ProbeSet"], probe_set_order):
        sub_probe = valid[valid["ProbeSet"] == probe_set].copy()
        if len(sub_probe) == 0:
            continue

        wt_strains_present = [s for s in sorted(sub_probe["Strain"].unique(), key=strain_sort_key) if s in wt_strain_names]
        if len(wt_strains_present) == 0:
            continue

        wt_vals = sub_probe.loc[sub_probe["Strain"].isin(wt_strains_present), metric].dropna().values
        wt_label = "; ".join(format_strain_label(s) for s in wt_strains_present)

        for ko_strain in sorted(sub_probe["Strain"].unique(), key=strain_sort_key):
            if ko_strain in wt_strain_names:
                continue

            ko_vals = sub_probe.loc[sub_probe["Strain"] == ko_strain, metric].dropna().values
            row = {
                "ProbeSet": probe_set,
                "WTStrain": "; ".join(wt_strains_present),
                "WTLabel": wt_label,
                "KOStrain": ko_strain,
                "KOLabel": format_strain_label(ko_strain),
                "comparison": f"{wt_label} vs {format_strain_label(ko_strain)}",
                "metric": metric,
                "test": "Mann-Whitney U on series-level metric",
                "n_wt": int(len(wt_vals)),
                "n_ko": int(len(ko_vals)),
                "median_wt": float(np.median(wt_vals)) if len(wt_vals) else np.nan,
                "median_ko": float(np.median(ko_vals)) if len(ko_vals) else np.nan,
                "difference_median_ko_minus_wt": (
                    float(np.median(ko_vals) - np.median(wt_vals))
                    if len(wt_vals) and len(ko_vals) else np.nan
                ),
                "pvalue_raw": np.nan,
                "pvalue_adjusted": np.nan,
                "significance": "",
                "note": "",
            }

            if len(wt_vals) < minimum_replicates_for_stats or len(ko_vals) < minimum_replicates_for_stats:
                row["note"] = f"Skipped: needs at least {minimum_replicates_for_stats} series per group"
            else:
                try:
                    test = mannwhitneyu(wt_vals, ko_vals, alternative="two-sided")
                    row["statistic"] = float(test.statistic)
                    row["pvalue_raw"] = float(test.pvalue)
                except Exception as e:
                    row["note"] = f"Mann-Whitney failed: {e}"

            rows.append(row)

    stats = pd.DataFrame(rows)
    if stats.empty:
        return stats

    if wt_vs_ko_pvalue_correction == "none":
        stats["pvalue_adjusted"] = stats["pvalue_raw"]
    elif wt_vs_ko_pvalue_correction == "fdr_bh_global":
        stats["pvalue_adjusted"] = benjamini_hochberg_adjust(stats["pvalue_raw"].values)
    else:
        for probe_set, idx in stats.groupby("ProbeSet").groups.items():
            stats.loc[idx, "pvalue_adjusted"] = benjamini_hochberg_adjust(stats.loc[idx, "pvalue_raw"].values)

    stats["significance"] = stats["pvalue_adjusted"].apply(pvalue_to_stars)
    stats["pvalue_correction"] = wt_vs_ko_pvalue_correction
    return stats


def add_wt_vs_ko_annotations(ax, probe_set, strains, values_by_strain, stats_df):
    """Add significance brackets comparing WT to each KO on a single-probe plot."""
    if stats_df is None or stats_df.empty:
        return

    sub_stats = stats_df[stats_df["ProbeSet"] == probe_set].copy()
    if sub_stats.empty:
        return

    strain_to_x = {strain: i for i, strain in enumerate(strains)}
    wt_positions = [strain_to_x[s] for s in strains if s in wt_strain_names]
    if not wt_positions:
        return
    wt_x = wt_positions[0]

    all_vals = np.concatenate([np.asarray(v, dtype=float) for v in values_by_strain if len(v) > 0])
    if len(all_vals) == 0:
        return

    y_min = float(np.nanmin(all_vals))
    y_max = float(np.nanmax(all_vals))
    y_range = y_max - y_min
    if not np.isfinite(y_range) or y_range == 0:
        y_range = 1.0

    bracket_height = 0.025 * y_range
    row_step = 0.085 * y_range
    next_y = y_max + 0.08 * y_range
    drawn = 0

    for _, row in sub_stats.iterrows():
        ko_strain = row.get("KOStrain", "")
        if ko_strain not in strain_to_x:
            continue

        star = row.get("significance", "")
        if star == "ns" and not show_nonsignificant_annotations:
            continue
        if star == "":
            continue

        ko_x = strain_to_x[ko_strain]
        x1, x2 = sorted([wt_x, ko_x])
        y = next_y + drawn * row_step

        ax.plot([x1, x1, x2, x2], [y, y + bracket_height, y + bracket_height, y],
                color="black", linewidth=1.1, clip_on=False)
        ax.text((x1 + x2) / 2, y + bracket_height + 0.01 * y_range, star,
                ha="center", va="bottom", fontsize=11, color="black", clip_on=False)
        drawn += 1

    if drawn > 0:
        ax.set_ylim(top=next_y + drawn * row_step + 0.08 * y_range)


def save_current_figure(output_path):
    plt.savefig(output_path, dpi=300, facecolor="white", bbox_inches="tight")
    if save_svg_copies:
        svg_path = output_path.with_suffix(".svg")
        plt.savefig(svg_path, facecolor="white", bbox_inches="tight")


def color_map_for(labels, mode="generic"):
    labels = list(labels)
    if mode == "strain":
        wt_color = "#4D4D4D"
        mutant_labels = [lab for lab in labels if str(lab).strip() not in wt_strain_names]
        cmap = plt.get_cmap("tab10")
        out = {}
        for lab in labels:
            if str(lab).strip() in wt_strain_names:
                out[lab] = wt_color
        for i, lab in enumerate(mutant_labels):
            out[lab] = cmap(i % cmap.N)
        return out

    if mode == "probe":
        fixed = {
            "ATP6/8": "#7A5195",
            "ATP2": "#EF5675",
            "ATP3": "#3690C0",
            "TIM50": "#FFA600",
        }
        cmap = plt.get_cmap("Set2")
        out = {}
        extra = [lab for lab in labels if lab not in fixed]
        for lab in labels:
            if lab in fixed:
                out[lab] = fixed[lab]
        for i, lab in enumerate(extra):
            out[lab] = cmap(i % cmap.N)
        return out

    cmap = plt.get_cmap("tab20")
    return {label: cmap(i % cmap.N) for i, label in enumerate(labels)}


def apply_plot_style():
    plt.rcParams.update({
        "font.size": base_font_size,
        "axes.titlesize": title_font_size,
        "axes.labelsize": axis_label_font_size,
        "xtick.labelsize": tick_font_size,
        "ytick.labelsize": tick_font_size,
        "legend.fontsize": legend_font_size,
        "legend.title_fontsize": legend_font_size,
    })


def style_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.25)



def plot_grouped_boxes_by_probe(series_metrics, metric, output_folder):
    """
    Main plot: x-axis = probe set, boxes = strains, dots = individual series.
    This is optimized for asking whether strains differ within each probe set.
    """
    valid = series_metrics[series_metrics["Status"] == "processed"].copy()
    if metric not in valid.columns:
        print(f"Metric not found for plotting: {metric}")
        return None

    valid[metric] = pd.to_numeric(valid[metric], errors="coerce")
    valid = valid.dropna(subset=[metric])

    if len(valid) == 0:
        print(f"No processed rows with numeric {metric}; plot not made.")
        return None

    probe_sets = ordered_unique(valid["ProbeSet"], probe_set_order)
    strains = sorted(valid["Strain"].unique(), key=strain_sort_key)
    colors = color_map_for(strains, mode="strain")

    n_strains = max(len(strains), 1)
    group_width = 0.82
    box_width = min(0.11, group_width / (n_strains + 1.5))
    offsets = np.linspace(-group_width / 2, group_width / 2, n_strains)
    rng = np.random.default_rng(random_jitter_seed)

    fig_width = max(10, len(probe_sets) * max(2.4, 0.55 * n_strains))
    apply_plot_style()
    apply_plot_style()
    plt.figure(figsize=(fig_width, 7), facecolor="white")
    ax = plt.gca()

    for strain_index, strain in enumerate(strains):
        positions = []
        box_data = []
        for probe_index, probe_set in enumerate(probe_sets):
            vals = valid.loc[
                (valid["ProbeSet"] == probe_set) & (valid["Strain"] == strain),
                metric,
            ].dropna().values
            if len(vals) == 0:
                continue
            pos = probe_index + offsets[strain_index]
            positions.append(pos)
            box_data.append(vals)

            jitter = rng.normal(0, box_width * 0.20, size=len(vals))
            ax.scatter(
                np.full(len(vals), pos) + jitter,
                vals,
                s=dot_size,
                alpha=0.85,
                color=colors[strain],
                edgecolor="black",
                linewidth=0.4,
                zorder=3,
            )

        if len(box_data) == 0:
            continue

        bp = ax.boxplot(
            box_data,
            positions=positions,
            widths=box_width,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "black", "linewidth": 1.5},
            whiskerprops={"color": "black", "linewidth": 1.0},
            capprops={"color": "black", "linewidth": 1.0},
            boxprops={"color": "black", "linewidth": 1.0},
        )
        for patch in bp["boxes"]:
            patch.set_facecolor(colors[strain])
            patch.set_alpha(0.35)

    ax.axhline(0, color="black", linewidth=1, linestyle="--", alpha=0.8)
    ax.set_xticks(range(len(probe_sets)))
    ax.set_xticklabels(probe_sets, rotation=main_xtick_rotation)
    ax.set_ylabel(metric_axis_label(metric))
    ax.set_xlabel("Probe set")
    ax.set_title("All probes")

    legend_handles = [Patch(facecolor=colors[s], edgecolor="black", alpha=0.35, label=format_strain_label(s)) for s in strains]
    ax.legend(handles=legend_handles, title="Strain", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)

    style_axes(ax)
    plt.tight_layout()

    output_path = output_folder / f"strain_probe_{metric}_grouped_boxplot.png"
    save_current_figure(output_path)
    plt.close()
    return output_path


def plot_grouped_boxes_by_strain(series_metrics, metric, output_folder):
    """
    Alternate plot: x-axis = strain, boxes = probe sets.
    This is optimized for seeing a perturbation profile within each strain.
    """
    valid = series_metrics[series_metrics["Status"] == "processed"].copy()
    if metric not in valid.columns:
        return None

    valid[metric] = pd.to_numeric(valid[metric], errors="coerce")
    valid = valid.dropna(subset=[metric])
    if len(valid) == 0:
        return None

    strains = sorted(valid["Strain"].unique(), key=strain_sort_key)
    probe_sets = ordered_unique(valid["ProbeSet"], probe_set_order)
    colors = color_map_for(probe_sets, mode="probe")

    n_probes = max(len(probe_sets), 1)
    group_width = 0.82
    box_width = min(0.13, group_width / (n_probes + 1.5))
    offsets = np.linspace(-group_width / 2, group_width / 2, n_probes)
    rng = np.random.default_rng(random_jitter_seed)

    fig_width = max(11, len(strains) * max(1.8, 0.50 * n_probes))
    plt.figure(figsize=(fig_width, 7), facecolor="white")
    ax = plt.gca()

    for probe_index, probe_set in enumerate(probe_sets):
        positions = []
        box_data = []
        for strain_index, strain in enumerate(strains):
            vals = valid.loc[
                (valid["ProbeSet"] == probe_set) & (valid["Strain"] == strain),
                metric,
            ].dropna().values
            if len(vals) == 0:
                continue
            pos = strain_index + offsets[probe_index]
            positions.append(pos)
            box_data.append(vals)

            jitter = rng.normal(0, box_width * 0.20, size=len(vals))
            ax.scatter(
                np.full(len(vals), pos) + jitter,
                vals,
                s=dot_size,
                alpha=0.85,
                color=colors[probe_set],
                edgecolor="black",
                linewidth=0.4,
                zorder=3,
            )

        if len(box_data) == 0:
            continue

        bp = ax.boxplot(
            box_data,
            positions=positions,
            widths=box_width,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "black", "linewidth": 1.5},
            whiskerprops={"color": "black", "linewidth": 1.0},
            capprops={"color": "black", "linewidth": 1.0},
            boxprops={"color": "black", "linewidth": 1.0},
        )
        for patch in bp["boxes"]:
            patch.set_facecolor(colors[probe_set])
            patch.set_alpha(0.35)

    ax.axhline(0, color="black", linewidth=1, linestyle="--", alpha=0.8)
    ax.set_xticks(range(len(strains)))
    ax.set_xticklabels([format_strain_label(s) for s in strains], rotation=strain_xtick_rotation, ha="right")
    ax.set_ylabel(metric_axis_label(metric))
    ax.set_xlabel("Strain")
    ax.set_title("All probes")

    legend_handles = [Patch(facecolor=colors[p], edgecolor="black", alpha=0.35, label=p) for p in probe_sets]
    ax.legend(handles=legend_handles, title="Probe set", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)

    style_axes(ax)
    plt.tight_layout()

    output_path = output_folder / f"strain_probe_{metric}_grouped_by_strain_boxplot.png"
    save_current_figure(output_path)
    plt.close()
    return output_path


def plot_individual_probe_plots(series_metrics, metric, output_folder, wt_stats=None):
    valid = series_metrics[series_metrics["Status"] == "processed"].copy()
    if metric not in valid.columns:
        return []

    valid[metric] = pd.to_numeric(valid[metric], errors="coerce")
    valid = valid.dropna(subset=[metric])
    if len(valid) == 0:
        return []

    output_paths = []
    rng = np.random.default_rng(random_jitter_seed)

    for probe_set in ordered_unique(valid["ProbeSet"], probe_set_order):
        sub = valid[valid["ProbeSet"] == probe_set].copy()
        if len(sub) == 0:
            continue

        strains = sorted(sub["Strain"].unique(), key=strain_sort_key)
        data = [sub.loc[sub["Strain"] == strain, metric].dropna().values for strain in strains]
        positions = np.arange(len(strains))

        apply_plot_style()
        plt.figure(figsize=(max(8, len(strains) * 1.25), 6), facecolor="white")
        ax = plt.gca()

        bp = ax.boxplot(
            data,
            positions=positions,
            widths=0.55,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "black", "linewidth": 1.5},
            whiskerprops={"color": "black", "linewidth": 1.0},
            capprops={"color": "black", "linewidth": 1.0},
            boxprops={"color": "black", "linewidth": 1.0},
        )
        for patch in bp["boxes"]:
            patch.set_alpha(0.25)

        for idx, vals in enumerate(data):
            jitter = rng.normal(0, 0.055, size=len(vals))
            ax.scatter(
                np.full(len(vals), idx) + jitter,
                vals,
                s=dot_size,
                alpha=0.85,
                edgecolor="black",
                linewidth=0.4,
                zorder=3,
            )

        ax.axhline(0, color="black", linewidth=1, linestyle="--", alpha=0.8)
        ax.set_xticks(positions)
        ax.set_xticklabels([format_strain_label(s) for s in strains], rotation=strain_xtick_rotation, ha="right")
        ax.set_ylabel(metric_axis_label(metric))
        ax.set_title(str(probe_set))
        style_axes(ax)
        if annotate_wt_vs_ko:
            add_wt_vs_ko_annotations(ax, probe_set, strains, data, wt_stats)
        plt.tight_layout()

        output_path = output_folder / f"{sanitize_name(probe_set)}_{metric}_by_strain_boxplot.png"
        save_current_figure(output_path)
        plt.close()
        output_paths.append(output_path)

    return output_paths


def plot_summary_heatmap(summary_df, metric, output_folder):
    if summary_df.empty:
        return None

    value_col = f"{metric}_median_across_series"
    if value_col not in summary_df.columns:
        return None

    strains = sorted(summary_df["Strain"].unique(), key=strain_sort_key)
    probe_sets = ordered_unique(summary_df["ProbeSet"], probe_set_order)
    mat = np.full((len(probe_sets), len(strains)), np.nan)

    for i, probe_set in enumerate(probe_sets):
        for j, strain in enumerate(strains):
            vals = summary_df.loc[
                (summary_df["ProbeSet"] == probe_set) & (summary_df["Strain"] == strain),
                value_col,
            ].values
            if len(vals):
                mat[i, j] = vals[0]

    apply_plot_style()
    plt.figure(figsize=(max(8, len(strains) * 1.1), max(4.5, len(probe_sets) * 0.8)), facecolor="white")
    ax = plt.gca()
    im = ax.imshow(mat, aspect="auto")
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(metric_axis_label(metric))

    ax.set_xticks(np.arange(len(strains)))
    ax.set_xticklabels([format_strain_label(s) for s in strains], rotation=strain_xtick_rotation, ha="right")
    ax.set_yticks(np.arange(len(probe_sets)))
    ax.set_yticklabels(probe_sets)
    ax.set_title("Summary")

    for i in range(len(probe_sets)):
        for j in range(len(strains)):
            if np.isfinite(mat[i, j]):
                ax.text(j, i, f"{mat[i, j]:.3g}", ha="center", va="center", fontsize=8)

    plt.tight_layout()
    output_path = output_folder / f"strain_probe_{metric}_summary_heatmap.png"
    save_current_figure(output_path)
    plt.close()
    return output_path


# ==========================================================
# Main
# ==========================================================
def main():
    if str(parent_folder).strip() == "." or str(parent_folder).strip() == "":
        print("No parent folder selected. Exiting.")
        return

    if not parent_folder.exists():
        raise FileNotFoundError(f"Parent folder does not exist: {parent_folder}")

    output_folder = parent_folder / output_folder_name
    output_folder.mkdir(parents=True, exist_ok=True)

    condition_folders = find_condition_folders(parent_folder)

    print(f"\nParent folder: {parent_folder}")
    print(f"Found {len(condition_folders)} folders with pooled outputs")

    if len(condition_folders) == 0:
        print("\nNo folders found.")
        print("Expected each probe folder to contain:")
        print("  pooled_real_vs_random_comparisons/pooled_input_file_pairing_summary.csv")
        return

    print("Indexing .npy files under the selected parent folder for robust path recovery...")
    file_index = build_npy_file_index(parent_folder)
    print(f"Indexed {sum(len(v) for v in file_index.values())} .npy files")

    all_rows = []
    folder_map_rows = []

    for condition_folder in condition_folders:
        strain, probe_folder, probe_set = infer_strain_and_probe(condition_folder, parent_folder)
        folder_map_rows.append({
            "ConditionFolder": str(condition_folder),
            "Strain": strain,
            "ProbeFolder": probe_folder,
            "ProbeSet": probe_set,
        })

        print(f"\nStrain:    {strain}")
        print(f"Probe dir: {probe_folder}")
        print(f"Probe set: {probe_set}")
        print(f"Folder:    {condition_folder}")

        pairing = read_pairing_summary(condition_folder)
        if pairing is None:
            print("  Missing or invalid pooled_input_file_pairing_summary.csv; skipping.")
            continue

        pooled_rows = pairing[pairing["Status"].astype(str).str.lower() == "pooled"].copy()
        print(f"  Pooled real/random file pairs: {len(pooled_rows)}")

        for _, row in pooled_rows.iterrows():
            metrics = compute_metrics_for_pair(
                strain=strain,
                probe_folder=probe_folder,
                folder_probe_set=probe_set,
                condition_folder=condition_folder,
                row=row,
                experiment_root=parent_folder,
                file_index=file_index,
            )
            all_rows.append(metrics)

    folder_map = pd.DataFrame(folder_map_rows)
    folder_map_path = output_folder / "strain_probe_folder_mapping.csv"
    folder_map.to_csv(folder_map_path, index=False)

    if len(all_rows) == 0:
        print("\nNo real/random pairs were found.")
        print(f"Folder mapping written to: {folder_map_path}")
        return

    series_metrics = pd.DataFrame(all_rows)
    series_metrics_path = output_folder / "strain_probe_series_level_metrics.csv"
    series_metrics.to_csv(series_metrics_path, index=False)

    processed = series_metrics[series_metrics["Status"] == "processed"].copy()
    print(f"\nProcessed series-level rows: {len(processed)}")
    status_counts = series_metrics["Status"].value_counts(dropna=False)
    print("Status counts:")
    for status, count in status_counts.items():
        print(f"  {status}: {count}")

    if len(processed) > 0:
        print("Processed rows by probe set:")
        for probe_set, count in processed["ProbeSet"].value_counts().sort_index().items():
            print(f"  {probe_set}: {count}")

        print("Probe assignment sources:")
        for source, count in processed["ProbeAssignmentSource"].value_counts().items():
            print(f"  {source}: {count}")

        mrna_probe_counts = (
            processed.groupby(["ProbeSet", "mRNA"])
            .size()
            .reset_index(name="n_rows")
            .sort_values(["ProbeSet", "mRNA"], key=lambda col: col.map(str))
        )
        mrna_probe_counts_path = output_folder / "probe_set_mRNA_assignment_counts.csv"
        mrna_probe_counts.to_csv(mrna_probe_counts_path, index=False)
        print(f"Probe/mRNA assignment counts: {mrna_probe_counts_path}")

    if len(processed) == 0:
        print("No valid processed rows. Check missing-file statuses in:")
        print(f"  {series_metrics_path}")
        print("Common causes: stale absolute NPY paths in the pairing CSV, renamed folders, or missing NPY files.")
        return

    summary_df = summarize_by_strain_probe(series_metrics, primary_metric)
    summary_path = output_folder / "strain_probe_summary.csv"
    summary_df.to_csv(summary_path, index=False)

    stats_by_probe = statistics_by_probe_across_strains(series_metrics, primary_metric)
    stats_by_probe_path = output_folder / f"statistics_by_probe_across_strains_{primary_metric}.csv"
    stats_by_probe.to_csv(stats_by_probe_path, index=False)

    stats_by_strain = statistics_by_strain_across_probes(series_metrics, primary_metric)
    stats_by_strain_path = output_folder / f"statistics_by_strain_across_probes_{primary_metric}.csv"
    stats_by_strain.to_csv(stats_by_strain_path, index=False)

    wt_stats = wt_vs_ko_statistics(series_metrics, primary_metric)
    wt_stats_path = output_folder / f"wt_vs_ko_statistics_{primary_metric}.csv"
    wt_stats.to_csv(wt_stats_path, index=False)

    plot_paths = []
    main_plot = plot_grouped_boxes_by_probe(series_metrics, primary_metric, output_folder)
    if main_plot is not None:
        plot_paths.append(main_plot)

    if make_alternate_strain_grouped_plot:
        alt_plot = plot_grouped_boxes_by_strain(series_metrics, primary_metric, output_folder)
        if alt_plot is not None:
            plot_paths.append(alt_plot)

    if make_individual_probe_plots:
        plot_paths.extend(plot_individual_probe_plots(series_metrics, primary_metric, output_folder, wt_stats=wt_stats))

    if make_summary_heatmap:
        heatmap = plot_summary_heatmap(summary_df, primary_metric, output_folder)
        if heatmap is not None:
            plot_paths.append(heatmap)

    print("\nDone.")
    print(f"Folder mapping:       {folder_map_path}")
    print(f"Series-level metrics: {series_metrics_path}")
    print(f"Strain/probe summary: {summary_path}")
    print(f"Stats by probe:       {stats_by_probe_path}")
    print(f"Stats by strain:      {stats_by_strain_path}")
    print(f"WT vs KO stats:       {wt_stats_path}")
    print("Plots:")
    for path in plot_paths:
        print(f"  {path}")
        if save_svg_copies:
            print(f"  {path.with_suffix('.svg')}")

    print("\nInterpretation note:")
    print("Each dot is one pooled series-level real-vs-random comparison.")
    print("The box summarizes those series-level values; distance-level values are not treated as independent biological replicates.")
    print("ATP6, ATP8, and row-level MS2 labels are combined into the ATP6/8 probe set before plotting and statistics.")
    print("WT-vs-KO stars on individual probe plots use Mann-Whitney U tests on series-level values with the configured p-value correction.")
    print("For mixed folders such as MS2 (ATP2), each pairing-summary row is assigned from its own mRNA/file label first; the folder name is only a fallback.")


if __name__ == "__main__":
    main()
