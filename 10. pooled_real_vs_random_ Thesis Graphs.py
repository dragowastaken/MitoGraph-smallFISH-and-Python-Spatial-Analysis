import re
import tkinter as tk
from pathlib import Path
from tkinter import filedialog

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from scipy.stats import ks_2samp, mannwhitneyu, wilcoxon

    HAS_SCIPY = True
except Exception:
    HAS_SCIPY = False


# ==========================================================
# Parent folder representing ONE experimental condition.
# All Series folders below this parent will be pooled.
# ==========================================================
root = tk.Tk()
root.withdraw()

parent_folder = Path(filedialog.askdirectory())

print("User selected:", parent_folder)


# ==========================================================
# Plot/settings
# ==========================================================
BIN_WIDTH_UM = 0.05
MAX_DISTANCE_UM = 2.0

# Diffraction limit floor:
# Any real/random NN distance below this value is pooled into this value.
DIFFRACTION_LIMIT_UM = 0.20

# Practical effect threshold:
# This does not change the p-value. It helps interpret whether a statistically
# significant result is large enough to be biologically meaningful.
# Set to 0.10 um so small-but-statistically-significant shifts do not get stars.
PRACTICAL_EFFECT_THRESHOLD_UM = 0.10

# Strict histogram inference gate:
# The histogram remains the main figure, but stars/brackets are only drawn when
# ALL of the following pass. This prevents large pooled spot-level n from driving
# overconfident visual conclusions.
STRICT_ALPHA = 0.001
STRICT_EFFECT_THRESHOLD_UM = PRACTICAL_EFFECT_THRESHOLD_UM
STRICT_MIN_REPLICATE_PAIRS = 3
STRICT_MIN_DIRECTION_AGREEMENT = 0.75
STRICT_REQUIRE_CI_EXCLUDES_ZERO = True

# If True, a positive effect requires the 95% CI lower bound to be >= +threshold,
# and a negative effect requires the 95% CI upper bound to be <= -threshold.
# This is ON by default so the entire confidence interval must clear the effect threshold.
STRICT_REQUIRE_CI_EXCLUDES_EFFECT_THRESHOLD = True

# Main statistic shown on plot:
# Use "replicate_median_permutation" to avoid over-sensitivity from pooled spot-level n.
MAIN_PLOT_TEST = "replicate_median_permutation"

# Number of random sign-flips for paired permutation test.
N_PERMUTATIONS = 20000
RANDOM_SEED = 12345

REAL_COLORS = {
    "MS2": "red",
    "ATP2": "orange",  # dark orange
    "ATP3": "green",  # dark green
    "TIM50": "purple",  # dark magenta / purple
}

RANDOM_COLOR = "#56B4E9"  # sky blue

# Display labels only. Internal filenames/summaries can still use MS2.
DISPLAY_NAME_OVERRIDES = {
    "MS2": "ATP6/8",
}


# ==========================================================
# Output folder, created inside parent_folder
# ==========================================================
POOLED_OUTPUT_FOLDER_NAME = "pooled_real_vs_random_comparisons"

REAL_ALPHA = 0.65
RANDOM_LINEWIDTH = 3.0

# Descriptive plot annotations.
# These show the real-vs-random effect directly on every histogram without
# depending on a binary significance threshold.
SHOW_MEDIAN_LINES = True
SHOW_DELTA_MEDIAN_ARROW = True
SHOW_SHIFT_DIRECTION_LABEL = False

# Keep strict stars optional. The strict call is still saved in the CSV.
# Set to True only if you want brackets/stars on plots that pass the strict gate.
SHOW_STRICT_SIGNIFICANCE_STARS = False

# Regular dashed median lines.
MEDIAN_LINESTYLE = "--"
MEDIAN_LINEWIDTH = 2.5
DELTA_ARROW_FONTSIZE = 16


def sanitize_name(name):
    name = str(name).strip()
    if name == "":
        name = "mRNA"
    for bad in ["\\", "/", ":", "*", "?", '"', "<", ">", "|"]:
        name = name.replace(bad, "_")
    name = re.sub(r"\s+", "_", name)
    return name


def move_column_to_end(df, column, value=None):
    """Add/update a column and move it to the end of a DataFrame."""
    df = df.copy()
    if value is not None:
        df[column] = value
    if column not in df.columns:
        return df
    return df[[c for c in df.columns if c != column] + [column]]


def output_csv_name(base_stem, strain_name):
    """Append the strain/condition name to CSV filenames."""
    strain_suffix = sanitize_name(strain_name)
    return f"{base_stem}_{strain_suffix}.csv"


def display_mrna_name(mrna):
    raw = str(mrna).strip()
    upper = raw.upper()

    if upper in DISPLAY_NAME_OVERRIDES:
        return DISPLAY_NAME_OVERRIDES[upper]

    return raw


def extract_series_number_from_name(name):
    m = re.search(r"series\s*0*(\d+)", str(name), flags=re.IGNORECASE)
    if m:
        return int(m.group(1))
    return None


def make_series_range_label(file_records):
    nums = []

    for rec in file_records:
        series_folder = rec.get("SeriesFolder", "")
        if series_folder:
            n = extract_series_number_from_name(Path(series_folder).name)
            if n is not None:
                nums.append(n)

    nums = sorted(set(nums))

    if len(nums) == 0:
        return ""

    if len(nums) == 1:
        return f"Series {nums[0]}"

    expected = list(range(nums[0], nums[-1] + 1))

    if nums == expected:
        return f"Series {nums[0]}-{nums[-1]}"

    return "Series " + ", ".join(str(n) for n in nums)


def make_plot_title(strain_name, mrna, file_records):
    display_name = display_mrna_name(mrna)
    series_label = make_series_range_label(file_records)

    if series_label:
        return f"{display_name}"

    return f"{display_name}"


def infer_mrna_name(path):
    """
    Infer MS2 / ATP2 / ATP3 / TIM50 from filename or path.
    """
    known = ["MS2", "ATP2", "ATP3", "TIM50"]
    p = Path(path)

    pieces = list(p.parts)
    pieces.append(p.name)

    for part in reversed(pieces):
        upper = part.upper()
        for name in known:
            if re.search(rf"(^|[_\-\s]){re.escape(name)}([_\-\s]|$)", upper):
                return name

    # fallback from filename like ABC_NN_distance_um.npy
    stem = p.stem
    stem = re.sub(r"_NN_distance_um$", "", stem, flags=re.IGNORECASE)
    stem = re.sub(r"^random_", "", stem, flags=re.IGNORECASE)
    stem = re.sub(r"_nn$", "", stem, flags=re.IGNORECASE)
    return sanitize_name(stem)


def find_series_folder(path):
    p = Path(path)
    for ancestor in p.parents:
        if ancestor.name.lower().startswith("series"):
            return ancestor
    return None


def is_random_or_generated_file(path):
    p = Path(path)
    lower_parts = [part.lower() for part in p.parts]
    name = p.name.lower()

    if name.startswith("random_"):
        return True
    if any(
        part.startswith("random_") and part.endswith("_output") for part in lower_parts
    ):
        return True
    if POOLED_OUTPUT_FOLDER_NAME.lower() in lower_parts:
        return True

    return False


def find_real_distance_files(parent):
    files = []
    for p in parent.rglob("*_NN_distance_um.npy"):
        if is_random_or_generated_file(p):
            continue
        files.append(p)
    return sorted(files)


def find_random_file_for_real(real_file, mrna):
    real_file = Path(real_file)
    spots_folder = real_file.parent
    safe_mrna = sanitize_name(mrna)

    candidates = [
        spots_folder / f"random_{safe_mrna}_output" / f"random_{safe_mrna}_nn.npy",
        spots_folder
        / f"random_{safe_mrna.upper()}_output"
        / f"random_{safe_mrna.upper()}_nn.npy",
        spots_folder
        / f"random_{safe_mrna.lower()}_output"
        / f"random_{safe_mrna.lower()}_nn.npy",
    ]

    for c in candidates:
        if c.exists():
            return c

    target = f"random_{safe_mrna}_nn.npy".lower()
    for p in spots_folder.rglob("*.npy"):
        if p.name.lower() == target:
            return p

    return None


def load_clean_npy(path):
    arr = np.load(path)
    arr = np.asarray(arr, dtype=float).ravel()
    arr = arr[np.isfinite(arr)]
    return arr


def apply_diffraction_floor(arr):
    arr = np.asarray(arr, dtype=float).copy()
    arr[arr < DIFFRACTION_LIMIT_UM] = DIFFRACTION_LIMIT_UM
    return arr


def format_p_value(p):
    if pd.isna(p):
        return "p = NA"
    if p < 1e-4:
        return "p < 1e-4"
    return f"p = {p:.4g}"


def significance_stars(p):
    if pd.isna(p):
        return "n/a"
    if p < 0.0001:
        return "****"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def is_significant(p):
    return (not pd.isna(p)) and p < 0.05


def passes_p_value_gate(p, alpha=STRICT_ALPHA):
    return (not pd.isna(p)) and p < alpha


def direction_agreement_fraction(differences, effect):
    """Fraction of replicate-level differences pointing in the effect direction."""
    differences = np.asarray(differences, dtype=float)
    differences = differences[np.isfinite(differences)]

    if len(differences) == 0 or pd.isna(effect) or effect == 0:
        return np.nan

    if effect > 0:
        return float(np.mean(differences > 0))
    return float(np.mean(differences < 0))


def strict_histogram_inference_call(p_value, effect, ci_low, ci_high, differences):
    """
    Decide whether the histogram should display a significance bracket/stars.

    This is intentionally stricter than a p-value alone. It requires:
      1. replicate-level p < STRICT_ALPHA,
      2. effect magnitude >= STRICT_EFFECT_THRESHOLD_UM,
      3. enough replicate pairs,
      4. enough replicate-level direction agreement,
      5. 95% CI excluding zero, when enabled,
      6. 95% CI clearing the strict effect threshold, when enabled.
    """
    differences = np.asarray(differences, dtype=float)
    differences = differences[np.isfinite(differences)]

    agreement = direction_agreement_fraction(differences, effect)

    p_pass = passes_p_value_gate(p_value)
    effect_pass = (not pd.isna(effect)) and abs(effect) >= STRICT_EFFECT_THRESHOLD_UM
    n_pass = len(differences) >= STRICT_MIN_REPLICATE_PAIRS
    direction_pass = (
        not pd.isna(agreement)
    ) and agreement >= STRICT_MIN_DIRECTION_AGREEMENT

    if not STRICT_REQUIRE_CI_EXCLUDES_ZERO:
        ci_zero_pass = True
    elif pd.isna(ci_low) or pd.isna(ci_high):
        ci_zero_pass = False
    else:
        ci_zero_pass = (ci_low > 0) or (ci_high < 0)

    if not STRICT_REQUIRE_CI_EXCLUDES_EFFECT_THRESHOLD:
        ci_threshold_pass = True
    elif pd.isna(ci_low) or pd.isna(ci_high) or pd.isna(effect):
        ci_threshold_pass = False
    elif effect > 0:
        ci_threshold_pass = ci_low >= STRICT_EFFECT_THRESHOLD_UM
    elif effect < 0:
        ci_threshold_pass = ci_high <= -STRICT_EFFECT_THRESHOLD_UM
    else:
        ci_threshold_pass = False

    passed = all(
        [p_pass, effect_pass, n_pass, direction_pass, ci_zero_pass, ci_threshold_pass]
    )

    reasons = []
    if not p_pass:
        reasons.append(f"p >= {STRICT_ALPHA}")
    if not effect_pass:
        reasons.append(f"abs(effect) < {STRICT_EFFECT_THRESHOLD_UM} um")
    if not n_pass:
        reasons.append(f"replicate pairs < {STRICT_MIN_REPLICATE_PAIRS}")
    if not direction_pass:
        reasons.append(f"direction agreement < {STRICT_MIN_DIRECTION_AGREEMENT:.0%}")
    if not ci_zero_pass:
        reasons.append("95% CI overlaps 0")
    if not ci_threshold_pass:
        reasons.append("95% CI does not clear effect threshold")

    if passed:
        call = "STRICT_PASS"
        reason = "passes strict histogram inference gate"
    else:
        call = "STRICT_FAIL"
        reason = "; ".join(reasons) if reasons else "strict gate failed"

    return {
        "strict_histogram_call": call,
        "strict_histogram_call_pass": bool(passed),
        "strict_histogram_call_reason": reason,
        "strict_histogram_direction_agreement": agreement,
        "strict_histogram_n_pairs_used": int(len(differences)),
        "strict_alpha": STRICT_ALPHA,
        "strict_effect_threshold_um": STRICT_EFFECT_THRESHOLD_UM,
        "strict_min_replicate_pairs": STRICT_MIN_REPLICATE_PAIRS,
        "strict_min_direction_agreement": STRICT_MIN_DIRECTION_AGREEMENT,
        "strict_require_ci_excludes_zero": STRICT_REQUIRE_CI_EXCLUDES_ZERO,
        "strict_require_ci_excludes_effect_threshold": STRICT_REQUIRE_CI_EXCLUDES_EFFECT_THRESHOLD,
    }


def rank_biserial_from_mannwhitney_u(u_stat, n1, n2):
    if n1 == 0 or n2 == 0 or pd.isna(u_stat):
        return np.nan
    return (2.0 * u_stat / (n1 * n2)) - 1.0


def paired_effect_size_dz(differences):
    differences = np.asarray(differences, dtype=float)
    differences = differences[np.isfinite(differences)]

    if len(differences) < 2:
        return np.nan

    sd = np.std(differences, ddof=1)

    if sd == 0:
        return np.nan

    return float(np.mean(differences) / sd)


def paired_sign_flip_permutation_test(
    differences, n_permutations=N_PERMUTATIONS, seed=RANDOM_SEED
):
    """
    Two-sided paired sign-flip permutation test on replicate-level differences.

    Null hypothesis:
      Real-random differences are symmetrically distributed around zero.

    This is more conservative than treating each mRNA distance as independent.
    """
    differences = np.asarray(differences, dtype=float)
    differences = differences[np.isfinite(differences)]
    differences = differences[differences != 0]

    if len(differences) == 0:
        return np.nan, np.nan, 0

    observed = float(np.mean(differences))

    rng = np.random.default_rng(seed)
    signs = rng.choice([-1.0, 1.0], size=(int(n_permutations), len(differences)))
    permuted = np.mean(signs * differences, axis=1)

    # Add one to numerator/denominator for a stable permutation p-value.
    p = (np.sum(np.abs(permuted) >= abs(observed)) + 1.0) / (len(permuted) + 1.0)

    return observed, float(p), int(len(differences))


def bootstrap_ci_mean(differences, n_boot=10000, seed=RANDOM_SEED):
    differences = np.asarray(differences, dtype=float)
    differences = differences[np.isfinite(differences)]

    if len(differences) == 0:
        return np.nan, np.nan

    rng = np.random.default_rng(seed)
    boots = []

    for _ in range(int(n_boot)):
        sample = rng.choice(differences, size=len(differences), replace=True)
        boots.append(np.mean(sample))

    low, high = np.percentile(boots, [2.5, 97.5])
    return float(low), float(high)


def run_point_level_stats(real_for_stats, random_for_stats):
    """
    Point-level tests. These are kept in the CSV but are not the main plot
    significance call, because pooled spot-level n can make tiny differences
    look highly significant.
    """
    stats = {
        "point_mannwhitney_u": np.nan,
        "point_mannwhitney_p": np.nan,
        "point_mannwhitney_rank_biserial": np.nan,
        "point_ks_statistic": np.nan,
        "point_ks_p": np.nan,
        "point_stats_note": "",
    }

    if len(real_for_stats) == 0 or len(random_for_stats) == 0:
        stats["point_stats_note"] = "missing_real_or_random_values"
        return stats

    if not HAS_SCIPY:
        stats["point_stats_note"] = "scipy_not_available_point_stats_not_run"
        return stats

    try:
        mw = mannwhitneyu(real_for_stats, random_for_stats, alternative="two-sided")
        stats["point_mannwhitney_u"] = float(mw.statistic)
        stats["point_mannwhitney_p"] = float(mw.pvalue)
        stats["point_mannwhitney_rank_biserial"] = float(
            rank_biserial_from_mannwhitney_u(
                stats["point_mannwhitney_u"],
                len(real_for_stats),
                len(random_for_stats),
            )
        )
    except Exception as e:
        stats["point_stats_note"] += f"mannwhitney_error: {e}; "

    try:
        ks = ks_2samp(real_for_stats, random_for_stats, alternative="two-sided")
        stats["point_ks_statistic"] = float(ks.statistic)
        stats["point_ks_p"] = float(ks.pvalue)
    except Exception as e:
        stats["point_stats_note"] += f"ks_error: {e}; "

    return stats


def run_replicate_level_stats(replicate_df):
    """
    Main conservative statistics.
    Each real/random file pair contributes one biological/series-level value.

    Metrics tested:
      - median difference: real median - random median
      - mean difference: real mean - random mean

    The plot uses the replicate median permutation p-value by default.
    """
    stats = {
        "replicate_n_pairs": 0,
        "replicate_median_difference_mean_um": np.nan,
        "replicate_median_difference_median_um": np.nan,
        "replicate_median_difference_ci95_low_um": np.nan,
        "replicate_median_difference_ci95_high_um": np.nan,
        "replicate_median_permutation_p": np.nan,
        "replicate_median_wilcoxon_p": np.nan,
        "replicate_median_effect_dz": np.nan,
        "replicate_mean_difference_mean_um": np.nan,
        "replicate_mean_difference_median_um": np.nan,
        "replicate_mean_difference_ci95_low_um": np.nan,
        "replicate_mean_difference_ci95_high_um": np.nan,
        "replicate_mean_permutation_p": np.nan,
        "replicate_mean_wilcoxon_p": np.nan,
        "replicate_mean_effect_dz": np.nan,
        "main_plot_test": MAIN_PLOT_TEST,
        "main_plot_p": np.nan,
        "main_plot_effect_um": np.nan,
        "main_plot_ci95_low_um": np.nan,
        "main_plot_ci95_high_um": np.nan,
        "main_plot_significance": "n/a",
        "main_plot_practical_effect_call": "n/a",
        "strict_histogram_call": "n/a",
        "strict_histogram_call_pass": False,
        "strict_histogram_call_reason": "n/a",
        "strict_histogram_direction_agreement": np.nan,
        "strict_histogram_n_pairs_used": 0,
        "strict_alpha": STRICT_ALPHA,
        "strict_effect_threshold_um": STRICT_EFFECT_THRESHOLD_UM,
        "strict_min_replicate_pairs": STRICT_MIN_REPLICATE_PAIRS,
        "strict_min_direction_agreement": STRICT_MIN_DIRECTION_AGREEMENT,
        "strict_require_ci_excludes_zero": STRICT_REQUIRE_CI_EXCLUDES_ZERO,
        "strict_require_ci_excludes_effect_threshold": STRICT_REQUIRE_CI_EXCLUDES_EFFECT_THRESHOLD,
    }

    if len(replicate_df) == 0:
        return stats

    median_diff = (
        replicate_df["real_median_diffraction_limited"]
        - replicate_df["random_median_diffraction_limited"]
    )
    mean_diff = (
        replicate_df["real_mean_diffraction_limited"]
        - replicate_df["random_mean_diffraction_limited"]
    )

    median_diff = (
        median_diff.replace([np.inf, -np.inf], np.nan).dropna().to_numpy(dtype=float)
    )
    mean_diff = (
        mean_diff.replace([np.inf, -np.inf], np.nan).dropna().to_numpy(dtype=float)
    )

    stats["replicate_n_pairs"] = int(len(replicate_df))

    if len(median_diff) > 0:
        obs, p_perm, n_used = paired_sign_flip_permutation_test(median_diff)
        ci_low, ci_high = bootstrap_ci_mean(median_diff)

        stats["replicate_median_difference_mean_um"] = float(np.mean(median_diff))
        stats["replicate_median_difference_median_um"] = float(np.median(median_diff))
        stats["replicate_median_difference_ci95_low_um"] = ci_low
        stats["replicate_median_difference_ci95_high_um"] = ci_high
        stats["replicate_median_permutation_p"] = p_perm
        stats["replicate_median_effect_dz"] = paired_effect_size_dz(median_diff)

        if HAS_SCIPY and len(median_diff[median_diff != 0]) > 0:
            try:
                w = wilcoxon(median_diff, alternative="two-sided", zero_method="wilcox")
                stats["replicate_median_wilcoxon_p"] = float(w.pvalue)
            except Exception:
                pass

    if len(mean_diff) > 0:
        obs, p_perm, n_used = paired_sign_flip_permutation_test(mean_diff)
        ci_low, ci_high = bootstrap_ci_mean(mean_diff)

        stats["replicate_mean_difference_mean_um"] = float(np.mean(mean_diff))
        stats["replicate_mean_difference_median_um"] = float(np.median(mean_diff))
        stats["replicate_mean_difference_ci95_low_um"] = ci_low
        stats["replicate_mean_difference_ci95_high_um"] = ci_high
        stats["replicate_mean_permutation_p"] = p_perm
        stats["replicate_mean_effect_dz"] = paired_effect_size_dz(mean_diff)

        if HAS_SCIPY and len(mean_diff[mean_diff != 0]) > 0:
            try:
                w = wilcoxon(mean_diff, alternative="two-sided", zero_method="wilcox")
                stats["replicate_mean_wilcoxon_p"] = float(w.pvalue)
            except Exception:
                pass

    if MAIN_PLOT_TEST == "replicate_mean_permutation":
        main_p = stats["replicate_mean_permutation_p"]
        main_effect = stats["replicate_mean_difference_mean_um"]
        main_ci_low = stats["replicate_mean_difference_ci95_low_um"]
        main_ci_high = stats["replicate_mean_difference_ci95_high_um"]
        main_differences = mean_diff
    else:
        main_p = stats["replicate_median_permutation_p"]
        main_effect = stats["replicate_median_difference_mean_um"]
        main_ci_low = stats["replicate_median_difference_ci95_low_um"]
        main_ci_high = stats["replicate_median_difference_ci95_high_um"]
        main_differences = median_diff

    stats["main_plot_p"] = main_p
    stats["main_plot_effect_um"] = main_effect
    stats["main_plot_ci95_low_um"] = main_ci_low
    stats["main_plot_ci95_high_um"] = main_ci_high
    stats["main_plot_significance"] = significance_stars(main_p)

    strict_call = strict_histogram_inference_call(
        p_value=main_p,
        effect=main_effect,
        ci_low=main_ci_low,
        ci_high=main_ci_high,
        differences=main_differences,
    )
    stats.update(strict_call)

    if pd.isna(main_effect):
        stats["main_plot_practical_effect_call"] = "n/a"
    elif abs(main_effect) >= PRACTICAL_EFFECT_THRESHOLD_UM:
        stats["main_plot_practical_effect_call"] = (
            f"effect >= {PRACTICAL_EFFECT_THRESHOLD_UM} um"
        )
    else:
        stats["main_plot_practical_effect_call"] = (
            f"effect < {PRACTICAL_EFFECT_THRESHOLD_UM} um"
        )

    return stats


def summarize_replicates_for_mrna(mrna, real_arrays, random_arrays, file_records):
    rows = []

    for real_arr, random_arr, rec in zip(real_arrays, random_arrays, file_records):
        real = apply_diffraction_floor(real_arr)
        random = apply_diffraction_floor(random_arr)

        row = dict(rec)
        row.update(
            {
                "mRNA": mrna,
                "mRNA_display_name": display_mrna_name(mrna),
                "real_n": int(len(real)),
                "random_n": int(len(random)),
                "real_mean_raw": float(np.mean(real_arr)) if len(real_arr) else np.nan,
                "random_mean_raw": float(np.mean(random_arr))
                if len(random_arr)
                else np.nan,
                "real_median_raw": float(np.median(real_arr))
                if len(real_arr)
                else np.nan,
                "random_median_raw": float(np.median(random_arr))
                if len(random_arr)
                else np.nan,
                "real_mean_diffraction_limited": float(np.mean(real))
                if len(real)
                else np.nan,
                "random_mean_diffraction_limited": float(np.mean(random))
                if len(random)
                else np.nan,
                "real_median_diffraction_limited": float(np.median(real))
                if len(real)
                else np.nan,
                "random_median_diffraction_limited": float(np.median(random))
                if len(random)
                else np.nan,
                "mean_difference_real_minus_random_um": float(
                    np.mean(real) - np.mean(random)
                )
                if len(real) and len(random)
                else np.nan,
                "median_difference_real_minus_random_um": float(
                    np.median(real) - np.median(random)
                )
                if len(real) and len(random)
                else np.nan,
                "real_fraction_at_diffraction_limit": float(
                    np.mean(real == DIFFRACTION_LIMIT_UM)
                )
                if len(real)
                else np.nan,
                "random_fraction_at_diffraction_limit": float(
                    np.mean(random == DIFFRACTION_LIMIT_UM)
                )
                if len(random)
                else np.nan,
            }
        )

        rows.append(row)

    return pd.DataFrame(rows)


def add_significance_stars_between_means(
    ax, real_mean, random_mean, p_value, y_top, strict_pass
):
    """
    Add significance stars centered between the two mean lines.
    Stars/brackets are drawn only when the strict histogram inference gate passes.
    """
    if not strict_pass:
        return

    stars = significance_stars(p_value)

    x1, x2 = sorted([real_mean, random_mean])

    min_width = 0.06
    if abs(x2 - x1) < min_width:
        midpoint = (x1 + x2) / 2.0
        x1 = midpoint - min_width / 2.0
        x2 = midpoint + min_width / 2.0

    y = y_top * 1.05 if y_top > 0 else 0.05
    h = y_top * 0.04 if y_top > 0 else 0.02

    ax.plot(
        [x1, x1, x2, x2],
        [y, y + h, y + h, y],
        color="black",
        linewidth=1.8,
        clip_on=False,
    )

    ax.text(
        (x1 + x2) / 2.0,
        y + h * 1.10,
        stars,
        ha="center",
        va="bottom",
        fontsize=18,
        fontweight="bold",
        clip_on=False,
    )

    ax.set_ylim(0, y + h * 3.5)


def add_delta_arrow_between_medians(ax, real_median, random_median, y_top):
    """
    Add a horizontal arrow showing the descriptive median real-minus-random shift.
    This annotation is independent of p-values and practical-effect thresholds.
    """
    if not SHOW_DELTA_MEDIAN_ARROW:
        return

    if pd.isna(real_median) or pd.isna(random_median):
        return

    delta = float(real_median - random_median)
    x1 = float(random_median)
    x2 = float(real_median)

    # If the medians are very close, use a minimum arrow width so the label is visible.
    min_width = 0.045
    arrow_x1, arrow_x2 = x1, x2
    if abs(arrow_x2 - arrow_x1) < min_width:
        midpoint = (arrow_x1 + arrow_x2) / 2.0
        arrow_x1 = midpoint - min_width / 2.0
        arrow_x2 = midpoint + min_width / 2.0

    y = y_top * 0.86 if y_top > 0 else 0.05
    label_y = y + (y_top * 0.035 if y_top > 0 else 0.015)

    ax.annotate(
        "",
        xy=(arrow_x2, y),
        xytext=(arrow_x1, y),
        arrowprops={
            "arrowstyle": "<->",
            "linewidth": 2.2,
            "color": "black",
            "shrinkA": 0,
            "shrinkB": 0,
        },
        clip_on=False,
    )

    ax.text(
        (arrow_x1 + arrow_x2) / 2.0,
        label_y,
        f"Δ median = {delta:+.3f} µm",
        ha="center",
        va="bottom",
        fontsize=DELTA_ARROW_FONTSIZE,
        color="black",
        bbox={
            "boxstyle": "round,pad=0.25",
            "facecolor": "white",
            "edgecolor": "none",
            "alpha": 0.85,
        },
        clip_on=False,
    )


def pooled_plot(mrna, real_arrays, random_arrays, file_records, output_folder):
    safe_mrna = sanitize_name(mrna)
    display_name = display_mrna_name(safe_mrna)
    strain_name = parent_folder.name

    real_raw = np.concatenate(real_arrays) if len(real_arrays) > 0 else np.array([])
    random_raw = (
        np.concatenate(random_arrays) if len(random_arrays) > 0 else np.array([])
    )

    real = apply_diffraction_floor(real_raw)
    random = apply_diffraction_floor(random_raw)

    output_folder.mkdir(parents=True, exist_ok=True)

    pooled_real_raw_npy = (
        output_folder / f"{safe_mrna}_POOLED_real_NN_distance_um_RAW.npy"
    )
    pooled_random_raw_npy = (
        output_folder / f"{safe_mrna}_POOLED_random_NN_distance_um_RAW.npy"
    )
    pooled_real_floor_npy = (
        output_folder
        / f"{safe_mrna}_POOLED_real_NN_distance_um_DIFFRACTION_FLOORED.npy"
    )
    pooled_random_floor_npy = (
        output_folder
        / f"{safe_mrna}_POOLED_random_NN_distance_um_DIFFRACTION_FLOORED.npy"
    )

    np.save(pooled_real_raw_npy, real_raw)
    np.save(pooled_random_raw_npy, random_raw)
    np.save(pooled_real_floor_npy, real)
    np.save(pooled_random_floor_npy, random)

    pooled_values_csv = output_folder / output_csv_name(
        f"{safe_mrna}_POOLED_real_vs_random_values", strain_name
    )

    pooled_df = pd.concat(
        [
            pd.DataFrame(
                {
                    "mRNA": safe_mrna,
                    "mRNA_display_name": display_name,
                    "group": "real",
                    "distance_um_raw": real_raw,
                    "distance_um_diffraction_limited": real,
                    "was_below_diffraction_limit": real_raw < DIFFRACTION_LIMIT_UM,
                }
            ),
            pd.DataFrame(
                {
                    "mRNA": safe_mrna,
                    "mRNA_display_name": display_name,
                    "group": "random",
                    "distance_um_raw": random_raw,
                    "distance_um_diffraction_limited": random,
                    "was_below_diffraction_limit": random_raw < DIFFRACTION_LIMIT_UM,
                }
            ),
        ],
        ignore_index=True,
    )
    pooled_df = move_column_to_end(pooled_df, "strain", strain_name)
    pooled_df.to_csv(pooled_values_csv, index=False)

    replicate_df = summarize_replicates_for_mrna(
        safe_mrna, real_arrays, random_arrays, file_records
    )
    replicate_values_csv = output_folder / output_csv_name(
        f"{safe_mrna}_REPLICATE_LEVEL_real_vs_random_summary", strain_name
    )
    replicate_df = move_column_to_end(replicate_df, "strain", strain_name)
    replicate_df.to_csv(replicate_values_csv, index=False)

    point_stats = run_point_level_stats(real, random)
    replicate_stats = run_replicate_level_stats(replicate_df)

    summary = {
        "mRNA": safe_mrna,
        "mRNA_display_name": display_name,
        "strain": strain_name,
        "series_range_label": make_series_range_label(file_records),
        "diffraction_limit_um": DIFFRACTION_LIMIT_UM,
        "practical_effect_threshold_um": PRACTICAL_EFFECT_THRESHOLD_UM,
        "strict_alpha": STRICT_ALPHA,
        "strict_effect_threshold_um": STRICT_EFFECT_THRESHOLD_UM,
        "strict_min_replicate_pairs": STRICT_MIN_REPLICATE_PAIRS,
        "strict_min_direction_agreement": STRICT_MIN_DIRECTION_AGREEMENT,
        "strict_require_ci_excludes_zero": STRICT_REQUIRE_CI_EXCLUDES_ZERO,
        "strict_require_ci_excludes_effect_threshold": STRICT_REQUIRE_CI_EXCLUDES_EFFECT_THRESHOLD,
        "real_n_distances": len(real),
        "random_n_distances": len(random),
        "real_n_below_diffraction_limit": int(np.sum(real_raw < DIFFRACTION_LIMIT_UM))
        if len(real_raw)
        else 0,
        "random_n_below_diffraction_limit": int(
            np.sum(random_raw < DIFFRACTION_LIMIT_UM)
        )
        if len(random_raw)
        else 0,
        "real_fraction_below_diffraction_limit": float(
            np.mean(real_raw < DIFFRACTION_LIMIT_UM)
        )
        if len(real_raw)
        else np.nan,
        "random_fraction_below_diffraction_limit": float(
            np.mean(random_raw < DIFFRACTION_LIMIT_UM)
        )
        if len(random_raw)
        else np.nan,
        "real_raw_mean": float(np.mean(real_raw)) if len(real_raw) else np.nan,
        "random_raw_mean": float(np.mean(random_raw)) if len(random_raw) else np.nan,
        "real_raw_median": float(np.median(real_raw)) if len(real_raw) else np.nan,
        "random_raw_median": float(np.median(random_raw))
        if len(random_raw)
        else np.nan,
        "real_diffraction_limited_mean": float(np.mean(real)) if len(real) else np.nan,
        "random_diffraction_limited_mean": float(np.mean(random))
        if len(random)
        else np.nan,
        "real_diffraction_limited_median": float(np.median(real))
        if len(real)
        else np.nan,
        "random_diffraction_limited_median": float(np.median(random))
        if len(random)
        else np.nan,
        "real_diffraction_limited_std": float(np.std(real, ddof=1))
        if len(real) > 1
        else np.nan,
        "random_diffraction_limited_std": float(np.std(random, ddof=1))
        if len(random) > 1
        else np.nan,
        "real_files_pooled": len(real_arrays),
        "random_files_pooled": len(random_arrays),
        "pooled_real_raw_npy": str(pooled_real_raw_npy),
        "pooled_random_raw_npy": str(pooled_random_raw_npy),
        "pooled_real_diffraction_limited_npy": str(pooled_real_floor_npy),
        "pooled_random_diffraction_limited_npy": str(pooled_random_floor_npy),
        "pooled_values_csv": str(pooled_values_csv),
        "replicate_values_csv": str(replicate_values_csv),
    }

    summary.update(point_stats)
    summary.update(replicate_stats)

    if len(real) == 0 or len(random) == 0:
        note_path = output_folder / f"{safe_mrna}_POOLED_real_vs_random_NOT_CREATED.txt"
        note_path.write_text(
            f"Could not create pooled plot for {safe_mrna}.\n"
            f"Real distances: {len(real)}\n"
            f"Random distances: {len(random)}\n",
            encoding="utf-8",
        )
        summary["status"] = "skipped_empty_real_or_random"
        summary["plot"] = ""
        summary["note"] = str(note_path)
        return summary

    bins = np.arange(
        DIFFRACTION_LIMIT_UM,
        MAX_DISTANCE_UM + BIN_WIDTH_UM,
        BIN_WIDTH_UM,
    )

    if bins[-1] < MAX_DISTANCE_UM:
        bins = np.append(bins, MAX_DISTANCE_UM)

    output_png = (
        output_folder
        / f"{safe_mrna}_POOLED_real_vs_random_REPLICATE_STATS_{strain_name}.png"
    )

    fig, ax = plt.subplots(figsize=(12, 8), facecolor="white")

    real_weights = np.ones_like(real, dtype=float) / len(real)
    random_weights = np.ones_like(random, dtype=float) / len(random)

    real_heights, _, _ = ax.hist(
        real,
        bins=bins,
        weights=real_weights,
        color=REAL_COLORS.get(safe_mrna.upper(), "gray"),
        edgecolor="black",
        alpha=REAL_ALPHA,
        label=f"Real {display_name}",
    )

    random_heights, _, _ = ax.hist(
        random,
        bins=bins,
        weights=random_weights,
        histtype="step",
        color=RANDOM_COLOR,
        linewidth=RANDOM_LINEWIDTH,
        alpha=1.0,
        label=f"Random {display_name}",
    )

    real_median = summary["real_diffraction_limited_median"]
    random_median = summary["random_diffraction_limited_median"]

    if SHOW_MEDIAN_LINES:
        ax.axvline(
            real_median,
            linestyle=MEDIAN_LINESTYLE,
            linewidth=MEDIAN_LINEWIDTH,
            color=REAL_COLORS.get(safe_mrna.upper(), "gray"),
            label=f"Real median = {real_median:.3f} µm",
        )

        ax.axvline(
            random_median,
            linestyle=MEDIAN_LINESTYLE,
            linewidth=MEDIAN_LINEWIDTH,
            color=RANDOM_COLOR,
            label=f"Random median = {random_median:.3f} µm",
        )

    main_p = summary["main_plot_p"]
    y_top = max(
        float(np.max(real_heights)) if len(real_heights) else 0.0,
        float(np.max(random_heights)) if len(random_heights) else 0.0,
    )

    add_delta_arrow_between_medians(
        ax,
        real_median=real_median,
        random_median=random_median,
        y_top=y_top,
    )

    if SHOW_STRICT_SIGNIFICANCE_STARS:
        add_significance_stars_between_means(
            ax,
            real_median,
            random_median,
            main_p,
            y_top,
            strict_pass=bool(summary.get("strict_histogram_call_pass", False)),
        )

    # Keep the final plot free of a statistical summary box.
    # P-values, CIs, and strict-call details are saved in the summary CSV.

    ax.set_xlabel("Nearest-neighbor distance (µm)", fontsize=18)
    ax.set_ylabel("Proportion of distances", fontsize=18)
    ax.set_title(make_plot_title(strain_name, safe_mrna, file_records), fontsize=18)
    ax.legend(fontsize=18, loc="upper right", bbox_to_anchor=(1.0, 0.72), frameon=True)

    ax.set_xticks(np.arange(DIFFRACTION_LIMIT_UM, MAX_DISTANCE_UM + 0.1, 0.2))
    ax.tick_params(axis="x", labelsize=18)
    ax.tick_params(axis="y", labelsize=18)
    ax.set_xlim(DIFFRACTION_LIMIT_UM, MAX_DISTANCE_UM)

    ax.grid(axis="y", alpha=0.25)

    plt.tight_layout()
    plt.savefig(output_png, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)

    summary["status"] = "plotted"
    summary["plot"] = str(output_png)

    return summary


def main():
    if not parent_folder.exists():
        raise FileNotFoundError(f"Parent folder does not exist: {parent_folder}")

    output_folder = parent_folder / POOLED_OUTPUT_FOLDER_NAME
    output_folder.mkdir(parents=True, exist_ok=True)

    real_files = find_real_distance_files(parent_folder)

    print(f"Parent folder: {parent_folder}")
    print(f"Strain/condition title: {parent_folder.name}")
    print(f"Found {len(real_files)} real NN distance files")
    print(f"Diffraction limit floor: {DIFFRACTION_LIMIT_UM} um")
    print(f"SciPy available for point-level statistics: {HAS_SCIPY}")
    print(f"Main plot test: {MAIN_PLOT_TEST}")
    print("Strict histogram gate:")
    print(f"  p < {STRICT_ALPHA}")
    print(f"  abs(effect) >= {STRICT_EFFECT_THRESHOLD_UM} um")
    print(f"  replicate pairs >= {STRICT_MIN_REPLICATE_PAIRS}")
    print(f"  direction agreement >= {STRICT_MIN_DIRECTION_AGREEMENT:.0%}")
    print(f"  CI excludes zero: {STRICT_REQUIRE_CI_EXCLUDES_ZERO}")
    print(
        f"  CI clears effect threshold: {STRICT_REQUIRE_CI_EXCLUDES_EFFECT_THRESHOLD}"
    )

    pooled = {}
    records = []

    for real_file in real_files:
        mrna = sanitize_name(infer_mrna_name(real_file))
        random_file = find_random_file_for_real(real_file, mrna)
        series_folder = find_series_folder(real_file)
        series_number = (
            extract_series_number_from_name(series_folder.name)
            if series_folder
            else None
        )

        record = {
            "mRNA": mrna,
            "mRNA_display_name": display_mrna_name(mrna),
            "SeriesFolder": str(series_folder) if series_folder else "",
            "SeriesNumber": series_number if series_number is not None else "",
            "RealFile": str(real_file),
            "RandomFile": str(random_file) if random_file else "",
            "Status": "",
            "RealN": "",
            "RandomN": "",
            "strain": parent_folder.name,
        }

        print(f"\nReal file: {real_file}")
        print(f"  mRNA/channel name: {mrna}")
        print(f"  Plot display name: {display_mrna_name(mrna)}")

        if random_file is None:
            print("  No matching random file found; skipping this real file.")
            record["Status"] = "missing_random_file"
            records.append(record)
            continue

        print(f"  Random file: {random_file}")

        real = load_clean_npy(real_file)
        random = load_clean_npy(random_file)

        record["RealN"] = len(real)
        record["RandomN"] = len(random)
        record["RealNBelowDiffractionLimit"] = int(np.sum(real < DIFFRACTION_LIMIT_UM))
        record["RandomNBelowDiffractionLimit"] = int(
            np.sum(random < DIFFRACTION_LIMIT_UM)
        )

        if len(real) == 0 or len(random) == 0:
            print(
                f"  Empty real or random array; real={len(real)}, random={len(random)}"
            )
            record["Status"] = "empty_real_or_random"
            records.append(record)
            continue

        if mrna not in pooled:
            pooled[mrna] = {
                "real_arrays": [],
                "random_arrays": [],
                "file_records": [],
            }

        pooled[mrna]["real_arrays"].append(real)
        pooled[mrna]["random_arrays"].append(random)
        pooled[mrna]["file_records"].append(record)

        record["Status"] = "pooled"
        records.append(record)

    strain_name = parent_folder.name
    pair_summary_path = output_folder / output_csv_name(
        "pooled_input_file_pairing_summary", strain_name
    )
    records_df = move_column_to_end(pd.DataFrame(records), "strain", strain_name)
    records_df.to_csv(pair_summary_path, index=False)

    pooled_summaries = []

    for mrna, data in sorted(pooled.items()):
        print(f"\nCreating pooled plot for {mrna} as {display_mrna_name(mrna)}")
        print(f"  Real files pooled: {len(data['real_arrays'])}")
        print(f"  Random files pooled: {len(data['random_arrays'])}")
        print(f"  Series label: {make_series_range_label(data['file_records'])}")

        summary = pooled_plot(
            mrna=mrna,
            real_arrays=data["real_arrays"],
            random_arrays=data["random_arrays"],
            file_records=data["file_records"],
            output_folder=output_folder,
        )

        pooled_summaries.append(summary)

        if summary["status"] == "plotted":
            print(f"  Saved plot: {summary['plot']}")
            print(
                f"  Main replicate-level p={summary['main_plot_p']}, "
                f"strict_call={summary['strict_histogram_call']}, "
                f"reason={summary['strict_histogram_call_reason']}"
            )
            print(f"  Point-level MW p={summary['point_mannwhitney_p']}")
        else:
            print(f"  Skipped plot: {summary['status']}")

    pooled_summary_path = output_folder / output_csv_name(
        "pooled_real_vs_random_summary", strain_name
    )
    pooled_summary_df = move_column_to_end(
        pd.DataFrame(pooled_summaries), "strain", strain_name
    )
    pooled_summary_df.to_csv(pooled_summary_path, index=False)

    print("\nDone.")
    print(f"Pooled output folder: {output_folder}")
    print(f"Input pairing summary: {pair_summary_path}")
    print(f"Pooled summary: {pooled_summary_path}")
    print("\nMain interpretation:")
    print(
        "  Use the histogram as the main visual inference, but only treat starred/bracketed plots as strict positives."
    )
    print(
        "  Strict positives require replicate-level p, effect-size, CI-threshold, n, and direction-agreement gates to pass."
    )
    print(
        "  The point-level Mann-Whitney p-value is saved, but it can be too sensitive when many spots are pooled."
    )
    print("\nFinal plots use diffraction-limited distances:")
    print(
        f"  all values < {DIFFRACTION_LIMIT_UM} um are plotted/statistically tested as {DIFFRACTION_LIMIT_UM} um"
    )


if __name__ == "__main__":
    main()
