import re
import tkinter as tk
from pathlib import Path
from tkinter import filedialog

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ==========================================================
# Parent folder containing all Series folders
# Change this to your experiment parent folder.
# ==========================================================
root = tk.Tk()
root.withdraw()

parent_folder = Path(filedialog.askdirectory())

print("User selected:", parent_folder)
# ==========================================================
# Plot settings
# ==========================================================
BIN_WIDTH_UM = 0.05
MAX_DISTANCE_UM = 2.0

REAL_COLORS = {
    "MS2": "red",
    "ATP2": "yellow",
    "ATP3": "yellow",
    "TIM50": "yellow",
}

RANDOM_COLOR = "cornflowerblue"

# ==========================================================
# Names produced by previous recursive scripts
# ==========================================================
REAL_PATTERN = "*_NN_distance_um.npy"
RANDOM_FOLDER_SUFFIX = "_output"


def sanitize_name(name):
    name = str(name).strip()
    if name == "":
        name = "mRNA"
    for bad in ["\\", "/", ":", "*", "?", '"', "<", ">", "|"]:
        name = name.replace(bad, "_")
    name = re.sub(r"\s+", "_", name)
    return name


def infer_mrna_name(path):
    """
    Infer MS2 / ATP2 from filename or path.
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
    """
    Finds nearest ancestor whose name starts with Series.
    """
    p = Path(path)
    for ancestor in p.parents:
        if ancestor.name.lower().startswith("series"):
            return ancestor
    return None


def is_generated_or_random_real_file(path):
    """
    Prevent accidental treatment of random output npy files as real files.
    """
    p = Path(path)
    lower_parts = [part.lower() for part in p.parts]
    name = p.name.lower()

    if "random_" in name:
        return True
    if any(
        part.startswith("random_") and part.endswith("_output") for part in lower_parts
    ):
        return True

    return False


def find_real_distance_files(parent):
    """
    Find real nearest-neighbor distance files from the previous real mRNA NN script.
    Expected examples:
      ATP2_NN_distance_um.npy
      MS2_NN_distance_um.npy
      ATP3_NN_distance_um.npy
      TIM50_NN_distance_um.npy
    """
    files = []
    for p in parent.rglob(REAL_PATTERN):
        if is_generated_or_random_real_file(p):
            continue
        files.append(p)

    return sorted(files)


def find_random_file_for_real(real_file, mrna):
    """
    Find matching random null distribution for a real NN file.

    Expected location after random_mrna_null_recursive_shared_skeleton_SHORT_PATHS.py:
      same spots_extraction folder/random_ATP2_output/random_ATP2_nn.npy
      same spots_extraction folder/random_MS2_output/random_MS2_nn.npy

    Also includes fallback recursive search within the same Series folder.
    """
    real_file = Path(real_file)
    spots_folder = real_file.parent
    safe_mrna = sanitize_name(mrna)

    direct_candidates = [
        spots_folder / f"random_{safe_mrna}_output" / f"random_{safe_mrna}_nn.npy",
        spots_folder
        / f"random_{safe_mrna.upper()}_output"
        / f"random_{safe_mrna.upper()}_nn.npy",
        spots_folder
        / f"random_{safe_mrna.lower()}_output"
        / f"random_{safe_mrna.lower()}_nn.npy",
    ]

    for candidate in direct_candidates:
        if candidate.exists():
            return candidate

    # Case-insensitive search directly under the spots folder.
    target_name = f"random_{safe_mrna}_nn.npy".lower()
    for p in spots_folder.rglob("*.npy"):
        if p.name.lower() == target_name:
            return p

    # Fallback: search within same mRNA branch, then same series.
    series_folder = find_series_folder(real_file)

    search_roots = [spots_folder]
    if series_folder is not None:
        # Prefer same mRNA branch under the series.
        for child in series_folder.iterdir():
            if child.is_dir() and child.name.upper() == safe_mrna.upper():
                search_roots.append(child)
        search_roots.append(series_folder)

    seen = set()
    for root in search_roots:
        root = Path(root)
        if root in seen or not root.exists():
            continue
        seen.add(root)

        for p in root.rglob("*.npy"):
            name_lower = p.name.lower()
            parts_lower = [part.lower() for part in p.parts]

            if name_lower == target_name:
                return p

            if (
                name_lower.endswith("_nn.npy")
                and name_lower.startswith(f"random_{safe_mrna.lower()}")
                and any(
                    part == f"random_{safe_mrna.lower()}_output" for part in parts_lower
                )
            ):
                return p

    return None


def load_clean_npy(path):
    arr = np.load(path)
    arr = np.asarray(arr, dtype=float).ravel()
    arr = arr[np.isfinite(arr)]
    return arr


def make_comparison_plot(real_file, random_file, mrna):
    """
    Load real and random NN distributions, plot density-normalized histograms,
    and save comparison outputs next to the real file.
    """
    real = load_clean_npy(real_file)
    random = load_clean_npy(random_file)

    output_folder = Path(real_file).parent / "plots"
    output_folder.mkdir(parents=True, exist_ok=True)

    safe_mrna = sanitize_name(mrna)
    output_png = output_folder / f"{safe_mrna}_real_vs_random.png"

    if len(real) == 0 or len(random) == 0:
        note = output_folder / f"{safe_mrna}_real_vs_random_NOT_CREATED.txt"
        note.write_text(
            f"Could not create plot for {safe_mrna}.\n"
            f"Real file: {real_file}\n"
            f"Random file: {random_file}\n"
            f"Real distances: {len(real)}\n"
            f"Random distances: {len(random)}\n",
            encoding="utf-8",
        )

        return {
            "mrna": safe_mrna,
            "real_file": str(real_file),
            "random_file": str(random_file),
            "status": "skipped_empty_real_or_random",
            "real_n": len(real),
            "random_n": len(random),
            "plot": "",
            "note": str(note),
        }

    bins = np.arange(0, MAX_DISTANCE_UM + BIN_WIDTH_UM, BIN_WIDTH_UM)

    plt.figure(figsize=(12, 8), facecolor="white")

    plt.hist(
        real,
        bins=bins,
        density=True,
        color=REAL_COLORS.get(safe_mrna.upper(), "gray"),
        edgecolor="black",
        alpha=0.7,
        label=f"nn_real_{safe_mrna}_distance",
    )

    plt.hist(
        random,
        bins=bins,
        density=True,
        color=RANDOM_COLOR,
        edgecolor="black",
        alpha=0.7,
        label=f"nn_random_{safe_mrna}_distance",
    )

    plt.xlabel("Distance (µm)", fontsize=18)
    plt.ylabel("Proportion", fontsize=18)
    plt.title(f"{safe_mrna}: real vs random", fontsize=18)
    plt.legend(fontsize=14)

    plt.xticks(
        np.arange(0, MAX_DISTANCE_UM + 0.1, 0.5),
        fontsize=14,
    )
    plt.yticks(fontsize=14)
    plt.xlim(0, MAX_DISTANCE_UM)

    plt.tight_layout()
    plt.savefig(output_png, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close()

    # Save basic numerical summary for each comparison.
    summary_csv = output_folder / f"{safe_mrna}_real_vs_random_summary.csv"

    summary_df = pd.DataFrame(
        [
            {
                "mRNA": safe_mrna,
                "real_file": str(real_file),
                "random_file": str(random_file),
                "real_n": len(real),
                "random_n": len(random),
                "real_mean": float(np.mean(real)),
                "random_mean": float(np.mean(random)),
                "real_median": float(np.median(real)),
                "random_median": float(np.median(random)),
                "real_std": float(np.std(real, ddof=1)) if len(real) > 1 else np.nan,
                "random_std": float(np.std(random, ddof=1))
                if len(random) > 1
                else np.nan,
            }
        ]
    )
    summary_df.to_csv(summary_csv, index=False)

    return {
        "mrna": safe_mrna,
        "real_file": str(real_file),
        "random_file": str(random_file),
        "status": "plotted",
        "real_n": len(real),
        "random_n": len(random),
        "plot": str(output_png),
        "summary_csv": str(summary_csv),
    }


def main():
    if not parent_folder.exists():
        raise FileNotFoundError(f"Parent folder does not exist: {parent_folder}")

    real_files = find_real_distance_files(parent_folder)

    print(f"Parent folder: {parent_folder}")
    print(f"Found {len(real_files)} real NN distance files")

    if len(real_files) == 0:
        print("No real distance files found.")
        print("Expected files like MS2_NN_distance_um.npy or ATP2_NN_distance_um.npy")
        return

    run_summary = []

    for real_file in real_files:
        mrna = infer_mrna_name(real_file)
        random_file = find_random_file_for_real(real_file, mrna)

        print(f"\nReal file: {real_file}")
        print(f"  mRNA/channel name: {mrna}")

        if random_file is None:
            print("  No matching random NN file found.")
            output_folder = real_file.parent / "plots"
            output_folder.mkdir(parents=True, exist_ok=True)

            note = output_folder / f"{sanitize_name(mrna)}_NO_RANDOM_FILE_FOUND.txt"
            note.write_text(
                f"No matching random NN file was found for:\n{real_file}\n",
                encoding="utf-8",
            )

            run_summary.append(
                {
                    "mrna": sanitize_name(mrna),
                    "real_file": str(real_file),
                    "random_file": "",
                    "status": "missing_random_file",
                    "real_n": "",
                    "random_n": "",
                    "plot": "",
                    "note": str(note),
                }
            )
            continue

        print(f"  Random file: {random_file}")

        result = make_comparison_plot(real_file, random_file, mrna)
        run_summary.append(result)

        if result["status"] == "plotted":
            print(f"  Saved plot: {result['plot']}")
        else:
            print(f"  Skipped plot: {result['status']}")

    run_summary_path = parent_folder / "real_vs_random_recursive_comparison_summary.csv"
    pd.DataFrame(run_summary).to_csv(run_summary_path, index=False)

    print("\nDone.")
    print(f"Comparisons checked: {len(run_summary)}")
    print(f"Overall summary: {run_summary_path}")


if __name__ == "__main__":
    main()
