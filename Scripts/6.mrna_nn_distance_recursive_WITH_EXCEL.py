import ast
import os
import re
import tkinter as tk
from pathlib import Path
from tkinter import filedialog

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

# ==========================================================
# Parent folder containing all series folders
# Change this to your experiment parent folder.
# ==========================================================
root = tk.Tk()
root.withdraw()

parent_folder = Path(filedialog.askdirectory())

print("User selected:", parent_folder)

# ==========================================================
# Voxel size / pixel-to-micron conversion
# Coordinates are assumed to be stored as [z, y, x] in pixels
# when a single coordinates column is present.
# ==========================================================
X_SCALE = 0.0645
Y_SCALE = 0.0645
Z_SCALE = 0.2

# ==========================================================
# Plot settings
# ==========================================================
MAX_DISTANCE_UM = 2.5
BIN_WIDTH_UM = 0.05

MRNA_COLORS = {
    "MS2": "purple",
    "ATP2": "red",
    "ATP3": "red",
    "TIM50": "red",
}

COORD_FOLDER_NAME = "converted_coordinates"
SUMMARY_SUFFIX = "_spot_counts.csv"
SUMMARY_EXCEL_SUFFIX = "_spot_counts.xlsx"
NN_SUFFIX = "_NN_distance_um.npy"
NN_EXCEL_SUFFIX = "_NN_distance_um.xlsx"
PLOT_SUFFIX = "_NN_distribution.png"

# Excel has a hard row limit per sheet. If a distance table is very large,
# the script automatically splits it across multiple sheets in one workbook.
MAX_EXCEL_ROWS_PER_SHEET = 1_048_576

# If True, prints the columns from skipped files so you can diagnose format issues.
PRINT_SKIPPED_COLUMNS = True
MAX_COLUMN_EXAMPLES_PER_FOLDER = 3


def sanitize_name(name):
    name = str(name).strip()
    if name == "":
        name = "mRNA"
    for bad in ["\\", "/", ":", "*", "?", '"', "<", ">", "|"]:
        name = name.replace(bad, "_")
    name = re.sub(r"\s+", "_", name)
    return name


def safe_to_csv(df, path, index=False):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    try:
        df.to_csv(path, index=index)
        return path
    except PermissionError:
        from datetime import datetime

        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        fallback = path.with_name(f"{path.stem}_{stamp}{path.suffix}")
        df.to_csv(fallback, index=index)
        print(f"  Permission denied writing {path.name}; saved as {fallback.name}")
        return fallback


def safe_to_excel(df, path, index=False, sheet_base_name="Sheet"):
    """
    Save a dataframe as an Excel workbook.

    If the table is larger than Excel's per-sheet row limit, split it across
    multiple sheets in the same workbook. This keeps the .npy file as the
    downstream analysis file while making the values easy to inspect manually.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    def _write_excel(target_path):
        # One row is used by the header when index=False, so leave room for it.
        max_data_rows = MAX_EXCEL_ROWS_PER_SHEET - 1

        with pd.ExcelWriter(target_path, engine="openpyxl") as writer:
            if len(df) <= max_data_rows:
                df.to_excel(writer, sheet_name=sheet_base_name[:31], index=index)
            else:
                for sheet_i, start in enumerate(range(0, len(df), max_data_rows), start=1):
                    stop = start + max_data_rows
                    sheet_name = f"{sheet_base_name}_{sheet_i}"[:31]
                    df.iloc[start:stop].to_excel(writer, sheet_name=sheet_name, index=index)

    try:
        _write_excel(path)
        return path
    except PermissionError:
        from datetime import datetime

        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        fallback = path.with_name(f"{path.stem}_{stamp}{path.suffix}")
        _write_excel(fallback)
        print(f"  Permission denied writing {path.name}; saved as {fallback.name}")
        return fallback


def is_inside_generated_folder(path):
    lower_parts = [p.lower() for p in Path(path).parts]
    generated = {COORD_FOLDER_NAME.lower()}
    return any(part in generated for part in lower_parts)


def find_spots_extraction_folders(parent):
    return [
        p
        for p in parent.rglob("*")
        if p.is_dir() and p.name.lower() == "spots_extraction"
    ]


def infer_mrna_name(spots_folder):
    known_names = ["MS2", "ATP2", "ATP3", "TIM50"]

    for part in reversed(spots_folder.parts):
        upper = part.upper()
        for name in known_names:
            if re.search(rf"(^|[_\-\s]){re.escape(name)}([_\-\s]|$)", upper):
                return name

    try:
        return sanitize_name(spots_folder.parents[1].name)
    except Exception:
        return "mRNA"


def convert_coord_value(coord):
    """
    Convert a coordinate string/list in [z, y, x] pixel order to [x, y, z] microns.

    Handles:
      "[12, 345, 678]"
      "(12, 345, 678)"
      "12,345,678"
      "12 345 678"
    """
    try:
        if pd.isna(coord):
            return None

        if isinstance(coord, str):
            s = coord.strip()

            if s == "":
                return None

            try:
                parsed = ast.literal_eval(s)
            except Exception:
                s_clean = (
                    s.replace("[", "")
                    .replace("]", "")
                    .replace("(", "")
                    .replace(")", "")
                )
                if "," in s_clean:
                    parsed = [v.strip() for v in s_clean.split(",")]
                else:
                    parsed = s_clean.split()

            z, y, x = parsed
        else:
            z, y, x = coord

        x = float(x) * X_SCALE
        y = float(y) * Y_SCALE
        z = float(z) * Z_SCALE

        return [x, y, z]

    except Exception:
        return None


def standardize_columns(df):
    """
    Strip whitespace and normalize unnamed columns.
    """
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]
    return df


def read_spot_table_any_format(file_path):
    """
    Read spot extraction files robustly.

    Supports CSV and Excel. For CSV, tries multiple headers and separators.
    Returns a dataframe, even if coordinate columns still need to be detected later.
    """
    suffix = file_path.suffix.lower()

    if suffix == ".csv":
        attempts = []

        # Common case: normal CSV with first row as header.
        for header in [0, 1, 2, 3, 4, 5]:
            attempts.extend(
                [
                    {"header": header, "sep": None, "engine": "python"},
                    {"header": header, "sep": ","},
                    {"header": header, "sep": "\t"},
                    {"header": header, "sep": ";"},
                ]
            )

        # No header fallback.
        attempts.extend(
            [
                {"header": None, "sep": None, "engine": "python"},
                {"header": None, "sep": ","},
                {"header": None, "sep": "\t"},
                {"header": None, "sep": ";"},
            ]
        )

        best = None
        best_cols = -1

        for kwargs in attempts:
            try:
                df = pd.read_csv(file_path, **kwargs)
                df = standardize_columns(df)

                # Prefer the parse with the most columns.
                if df.shape[1] > best_cols:
                    best = df
                    best_cols = df.shape[1]

                # Return immediately if it already has recognizable coordinate info.
                coords = extract_coordinates_from_df(df)
                if coords is not None and len(coords) > 0:
                    return df

            except Exception:
                pass

        return best

    if suffix in [".xlsx", ".xlsm", ".xls"]:
        best = None
        for header in [0, 1, 2, 3, 4, 5]:
            try:
                df = pd.read_excel(file_path, header=header, engine="openpyxl")
                df = standardize_columns(df)
                best = df

                coords = extract_coordinates_from_df(df)
                if coords is not None and len(coords) > 0:
                    return df
            except Exception:
                pass

        return best

    return None


def find_col_case_insensitive(df, candidates):
    """
    Find a column whose normalized name exactly matches one of candidates.
    """
    norm_to_col = {}
    for col in df.columns:
        norm = str(col).strip().lower()
        norm_to_col[norm] = col

    for c in candidates:
        if c.lower() in norm_to_col:
            return norm_to_col[c.lower()]

    return None


def find_coordinate_triplet_columns(df):
    """
    Try to identify separate z/y/x columns.

    Supports common column names:
      z,y,x
      z_px,y_px,x_px
      spot_z,spot_y,spot_x
      axis-0,axis-1,axis-2
      axis_0,axis_1,axis_2
      coordinate-0,coordinate-1,coordinate-2

    Returns columns in z,y,x order.
    """
    columns_lower = {str(c).strip().lower(): c for c in df.columns}

    candidate_sets = [
        (["z"], ["y"], ["x"]),
        (
            ["z_px", "z pixel", "z_pixels", "spot_z", "spot z"],
            ["y_px", "y pixel", "y_pixels", "spot_y", "spot y"],
            ["x_px", "x pixel", "x_pixels", "spot_x", "spot x"],
        ),
        (
            ["axis-0", "axis_0", "axis 0", "dim-0", "dim_0"],
            ["axis-1", "axis_1", "axis 1", "dim-1", "dim_1"],
            ["axis-2", "axis_2", "axis 2", "dim-2", "dim_2"],
        ),
        (
            ["coordinate-0", "coordinate_0", "coordinate 0"],
            ["coordinate-1", "coordinate_1", "coordinate 1"],
            ["coordinate-2", "coordinate_2", "coordinate 2"],
        ),
        (["0"], ["1"], ["2"]),
    ]

    for z_names, y_names, x_names in candidate_sets:
        z_col = find_col_case_insensitive(df, z_names)
        y_col = find_col_case_insensitive(df, y_names)
        x_col = find_col_case_insensitive(df, x_names)

        if z_col is not None and y_col is not None and x_col is not None:
            return z_col, y_col, x_col

    # Fallback: look for columns containing z/y/x plus coordinate/spot.
    z_col = y_col = x_col = None
    for col in df.columns:
        c = str(col).strip().lower()
        if z_col is None and "z" in c and ("coord" in c or "spot" in c or "pixel" in c):
            z_col = col
        if y_col is None and "y" in c and ("coord" in c or "spot" in c or "pixel" in c):
            y_col = col
        if x_col is None and "x" in c and ("coord" in c or "spot" in c or "pixel" in c):
            x_col = col

    if z_col is not None and y_col is not None and x_col is not None:
        return z_col, y_col, x_col

    return None


def extract_coordinates_from_df(df):
    """
    Extract coordinates from a dataframe and return Nx3 array in x_um, y_um, z_um.

    Supports:
      1) single 'coordinates' column storing [z,y,x]
      2) separate z/y/x columns
      3) unnamed three-column CSV fallback
    """
    if df is None or df.shape[0] == 0:
        return None

    df = standardize_columns(df)

    # Case 1: single coordinates column.
    coord_col = find_col_case_insensitive(
        df,
        [
            "coordinates",
            "coordinate",
            "coords",
            "coord",
            "spot_coordinates",
            "spot coordinates",
            "spot_coordinate",
        ],
    )

    if coord_col is not None:
        converted = df[coord_col].dropna().apply(convert_coord_value).dropna()

        if len(converted) > 0:
            return np.array(converted.tolist(), dtype=float)

    # Case 2: separate z/y/x columns.
    triplet = find_coordinate_triplet_columns(df)
    if triplet is not None:
        z_col, y_col, x_col = triplet

        tmp = df[[z_col, y_col, x_col]].copy()
        for col in [z_col, y_col, x_col]:
            tmp[col] = pd.to_numeric(tmp[col], errors="coerce")

        tmp = tmp.dropna()

        if len(tmp) > 0:
            z = tmp[z_col].to_numpy(dtype=float) * Z_SCALE
            y = tmp[y_col].to_numpy(dtype=float) * Y_SCALE
            x = tmp[x_col].to_numpy(dtype=float) * X_SCALE

            return np.column_stack([x, y, z])

    # Case 3: headerless three-column fallback.
    if df.shape[1] == 3:
        tmp = df.copy()
        for col in tmp.columns:
            tmp[col] = pd.to_numeric(tmp[col], errors="coerce")
        tmp = tmp.dropna()

        if len(tmp) > 0:
            # Assume first three columns are z,y,x.
            z = tmp.iloc[:, 0].to_numpy(dtype=float) * Z_SCALE
            y = tmp.iloc[:, 1].to_numpy(dtype=float) * Y_SCALE
            x = tmp.iloc[:, 2].to_numpy(dtype=float) * X_SCALE

            return np.column_stack([x, y, z])

    return None


def find_spot_files(spots_folder):
    """
    Finds CSV and Excel spot files inside spots_extraction.
    Excludes files generated by this script.
    """
    files = []

    for ext in ["*.csv", "*.xlsx", "*.xlsm", "*.xls"]:
        files.extend(spots_folder.rglob(ext))

    filtered = []
    for p in files:
        name = p.name
        lower_name = name.lower()

        if name.startswith("~$"):
            continue
        if is_inside_generated_folder(p):
            continue
        if lower_name.endswith(SUMMARY_SUFFIX.lower()):
            continue
        if lower_name.endswith(SUMMARY_EXCEL_SUFFIX.lower()):
            continue
        if lower_name.endswith(NN_EXCEL_SUFFIX.lower()):
            continue
        if lower_name.endswith("_xyz_um.csv"):
            continue
        if lower_name.endswith("_xyz_um.xlsx"):
            continue
        if lower_name == "nn_distance_recursive_run_summary.csv":
            continue
        if lower_name == "nn_distance_recursive_run_summary.xlsx":
            continue
        if "recursive_run_summary" in lower_name:
            continue
        if "no_nn_distances" in lower_name:
            continue

        filtered.append(p)

    return sorted(filtered)


def report_file_extensions(spots_folder):
    counts = {}
    for p in spots_folder.rglob("*"):
        if p.is_file() and not is_inside_generated_folder(p):
            suffix = p.suffix.lower() if p.suffix else "[no extension]"
            counts[suffix] = counts.get(suffix, 0) + 1

    if not counts:
        print("  No files found inside this spots_extraction folder.")
        return

    print("  File extensions found inside this spots_extraction folder:")
    for ext, count in sorted(counts.items()):
        print(f"    {ext}: {count}")


def extract_cell_index_from_name(name):
    """
    Extract trailing cell index from names like:
      spots_extractions_..._000.csv
      spots_extractions_..._017.xlsx
    """
    stem = Path(str(name)).stem
    m = re.search(r"[_\-](\d{1,4})$", stem)
    if m:
        return int(m.group(1))
    return None


def make_short_output_name(file_path, spots_folder):
    """
    Create a short output base name to avoid Windows long-path failures.

    The previous script used the entire original filename, which produced paths like:
      converted_coordinates/spots_extractions_MS2_+_ATP2_..._000_xyz_um.csv

    That can exceed Windows path-length limits in deeply nested smallFISH folders.
    """
    cell_i = extract_cell_index_from_name(file_path.name)

    if cell_i is not None:
        return f"cell_{cell_i:03d}"

    # Fallback: use a short sanitized stem.
    stem = Path(file_path).stem
    safe = str(stem)
    for bad in ["\\", "/", ":", "*", "?", '"', "<", ">", "|"]:
        safe = safe.replace(bad, "_")
    safe = re.sub(r"\s+", "_", safe)

    if len(safe) > 40:
        safe = safe[-40:]

    return safe


def analyze_spots_folder(spots_folder):
    mrna_name = infer_mrna_name(spots_folder)
    safe_mrna = sanitize_name(mrna_name)

    coord_folder = spots_folder / COORD_FOLDER_NAME
    coord_folder.mkdir(parents=True, exist_ok=True)

    spot_files = find_spot_files(spots_folder)

    print(f"\nProcessing: {spots_folder}")
    print(f"  mRNA/channel name: {safe_mrna}")
    print(f"  Found {len(spot_files)} spot table files (.csv/.xlsx/.xlsm/.xls)")

    if len(spot_files) == 0:
        report_file_extensions(spots_folder)

    all_nn = []
    summary = []
    processed_files = 0
    skipped_files = 0
    column_examples_printed = 0

    for file_path in spot_files:
        filename = file_path.name

        try:
            df = read_spot_table_any_format(file_path)
            coords = extract_coordinates_from_df(df)

            if coords is None or len(coords) == 0:
                skipped_files += 1
                print(f"  Skipping {filename}: no usable coordinate columns found")

                if (
                    PRINT_SKIPPED_COLUMNS
                    and column_examples_printed < MAX_COLUMN_EXAMPLES_PER_FOLDER
                ):
                    if df is None:
                        print("    Could not read file into a table.")
                    else:
                        print(f"    Columns found: {list(df.columns)}")
                        print(f"    Shape: {df.shape}")
                        if df.shape[0] > 0:
                            print(f"    First row: {df.iloc[0].to_dict()}")
                    column_examples_printed += 1

                continue

            short_name = make_short_output_name(file_path, spots_folder)
            output_csv = coord_folder / f"{short_name}_xyz_um.csv"
            output_csv.parent.mkdir(parents=True, exist_ok=True)

            coord_df = pd.DataFrame(coords, columns=["x_um", "y_um", "z_um"])
            safe_to_csv(coord_df, output_csv, index=False)

            n_spots = len(coords)

            summary.append(
                {
                    "SourceFolder": str(spots_folder),
                    "mRNA": safe_mrna,
                    "Filename": filename,
                    "RelativeInputPath": str(file_path.relative_to(spots_folder)),
                    "SpotCount": n_spots,
                    "ConvertedCoordinatesFile": str(output_csv),
                }
            )

            processed_files += 1

            if n_spots < 2:
                print(f"  {filename}: {n_spots} spots; skipped NN distance")
                continue

            tree = cKDTree(coords)
            dist, _ = tree.query(coords, k=2)
            nn = dist[:, 1]

            all_nn.extend(nn)

            print(f"  {filename}: {n_spots} spots")

        except Exception as e:
            skipped_files += 1
            print(f"  Error: {filename}")
            print(f"    {e}")

    all_nn = np.array(all_nn, dtype=float)

    npy_file = spots_folder / f"{safe_mrna}{NN_SUFFIX}"
    np.save(npy_file, all_nn)

    nn_excel_df = pd.DataFrame(
        {
            "DistanceIndex": np.arange(1, len(all_nn) + 1, dtype=int),
            f"{safe_mrna}_NN_distance_um": all_nn,
        }
    )
    nn_excel_file = spots_folder / f"{safe_mrna}{NN_EXCEL_SUFFIX}"
    actual_nn_excel_file = safe_to_excel(
        nn_excel_df,
        nn_excel_file,
        index=False,
        sheet_base_name=f"{safe_mrna}_NN_um",
    )

    summary_df = pd.DataFrame(summary)
    summary_file = spots_folder / f"{safe_mrna}{SUMMARY_SUFFIX}"
    actual_summary_file = safe_to_csv(summary_df, summary_file, index=False)
    summary_excel_file = spots_folder / f"{safe_mrna}{SUMMARY_EXCEL_SUFFIX}"
    actual_summary_excel_file = safe_to_excel(
        summary_df,
        summary_excel_file,
        index=False,
        sheet_base_name=f"{safe_mrna}_spot_counts",
    )

    plot_file = spots_folder / f"{safe_mrna}{PLOT_SUFFIX}"

    if len(all_nn) > 0:
        plt.figure(figsize=(12, 8), facecolor="white")

        bin_edges = np.arange(0, MAX_DISTANCE_UM + BIN_WIDTH_UM, BIN_WIDTH_UM)

        color = MRNA_COLORS.get(safe_mrna.upper(), "gray")

        plt.hist(all_nn, bins=bin_edges, color=color, edgecolor="black", linewidth=1)

        plt.xlabel("Nearest-neighbor distance (µm)", fontsize=18)
        plt.ylabel("Frequency", fontsize=18)
        plt.title(f"{safe_mrna} nearest-neighbor distances", fontsize=18)

        plt.xticks(np.arange(0, MAX_DISTANCE_UM + 0.1, 0.5), fontsize=14)
        plt.yticks(fontsize=14)
        plt.xlim(0, MAX_DISTANCE_UM)

        plt.tight_layout()
        plt.savefig(plot_file, dpi=300, facecolor="white", bbox_inches="tight")
        plt.close()

    else:
        note_file = spots_folder / f"{safe_mrna}_NO_NN_DISTANCES.txt"
        note_file.write_text(
            "No nearest-neighbor distances were calculated. "
            "This usually means all files had fewer than two valid spots "
            "or no valid spot CSV/Excel files were found.\n",
            encoding="utf-8",
        )

    print(f"  Saved distances .npy:  {npy_file}")
    print(f"  Saved distances Excel: {actual_nn_excel_file}")
    print(f"  Saved summary CSV:     {actual_summary_file}")
    print(f"  Saved summary Excel:   {actual_summary_excel_file}")
    if len(all_nn) > 0:
        print(f"  Saved plot:            {plot_file}")
    else:
        print("  No histogram created because there were no NN distances.")

    return {
        "spots_folder": str(spots_folder),
        "mrna": safe_mrna,
        "spot_files_found": len(spot_files),
        "processed_files": processed_files,
        "skipped_files": skipped_files,
        "files_with_spot_counts": len(summary),
        "pooled_nn_distances": len(all_nn),
        "summary_file": str(actual_summary_file),
        "summary_excel_file": str(actual_summary_excel_file),
        "distance_file": str(npy_file),
        "distance_excel_file": str(actual_nn_excel_file),
        "plot_file": str(plot_file) if len(all_nn) > 0 else "",
    }


def main():
    if not parent_folder.exists():
        raise FileNotFoundError(f"Parent folder does not exist: {parent_folder}")

    spots_folders = find_spots_extraction_folders(parent_folder)

    print(f"Parent folder: {parent_folder}")
    print(f"Found {len(spots_folders)} spots_extraction folders")

    if len(spots_folders) == 0:
        print("No spots_extraction folders found.")
        return

    run_summary = []

    for spots_folder in spots_folders:
        result = analyze_spots_folder(spots_folder)
        run_summary.append(result)

    run_summary_df = pd.DataFrame(run_summary)
    run_summary_path = parent_folder / "NN_distance_recursive_run_summary.csv"
    actual_run_summary = safe_to_csv(run_summary_df, run_summary_path, index=False)

    run_summary_excel_path = parent_folder / "NN_distance_recursive_run_summary.xlsx"
    actual_run_summary_excel = safe_to_excel(
        run_summary_df,
        run_summary_excel_path,
        index=False,
        sheet_base_name="run_summary",
    )

    print("\nDone.")
    print(f"spots_extraction folders processed: {len(spots_folders)}")
    print(f"Overall run summary CSV:   {actual_run_summary}")
    print(f"Overall run summary Excel: {actual_run_summary_excel}")


if __name__ == "__main__":
    main()
