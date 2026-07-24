import math
import tkinter as tk
from pathlib import Path
from tkinter import filedialog

import numpy as np
import pandas as pd


# ==========================================================
# Recursive .npy to .xlsx converter
# ==========================================================
# Purpose
# -------
# Select one parent folder, recursively find every .npy file inside it,
# and save a readable .xlsx copy beside each .npy file.
#
# Example:
#   ATP2_NN_distance_um.npy
#
# becomes:
#   ATP2_NN_distance_um.xlsx
#
# The original .npy file is never deleted or modified.
# ==========================================================


# ==========================================================
# Folder picker
# ==========================================================
root = tk.Tk()
root.withdraw()

parent_folder = Path(
    filedialog.askdirectory(
        title="Select parent folder containing .npy files"
    )
)

print("User selected:", parent_folder)


# ==========================================================
# User settings
# ==========================================================
# If False, existing .xlsx outputs are skipped.
# If True, existing .xlsx outputs are replaced.
OVERWRITE_EXISTING = False

# Safer default: do not load pickled object arrays.
# Most numeric analysis outputs, including NN-distance arrays, do not need pickle.
ALLOW_PICKLE = False

# Write a metadata sheet in every workbook.
WRITE_METADATA_SHEET = True

# Write a parent-level conversion summary.
WRITE_PARENT_SUMMARY = True

# Excel limits.
EXCEL_MAX_ROWS = 1_048_576
EXCEL_MAX_COLUMNS = 16_384

# Keep each data sheet below Excel's hard row limit.
# One row is reserved for the header, so data rows per sheet are one less.
MAX_DATA_ROWS_PER_SHEET = EXCEL_MAX_ROWS - 1

# Skip these folders while searching.
SKIP_FOLDER_KEYWORDS = [
    "__pycache__",
    ".ipynb_checkpoints",
    "npy_to_xlsx_conversion_summary",
]


# ==========================================================
# Helper functions
# ==========================================================
def is_inside_skipped_folder(path):
    lower_parts = [str(part).lower() for part in Path(path).parts]
    joined = "/".join(lower_parts)

    return any(keyword.lower() in joined for keyword in SKIP_FOLDER_KEYWORDS)


def find_npy_files(parent):
    return sorted(
        [
            p
            for p in Path(parent).rglob("*.npy")
            if p.is_file()
            and not p.name.startswith("~$")
            and not is_inside_skipped_folder(p)
        ]
    )


def make_output_path(npy_path):
    return Path(npy_path).with_suffix(".xlsx")


def safe_sheet_name(name):
    bad_chars = ["\\", "/", "?", "*", "[", "]", ":"]
    out = str(name)

    for bad in bad_chars:
        out = out.replace(bad, "_")

    out = out.strip()

    if out == "":
        out = "Sheet"

    return out[:31]


def flatten_array_to_dataframe(arr):
    """
    Convert any-dimensional array to a long-form table.

    Columns:
      flat_index
      dim_0
      dim_1
      ...
      value
    """
    arr = np.asarray(arr)
    flat = arr.ravel()

    if arr.ndim == 0:
        return pd.DataFrame(
            {
                "flat_index": [0],
                "value": [arr.item()],
            }
        )

    unravel = np.array(np.unravel_index(np.arange(flat.size), arr.shape)).T

    data = {
        "flat_index": np.arange(flat.size, dtype=np.int64),
    }

    for dim_i in range(arr.ndim):
        data[f"dim_{dim_i}"] = unravel[:, dim_i]

    data["value"] = flat

    return pd.DataFrame(data)


def array_to_dataframes(arr):
    """
    Return one or more named dataframes for writing to Excel.

    1D arrays:
      index, value

    2D arrays:
      matrix format if it fits Excel columns/rows.
      otherwise flattened long format.

    0D or >=3D arrays:
      flattened long format with dimension indices.
    """
    arr = np.asarray(arr)

    if arr.ndim == 0:
        return [("data", pd.DataFrame({"value": [arr.item()]}))]

    if arr.ndim == 1:
        df = pd.DataFrame(
            {
                "index": np.arange(arr.shape[0], dtype=np.int64),
                "value": arr,
            }
        )
        return [("data", df)]

    if arr.ndim == 2:
        n_rows, n_cols = arr.shape

        # Excel sheet also needs the index column if included.
        if n_rows <= MAX_DATA_ROWS_PER_SHEET and n_cols + 1 <= EXCEL_MAX_COLUMNS:
            columns = [f"col_{i:04d}" for i in range(n_cols)]
            df = pd.DataFrame(arr, columns=columns)
            df.insert(0, "row_index", np.arange(n_rows, dtype=np.int64))
            return [("data_matrix", df)]

        flat_df = flatten_array_to_dataframe(arr)
        return split_large_dataframe(flat_df, base_sheet_name="data_flat")

    flat_df = flatten_array_to_dataframe(arr)
    return split_large_dataframe(flat_df, base_sheet_name="data_flat")


def split_large_dataframe(df, base_sheet_name="data"):
    """
    Split a dataframe across multiple sheets if it exceeds Excel's row limit.
    """
    if len(df) <= MAX_DATA_ROWS_PER_SHEET:
        return [(base_sheet_name, df)]

    out = []
    n_sheets = int(math.ceil(len(df) / MAX_DATA_ROWS_PER_SHEET))

    for i in range(n_sheets):
        start = i * MAX_DATA_ROWS_PER_SHEET
        end = min((i + 1) * MAX_DATA_ROWS_PER_SHEET, len(df))
        sheet = f"{base_sheet_name}_{i + 1:03d}"
        out.append((sheet, df.iloc[start:end].copy()))

    return out


def make_metadata_dataframe(npy_path, arr, output_path):
    arr = np.asarray(arr)

    return pd.DataFrame(
        [
            {
                "npy_file": str(npy_path),
                "xlsx_file": str(output_path),
                "array_shape": str(tuple(arr.shape)),
                "array_ndim": int(arr.ndim),
                "array_dtype": str(arr.dtype),
                "array_size": int(arr.size),
                "overwrite_existing": OVERWRITE_EXISTING,
                "allow_pickle": ALLOW_PICKLE,
            }
        ]
    )


def write_array_to_xlsx(npy_path, arr, output_path):
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    dataframes = array_to_dataframes(arr)

    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        if WRITE_METADATA_SHEET:
            metadata_df = make_metadata_dataframe(npy_path, arr, output_path)
            metadata_df.to_excel(
                writer,
                sheet_name="metadata",
                index=False,
            )

        for sheet_name, df in dataframes:
            df.to_excel(
                writer,
                sheet_name=safe_sheet_name(sheet_name),
                index=False,
            )

    return output_path


def save_summary(summary_rows, parent):
    summary_df = pd.DataFrame(summary_rows)

    csv_path = Path(parent) / "npy_to_xlsx_conversion_summary.csv"
    xlsx_path = Path(parent) / "npy_to_xlsx_conversion_summary.xlsx"

    summary_df.to_csv(csv_path, index=False)

    try:
        summary_df.to_excel(xlsx_path, index=False, sheet_name="summary")
    except Exception as e:
        print(f"Could not write Excel summary: {e}")

    return csv_path, xlsx_path


# ==========================================================
# Main
# ==========================================================
def main():
    if not parent_folder.exists():
        raise FileNotFoundError(f"Selected parent folder does not exist: {parent_folder}")

    npy_files = find_npy_files(parent_folder)

    print(f"\nFound {len(npy_files)} .npy file(s).")

    summary_rows = []

    for i, npy_path in enumerate(npy_files, start=1):
        output_path = make_output_path(npy_path)

        print(f"\n[{i}/{len(npy_files)}] {npy_path}")

        if output_path.exists() and not OVERWRITE_EXISTING:
            print(f"  Skipped because output already exists: {output_path.name}")
            summary_rows.append(
                {
                    "npy_file": str(npy_path),
                    "xlsx_file": str(output_path),
                    "status": "skipped_xlsx_already_exists",
                    "array_shape": "",
                    "array_ndim": "",
                    "array_dtype": "",
                    "array_size": "",
                    "message": "Set OVERWRITE_EXISTING = True to replace existing output.",
                }
            )
            continue

        try:
            arr = np.load(npy_path, allow_pickle=ALLOW_PICKLE)
            arr_np = np.asarray(arr)

            write_array_to_xlsx(npy_path, arr_np, output_path)

            print(f"  Saved: {output_path.name}")
            print(f"  Shape: {arr_np.shape}; dtype: {arr_np.dtype}")

            summary_rows.append(
                {
                    "npy_file": str(npy_path),
                    "xlsx_file": str(output_path),
                    "status": "converted",
                    "array_shape": str(tuple(arr_np.shape)),
                    "array_ndim": int(arr_np.ndim),
                    "array_dtype": str(arr_np.dtype),
                    "array_size": int(arr_np.size),
                    "message": "",
                }
            )

        except Exception as e:
            print(f"  ERROR: {e}")
            summary_rows.append(
                {
                    "npy_file": str(npy_path),
                    "xlsx_file": str(output_path),
                    "status": "error",
                    "array_shape": "",
                    "array_ndim": "",
                    "array_dtype": "",
                    "array_size": "",
                    "message": str(e),
                }
            )

    if WRITE_PARENT_SUMMARY:
        csv_path, xlsx_path = save_summary(summary_rows, parent_folder)
        print("\nSaved conversion summary:")
        print(f"  {csv_path}")
        print(f"  {xlsx_path}")

    converted = sum(1 for row in summary_rows if row["status"] == "converted")
    skipped = sum(1 for row in summary_rows if row["status"].startswith("skipped"))
    errors = sum(1 for row in summary_rows if row["status"] == "error")

    print("\nDone.")
    print(f"Converted: {converted}")
    print(f"Skipped: {skipped}")
    print(f"Errors: {errors}")


if __name__ == "__main__":
    main()
