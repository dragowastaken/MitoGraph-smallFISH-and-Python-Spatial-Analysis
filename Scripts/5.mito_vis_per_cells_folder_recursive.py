import glob
import os
import sys
import tkinter as tk
from pathlib import Path
from tkinter import filedialog

import numpy as np
import pandas as pd
import plotly.graph_objects as go

print(sys.executable)
print(sys.version)

# ==========================================
# Parent folder to search recursively
# ==========================================
root = tk.Tk()
root.withdraw()

parent_folder = Path(filedialog.askdirectory())

print("User selected:", parent_folder)

# ==========================================
# Output folder name created inside each cells folder
# ==========================================
output_folder_name = "mito_visualization_interactive"

# ==========================================
# Whether to search inside subfolders of each cells folder
# This is useful because MitoGraph may save txt output in nested folders.
# ==========================================
search_recursively_inside_cells = True

# ==========================================
# Required columns in MitoGraph txt files
# ==========================================
required_columns = {"x", "y", "z", "width_(um)"}


def is_cells_folder(path: Path) -> bool:
    return path.is_dir() and path.name.lower() == "cells"


def is_inside_generated_output(path: Path) -> bool:
    parts_lower = [p.lower() for p in path.parts]
    return output_folder_name.lower() in parts_lower


def find_cells_folders(parent: Path):
    return [
        p
        for p in parent.rglob("*")
        if is_cells_folder(p) and not is_inside_generated_output(p)
    ]


def find_txt_files(cells_folder: Path):
    if search_recursively_inside_cells:
        files = [
            p for p in cells_folder.rglob("*.txt") if not is_inside_generated_output(p)
        ]
    else:
        files = [p for p in cells_folder.glob("*.txt")]
    return files


def report_file_extensions(folder: Path):
    counts = {}
    for p in folder.rglob("*"):
        if p.is_file() and not is_inside_generated_output(p):
            suffix = p.suffix.lower() if p.suffix else "[no extension]"
            counts[suffix] = counts.get(suffix, 0) + 1

    if not counts:
        print("  No files found anywhere inside this cells folder.")
        return

    print("  File extensions found inside this cells folder:")
    for ext, count in sorted(counts.items()):
        print(f"    {ext}: {count}")


def make_safe_relative_name(file_path: Path, base_folder: Path) -> str:
    """
    If txt files are inside nested folders, this prevents duplicate HTML names
    by including the relative folder path in the output filename.
    """
    rel = file_path.relative_to(base_folder)
    no_suffix = rel.with_suffix("")
    safe = str(no_suffix)
    for bad in ["\\", "/", ":", "*", "?", '"', "<", ">", "|"]:
        safe = safe.replace(bad, "_")
    return safe


cells_folders = find_cells_folders(parent_folder)

print(f"Found {len(cells_folders)} cells folders")

processed_folders = 0
created_html = 0
skipped_or_errors = 0

for cells_folder in cells_folders:
    processed_folders += 1
    print(f"\nProcessing cells folder: {cells_folder}")

    output_folder = cells_folder / output_folder_name
    output_folder.mkdir(exist_ok=True)

    txt_files = find_txt_files(cells_folder)
    print(f"  Found {len(txt_files)} txt files")

    if len(txt_files) == 0:
        report_file_extensions(cells_folder)
        continue

    for file_path in txt_files:
        try:
            # ----------------------------------
            # Read tab-delimited MitoGraph output
            # ----------------------------------
            df = pd.read_csv(file_path, sep="\t")

            missing = required_columns.difference(df.columns)
            if missing:
                print(f"  Skipping {file_path.name}: missing columns {sorted(missing)}")
                print(f"    Columns found: {list(df.columns)}")
                skipped_or_errors += 1
                continue

            # ----------------------------------
            # Coordinates and width
            # ----------------------------------
            x = df["x"].values
            y = df["y"].values
            z = df["z"].values
            width = df["width_(um)"].values

            # Avoid zero/negative marker sizes
            marker_size = np.maximum(width * 100, 1)

            # ----------------------------------
            # File name
            # ----------------------------------
            file_name = make_safe_relative_name(file_path, cells_folder)

            # ----------------------------------
            # Interactive 3D plot
            # ----------------------------------
            fig = go.Figure()

            fig.add_trace(
                go.Scatter3d(
                    x=x,
                    y=y,
                    z=z,
                    mode="markers",
                    marker=dict(
                        size=marker_size,
                        color=width,
                        colorscale="Viridis_r",
                        opacity=0.8,
                        colorbar=dict(title="Width (um)"),
                    ),
                )
            )

            fig.update_layout(
                title=file_name,
                scene=dict(
                    xaxis_title="X",
                    yaxis_title="Y",
                    zaxis_title="Z",
                    aspectmode="data",
                ),
                width=900,
                height=800,
            )

            # ----------------------------------
            # Save interactive HTML
            # ----------------------------------
            output_path = output_folder / f"{file_name}.html"
            fig.write_html(str(output_path))

            created_html += 1
            print(f"  Saved: {output_path}")

        except Exception as e:
            skipped_or_errors += 1
            print(f"  Error processing {file_path}")
            print(f"    {e}")

print("\nDone.")
print(f"Cells folders processed: {processed_folders}")
print(f"HTML files created: {created_html}")
print(f"Files skipped/errors: {skipped_or_errors}")
