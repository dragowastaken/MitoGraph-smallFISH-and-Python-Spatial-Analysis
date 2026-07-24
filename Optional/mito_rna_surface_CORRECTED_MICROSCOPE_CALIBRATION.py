import ast
import re
import tkinter as tk
from pathlib import Path
from tkinter import filedialog

import numpy as np
import pandas as pd
import plotly.graph_objects as go

try:
    from scipy.spatial import cKDTree

    HAS_SCIPY = True
except Exception:
    HAS_SCIPY = False


# ==========================================
# Parent folder to search recursively
# ==========================================
root = tk.Tk()
root.withdraw()

parent_folder = Path(filedialog.askdirectory(title="Select parent folder"))

print("User selected:", parent_folder)


# ==========================================================
# CORRECTED MICROSCOPE / MITOGRAPH CALIBRATION
# ==========================================================
# This version is for datasets where MitoGraph was processed with the
# corrected microscope XY calibration.
#
# IMPORTANT:
# The RNA coordinates and MitoGraph surface coordinates must be in the same
# physical coordinate system. If MitoGraph was run with:
#
#   MitoGraph.exe -xy 0.0645 -z 0.2 ...
#
# then the raw smallFISH pixel coordinates should also be converted using:
#
#   X_SCALE = 0.0645 um/pixel
#   Y_SCALE = 0.0645 um/pixel
#
# NOTE:
# You mentioned ".00645" in the request. For this microscope workflow, the
# value used elsewhere in the procedure is 0.0645 um/pixel. If your actual
# confirmed microscope calibration is truly 0.0065 um/pixel, change only the
# CORRECTED_XY_SCALE_UM_PER_PIXEL constant below.
#
# This corrected version removes the old MS2-positive-control scale/offset
# calibration:
#
#   old X/Y scale: 0.05805 um/pixel
#   old dx/dy offset: -0.125, -0.250 um
#
# and instead uses the microscope/MitoGraph calibration directly.
# ==========================================================

FIXED_TRANSFORM = "flip_y"

X_SCALE = 0.0645
Y_SCALE = 0.0645
Z_SCALE = 0.2

# Set these to nonzero only if a new post-correction alignment check shows
# that a global x/y/z offset is still needed.
GLOBAL_DX_UM = 0.000
GLOBAL_DY_UM = 0.000
GLOBAL_DZ_UM = 0.000


# ==========================================================
# Crop size used for flip_y
# ==========================================================
# Must match the CropCells crop size.
CROP_WIDTH_PIXELS = 200
CROP_HEIGHT_PIXELS = 200


# ==========================================================
# Channels to analyze
# ==========================================================
RNA_CHANNELS = ["MS2", "ATP2", "ATP3", "Tim50"]


# ==========================================================
# Colocalization thresholds
# ==========================================================
COLOCALIZATION_THRESHOLDS_UM = [0.25, 0.5, 0.75, 1.0]


# ==========================================================
# Output
# ==========================================================
OUTPUT_FOLDER_NAME = "mito_rna_surface_CORRECTED_XY_0p065"

MAKE_HTML_OVERLAYS = True
SAVE_TRANSFORMED_COORDINATES = True
SHOW_NEAREST_SURFACE_LINKS = True


# ==========================================================
# Strict input filtering
# ==========================================================
# Only raw smallFISH spots_extractions files should be read.
RAW_SPOT_FILE_REQUIRED_SUBSTRING = "spots_extractions"

SKIP_FOLDER_KEYWORDS = [
    "random_ms2_output",
    "random_atp2_output",
    "random_atp3_output",
    "random_tim50_output",
    "random_output",
    "converted_coordinates",
    "mito_rna_3d_overlay",
    "mito_rna_colocalization",
    "mito_visualization_interactive",
    "mito_rna_surface_overlay",
    "mito_rna_surface_scale_offset",
    "mito_rna_surface_global_calibration",
    "mito_rna_surface_corrected_xy_0p065",
    OUTPUT_FOLDER_NAME.lower(),
]


def sanitize_name(name):
    name = str(name).strip()
    if name == "":
        name = "unknown"

    for bad in ["\\", "/", ":", "*", "?", '"', "<", ">", "|"]:
        name = name.replace(bad, "_")

    name = re.sub(r"\s+", "_", name)
    return name


def is_inside_skipped_folder(path):
    parts = [str(p).lower() for p in Path(path).parts]
    joined = "/".join(parts)

    for key in SKIP_FOLDER_KEYWORDS:
        if key.lower() in joined:
            return True

    return False


def extract_cell_index(path_or_name):
    """
    Extract true cell index from names like:
      spots_extractions_..._000.csv
      ..._000_mitosurface.vtk
      ..._000.txt
      cell_000_random_MS2.csv
    """
    stem = Path(str(path_or_name)).stem

    m = re.search(
        r"(?:^|[_\-])cell[_\-](\d{1,4})(?:[_\-]|$)", stem, flags=re.IGNORECASE
    )
    if m:
        return int(m.group(1))

    m = re.search(
        r"[_\-](\d{1,4})(?:[_\-](?:mitosurface|skeleton|nodes)|$)",
        stem,
        flags=re.IGNORECASE,
    )
    if m:
        return int(m.group(1))

    m = re.search(
        r"[_\-](\d{1,4})(?:$|[_\-](?:random|MS2|ATP2|ATP3|Tim50))",
        stem,
        flags=re.IGNORECASE,
    )
    if m:
        return int(m.group(1))

    m = re.search(r"[_\-](\d{1,4})$", stem)
    if m:
        return int(m.group(1))

    return None


def find_series_folders(parent):
    parent = Path(parent)

    # Allow user to select a Series folder directly.
    if re.match(r"^series(\s+\d+|\s*$|\d+)", parent.name, flags=re.IGNORECASE):
        return [parent]

    out = []

    for p in parent.rglob("*"):
        if not p.is_dir():
            continue
        if is_inside_skipped_folder(p):
            continue

        if re.match(r"^series(\s+\d+|\s*$|\d+)", p.name, flags=re.IGNORECASE):
            out.append(p)

    out = sorted(set(out))

    # If no Series folders are found, treat the selected folder as one analysis unit.
    # This makes the script easier to troubleshoot on copied test folders.
    if len(out) == 0:
        return [parent]

    return out


def find_cells_folders(series_folder):
    return sorted(
        [
            p
            for p in series_folder.rglob("*")
            if p.is_dir()
            and p.name.lower() == "cells"
            and not is_inside_skipped_folder(p)
        ]
    )


def find_mitosurface_vtk_files(cells_folder):
    return sorted(
        [
            p
            for p in cells_folder.rglob("*_mitosurface.vtk")
            if p.is_file() and not is_inside_skipped_folder(p)
        ]
    )


def find_skeleton_txt_files(cells_folder):
    return sorted(
        [
            p
            for p in cells_folder.rglob("*.txt")
            if p.is_file() and not is_inside_skipped_folder(p)
        ]
    )


def read_mitograph_txt(path):
    df = pd.read_csv(path, sep="\t")
    df.columns = [str(c).strip() for c in df.columns]

    required = {"x", "y", "z"}

    if not required.issubset(df.columns):
        raise ValueError(f"MitoGraph file missing required columns {required}: {path}")

    if "width_(um)" not in df.columns:
        df["width_(um)"] = 0.1

    out = df[["x", "y", "z", "width_(um)"]].copy()

    for col in ["x", "y", "z", "width_(um)"]:
        out[col] = pd.to_numeric(out[col], errors="coerce")

    out = out.dropna(subset=["x", "y", "z"])
    return out


def path_part_matches_channel(part, channel):
    p = str(part).strip().upper()
    c = str(channel).strip().upper()

    if p == c:
        return True
    if p.startswith(c + "_"):
        return True
    if p.startswith(c + " "):
        return True
    if p.startswith(c + "-"):
        return True

    return False


def spots_folder_belongs_to_channel(spots_folder, series_folder, channel):
    """
    Strict channel matching based only on the path relative to the Series folder.
    This prevents parent folder names like 'MS2 + ATP2,ATP3,Tim50' from matching both channels.
    """
    try:
        rel = spots_folder.relative_to(series_folder)
        parts = rel.parts
    except Exception:
        parts = spots_folder.parts

    useful_parts = []

    for part in parts:
        useful_parts.append(part)

        if str(part).lower() == "spots_extraction":
            break

    return any(path_part_matches_channel(part, channel) for part in useful_parts)


def find_channel_spots_folders(series_folder, channel):
    folders = []

    for p in series_folder.rglob("spots_extraction"):
        if not p.is_dir():
            continue
        if is_inside_skipped_folder(p):
            continue

        if spots_folder_belongs_to_channel(p, series_folder, channel):
            folders.append(p)

    return sorted(folders)


def is_raw_spot_file(path):
    if is_inside_skipped_folder(path):
        return False

    if path.name.startswith("~$"):
        return False

    lower_name = path.name.lower()
    lower_path = str(path).lower()

    if RAW_SPOT_FILE_REQUIRED_SUBSTRING not in lower_name:
        return False

    generated_terms = [
        "random_",
        "_random_",
        "spot_counts",
        "summary",
        "nn_distance",
        "overlay_coordinates",
        "transform_scores",
        "pairing_transform",
        "colocalization",
        "mito_distance",
        "surface_coordinates",
        "scale_offset",
        "converted",
    ]

    if any(term in lower_name for term in generated_terms):
        return False

    if any(
        term in lower_path
        for term in [
            "random_ms2_output",
            "random_atp2_output",
            "random_atp3_output",
            "random_tim50_output",
            "random_output",
        ]
    ):
        return False

    return True


def build_spot_index(series_folder):
    """
    Returns:
      index = {channel: {cell_index: raw_spot_file_path}}
      accepted_df = accepted raw spot files
      rejected_df = rejected candidate files
    """
    index = {channel: {} for channel in RNA_CHANNELS}
    accepted_rows = []
    rejected_rows = []

    for channel in RNA_CHANNELS:
        spots_folders = find_channel_spots_folders(series_folder, channel)

        for spots_folder in spots_folders:
            all_candidates = []

            for ext in ["*.csv", "*.xlsx", "*.xlsm", "*.xls"]:
                all_candidates.extend(spots_folder.rglob(ext))

            for file_path in sorted(all_candidates):
                if is_raw_spot_file(file_path):
                    cell_i = extract_cell_index(file_path.name)

                    accepted_rows.append(
                        {
                            "SeriesFolder": str(series_folder),
                            "Channel": channel,
                            "SpotsFolder": str(spots_folder),
                            "SpotFile": str(file_path),
                            "CellIndex": cell_i,
                        }
                    )

                    if cell_i is None:
                        continue

                    if cell_i in index[channel]:
                        print(
                            f"WARNING: duplicate raw {channel} cell index {cell_i}; keeping first file:"
                        )
                        print(f"  kept:    {index[channel][cell_i]}")
                        print(f"  ignored: {file_path}")
                    else:
                        index[channel][cell_i] = file_path
                else:
                    rejected_rows.append(
                        {
                            "SeriesFolder": str(series_folder),
                            "Channel": channel,
                            "SpotsFolder": str(spots_folder),
                            "RejectedFile": str(file_path),
                            "Reason": "not_raw_spots_extractions_file_or_generated_output",
                        }
                    )

    return index, pd.DataFrame(accepted_rows), pd.DataFrame(rejected_rows)


def read_spot_table(path):
    suffix = path.suffix.lower()

    if suffix == ".csv":
        for sep in [";", ",", "\t"]:
            try:
                df = pd.read_csv(path, sep=sep)
                if df.shape[1] > 1:
                    df.columns = [str(c).strip() for c in df.columns]
                    return df
            except Exception:
                pass

        df = pd.read_csv(path, sep=None, engine="python")
        df.columns = [str(c).strip() for c in df.columns]
        return df

    if suffix in [".xlsx", ".xlsm", ".xls"]:
        df = pd.read_excel(path)
        df.columns = [str(c).strip() for c in df.columns]
        return df

    raise ValueError(f"Unsupported spot file type: {path}")


def parse_coordinate_triplet(value):
    """
    Parse smallFISH coordinate entry in (z, y, x) pixel order.
    """
    if pd.isna(value):
        return None

    if isinstance(value, str):
        s = value.strip()

        if s == "":
            return None

        try:
            parsed = ast.literal_eval(s)
        except Exception:
            s = s.replace("[", "").replace("]", "")
            s = s.replace("(", "").replace(")", "")

            if "," in s:
                parsed = [v.strip() for v in s.split(",")]
            else:
                parsed = s.split()

        if len(parsed) < 3:
            return None

        z, y, x = parsed[:3]
        return float(z), float(y), float(x)

    try:
        z, y, x = value
        return float(z), float(y), float(x)
    except Exception:
        return None


def find_column(df, candidates):
    lower_to_col = {str(c).strip().lower(): c for c in df.columns}

    for c in candidates:
        if c.lower() in lower_to_col:
            return lower_to_col[c.lower()]

    return None


def read_rna_spots_px(path):
    """
    Reads raw smallFISH spot coordinates in pixels as z_px, y_px, x_px.
    """
    df = read_spot_table(path)
    df.columns = [str(c).strip() for c in df.columns]

    coord_col = find_column(
        df,
        [
            "coordinates",
            "coordinate",
            "coords",
            "coord",
            "spot_coordinates",
            "spot coordinates",
        ],
    )

    if coord_col is None:
        raise ValueError(f"No raw coordinates column found in: {path}")

    coords = []

    for value in df[coord_col].dropna():
        parsed = parse_coordinate_triplet(value)

        if parsed is None:
            continue

        z_px, y_px, x_px = parsed

        coords.append(
            {
                "z_px": z_px,
                "y_px": y_px,
                "x_px": x_px,
            }
        )

    return pd.DataFrame(coords)


def transform_pixels_to_um(px_df):
    """
    Convert raw smallFISH pixels to microns using the corrected microscope calibration.

    Raw smallFISH coordinates are interpreted as:
      z_px, y_px, x_px

    Converted coordinates are:
      x_um, y_um, z_um

    Default transform:
      FIXED_TRANSFORM = "flip_y"

    Supported transform choices:
      "none"    -> x = x_px * X_SCALE;               y = y_px * Y_SCALE
      "flip_y"  -> x = x_px * X_SCALE;               y = crop_height_um - y_px * Y_SCALE
      "flip_x"  -> x = crop_width_um - x_px*X_SCALE; y = y_px * Y_SCALE
      "flip_xy" -> x = crop_width_um - x_px*X_SCALE; y = crop_height_um - y_px * Y_SCALE

    The default offset is now zero because MitoGraph and RNA coordinates should
    already be in the same corrected microscope coordinate system.
    """
    x_unflipped = px_df["x_px"].to_numpy(dtype=float) * X_SCALE
    y_unflipped = px_df["y_px"].to_numpy(dtype=float) * Y_SCALE
    z = px_df["z_px"].to_numpy(dtype=float) * Z_SCALE

    crop_w_um = CROP_WIDTH_PIXELS * X_SCALE
    crop_h_um = CROP_HEIGHT_PIXELS * Y_SCALE

    transform = str(FIXED_TRANSFORM).strip().lower()

    if transform == "none":
        x = x_unflipped
        y = y_unflipped
    elif transform == "flip_y":
        x = x_unflipped
        y = crop_h_um - y_unflipped
    elif transform == "flip_x":
        x = crop_w_um - x_unflipped
        y = y_unflipped
    elif transform == "flip_xy":
        x = crop_w_um - x_unflipped
        y = crop_h_um - y_unflipped
    else:
        raise ValueError(
            "Unsupported FIXED_TRANSFORM. Use one of: 'none', 'flip_y', 'flip_x', 'flip_xy'."
        )

    out = pd.DataFrame(
        {
            "x_um": x + GLOBAL_DX_UM,
            "y_um": y + GLOBAL_DY_UM,
            "z_um": z + GLOBAL_DZ_UM,
            "x_um_before_offset": x,
            "y_um_before_offset": y,
            "z_um_before_offset": z,
            "x_px": px_df["x_px"].to_numpy(dtype=float),
            "y_px": px_df["y_px"].to_numpy(dtype=float),
            "z_px": px_df["z_px"].to_numpy(dtype=float),
            "x_scale_um_per_px": X_SCALE,
            "y_scale_um_per_px": Y_SCALE,
            "z_scale_um_per_px": Z_SCALE,
            "crop_width_pixels": CROP_WIDTH_PIXELS,
            "crop_height_pixels": CROP_HEIGHT_PIXELS,
            "crop_width_um": crop_w_um,
            "crop_height_um": crop_h_um,
            "global_dx_um": GLOBAL_DX_UM,
            "global_dy_um": GLOBAL_DY_UM,
            "global_dz_um": GLOBAL_DZ_UM,
            "transform": FIXED_TRANSFORM,
            "calibration_version": "corrected_microscope_xy_0p065",
        }
    )

    return out


def _parse_binary_legacy_vtk_polydata(vtk_path):
    """
    Minimal parser for MitoGraph legacy binary VTK POLYDATA files.
    Reads:
      POINTS
      POLYGONS
    """
    data = Path(vtk_path).read_bytes()

    pts_match = re.search(rb"POINTS\s+(\d+)\s+(\w+)\n", data)
    if pts_match is None:
        raise ValueError(f"Could not find POINTS section in {vtk_path}")

    n_points = int(pts_match.group(1))
    dtype_name = pts_match.group(2).decode("ascii").lower()

    if dtype_name not in {"float", "double"}:
        raise ValueError(f"Unsupported VTK point dtype '{dtype_name}' in {vtk_path}")

    pts_start = pts_match.end()
    itemsize = 4 if dtype_name == "float" else 8
    np_dtype = ">f4" if dtype_name == "float" else ">f8"

    pts_n_values = n_points * 3
    pts_n_bytes = pts_n_values * itemsize

    points = np.frombuffer(
        data[pts_start : pts_start + pts_n_bytes], dtype=np_dtype
    ).astype(float)
    points = points.reshape(n_points, 3)

    remainder = data[pts_start + pts_n_bytes :]
    poly_match = re.search(rb"POLYGONS\s+(\d+)\s+(\d+)\n", remainder)

    if poly_match is None:
        return points, np.empty((0, 3), dtype=int)

    n_polys = int(poly_match.group(1))
    poly_size = int(poly_match.group(2))
    poly_start = pts_start + pts_n_bytes + poly_match.end()

    poly_raw = np.frombuffer(
        data[poly_start : poly_start + poly_size * 4], dtype=">i4"
    ).astype(int)

    faces = []
    idx = 0

    for _ in range(n_polys):
        nverts = poly_raw[idx]
        idx += 1

        verts = poly_raw[idx : idx + nverts]
        idx += nverts

        if nverts == 3:
            faces.append(verts.tolist())
        elif nverts > 3:
            v0 = int(verts[0])
            for j in range(1, nverts - 1):
                faces.append([v0, int(verts[j]), int(verts[j + 1])])

    return points, np.asarray(faces, dtype=int)


def read_mitosurface_vtk(vtk_path):
    return _parse_binary_legacy_vtk_polydata(vtk_path)


def nearest_surface_vertex_info(surface_points, query_points):
    """
    Approximate RNA-to-mitosurface distance using nearest surface vertex.
    """
    if len(query_points) == 0 or len(surface_points) == 0:
        return np.array([]), np.empty((0, 3), dtype=float)

    if HAS_SCIPY:
        tree = cKDTree(surface_points)
        dist, idx = tree.query(query_points, k=1)
        nearest = surface_points[idx]
        return np.asarray(dist, dtype=float), np.asarray(nearest, dtype=float)

    distances = []
    nearest = []

    for p in query_points:
        d = np.sqrt(np.sum((surface_points - p) ** 2, axis=1))
        j = int(np.argmin(d))
        distances.append(float(d[j]))
        nearest.append(surface_points[j])

    return np.asarray(distances, dtype=float), np.asarray(nearest, dtype=float)


def summarize_distances(distances):
    distances = np.asarray(distances, dtype=float)
    distances = distances[np.isfinite(distances)]

    out = {
        "Nspots": int(len(distances)),
        "MeanDistanceToSurface_um": np.nan,
        "MedianDistanceToSurface_um": np.nan,
        "MinDistanceToSurface_um": np.nan,
        "MaxDistanceToSurface_um": np.nan,
    }

    if len(distances) > 0:
        out["MeanDistanceToSurface_um"] = float(np.mean(distances))
        out["MedianDistanceToSurface_um"] = float(np.median(distances))
        out["MinDistanceToSurface_um"] = float(np.min(distances))
        out["MaxDistanceToSurface_um"] = float(np.max(distances))

    for threshold in COLOCALIZATION_THRESHOLDS_UM:
        key = str(threshold).replace(".", "p")

        if len(distances) > 0:
            out[f"FractionWithin{key}um"] = float(np.mean(distances <= threshold))
            out[f"SpotsWithin{key}um"] = int(np.sum(distances <= threshold))
        else:
            out[f"FractionWithin{key}um"] = np.nan
            out[f"SpotsWithin{key}um"] = 0

    return out


def make_surface_overlay_figure(
    surface_points, surface_faces, rna_traces, skeleton_df=None, title=""
):
    fig = go.Figure()

    if len(surface_faces) > 0:
        fig.add_trace(
            go.Mesh3d(
                x=surface_points[:, 0],
                y=surface_points[:, 1],
                z=surface_points[:, 2],
                i=surface_faces[:, 0],
                j=surface_faces[:, 1],
                k=surface_faces[:, 2],
                color="lightgray",
                opacity=0.28,
                name="Mito surface",
                hoverinfo="skip",
                lighting=dict(ambient=0.6, diffuse=0.7, roughness=0.9, specular=0.1),
            )
        )
    else:
        fig.add_trace(
            go.Scatter3d(
                x=surface_points[:, 0],
                y=surface_points[:, 1],
                z=surface_points[:, 2],
                mode="markers",
                marker=dict(size=1.5, color="lightgray", opacity=0.55),
                name="Mito surface vertices",
            )
        )

    if skeleton_df is not None and len(skeleton_df) > 0:
        fig.add_trace(
            go.Scatter3d(
                x=skeleton_df["x"],
                y=skeleton_df["y"],
                z=skeleton_df["z"],
                mode="markers",
                marker=dict(size=2.0, color="black", opacity=0.50),
                name="Mito skeleton",
            )
        )

    channel_colors = {
        "MS2": "magenta",
        "ATP2": "green",
        "ATP3": "green",
        "Tim50": "green",
    }

    for trace in rna_traces:
        channel = trace["channel"]
        coords = trace["coords"]
        distances = trace["distances"]
        nearest_pts = trace["nearest_pts"]
        rna_file = trace["rna_file"]

        hover = []

        for i in range(len(coords)):
            if i < len(distances):
                d_txt = f"{distances[i]:.3f} um"
            else:
                d_txt = "NA"

            hover.append(
                f"{channel} spot {i}<br>"
                f"distance to mito surface = {d_txt}<br>"
                f"x={coords['x_um'].iloc[i]:.3f}, y={coords['y_um'].iloc[i]:.3f}, z={coords['z_um'].iloc[i]:.3f}<br>"
                f"file={Path(rna_file).name}"
            )

        fig.add_trace(
            go.Scatter3d(
                x=coords["x_um"],
                y=coords["y_um"],
                z=coords["z_um"],
                mode="markers",
                marker=dict(
                    size=5,
                    color=channel_colors.get(channel.upper(), "green"),
                    opacity=0.95,
                ),
                name=f"{channel} spots",
                text=hover,
                hoverinfo="text",
            )
        )

        if SHOW_NEAREST_SURFACE_LINKS and len(nearest_pts) == len(coords):
            xs = []
            ys = []
            zs = []

            for idx_row in range(len(coords)):
                xs.extend([coords["x_um"].iloc[idx_row], nearest_pts[idx_row, 0], None])
                ys.extend([coords["y_um"].iloc[idx_row], nearest_pts[idx_row, 1], None])
                zs.extend([coords["z_um"].iloc[idx_row], nearest_pts[idx_row, 2], None])

            fig.add_trace(
                go.Scatter3d(
                    x=xs,
                    y=ys,
                    z=zs,
                    mode="lines",
                    line=dict(
                        width=2, color=channel_colors.get(channel.upper(), "green")
                    ),
                    opacity=0.30,
                    name=f"{channel} nearest-surface links",
                    hoverinfo="skip",
                    showlegend=False,
                )
            )

    fig.update_layout(
        title=title,
        scene=dict(
            xaxis_title="X (um)",
            yaxis_title="Y (um)",
            zaxis_title="Z (um)",
            aspectmode="data",
        ),
        width=1200,
        height=980,
    )

    return fig


def build_surface_file_index(cells_folder):
    surface_files = find_mitosurface_vtk_files(cells_folder)

    surface_by_cell = {}

    for sf in surface_files:
        cell_i = extract_cell_index(sf.name)

        if cell_i is not None and cell_i not in surface_by_cell:
            surface_by_cell[cell_i] = sf

    return surface_by_cell


def build_skeleton_file_index(cells_folder):
    skeleton_files = find_skeleton_txt_files(cells_folder)

    skeleton_by_cell = {}

    for sf in skeleton_files:
        cell_i = extract_cell_index(sf.name)

        if cell_i is not None and cell_i not in skeleton_by_cell:
            skeleton_by_cell[cell_i] = sf

    return skeleton_by_cell


def process_series(series_folder):
    print(f"\nSeries: {series_folder}")

    spot_index, accepted_df, rejected_df = build_spot_index(series_folder)
    cells_folders = find_cells_folders(series_folder)

    print(f"  cells folders: {len(cells_folders)}")

    for channel in RNA_CHANNELS:
        print(f"  accepted raw {channel} files: {len(spot_index[channel])}")

    series_summary_rows = []
    series_spot_rows = []

    for cells_folder in cells_folders:
        output_folder = cells_folder / OUTPUT_FOLDER_NAME
        output_folder.mkdir(parents=True, exist_ok=True)

        if len(accepted_df) > 0:
            accepted_df.to_csv(
                output_folder / "CORRECTED_XY_0p065_ACCEPTED_raw_spot_files.csv",
                index=False,
            )

        if len(rejected_df) > 0:
            rejected_df.to_csv(
                output_folder / "CORRECTED_XY_0p065_REJECTED_nonraw_files.csv",
                index=False,
            )

        surface_by_cell = build_surface_file_index(cells_folder)
        skeleton_by_cell = build_skeleton_file_index(cells_folder)

        print(f"  cells folder: {cells_folder}")
        print(f"    mitosurface VTK files: {len(surface_by_cell)}")
        print(f"    skeleton TXT files: {len(skeleton_by_cell)}")

        all_cell_indices = sorted(set(surface_by_cell.keys()))

        for cell_i in all_cell_indices:
            surface_file = surface_by_cell.get(cell_i)
            skeleton_file = skeleton_by_cell.get(cell_i)

            if surface_file is None:
                series_summary_rows.append(
                    {
                        "SeriesFolder": str(series_folder),
                        "CellsFolder": str(cells_folder),
                        "CellIndex": cell_i,
                        "Status": "missing_mitosurface_vtk",
                    }
                )
                continue

            try:
                surface_points, surface_faces = read_mitosurface_vtk(surface_file)
            except Exception as e:
                print(f"    Could not read mitosurface {surface_file.name}: {e}")
                series_summary_rows.append(
                    {
                        "SeriesFolder": str(series_folder),
                        "CellsFolder": str(cells_folder),
                        "CellIndex": cell_i,
                        "MitoSurfaceFile": str(surface_file),
                        "Status": f"surface_read_error: {e}",
                    }
                )
                continue

            skeleton_df = None

            if skeleton_file is not None:
                try:
                    skeleton_df = read_mitograph_txt(skeleton_file)
                except Exception:
                    skeleton_df = None

            rna_traces = []

            for channel in RNA_CHANNELS:
                rna_file = spot_index.get(channel, {}).get(cell_i)

                if rna_file is None:
                    series_summary_rows.append(
                        {
                            "SeriesFolder": str(series_folder),
                            "CellsFolder": str(cells_folder),
                            "CellIndex": cell_i,
                            "Channel": channel,
                            "MitoSurfaceFile": str(surface_file),
                            "RNAFile": "",
                            "Status": "no_same_index_raw_spot_file",
                        }
                    )
                    continue

                try:
                    px_df = read_rna_spots_px(rna_file)
                    coords = transform_pixels_to_um(px_df)

                    query_xyz = coords[["x_um", "y_um", "z_um"]].to_numpy(dtype=float)
                    distances, nearest_pts = nearest_surface_vertex_info(
                        surface_points, query_xyz
                    )

                    per_spot = coords.copy()
                    per_spot["DistanceToMitoSurface_um"] = distances

                    if len(nearest_pts) == len(per_spot):
                        per_spot["nearest_surface_x_um"] = nearest_pts[:, 0]
                        per_spot["nearest_surface_y_um"] = nearest_pts[:, 1]
                        per_spot["nearest_surface_z_um"] = nearest_pts[:, 2]

                    per_spot["SeriesFolder"] = str(series_folder)
                    per_spot["CellsFolder"] = str(cells_folder)
                    per_spot["CellIndex"] = cell_i
                    per_spot["Channel"] = channel
                    per_spot["MitoSurfaceFile"] = str(surface_file)
                    per_spot["RNAFile"] = str(rna_file)

                    for threshold in COLOCALIZATION_THRESHOLDS_UM:
                        key = str(threshold).replace(".", "p")
                        per_spot[f"within_{key}_um"] = (
                            per_spot["DistanceToMitoSurface_um"] <= threshold
                        )

                    series_spot_rows.append(per_spot)

                    if SAVE_TRANSFORMED_COORDINATES:
                        coord_out = (
                            output_folder
                            / f"cell_{cell_i:03d}_{channel}_CORRECTED_XY_0p065_surface_distances.csv"
                        )
                        per_spot.to_csv(coord_out, index=False)

                    summary = {
                        "SeriesFolder": str(series_folder),
                        "CellsFolder": str(cells_folder),
                        "CellIndex": cell_i,
                        "Channel": channel,
                        "MitoSurfaceFile": str(surface_file),
                        "RNAFile": str(rna_file),
                        "Transform": FIXED_TRANSFORM,
                        "X_SCALE_um_per_px": X_SCALE,
                        "Y_SCALE_um_per_px": Y_SCALE,
                        "Z_SCALE_um_per_px": Z_SCALE,
                        "GlobalDX_um": GLOBAL_DX_UM,
                        "GlobalDY_um": GLOBAL_DY_UM,
                        "GlobalDZ_um": GLOBAL_DZ_UM,
                        "Status": "processed_corrected_microscope_calibration",
                    }
                    summary.update(summarize_distances(distances))
                    series_summary_rows.append(summary)

                    rna_traces.append(
                        {
                            "channel": channel,
                            "coords": coords,
                            "distances": distances,
                            "nearest_pts": nearest_pts,
                            "rna_file": rna_file,
                        }
                    )

                except Exception as e:
                    print(f"    Error processing {channel} cell {cell_i:03d}: {e}")
                    series_summary_rows.append(
                        {
                            "SeriesFolder": str(series_folder),
                            "CellsFolder": str(cells_folder),
                            "CellIndex": cell_i,
                            "Channel": channel,
                            "MitoSurfaceFile": str(surface_file),
                            "RNAFile": str(rna_file),
                            "Status": f"rna_error: {e}",
                        }
                    )

            if MAKE_HTML_OVERLAYS and len(rna_traces) > 0:
                title = (
                    f"{series_folder.name} cell {cell_i:03d}: mito surface + RNA | "
                    f"CORRECTED microscope calibration: {FIXED_TRANSFORM}, XY={X_SCALE:.5f}, Z={Z_SCALE:.3f}, "
                    f"dx={GLOBAL_DX_UM:.3f}, dy={GLOBAL_DY_UM:.3f}, dz={GLOBAL_DZ_UM:.3f}"
                )

                fig = make_surface_overlay_figure(
                    surface_points=surface_points,
                    surface_faces=surface_faces,
                    rna_traces=rna_traces,
                    skeleton_df=skeleton_df,
                    title=title,
                )

                html_out = (
                    output_folder
                    / f"cell_{cell_i:03d}_CORRECTED_XY_0p065_surface_overlay.html"
                )
                fig.write_html(html_out)
                print(f"    Saved overlay: {html_out.name}")

    summary_df = pd.DataFrame(series_summary_rows)

    if len(series_spot_rows) > 0:
        spot_df = pd.concat(series_spot_rows, ignore_index=True)
    else:
        spot_df = pd.DataFrame()

    return summary_df, spot_df


def summarize_parent_cell_level(parent_cell_df):
    if len(parent_cell_df) == 0:
        return pd.DataFrame()

    valid = parent_cell_df[
        parent_cell_df["Status"].astype(str).isin(["processed_global_calibration", "processed_corrected_microscope_calibration"])
    ].copy()

    if len(valid) == 0:
        return pd.DataFrame()

    metric_cols = [
        "Nspots",
        "MeanDistanceToSurface_um",
        "MedianDistanceToSurface_um",
        "MinDistanceToSurface_um",
        "MaxDistanceToSurface_um",
    ]

    for threshold in COLOCALIZATION_THRESHOLDS_UM:
        key = str(threshold).replace(".", "p")
        metric_cols.append(f"FractionWithin{key}um")
        metric_cols.append(f"SpotsWithin{key}um")

    rows = []

    for channel, sub in valid.groupby("Channel"):
        row = {
            "Channel": channel,
            "Ncells": int(sub["CellIndex"].nunique()),
            "TotalSpots": int(
                pd.to_numeric(sub["Nspots"], errors="coerce").fillna(0).sum()
            ),
        }

        for col in metric_cols:
            vals = pd.to_numeric(sub[col], errors="coerce").dropna()

            if len(vals) == 0:
                row[f"{col}_mean_across_cells"] = np.nan
                row[f"{col}_median_across_cells"] = np.nan
                row[f"{col}_sem_across_cells"] = np.nan
            else:
                row[f"{col}_mean_across_cells"] = float(vals.mean())
                row[f"{col}_median_across_cells"] = float(vals.median())
                row[f"{col}_sem_across_cells"] = (
                    float(vals.std(ddof=1) / np.sqrt(len(vals)))
                    if len(vals) > 1
                    else np.nan
                )

        rows.append(row)

    return pd.DataFrame(rows)


def main():
    if not parent_folder.exists():
        raise FileNotFoundError(
            f"Selected parent folder does not exist: {parent_folder}"
        )

    parent_output = parent_folder / OUTPUT_FOLDER_NAME
    parent_output.mkdir(parents=True, exist_ok=True)
    settings_df = pd.DataFrame(
        [
            {
                "calibration_version": "corrected_microscope_xy_0p065",
                "corrected_xy_scale_um_per_pixel": CORRECTED_XY_SCALE_UM_PER_PIXEL,
                "corrected_z_scale_um_per_pixel": CORRECTED_Z_SCALE_UM_PER_PIXEL,
                "fixed_transform": FIXED_TRANSFORM,
                "x_scale_um_per_pixel_used_for_rna": X_SCALE,
                "y_scale_um_per_pixel_used_for_rna": Y_SCALE,
                "z_scale_um_per_pixel_used_for_rna": Z_SCALE,
                "global_dx_um": GLOBAL_DX_UM,
                "global_dy_um": GLOBAL_DY_UM,
                "global_dz_um": GLOBAL_DZ_UM,
                "crop_width_pixels": CROP_WIDTH_PIXELS,
                "crop_height_pixels": CROP_HEIGHT_PIXELS,
                "old_scale_removed_from_prior_script": "0.05805 um/pixel",
                "old_offsets_removed_from_prior_script": "dx=-0.125, dy=-0.250, dz=0.000 um",
                "parent_folder": str(parent_folder),
                "output_folder": str(parent_output),
            }
        ]
    )
    settings_df.to_csv(
        parent_output / "CORRECTED_XY_0p065_RUN_SETTINGS.csv",
        index=False,
    )


    series_folders = find_series_folders(parent_folder)

    print(f"\nFound {len(series_folders)} Series folders")
    print(f"Output folder name: {OUTPUT_FOLDER_NAME}")
    print(f"SciPy available: {HAS_SCIPY}")
    print("\nCorrected microscope/MitoGraph calibration:")
    print(f"  Transform: {FIXED_TRANSFORM}")
    print(f"  X_SCALE: {X_SCALE} um/px")
    print(f"  Y_SCALE: {Y_SCALE} um/px")
    print(f"  Z_SCALE: {Z_SCALE} um/px")
    print(f"  dx: {GLOBAL_DX_UM} um")
    print(f"  dy: {GLOBAL_DY_UM} um")
    print(f"  dz: {GLOBAL_DZ_UM} um")

    all_cell_summary = []
    all_spot_distances = []

    for series_folder in series_folders:
        cell_summary_df, spot_df = process_series(series_folder)

        if len(cell_summary_df) > 0:
            all_cell_summary.append(cell_summary_df)

        if len(spot_df) > 0:
            all_spot_distances.append(spot_df)

    if len(all_cell_summary) > 0:
        parent_cell_summary = pd.concat(all_cell_summary, ignore_index=True)
    else:
        parent_cell_summary = pd.DataFrame()

    if len(all_spot_distances) > 0:
        parent_spot_distances = pd.concat(all_spot_distances, ignore_index=True)
    else:
        parent_spot_distances = pd.DataFrame()

    parent_output.mkdir(parents=True, exist_ok=True)

    if len(parent_cell_summary) > 0:
        parent_cell_summary.to_csv(
            parent_output
            / "PARENT_cell_level_surface_colocalization_CORRECTED_XY_0p065.csv",
            index=False,
        )

    if len(parent_spot_distances) > 0:
        parent_spot_distances.to_csv(
            parent_output
            / "PARENT_spot_level_surface_distances_CORRECTED_XY_0p065.csv",
            index=False,
        )

    channel_summary = summarize_parent_cell_level(parent_cell_summary)

    if len(channel_summary) > 0:
        channel_summary.to_csv(
            parent_output / "PARENT_channel_summary_CORRECTED_XY_0p065.csv", index=False
        )

    print("\nDone.")
    print(f"Parent output folder: {parent_output}")
    print("\nMain outputs:")
    print(
        f"  {parent_output / 'PARENT_cell_level_surface_colocalization_GLOBAL_CALIBRATION.csv'}"
    )
    print(
        f"  {parent_output / 'PARENT_spot_level_surface_distances_GLOBAL_CALIBRATION.csv'}"
    )
    print(f"  {parent_output / 'PARENT_channel_summary_GLOBAL_CALIBRATION.csv'}")
    print("\nPer-cell HTML overlays are saved inside each cells folder:")
    print(f"  cells/{OUTPUT_FOLDER_NAME}/")
    print("\nThis script uses the corrected microscope/MitoGraph calibration for all cells/channels.")
    print("It does not perform per-cell scale/offset optimization.")
    print("If overlays are still shifted after this correction, adjust only GLOBAL_DX_UM/GLOBAL_DY_UM/GLOBAL_DZ_UM after visual QC.")


if __name__ == "__main__":
    main()
