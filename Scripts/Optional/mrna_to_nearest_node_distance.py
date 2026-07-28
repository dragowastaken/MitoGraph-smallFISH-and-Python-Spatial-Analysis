import re
import tkinter as tk
from collections import defaultdict
from pathlib import Path
from tkinter import filedialog

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

# ==========================================================
# Folder picker
# ==========================================================
# Select the parent strain/experiment folder that contains Series folders.
# The script analyzes each Series folder independently, then pools all Series.
# ==========================================================
root = tk.Tk()
root.withdraw()

parent_folder = Path(
    filedialog.askdirectory(
        title="Select parent folder containing Series folders"
    )
)

print("User selected:", parent_folder)

# ==========================================================
# Main settings
# ==========================================================
# Input mode:
#   True  = only use converted mRNA coordinate files like:
#           .../converted_coordinates/cell_000_xyz_um.csv
#           This is recommended.
#
#   False = recursively consider any CSV/XLSX table with x/y/z columns,
#           while still excluding derived outputs.
PROCESS_ONLY_CONVERTED_COORDINATES = True

# Distances are measured in the coordinate units present in your files.
# If both mRNA coordinate CSVs and node VTK files are already in microns,
# leave all scale factors at 1.0.
MRNA_X_SCALE_UM = 1.0
MRNA_Y_SCALE_UM = 1.0
MRNA_Z_SCALE_UM = 1.0

NODE_X_SCALE_UM = 1.0
NODE_Y_SCALE_UM = 1.0
NODE_Z_SCALE_UM = 1.0

# The MitoGraph *_nodes.vtk file is usually a mesh representation of nodes,
# not one point per node. This mode groups connected mesh components and uses
# each component center as the node coordinate.
#
# Options:
#   "connected_component_centers"  recommended for MitoGraph *_nodes.vtk files
#   "all_points"                  use every VTK POINT as a node coordinate
NODE_POSITION_MODE = "connected_component_centers"

# Keep True unless you have unusual filenames.
# This prevents accidental matching between the wrong mRNA cell and node file.
REQUIRE_CELL_INDEX_MATCH = True

# Count mRNAs within these nearest-node distance thresholds.
DISTANCE_THRESHOLDS_UM = [0.1, 0.25, 0.5, 1.0, 2.0, 5.0]

# Plot settings.
HIST_BINS = 50

# Short output folder name to avoid Windows path-length failures.
OUTPUT_FOLDER_NAME = "node_mrna_dist_BY_SERIES"

KNOWN_MRNA_NAMES = ["MS2", "ATP2", "ATP3", "TIM50"]

# ==========================================================
# File discovery settings
# ==========================================================
NODE_VTK_PATTERNS = [
    "*_nodes.vtk",
    "*nodes.vtk",
]

COORDINATE_TABLE_PATTERNS = [
    "*.csv",
    "*.xlsx",
    "*.xlsm",
    "*.xls",
]

GENERATED_FOLDER_KEYWORDS = {
    "node_mrna_distance_output",
    "node_mrna_distance_output_by_series",
    "node_mrna_dist_by_series",
    "node_mrna_dist",
    "mito_visualization_interactive",
}

RANDOM_FOLDER_KEYWORDS = {
    "random_ms2_output",
    "random_atp2_output",
    "random_atp3_output",
    "random_tim50_output",
    "random_output",
    "randomized_mrna_output",
}

DERIVED_DISTANCE_FOLDER_KEYWORDS = {
    "mito_rna_surface_global_calibration",
    "mito_rna_surface",
    "surface_distance",
    "surface_distances",
}

EXCLUDE_COORDINATE_FILENAME_KEYWORDS = [
    "spot_counts",
    "_counts",
    "count_table",
    "nn_distance",
    "nearest_node",
    "nearest_mrna",
    "run_summary",
    "per_cell_summary",
    "summary",
    "distribution",
    "histogram",
    "cumulative",
    "node_occupancy",
    "surface_distances",
    "global_calibration_surface",
]


# ==========================================================
# General utility functions
# ==========================================================
def sanitize_name(name):
    name = str(name).strip()
    if name == "":
        name = "unnamed"
    for bad in ["\\", "/", ":", "*", "?", '"', "<", ">", "|"]:
        name = name.replace(bad, "_")
    name = re.sub(r"\s+", "_", name)
    return name


def normalize_stem(name):
    s = Path(str(name)).stem.lower()

    removable = [
        "_nodes",
        "_node",
        "_vtk",
        "_xyz_um",
        "_coordinates_um",
        "_coordinates",
        "_coords_um",
        "_coords",
        "_spots",
        "spots_extractions_",
        "spots_extraction_",
        "converted_coordinates_",
        "converted_coordinate_",
        "_random",
    ]

    for token in removable:
        s = s.replace(token, "")

    s = re.sub(r"[^a-z0-9]+", "", s)
    return s


def extract_cell_index(name):
    """
    Extract the true cell index from either:
      - cell_000_xyz_um.csv
      - cell_000_ATP3_...csv
      - ..._000_nodes.vtk

    Important fix:
    This does NOT use earlier numbers in the ND2 filename, such as "-002-1".
    It prioritizes the final "_000_nodes" part for node files.
    """
    stem = Path(str(name)).stem

    # mRNA coordinate files: cell_000_xyz_um.csv, cell_000_ATP3_...
    m = re.search(r"(?:^|[_\-])cell[_\-]?(\d{1,4})(?:[_\-]|$)", stem, flags=re.IGNORECASE)
    if m:
        return int(m.group(1))

    # MitoGraph node files: long_name_000_nodes.vtk
    m = re.search(r"[_\-](\d{1,4})[_\-]nodes?$", stem, flags=re.IGNORECASE)
    if m:
        return int(m.group(1))

    # Other converted-coordinate style names.
    m = re.search(
        r"[_\-](\d{1,4})(?:[_\-]xyz[_\-]um|[_\-]coords|[_\-]coordinates)?$",
        stem,
        flags=re.IGNORECASE,
    )
    if m:
        return int(m.group(1))

    return None


def make_cell_label(cell_index):
    if cell_index is None:
        return "cell_unknown"
    return f"cell_{int(cell_index):03d}"


def infer_mrna_name(path):
    pieces = list(Path(path).parts)
    pieces.append(Path(path).name)

    for part in reversed(pieces):
        upper = part.upper()
        for name in KNOWN_MRNA_NAMES:
            if re.search(rf"(^|[_\-\s]){re.escape(name)}([_\-\s]|$)", upper):
                return name

    stem = Path(path).stem
    stem = re.sub(r"(_xyz_um|_coordinates|_coords|_spots)$", "", stem, flags=re.IGNORECASE)
    return sanitize_name(stem)


def is_inside_any_keyword_folder(path, keywords):
    lower_parts = {p.lower() for p in Path(path).parts}
    return any(k.lower() in lower_parts for k in keywords)


def is_inside_generated_or_derived_folder(path):
    return (
        is_inside_any_keyword_folder(path, GENERATED_FOLDER_KEYWORDS)
        or is_inside_any_keyword_folder(path, RANDOM_FOLDER_KEYWORDS)
        or is_inside_any_keyword_folder(path, DERIVED_DISTANCE_FOLDER_KEYWORDS)
    )


def find_series_folder(path):
    p = Path(path)
    for ancestor in p.parents:
        if ancestor.name.lower().startswith("series"):
            return ancestor
    return p.parent


def find_series_folders(parent):
    """
    Find Series folders under the selected parent.

    Each Series is analyzed independently. mRNA files from one Series are never
    matched to node files from another Series.
    """
    parent = Path(parent)

    if parent.name.lower().startswith("series"):
        return [parent]

    series_folders = [
        p
        for p in parent.rglob("*")
        if p.is_dir()
        and p.name.lower().startswith("series")
        and not is_inside_generated_or_derived_folder(p)
    ]

    series_folders = sorted(set(series_folders))

    # Avoid nested Series folders causing duplicated analysis.
    non_nested = []
    for series in series_folders:
        if not any(other != series and other in series.parents for other in series_folders):
            non_nested.append(series)

    if len(non_nested) == 0:
        return [parent]

    return non_nested


def make_series_output_label(series_folder, parent_folder):
    series_folder = Path(series_folder)
    parent_folder = Path(parent_folder)

    try:
        rel = series_folder.relative_to(parent_folder)
        label = "__".join(rel.parts)
    except Exception:
        label = series_folder.name

    return sanitize_name(label)


def format_distance_value(value):
    text = f"{float(value):g}"
    text = text.replace("-", "neg").replace(".", "p")
    return text


def table_has_xyz_columns(df):
    try:
        find_xyz_columns(df)
        return True
    except Exception:
        return False


# ==========================================================
# VTK parsing
# ==========================================================
def _binary_dtype_for_vtk(vtk_type):
    vtk_type = str(vtk_type).lower()
    if vtk_type in ["float", "float32"]:
        return ">f4"
    if vtk_type in ["double", "float64"]:
        return ">f8"
    raise ValueError(f"Unsupported VTK point type: {vtk_type}")


def read_vtk_points_and_polygons(vtk_file):
    """
    Read POINTS and POLYGONS from a legacy VTK POLYDATA file.

    Supports common MitoGraph binary VTK files:
        BINARY
        DATASET POLYDATA
        POINTS N float
        [big-endian binary point data]
        POLYGONS M K
        [big-endian int32 polygon connectivity]
    """
    data = Path(vtk_file).read_bytes()

    header_preview = data[:500].decode("latin1", errors="ignore").upper()
    is_binary = "BINARY" in header_preview

    points_match = re.search(rb"(?:^|\n)POINTS\s+(\d+)\s+([A-Za-z0-9_]+)\s*\n", data)
    if points_match is None:
        raise ValueError(f"No POINTS block found in VTK file: {vtk_file}")

    n_points = int(points_match.group(1))
    vtk_point_type = points_match.group(2).decode("ascii", errors="ignore")
    points_start = points_match.end()

    if is_binary:
        dtype = np.dtype(_binary_dtype_for_vtk(vtk_point_type))
        n_values = n_points * 3
        n_bytes = n_values * dtype.itemsize

        raw = data[points_start : points_start + n_bytes]
        if len(raw) != n_bytes:
            raise ValueError(
                f"VTK file ended before all POINTS were read: {vtk_file}"
            )

        points = np.frombuffer(raw, dtype=dtype).astype(float).reshape(n_points, 3)
        after_points = points_start + n_bytes

        polygons = []
        poly_match = re.search(
            rb"\nPOLYGONS\s+(\d+)\s+(\d+)\s*\n",
            data[after_points : after_points + 10000],
        )

        if poly_match is not None:
            n_polygons = int(poly_match.group(1))
            total_ints = int(poly_match.group(2))
            ints_start = after_points + poly_match.end()

            int_raw = data[ints_start : ints_start + total_ints * 4]
            if len(int_raw) == total_ints * 4:
                ints = np.frombuffer(int_raw, dtype=">i4").astype(int)
                i = 0
                for _ in range(n_polygons):
                    if i >= len(ints):
                        break
                    k = int(ints[i])
                    i += 1
                    poly = ints[i : i + k].tolist()
                    i += k
                    if len(poly) == k:
                        polygons.append(poly)

        return points, polygons

    # ASCII fallback.
    text = data.decode("latin1", errors="replace")
    point_header = re.search(
        r"(?:^|\n)POINTS\s+(\d+)\s+([A-Za-z0-9_]+)\s*\n",
        text,
        flags=re.IGNORECASE,
    )
    if point_header is None:
        raise ValueError(f"No ASCII POINTS block found in VTK file: {vtk_file}")

    n_points = int(point_header.group(1))
    rest = text[point_header.end() :]
    tokens = rest.split()

    values = []
    token_i = 0
    while token_i < len(tokens) and len(values) < n_points * 3:
        try:
            values.append(float(tokens[token_i]))
            token_i += 1
        except ValueError:
            break

    if len(values) < n_points * 3:
        raise ValueError(f"Could not read enough ASCII VTK points from: {vtk_file}")

    points = np.array(values, dtype=float).reshape(n_points, 3)
    polygons = []

    return points, polygons


def connected_component_centers_from_polygons(points, polygons):
    """
    Use polygon connectivity to group the VTK node mesh into connected components.
    Each component center is treated as one node.
    """
    n = len(points)
    if n == 0:
        return np.empty((0, 3), dtype=float), []

    if len(polygons) == 0:
        return points.copy(), [[i] for i in range(n)]

    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra = find(a)
        rb = find(b)
        if ra != rb:
            parent[rb] = ra

    for poly in polygons:
        if len(poly) < 2:
            continue
        first = poly[0]
        for v in poly[1:]:
            if 0 <= first < n and 0 <= v < n:
                union(first, v)

    components = defaultdict(list)
    for i in range(n):
        components[find(i)].append(i)

    component_lists = list(components.values())
    centers = np.array([points[idxs].mean(axis=0) for idxs in component_lists], dtype=float)

    return centers, component_lists


def load_node_positions(vtk_file):
    points, polygons = read_vtk_points_and_polygons(vtk_file)

    points = points.copy()
    points[:, 0] *= NODE_X_SCALE_UM
    points[:, 1] *= NODE_Y_SCALE_UM
    points[:, 2] *= NODE_Z_SCALE_UM

    if NODE_POSITION_MODE == "all_points":
        node_positions = points
        component_lists = [[i] for i in range(len(points))]
    elif NODE_POSITION_MODE == "connected_component_centers":
        node_positions, component_lists = connected_component_centers_from_polygons(
            points,
            polygons,
        )
    else:
        raise ValueError(
            "NODE_POSITION_MODE must be 'connected_component_centers' or 'all_points'."
        )

    node_positions = node_positions[np.isfinite(node_positions).all(axis=1)]
    if len(node_positions) == 0:
        raise ValueError(f"No valid node coordinates found in: {vtk_file}")

    return node_positions, {
        "raw_point_count": len(points),
        "polygon_count": len(polygons),
        "component_count": len(component_lists),
    }


# ==========================================================
# mRNA coordinate parsing
# ==========================================================
def read_raw_table(path):
    path = Path(path)
    suffix = path.suffix.lower()

    if suffix == ".csv":
        read_attempts = [
            {},
            {"sep": ","},
            {"sep": ";"},
            {"sep": "\t"},
            {"sep": None, "engine": "python"},
        ]

        last_error = None

        for kwargs in read_attempts:
            try:
                df = pd.read_csv(path, **kwargs)
                df.columns = [str(c).strip() for c in df.columns]

                # Avoid incorrectly accepting a semicolon-delimited table as one giant column.
                if df.shape[1] >= 3 and table_has_xyz_columns(df):
                    return df

                # Keep trying if only one column was parsed.
                if df.shape[1] == 1:
                    continue

            except Exception as e:
                last_error = e

        raise ValueError(f"Could not read CSV with usable x/y/z columns: {path}. Last error: {last_error}")

    if suffix in [".xlsx", ".xlsm", ".xls"]:
        df = pd.read_excel(path)
        df.columns = [str(c).strip() for c in df.columns]
        return df

    raise ValueError(f"Unsupported coordinate table type: {path}")


def normalized_column_name(col):
    s = str(col).strip().lower()
    s = s.replace("µ", "u")
    s = s.replace("μ", "u")
    s = re.sub(r"[^a-z0-9]+", "", s)
    return s


def find_xyz_columns(df):
    """
    Detect x/y/z coordinate columns in common coordinate table formats.
    """
    normalized = {col: normalized_column_name(col) for col in df.columns}

    x_names = {
        "x",
        "xum",
        "xmicron",
        "xmicrons",
        "xposition",
        "positionx",
        "xpos",
        "posx",
        "centroidx",
        "centerx",
        "spotx",
        "globalx",
    }
    y_names = {
        "y",
        "yum",
        "ymicron",
        "ymicrons",
        "yposition",
        "positiony",
        "ypos",
        "posy",
        "centroidy",
        "centery",
        "spoty",
        "globaly",
    }
    z_names = {
        "z",
        "zum",
        "zmicron",
        "zmicrons",
        "zposition",
        "positionz",
        "zpos",
        "posz",
        "centroidz",
        "centerz",
        "spotz",
        "globalz",
        "slice",
        "plane",
    }

    x_col = next((col for col, norm in normalized.items() if norm in x_names), None)
    y_col = next((col for col, norm in normalized.items() if norm in y_names), None)
    z_col = next((col for col, norm in normalized.items() if norm in z_names), None)

    if x_col is not None and y_col is not None and z_col is not None:
        return x_col, y_col, z_col

    # Fallback: first three mostly numeric columns.
    numeric_cols = []
    for col in df.columns:
        series = df[col].dropna()
        if len(series) == 0:
            continue
        numeric = pd.to_numeric(series, errors="coerce")
        if numeric.notna().mean() > 0.8:
            numeric_cols.append(col)

    if len(numeric_cols) >= 3:
        return numeric_cols[0], numeric_cols[1], numeric_cols[2]

    raise ValueError(
        "Could not identify x/y/z coordinate columns. "
        f"Columns found: {list(df.columns)}"
    )


def load_mrna_coordinates(path):
    df = read_raw_table(path)
    x_col, y_col, z_col = find_xyz_columns(df)

    coords = df[[x_col, y_col, z_col]].copy()
    coords.columns = ["x", "y", "z"]

    coords["x"] = pd.to_numeric(coords["x"], errors="coerce")
    coords["y"] = pd.to_numeric(coords["y"], errors="coerce")
    coords["z"] = pd.to_numeric(coords["z"], errors="coerce")

    coords = coords.dropna(subset=["x", "y", "z"])
    coords = coords[np.isfinite(coords[["x", "y", "z"]]).all(axis=1)]

    xyz = coords[["x", "y", "z"]].to_numpy(dtype=float)

    xyz[:, 0] *= MRNA_X_SCALE_UM
    xyz[:, 1] *= MRNA_Y_SCALE_UM
    xyz[:, 2] *= MRNA_Z_SCALE_UM

    if len(xyz) == 0:
        raise ValueError(f"No valid mRNA x/y/z rows found in: {path}")

    return xyz, {
        "x_column": x_col,
        "y_column": y_col,
        "z_column": z_col,
        "input_rows": len(df),
        "valid_xyz_rows": len(xyz),
    }


def looks_like_coordinate_table(path):
    path = Path(path)
    name = path.name.lower()
    lower_parts = {p.lower() for p in path.parts}

    if path.name.startswith("~$"):
        return False

    if is_inside_generated_or_derived_folder(path):
        return False

    if path.suffix.lower() not in [".csv", ".xlsx", ".xlsm", ".xls"]:
        return False

    for keyword in EXCLUDE_COORDINATE_FILENAME_KEYWORDS:
        if keyword in name:
            return False

    if PROCESS_ONLY_CONVERTED_COORDINATES:
        if "converted_coordinates" not in lower_parts:
            return False
        if not re.search(r"^cell[_\-]?\d{1,4}.*xyz.*um", name):
            return False

    return True


def find_node_vtk_files(series_folder):
    found = []
    for pattern in NODE_VTK_PATTERNS:
        found.extend(Path(series_folder).rglob(pattern))

    filtered = []
    for p in found:
        if p.name.startswith("~$"):
            continue
        if is_inside_generated_or_derived_folder(p):
            continue
        filtered.append(p)

    return sorted(set(filtered))


def find_coordinate_files(series_folder):
    if PROCESS_ONLY_CONVERTED_COORDINATES:
        found = list(Path(series_folder).rglob("cell_*_xyz_um.csv"))
    else:
        found = []
        for pattern in COORDINATE_TABLE_PATTERNS:
            found.extend(Path(series_folder).rglob(pattern))

    candidates = [p for p in found if looks_like_coordinate_table(p)]

    valid = []
    rejected = []

    for p in sorted(set(candidates)):
        try:
            df = read_raw_table(p)
            find_xyz_columns(df)
            valid.append(p)
        except Exception as e:
            rejected.append((p, str(e)))

    return valid, rejected


# ==========================================================
# File matching
# ==========================================================
def build_node_index(node_files):
    entries = []
    for vtk in node_files:
        entries.append(
            {
                "path": vtk,
                "series_folder": find_series_folder(vtk),
                "cell_index": extract_cell_index(vtk.name),
                "norm": normalize_stem(vtk.name),
            }
        )
    return entries


def find_matching_node_file(coord_file, node_entries):
    """
    Match an mRNA coordinate file to a node VTK file within the same Series.

    Primary match is by cell index:
        cell_014_xyz_um.csv -> *_014_nodes.vtk
    """
    if len(node_entries) == 0:
        return None

    coord_index = extract_cell_index(coord_file.name)
    coord_norm = normalize_stem(coord_file.name)

    if coord_index is not None:
        same_index = [
            e
            for e in node_entries
            if e["cell_index"] is not None and e["cell_index"] == coord_index
        ]

        if len(same_index) == 1:
            return same_index[0]["path"]

        if len(same_index) > 1:
            # Prefer an exact/partial normalized-name relationship if duplicated.
            exactish = [
                e
                for e in same_index
                if coord_norm and (coord_norm in e["norm"] or e["norm"] in coord_norm)
            ]
            if len(exactish) == 1:
                return exactish[0]["path"]

            # Otherwise choose the shortest filename as a stable fallback.
            return sorted(same_index, key=lambda e: len(e["path"].name))[0]["path"]

    if REQUIRE_CELL_INDEX_MATCH:
        return None

    # Conservative fallback if explicitly allowed.
    exact = [e for e in node_entries if e["norm"] == coord_norm]
    if len(exact) == 1:
        return exact[0]["path"]

    if len(node_entries) == 1:
        return node_entries[0]["path"]

    return None


# ==========================================================
# Distance calculations and output
# ==========================================================
def calculate_mrna_to_nearest_node(mrna_xyz, node_xyz):
    tree = cKDTree(node_xyz)
    distances, nearest_node_idx = tree.query(mrna_xyz, k=1)
    nearest_node_xyz = node_xyz[nearest_node_idx]

    return distances.astype(float), nearest_node_idx.astype(int), nearest_node_xyz


def summarize_distances(distances):
    distances = np.asarray(distances, dtype=float)

    if len(distances) == 0:
        summary = {
            "mRNA_count": 0,
            "mean_distance_um": np.nan,
            "median_distance_um": np.nan,
            "min_distance_um": np.nan,
            "max_distance_um": np.nan,
            "std_distance_um": np.nan,
        }
    else:
        summary = {
            "mRNA_count": int(len(distances)),
            "mean_distance_um": float(np.mean(distances)),
            "median_distance_um": float(np.median(distances)),
            "min_distance_um": float(np.min(distances)),
            "max_distance_um": float(np.max(distances)),
            "std_distance_um": float(np.std(distances, ddof=1)) if len(distances) > 1 else 0.0,
        }

    for threshold in DISTANCE_THRESHOLDS_UM:
        count = int(np.sum(distances <= threshold))
        percent = 100.0 * count / len(distances) if len(distances) > 0 else np.nan
        tag = format_distance_value(threshold)
        summary[f"count_within_{tag}_um"] = count
        summary[f"percent_within_{tag}_um"] = percent

    return summary


def summarize_group_distances(df, group_columns):
    rows = []

    if len(df) == 0:
        return pd.DataFrame(rows)

    grouped = df.groupby(group_columns, dropna=False)

    for keys, group in grouped:
        if not isinstance(keys, tuple):
            keys = (keys,)

        distances = group["nearest_node_distance_um"].to_numpy(dtype=float)
        summary = summarize_distances(distances)

        row = {}
        for col, value in zip(group_columns, keys):
            row[col] = value

        row.update(summary)
        rows.append(row)

    return pd.DataFrame(rows)


def make_histogram_plot(distances, title, out_path):
    distances = np.asarray(distances, dtype=float)

    plt.figure(figsize=(8, 6), facecolor="white")
    plt.hist(distances, bins=HIST_BINS, edgecolor="black")
    plt.xlabel("Distance from mRNA to nearest node (um)")
    plt.ylabel("mRNA count")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close()


def make_cumulative_threshold_plot(distances, title, out_path):
    distances = np.asarray(distances, dtype=float)
    thresholds = np.array(DISTANCE_THRESHOLDS_UM, dtype=float)
    counts = np.array([np.sum(distances <= t) for t in thresholds], dtype=int)

    plt.figure(figsize=(8, 6), facecolor="white")
    plt.plot(thresholds, counts, marker="o")
    plt.xlabel("Distance threshold (um)")
    plt.ylabel("Number of mRNAs within threshold")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close()


def make_xy_distance_scatter_plot(mrna_xyz, node_xyz, distances, title, out_path):
    distances = np.asarray(distances, dtype=float)

    plt.figure(figsize=(8, 7), facecolor="white")

    scatter = plt.scatter(
        mrna_xyz[:, 0],
        mrna_xyz[:, 1],
        c=distances,
        s=20,
        alpha=0.8,
    )
    plt.scatter(
        node_xyz[:, 0],
        node_xyz[:, 1],
        marker="x",
        s=25,
        alpha=0.8,
    )

    cbar = plt.colorbar(scatter)
    cbar.set_label("mRNA-to-nearest-node distance (um)")

    plt.xlabel("x (um)")
    plt.ylabel("y (um)")
    plt.title(title)
    plt.axis("equal")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close()


def make_node_occupancy_plot(node_occupancy_df, title, out_path):
    plt.figure(figsize=(8, 7), facecolor="white")

    counts = node_occupancy_df["assigned_mrna_count"].to_numpy(dtype=float)
    sizes = 20 + 20 * counts

    scatter = plt.scatter(
        node_occupancy_df["node_x"],
        node_occupancy_df["node_y"],
        c=counts,
        s=sizes,
        alpha=0.8,
    )

    cbar = plt.colorbar(scatter)
    cbar.set_label("mRNAs assigned to node")

    plt.xlabel("node x (um)")
    plt.ylabel("node y (um)")
    plt.title(title)
    plt.axis("equal")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close()


def save_pooled_summary_tables(all_distances_df, output_root):
    saved_paths = {}

    if len(all_distances_df) == 0:
        return saved_paths

    if "series_folder" in all_distances_df.columns:
        by_series = summarize_group_distances(all_distances_df, ["series_label"])
        by_series_path = output_root / "ALL_pooled_by_series_summary.csv"
        by_series.to_csv(by_series_path, index=False)
        saved_paths["by_series"] = str(by_series_path)

    if "mRNA" in all_distances_df.columns:
        by_mrna = summarize_group_distances(all_distances_df, ["mRNA"])
        by_mrna_path = output_root / "ALL_pooled_by_mRNA_summary.csv"
        by_mrna.to_csv(by_mrna_path, index=False)
        saved_paths["by_mRNA"] = str(by_mrna_path)

    if {"series_label", "mRNA"}.issubset(all_distances_df.columns):
        by_series_mrna = summarize_group_distances(all_distances_df, ["series_label", "mRNA"])
        by_series_mrna_path = output_root / "ALL_pooled_by_series_and_mRNA_summary.csv"
        by_series_mrna.to_csv(by_series_mrna_path, index=False)
        saved_paths["by_series_and_mRNA"] = str(by_series_mrna_path)

    return saved_paths


def make_overall_plots(all_distances_df, output_root, prefix="ALL"):
    if len(all_distances_df) == 0:
        return {}

    distances = all_distances_df["nearest_node_distance_um"].to_numpy(dtype=float)

    hist_path = output_root / f"{prefix}_mRNA_nearest_node_distance_histogram.png"
    make_histogram_plot(
        distances,
        f"{prefix}: mRNA distance to nearest node",
        hist_path,
    )

    cumulative_path = output_root / f"{prefix}_mRNAs_within_node_distance_thresholds.png"
    make_cumulative_threshold_plot(
        distances,
        f"{prefix}: counts within node-distance thresholds",
        cumulative_path,
    )

    return {
        "histogram_png": str(hist_path),
        "cumulative_threshold_png": str(cumulative_path),
    }


def analyze_coordinate_file(coord_file, node_file, series_output_root, node_cache, series_folder, series_label):
    mrna_name = infer_mrna_name(coord_file)
    cell_index = extract_cell_index(coord_file.name)
    cell_label = make_cell_label(cell_index)

    short_label = f"{sanitize_name(mrna_name)}_{cell_label}"
    out_folder = series_output_root / short_label
    out_folder.mkdir(parents=True, exist_ok=True)

    print(f"\n  Processing {mrna_name} {cell_label}")
    print(f"    mRNA coordinates: {coord_file}")
    print(f"    Matched node VTK: {node_file}")
    print(f"    Output folder: {out_folder}")

    mrna_xyz, mrna_meta = load_mrna_coordinates(coord_file)

    if node_file in node_cache:
        node_xyz, node_meta = node_cache[node_file]
    else:
        node_xyz, node_meta = load_node_positions(node_file)
        node_cache[node_file] = (node_xyz, node_meta)

    # Save node centers using a short name.
    node_centers_dir = series_output_root / "node_centers"
    node_centers_dir.mkdir(parents=True, exist_ok=True)
    node_centers_path = node_centers_dir / f"{cell_label}_node_centers.csv"

    if not node_centers_path.exists():
        pd.DataFrame(node_xyz, columns=["node_x", "node_y", "node_z"]).assign(
            node_index=np.arange(len(node_xyz))
        )[["node_index", "node_x", "node_y", "node_z"]].to_csv(node_centers_path, index=False)

    distances, nearest_node_idx, nearest_node_xyz = calculate_mrna_to_nearest_node(
        mrna_xyz,
        node_xyz,
    )

    per_mrna_df = pd.DataFrame(
        {
            "series_label": series_label,
            "series_folder": str(series_folder),
            "mRNA": mrna_name,
            "cell_index": cell_index,
            "mRNA_file": str(coord_file),
            "node_file": str(node_file),
            "mrna_index": np.arange(len(mrna_xyz)),
            "mrna_x": mrna_xyz[:, 0],
            "mrna_y": mrna_xyz[:, 1],
            "mrna_z": mrna_xyz[:, 2],
            "nearest_node_index": nearest_node_idx,
            "nearest_node_x": nearest_node_xyz[:, 0],
            "nearest_node_y": nearest_node_xyz[:, 1],
            "nearest_node_z": nearest_node_xyz[:, 2],
            "nearest_node_distance_um": distances,
        }
    )

    per_mrna_path = out_folder / f"{short_label}_distances.csv"
    per_mrna_df.to_csv(per_mrna_path, index=False)

    # Node occupancy: how many mRNAs have each node as their nearest node.
    node_occupancy = pd.DataFrame(
        {
            "series_label": series_label,
            "series_folder": str(series_folder),
            "mRNA": mrna_name,
            "cell_index": cell_index,
            "mRNA_file": str(coord_file),
            "node_file": str(node_file),
            "node_index": np.arange(len(node_xyz)),
            "node_x": node_xyz[:, 0],
            "node_y": node_xyz[:, 1],
            "node_z": node_xyz[:, 2],
        }
    )

    assigned_counts = pd.Series(nearest_node_idx).value_counts().sort_index()
    closest_distances = (
        pd.DataFrame({"node_index": nearest_node_idx, "distance": distances})
        .groupby("node_index")["distance"]
        .min()
    )

    node_occupancy["assigned_mrna_count"] = (
        node_occupancy["node_index"].map(assigned_counts).fillna(0).astype(int)
    )
    node_occupancy["closest_assigned_mrna_distance_um"] = (
        node_occupancy["node_index"].map(closest_distances)
    )

    node_occupancy_path = out_folder / f"{short_label}_node_occupancy.csv"
    node_occupancy.to_csv(node_occupancy_path, index=False)

    summary = summarize_distances(distances)
    summary.update(
        {
            "series_label": series_label,
            "series_folder": str(series_folder),
            "mRNA": mrna_name,
            "cell_index": cell_index,
            "mRNA_file": str(coord_file),
            "node_file": str(node_file),
            "output_folder": str(out_folder),
            "node_count": int(len(node_xyz)),
            "node_raw_vtk_point_count": int(node_meta.get("raw_point_count", 0)),
            "node_vtk_polygon_count": int(node_meta.get("polygon_count", 0)),
            "node_component_count": int(node_meta.get("component_count", 0)),
            "mRNA_input_rows": int(mrna_meta.get("input_rows", 0)),
            "mRNA_valid_xyz_rows": int(mrna_meta.get("valid_xyz_rows", 0)),
            "mRNA_x_column": mrna_meta.get("x_column", ""),
            "mRNA_y_column": mrna_meta.get("y_column", ""),
            "mRNA_z_column": mrna_meta.get("z_column", ""),
            "per_mrna_distances_csv": str(per_mrna_path),
            "node_occupancy_csv": str(node_occupancy_path),
            "node_centers_csv": str(node_centers_path),
        }
    )

    summary_path = out_folder / f"{short_label}_summary.csv"
    pd.DataFrame([summary]).to_csv(summary_path, index=False)

    hist_path = out_folder / f"{short_label}_histogram.png"
    make_histogram_plot(
        distances,
        f"{series_label} {mrna_name} {cell_label}: mRNA to nearest node",
        hist_path,
    )

    cumulative_path = out_folder / f"{short_label}_threshold_counts.png"
    make_cumulative_threshold_plot(
        distances,
        f"{series_label} {mrna_name} {cell_label}: mRNAs near nodes",
        cumulative_path,
    )

    xy_scatter_path = out_folder / f"{short_label}_xy_distance.png"
    make_xy_distance_scatter_plot(
        mrna_xyz,
        node_xyz,
        distances,
        f"{series_label} {mrna_name} {cell_label}: mRNAs colored by node distance",
        xy_scatter_path,
    )

    occupancy_plot_path = out_folder / f"{short_label}_node_occupancy.png"
    make_node_occupancy_plot(
        node_occupancy,
        f"{series_label} {mrna_name} {cell_label}: node occupancy",
        occupancy_plot_path,
    )

    summary.update(
        {
            "histogram_png": str(hist_path),
            "cumulative_threshold_png": str(cumulative_path),
            "xy_distance_scatter_png": str(xy_scatter_path),
            "node_occupancy_png": str(occupancy_plot_path),
            "summary_csv": str(summary_path),
        }
    )

    print(f"    mRNA points analyzed: {len(mrna_xyz)}")
    print(f"    Node centers analyzed: {len(node_xyz)}")
    print(f"    Median distance: {summary['median_distance_um']:.4f} um")
    print(f"    Mean distance: {summary['mean_distance_um']:.4f} um")

    return per_mrna_df, summary, node_occupancy


# ==========================================================
# Main
# ==========================================================
def main():
    if not parent_folder.exists():
        raise FileNotFoundError(f"Parent folder does not exist: {parent_folder}")

    output_root = parent_folder / OUTPUT_FOLDER_NAME
    output_root.mkdir(parents=True, exist_ok=True)

    print(f"Parent folder: {parent_folder}")
    print(f"Output folder: {output_root}")
    print(f"Node position mode: {NODE_POSITION_MODE}")
    print(f"Coordinate input mode: {'converted_coordinates only' if PROCESS_ONLY_CONVERTED_COORDINATES else 'all x/y/z tables'}")
    print("\nAnalysis mode: Series-scoped recursion")
    print("  Node and mRNA files are matched only within the same Series folder.")
    print("  Results are then pooled across all Series folders.")

    series_folders = find_series_folders(parent_folder)

    print(f"\nFound {len(series_folders)} Series analysis folder(s):")
    for series in series_folders:
        print(f"  SERIES: {series}")

    all_distances = []
    all_summaries = []
    all_node_occupancy = []
    unmatched = []
    rejected_all = []

    for series_folder in series_folders:
        series_label = make_series_output_label(series_folder, parent_folder)
        series_output_root = output_root / "per_series" / series_label
        series_output_root.mkdir(parents=True, exist_ok=True)

        print("\n" + "=" * 80)
        print(f"Processing Series folder: {series_folder}")
        print(f"Series output folder: {series_output_root}")

        node_files = find_node_vtk_files(series_folder)
        coordinate_files, rejected_coordinate_files = find_coordinate_files(series_folder)

        print(f"  Found {len(node_files)} node VTK file(s) in this Series.")
        print(f"  Found {len(coordinate_files)} mRNA coordinate file(s) in this Series.")

        if len(node_files) > 0:
            node_debug = [
                {
                    "node_file": str(p),
                    "cell_index": extract_cell_index(p.name),
                }
                for p in node_files
            ]
            pd.DataFrame(node_debug).to_csv(
                series_output_root / f"{series_label}_node_file_index.csv",
                index=False,
            )

        if len(coordinate_files) > 0:
            coord_debug = [
                {
                    "mRNA_file": str(p),
                    "mRNA": infer_mrna_name(p),
                    "cell_index": extract_cell_index(p.name),
                }
                for p in coordinate_files
            ]
            pd.DataFrame(coord_debug).to_csv(
                series_output_root / f"{series_label}_coordinate_file_index.csv",
                index=False,
            )

        for path, reason in rejected_coordinate_files:
            rejected_all.append(
                {
                    "series_label": series_label,
                    "series_folder": str(series_folder),
                    "file": str(path),
                    "reason": reason,
                }
            )

        if len(node_files) == 0:
            print("  No *_nodes.vtk files found in this Series. Skipping.")
            unmatched.append(
                {
                    "series_label": series_label,
                    "series_folder": str(series_folder),
                    "mRNA_file": "",
                    "reason": "no_node_vtk_files_in_series",
                }
            )
            continue

        if len(coordinate_files) == 0:
            print("  No mRNA coordinate files found in this Series. Skipping.")
            unmatched.append(
                {
                    "series_label": series_label,
                    "series_folder": str(series_folder),
                    "mRNA_file": "",
                    "reason": "no_coordinate_tables_in_series",
                }
            )
            continue

        node_entries = build_node_index(node_files)
        node_cache = {}

        series_distances = []
        series_summaries = []
        series_node_occupancy = []

        for coord_file in coordinate_files:
            node_file = find_matching_node_file(coord_file, node_entries)

            if node_file is None:
                print(f"\n  No matching node VTK found for: {coord_file}")
                unmatched.append(
                    {
                        "series_label": series_label,
                        "series_folder": str(series_folder),
                        "mRNA_file": str(coord_file),
                        "cell_index": extract_cell_index(coord_file.name),
                        "reason": "no_matching_node_vtk_with_same_cell_index",
                    }
                )
                continue

            try:
                per_mrna_df, summary, node_occupancy = analyze_coordinate_file(
                    coord_file=coord_file,
                    node_file=node_file,
                    series_output_root=series_output_root,
                    node_cache=node_cache,
                    series_folder=series_folder,
                    series_label=series_label,
                )

                series_distances.append(per_mrna_df)
                series_summaries.append(summary)
                series_node_occupancy.append(node_occupancy)

                all_distances.append(per_mrna_df)
                all_summaries.append(summary)
                all_node_occupancy.append(node_occupancy)

            except Exception as e:
                print(f"\n  Error processing {coord_file}")
                print(f"    Reason: {e}")
                unmatched.append(
                    {
                        "series_label": series_label,
                        "series_folder": str(series_folder),
                        "mRNA_file": str(coord_file),
                        "cell_index": extract_cell_index(coord_file.name),
                        "reason": str(e),
                    }
                )

        # Save pooled outputs for this Series only.
        if len(series_distances) > 0:
            series_distances_df = pd.concat(series_distances, ignore_index=True)
            series_distances_path = series_output_root / f"{series_label}_ALL_distances.csv"
            series_distances_df.to_csv(series_distances_path, index=False)

            series_summary_df = pd.DataFrame(series_summaries)
            series_summary_path = series_output_root / f"{series_label}_summary_by_file.csv"
            series_summary_df.to_csv(series_summary_path, index=False)

            series_node_occupancy_df = pd.concat(series_node_occupancy, ignore_index=True)
            series_node_occupancy_path = series_output_root / f"{series_label}_node_occupancy.csv"
            series_node_occupancy_df.to_csv(series_node_occupancy_path, index=False)

            make_overall_plots(series_distances_df, series_output_root, prefix=f"{series_label}_ALL")
            save_pooled_summary_tables(series_distances_df, series_output_root)

            print(f"\n  Series pooled distances: {series_distances_path}")
            print(f"  Series summary: {series_summary_path}")
            print(f"  Series node occupancy: {series_node_occupancy_path}")

    if len(rejected_all) > 0:
        rejected_path = output_root / "ALL_rejected_coordinate_tables.csv"
        pd.DataFrame(rejected_all).to_csv(rejected_path, index=False)
        print(f"\nRejected coordinate-table candidates saved to: {rejected_path}")

    if len(all_distances) > 0:
        all_distances_df = pd.concat(all_distances, ignore_index=True)
        all_distances_path = output_root / "ALL_SERIES_mRNA_to_nearest_node_distances.csv"
        all_distances_df.to_csv(all_distances_path, index=False)

        summary_df = pd.DataFrame(all_summaries)
        summary_path = output_root / "ALL_SERIES_summary_by_file.csv"
        summary_df.to_csv(summary_path, index=False)

        all_node_occupancy_df = pd.concat(all_node_occupancy, ignore_index=True)
        all_node_occupancy_path = output_root / "ALL_SERIES_node_occupancy.csv"
        all_node_occupancy_df.to_csv(all_node_occupancy_path, index=False)

        overall_plot_paths = make_overall_plots(all_distances_df, output_root, prefix="ALL_SERIES")
        pooled_summary_paths = save_pooled_summary_tables(all_distances_df, output_root)

        print("\n" + "=" * 80)
        print("Done.")
        print(f"All-Series pooled per-mRNA nearest-node distances: {all_distances_path}")
        print(f"All-Series per-file summary: {summary_path}")
        print(f"All-Series node occupancy: {all_node_occupancy_path}")

        if "by_series" in pooled_summary_paths:
            print(f"Pooled by Series summary: {pooled_summary_paths['by_series']}")
        if "by_mRNA" in pooled_summary_paths:
            print(f"Pooled by mRNA summary: {pooled_summary_paths['by_mRNA']}")
        if "by_series_and_mRNA" in pooled_summary_paths:
            print(f"Pooled by Series and mRNA summary: {pooled_summary_paths['by_series_and_mRNA']}")

        if "histogram_png" in overall_plot_paths:
            print(f"All-Series histogram: {overall_plot_paths['histogram_png']}")
        if "cumulative_threshold_png" in overall_plot_paths:
            print(f"All-Series cumulative plot: {overall_plot_paths['cumulative_threshold_png']}")

    else:
        print("\nNo mRNA coordinate files were successfully processed across any Series.")

    if len(unmatched) > 0:
        unmatched_path = output_root / "ALL_unmatched_or_failed_files.csv"
        pd.DataFrame(unmatched).to_csv(unmatched_path, index=False)
        print(f"Unmatched/failed files saved to: {unmatched_path}")


if __name__ == "__main__":
    main()
