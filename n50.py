"""Read and preview Kartverket N50 File Geodatabase vector data."""
from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pyogrio
from rasterio.features import rasterize
from shapely import constrained_delaunay_triangles
from shapely.geometry import MultiPolygon, Polygon, box
from shapely.ops import unary_union
import trimesh


WATER_AREA_PREFIXES = ("Innsj", "Havflat")


def list_layers(gdb_path):
    """Return the feature layers and their geometry types."""
    rows = pyogrio.list_layers(Path(gdb_path))
    return pd.DataFrame(rows, columns=["layer", "geometry_type"])


def read_features(gdb_path, bbox):
    """Read every N50 feature layer intersecting bbox=(xmin, ymin, xmax, ymax)."""
    result = {}
    for layer, _geometry_type in pyogrio.list_layers(Path(gdb_path)):
        # N50 contains year-1000 placeholder dates, outside pandas' datetime range.
        features = pyogrio.read_dataframe(
            gdb_path, layer=layer, bbox=bbox, datetime_as_string=True
        )
        if not features.empty:
            result[layer] = features
    if not result:
        raise ValueError("No N50 features intersect the requested bbox.")
    crs = {str(features.crs) for features in result.values()}
    if len(crs) != 1:
        raise ValueError(f"N50 layers use mixed CRSs: {crs}")
    return result


def water_areas(features):
    """Return N50 lake and sea polygons from the area-cover layer.

    N50 stores these as ``Innsjø`` and ``Havflate`` in the ``objtype`` field.
    Prefix matching also handles driver-level text-decoding differences.
    """
    candidates = [name for name in features if name.startswith("N50_Arealdekke_omrade")]
    if not candidates:
        raise ValueError("N50 area-cover polygon layer is missing.")
    areas = features[candidates[0]]
    if "objtype" not in areas:
        raise ValueError("N50 area-cover layer has no objtype field.")
    objtype = areas["objtype"].fillna("").astype(str)
    is_water = objtype.str.startswith(WATER_AREA_PREFIXES)
    return areas.loc[is_water].copy()


def water_summary(water):
    """Summarize water features without double-counting their clipped geometry."""
    return (water.assign(area_m2=water.geometry.area)
            .groupby("objtype", dropna=False)
            .agg(features=("objtype", "size"), mapped_area_m2=("area_m2", "sum"))
            .reset_index()
            .sort_values("objtype"))


def water_mask(water, dem):
    """Rasterize water polygons on precisely the DEM grid."""
    if str(water.crs) != str(dem["crs"]):
        raise ValueError(f"CRS mismatch: N50={water.crs}, DEM={dem['crs']}")
    return rasterize(
        ((geometry, 1) for geometry in water.geometry if geometry and not geometry.is_empty),
        out_shape=dem["elevation_m"].shape,
        transform=dem["transform"],
        fill=0,
        dtype="uint8",
        all_touched=False,
    ).astype(bool)


def plot_features(features, bbox, ax, title="N50 features"):
    """Draw an overview map, preserving the requested projected bounding box."""
    for layer, frame in features.items():
        kind = frame.geom_type.iloc[0] if not frame.empty else ""
        if "Polygon" in kind:
            frame.plot(ax=ax, facecolor="#d9d9d9", edgecolor="#666666", linewidth=0.25)
        elif "Line" in kind:
            frame.plot(ax=ax, color="#555555", linewidth=0.45)
        else:
            frame.plot(ax=ax, color="#222222", markersize=2)
    xmin, ymin, xmax, ymax = bbox
    ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), title=title,
           xlabel="UTM east (m)", ylabel="UTM north (m)")
    ax.set_aspect("equal")


def plot_water(water, bbox, ax, title="N50 lake and sea polygons", alpha=0.70):
    """Draw only the N50 sea and lake polygons."""
    water.plot(ax=ax, facecolor="#2b83ba", edgecolor="#0b3c5d", linewidth=0.5, alpha=alpha)
    xmin, ymin, xmax, ymax = bbox
    ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), title=title,
           xlabel="UTM east (m)", ylabel="UTM north (m)")
    ax.set_aspect("equal")


def water_overlay_mesh(water, dem, model_width_mm, vertical_exaggeration=1.0,
                       base_mm=3.0, thickness_mm=0.4, raise_mm=0.0):
    """Create a thin, watertight water part with N50-defined boundaries.

    The N50 polygons are triangulated as vectors, so their coastline is not
    snapped to DTM cells. Heights are sampled bilinearly from the DTM solely to
    drape this printable overlay on the terrain surface.
    """
    if thickness_mm <= 0 or model_width_mm <= 0 or vertical_exaggeration <= 0:
        raise ValueError("model width, vertical exaggeration, and thickness must be positive")
    if str(water.crs) != str(dem["crs"]):
        raise ValueError(f"CRS mismatch: N50={water.crs}, DEM={dem['crs']}")

    height = np.asarray(dem["elevation_m"], dtype="float64")
    transform = dem["transform"]
    resolution = float(dem["output_resolution_m"])
    rows, columns = height.shape
    xmin = float(transform.c)
    ymax = float(transform.f)
    xmax = xmin + columns * resolution
    ymin = ymax - rows * resolution
    crop = box(xmin, ymin, xmax, ymax)
    scale = model_width_mm / (columns * resolution)

    def terrain_height(x, y):
        # Coordinates refer to cell centres; bilinear sampling avoids a stair-step Z surface.
        col = np.clip((x - xmin) / resolution - 0.5, 0, columns - 1)
        row = np.clip((ymax - y) / resolution - 0.5, 0, rows - 1)
        c0, r0 = int(np.floor(col)), int(np.floor(row))
        c1, r1 = min(c0 + 1, columns - 1), min(r0 + 1, rows - 1)
        dc, dr = col - c0, row - r0
        return ((1 - dc) * (1 - dr) * height[r0, c0] + dc * (1 - dr) * height[r0, c1]
                + (1 - dc) * dr * height[r1, c0] + dc * dr * height[r1, c1])

    xy_to_index = {}
    xy = []

    def add_vertex(x, y):
        key = (round(float(x), 8), round(float(y), 8))
        if key not in xy_to_index:
            xy_to_index[key] = len(xy)
            xy.append((float(x), float(y)))
        return xy_to_index[key]

    # N50 commonly splits a continuous sea/lake surface into adjacent features.
    # Dissolving first removes their internal shared borders; otherwise each one
    # becomes a duplicate vertical wall and the combined mesh is non-manifold.
    merged_water = unary_union([geometry for geometry in water.geometry
                                if geometry is not None and not geometry.is_empty])
    if isinstance(merged_water, Polygon):
        source_polygons = (merged_water,)
    elif isinstance(merged_water, MultiPolygon):
        source_polygons = merged_water.geoms
    else:
        source_polygons = tuple(geometry for geometry in merged_water.geoms
                                if isinstance(geometry, Polygon))

    top_faces = []
    rings = []
    for geometry in source_polygons:
        clipped = geometry.intersection(crop)
        parts = clipped.geoms if isinstance(clipped, MultiPolygon) else (clipped,)
        for polygon in parts:
            if not isinstance(polygon, Polygon) or polygon.is_empty or polygon.area == 0:
                continue
            for triangle in constrained_delaunay_triangles(polygon).geoms:
                coords = list(triangle.exterior.coords)[:-1]
                if len(coords) == 3:
                    top_faces.append([add_vertex(x, y) for x, y in coords])
            rings.append(list(polygon.exterior.coords)[:-1])
            rings.extend(list(interior.coords)[:-1] for interior in polygon.interiors)
    if not top_faces:
        raise ValueError("No N50 water polygons remain inside the DEM crop.")
    # Preserve every original polygon-ring vertex, including colinear vertices
    # that triangulation may not need but which define the visible shoreline.
    for ring in rings:
        for x, y in ring:
            add_vertex(x, y)

    top = np.array([
        ((x - xmin) * scale, (y - ymin) * scale,
         base_mm + (terrain_height(x, y) - height.min()) * scale * vertical_exaggeration + raise_mm)
        for x, y in xy
    ])
    bottom = top.copy()
    bottom[:, 2] -= thickness_mm
    vertex_count = len(top)
    faces = [*top_faces, *[[a + vertex_count, c + vertex_count, b + vertex_count]
                           for a, b, c in top_faces]]
    for ring in rings:
        indices = [add_vertex(x, y) for x, y in ring]
        for current, following in zip(indices, indices[1:] + indices[:1]):
            faces.extend([[current, following, following + vertex_count],
                          [current, following + vertex_count, current + vertex_count]])
    mesh = trimesh.Trimesh(vertices=np.vstack([top, bottom]), faces=np.asarray(faces), process=False)
    trimesh.repair.fix_normals(mesh)
    if not mesh.is_watertight:
        raise RuntimeError("Generated water overlay is not watertight")
    return mesh, {
        "water_features": len(water), "vertices": len(mesh.vertices), "triangles": len(mesh.faces),
        "watertight": bool(mesh.is_watertight), "thickness_mm": thickness_mm,
        "raise_mm": raise_mm,
    }
