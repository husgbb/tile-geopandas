from typing import Any
from typeguard import typechecked
import numpy as np
import geopandas as gpd
import shapely as shp
import pandas as pd
from datetime import datetime


@typechecked
def tlog(msg: Any, verbose: bool):
    dt_str = datetime.now().isoformat(timespec="milliseconds")
    if verbose:
        print(f"{dt_str}: {msg}")

@typechecked
def fix_polygon_topology(gdf: gpd.GeoDataFrame, 
                         grid_size: float, 
                         gap_tol: float, 
                         simplify_tol: float, 
                         verbose: bool = False) -> gpd.GeoDataFrame:
    """Fixes the topology of a polygon GeoDataFrame.

    This function cleans polygon geometries by removing overlaps, filling gaps,
    and simplifying the boundary. It operates by first snapping all vertices
    to a grid of a given size, then rebuilding the polygons from their
    boundaries to enforce a 'planar' topology. Small gaps between the resulting
    polygons are identified and filled by assigning them to the adjacent
    polygon with the longest shared boundary. Finally, the cleaned topology
    is simplified.

    Args:
        gdf (gpd.GeoDataFrame): The input GeoDataFrame with Polygon or
            MultiPolygon geometries.
        grid_size (float): The grid size to snap vertices to. All vertices
            will be rounded to a multiple of this value. Must be > 0.
        gap_tol (float): The tolerance for identifying gaps. Gaps smaller
            than this tolerance (by negative buffering) will be ignored.
            Must be > 0.
        simplify_tol (float): The tolerance for simplifying the final
            geometries. This is a topology-preserving simplification.
            Must be > 0.
        verbose (bool, optional): If True, prints progress messages.
            Defaults to False.

    Raises:
        ValueError: If `grid_size`, `gap_tol`, or `simplify_tol` are not
            positive, or if the input GeoDataFrame contains non-polygonal
            geometries.

    Returns:
        gpd.GeoDataFrame: A new GeoDataFrame with the fixed polygon topology.
    """
    
    tlog("input check", verbose)
    if grid_size <= 0:
        raise ValueError("grid_size should be greater than 0")
    if gap_tol <= 0:
        raise ValueError("gap_tol should be greater than 0")
    if simplify_tol <= 0:
        raise ValueError("simplify_tol should be greater than 0")    
    geom_types = set(gdf.geometry.geom_type.unique().tolist())
    if geom_types - {"Polygon", "MultiPolygon"}:
        raise ValueError(f"input geometry should be polygonal geometry >>> {geom_types}")

    tlog("separate attribute", verbose)
    gdf["geom_id"] = range(len(gdf))
    attr_df = gdf.drop(columns="geometry")
    geom_df = gdf.loc[:, ["geom_id", "geometry"]]
    del gdf
    
    tlog("precision control", verbose)
    geom_df = (
        geom_df
        .set_index("geom_id")
        .explode()
        .set_precision(grid_size, mode="valid_output")
        .reset_index(name="geometry")
    )
    
    tlog("regenerate polygon", verbose)
    face_df = (
        geom_df
        .set_index("geom_id")
        .geometry.boundary
        .polygonize(node=True, full=False)
        .rename_axis("face_id")
        .reset_index(name="geometry")
    )
    rpoint_df = (
        face_df
        .set_index("face_id")
        .geometry.representative_point()
        .reset_index(name="geometry")
    )
    match_df = (
        rpoint_df
        .sjoin(geom_df, how="left", predicate="within")
        .drop_duplicates(subset="face_id")
        .loc[:, ["face_id", "geom_id"]]
    )
    face_df = (
        face_df
        .merge(match_df, how="left", on="face_id")
        .drop_duplicates(subset="face_id")
    )
    
    tlog("precision control", verbose)
    face_df.geometry = face_df.geometry.set_precision(grid_size)
    face_df = face_df[~face_df.geometry.is_empty]
    face_df["face_id"] = range(len(face_df))
    
    tlog("gap filling", verbose)
    gap_df = (
        face_df.loc[face_df["geom_id"].isna(), ["geometry"]]
        .pipe(lambda df: df[~df.geometry.buffer(-gap_tol).is_empty])
        .reset_index(drop=True).reset_index(names="gap_id")
    )
    face_df = (
        face_df.loc[face_df["geom_id"].notna(), ["geom_id", "geometry"]]
        .astype({"geom_id": int})
        .reset_index(drop=True).reset_index(names="face_id")
    )
    gap_topo_df = (
        gap_df
        .overlay(right=face_df, 
                 how='intersection', 
                 keep_geom_type=False, 
                 make_valid=False)
        .assign(touch_length=lambda df: df["geometry"].length)
        .pipe(lambda df: df[df["touch_length"] > 0])
    )
    
    tlog("setting precision...", verbose)
    # if multilinestring, its hole 'touches' the exterior
    gap_df["geom_id"] = -1
    for gap_id, gap_topo in gap_topo_df.groupby("gap_id"):
        idx = gap_topo["touch_length"].idxmax()
        alloc_geom_id = gap_topo.loc[idx, "geom_id"]
        gap_df.loc[gap_df["gap_id"]==gap_id, "geom_id"] = alloc_geom_id
    
    tlog("restore geometry...", verbose)
    face_df = pd.concat([
        face_df[["geom_id", "geometry"]], 
        gap_df[["geom_id", "geometry"]]
    ])
    geom_df = face_df.dissolve(by="geom_id", method="coverage", as_index=False)
    del gap_df, gap_topo_df, face_df
    
    tlog("simplifying...", verbose)
    stol_sqrt = np.sqrt(simpl_tol)
    geom_df.geometry = geom_df.geometry.simplify_coverage(stol_sqrt)
    
    tlog("result", verbose)
    gdf = geom_df.merge(attr_df, how="inner", on="geom_id")
    gdf.to_file("example_data/result.gpkg", index=False)
    
    tlog("done", verbose)
    
    return gdf

if __name__ == "__main__":
    file_path = "example_data/bnd_oa_00_2024_2Q.shp"
    save_path = "example_data/bnd_oa_00_2024_2Q_topo.gpkg"
    
    gdf = gpd.read_file(file_path)
    grid_size = 0.01 # must be set above... for topology error handlng
    gap_tol = 0.1
    simpl_tol = 10
    params = dict(
        gdf=gdf,
        grid_size=grid_size,
        gap_tol=gap_tol,
        simplify_tol=simpl_tol,
        verbose=True
    )
    gdf = fix_polygon_topology(**params)
    gdf.to_file(save_path, index=False)

