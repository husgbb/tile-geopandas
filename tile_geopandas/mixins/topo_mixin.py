from __future__ import annotations

import numpy as np
import geopandas as gpd
import shapely as shp
from shapely.geometry import Point, Polygon, LineString, MultiLineString
import pandas as pd
from itertools import combinations
import warnings


if __name__ == "__main__":
    
    from datetime import datetime
    def tlog(msg):
        print(f"{datetime.now()}: {msg}")
    
    file_path = "example_data/bnd_oa_00_2024_2Q.shp"
    save_path = "example_data/bnd_oa_00_2024_2Q_topo.gpkg"
    grid_size = 0.01 # must be set above... for topology error handlng
    gap_tol = 0.1
    simpl_tol = 10
    
    gdf = gpd.read_file(file_path)
    # gdf = gdf[gdf["ADM_CD"].str.startswith("1")]
    gdf.to_file("example_data/raw.gpkg", index=False)

    # =========================================================================
    tlog("separate attribute")
    gdf["geom_id"] = range(len(gdf))
    attr_df = gdf.drop(columns="geometry")
    geom_df = gdf.loc[:, ["geom_id", "geometry"]]
    del gdf
    geom_df.to_file("example_data/geom_df.gpkg", index=False)
    # =========================================================================
    tlog("precision control")
    geom_df = (
        geom_df
        .set_index("geom_id")
        .explode()
        .set_precision(grid_size, mode="valid_output")
        .reset_index(name="geometry")
    )
    geom_df.to_file("example_data/geom_df.gpkg", index=False)
    # =========================================================================
    tlog("regenerate polygon")
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
    rpoint_df.to_file("example_data/rpoint_df.gpkg", index=False)
    face_df.to_file("example_data/face_df.gpkg", index=False)
    del geom_df, rpoint_df, match_df
    # =========================================================================
    tlog("precision control")
    face_df.geometry = face_df.geometry.set_precision(grid_size)
    face_df = face_df[~face_df.geometry.is_empty]
    face_df["face_id"] = range(len(face_df))
    face_df.to_file("example_data/face_df.gpkg", index=False)
    # =========================================================================
    tlog("gap filling")
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
    gap_df.to_file("example_data/gap_df.gpkg", index=False)
    face_df.to_file("example_data/face_df.gpkg", index=False)
    # =========================================================================
    tlog("setting precision...")
    merge_lines = lambda ls: shp.line_merge(shp.disjoint_subset_union_all(ls)) 
    # if multilinestring, its hole 'touches' the exterior
    gap_df["geom_id"] = -1
    for gap_id, gap_topo in gap_topo_df.groupby("gap_id"):
        idx = gap_topo["touch_length"].idxmax()
        alloc_geom_id = gap_topo.loc[idx, "geom_id"]
        gap_df.loc[gap_df["gap_id"]==gap_id, "geom_id"] = alloc_geom_id
    # =========================================================================
    tlog("restore geometry...")
    face_df = pd.concat([
        face_df[["geom_id", "geometry"]], 
        gap_df[["geom_id", "geometry"]]
    ])
    geom_df = face_df.dissolve(by="geom_id", method="coverage", as_index=False)
    del gap_df, gap_topo_df, face_df
    # =========================================================================
    tlog("simplifying...")
    stol_sqrt = np.sqrt(simpl_tol)
    geom_df.geometry = geom_df.geometry.simplify_coverage(stol_sqrt)
    # =========================================================================
    tlog("result")
    gdf = geom_df.merge(attr_df, how="inner", on="geom_id")
    gdf.to_file("example_data/result.gpkg", index=False)
    # =========================================================================
    tlog("face table")
