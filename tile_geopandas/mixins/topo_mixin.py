from __future__ import annotations

import geopandas as gpd
import shapely as shp
from shapely.geometry import Point, Polygon, LineString
import pandas as pd
from itertools import combinations
import warnings


if __name__ == "__main__":
    
    from datetime import datetime
    def tlog(msg):
        print(f"{datetime.now()}: {msg}")
    
    file_path = "example_data/bnd_oa_00_2024_2Q.shp"
    save_path = "example_data/bnd_oa_00_2024_2Q_topo.gpkg"
    grid_size = 0.01 #must be set above... for topology error handlng
    simpl_tol = 10
    
    gdf = gpd.read_file(file_path)
    gdf = gdf[gdf["ADM_CD"].str.startswith("23")]
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
        .geometry
        .explode()
        .set_precision(grid_size, mode="valid_output")
        .reset_index(name="geometry")
    )
    geom_df.to_file("example_data/geom_df.gpkg", index=False)
    # =========================================================================
    tlog("boundary extraction")
    lines = []
    for idx, row in geom_df.iterrows():
        geom = row.geometry
        exterior = geom.exterior
        interiors = geom.interiors
        lines.append(exterior)
        lines.extend(interiors)
        # divide polygon to remove 'hole'
        for interior in interiors:
            is_split = False
            vertices = [Point(c) for c in interior.coords][1:]
            for v1, v2 in combinations(vertices, 2):
                l1 = shp.shortest_line(v1, exterior)
                l2 = shp.shortest_line(v2, exterior)
                conds = [l1.length > 0, 
                         l2.length > 0,
                         shp.intersects(l1, exterior),
                         shp.intersects(l2, exterior),
                         not shp.intersects(l1, l2)]
                if all(conds):
                    lines.extend([l1, l2])
                    is_split = True
                    break
            if not is_split:
                raise ValueError("hole not splitted")
    # =========================================================================
    tlog("to graph")
    lines = shp.disjoint_subset_union_all(lines)
    lines = shp.line_merge(lines)
    lines = shp.node(lines)
    lines = [l for l in lines.geoms]
    edge_df = (
        gpd.GeoDataFrame(geometry=lines, crs=geom_df.crs)
        .reset_index(drop=True)
        .reset_index(names="edge_id")
    )
    edge_df.to_file("example_data/edge_df.gpkg", index=False)
    # =========================================================================
    tlog("polygonize")
    faces = shp.polygonize(lines)
    faces = [f for f in faces.geoms if isinstance(f, Polygon)]
    face_df = (
        gpd.GeoSeries(faces, crs=geom_df.crs)
        .rename_axis("face_id")
        .reset_index(name="geometry")
    )
    face_df.to_file("example_data/face_df.gpkg", index=False)
    # =========================================================================
    tlog("build edge-face topology")
    edge_face_topo = (
        edge_df
        .overlay(face_df, how='intersection', keep_geom_type=True)
        .assign(length=lambda df: df["geometry"].length)
        .loc[:, ["edge_id", "face_id", "length"]]
        .sort_values(by=["edge_id"])
        .reset_index(drop=True)
    )
    # =========================================================================
    tlog("setting precision...")
    edge_df.geometry = edge_df.geometry.set_precision(grid_size)
    edge_df.to_file("example_data/edge_df.gpkg", index=False)
    # =========================================================================
    tlog("face matching")
    face_match = (
        face_df
        .geometry
        .representative_point()
        .rename_axis("face_id")
        .reset_index(name="geometry")
        .sjoin(geom_df, how="left", predicate="within")
        .drop_duplicates(subset="face_id")
        .loc[:, ["face_id", "geom_id"]]
    )
    face_df = face_df.merge(face_match, how="left", on="face_id")
    face_df.to_file("example_data/face_df.gpkg", index=False)
    # =========================================================================
    tlog("gap matching")
    gap_face_idx = (
        face_df
        .loc[face_df["geom_id"].isna(), "face_id"]
        .to_list()
    )
    for face_id in gap_face_idx:
        conn = (
            edge_face_topo
            .loc[edge_face_topo["face_id"]==face_id, ["edge_id"]]
            .merge(edge_face_topo, how="left", on="edge_id")
            .query("face_id != @face_id")
        )
        alloc_face_id = conn.at[conn["length"].idxmax(), "face_id"]
        alloc_geom_id = (
            face_df
            .loc[face_df["face_id"]==alloc_face_id, "geom_id"]
            .iloc[0]
        )
        face_df.loc[face_id, "geom_id"] = alloc_geom_id
    face_df.to_file("example_data/face_df.gpkg", index=False)
    face_df.drop(columns="geometry", inplace=True)
    # =========================================================================
    tlog("generate new polygon")
    new_geom = (
        edge_df
        .merge(edge_face_topo, how="inner", on="edge_id")
        .dissolve(by="face_id", method="disjoint_subset", as_index=False)
        .assign(geometry=lambda df: shp.line_merge(df.geometry))
        .pipe(lambda df: df[df.geometry.is_ring])
        .assign(geometry=lambda df: df.geometry.apply(Polygon))
        .merge(face_match, how="inner", on="face_id")
        .query("geom_id.notna()")
        .astype({"geom_id": int})
        .dissolve(by="geom_id", method="coverage", as_index=False)
    )
    new_geom.to_file("example_data/new_geom.gpkg", index=False)
    # for i, row in new_geom[new_geom["face_id"] == 3676].explode().iterrows():
    #     coords = row["geometry"].coords
    #     print(i, coords[0], coords[-1])
    # print(new_geom.shape)
    # print(face_df.shape)
    # raise Exception("STOP")
    # =========================================================================
    tlog("simplifying")
    new_geom.geometry = new_geom.geometry.simplify_coverage(simpl_tol)
    new_geom.to_file("example_data/new_geom.gpkg", index=False)
    # =========================================================================
    tlog("result")
    gdf = (
        new_geom
        .merge(attr_df, how="right", on="geom_id")
        .loc[:, list(attr_df.columns) + ["geometry"]]
        .drop(columns="geom_id")
        .reset_index(drop=True)
    )
    gdf.to_file("example_data/result.gpkg", index=False)
    

    # 이후 simplify
    # 마지막으로 dissolve
    tlog("face table")