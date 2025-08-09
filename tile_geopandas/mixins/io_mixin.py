"""
Mixin for io operations
"""
from __future__ import annotations

import geopandas as gpd
import pandas as pd
from typing import Callable
from geopandas import GeoDataFrame
import warnings

class IOMixin(GeoDataFrame):
    tiles: GeoDataFrame | None = None
    process: Callable
    self.process = lambda x: x

    def tile(self, tile_size: float, grid_size: float = 0) -> TileMixin:
        xmin, ymin, xmax, ymax = self.total_bounds
        # generate tile bounds
        x = np.arange(xmin - grid_size, xmax + grid_size, tile_size)
        y = np.arange(ymin - grid_size, ymax + grid_size, tile_size)
        # round to grid size
        if grid_size > 0:
            x = np.round(x / grid_size) * grid_size
            y = np.round(y / grid_size) * grid_size
        x = np.unique(x)
        y = np.unique(y)
        if len(x) == 1 and len(y) == 1:
            warnings.warn("only one tile generated. Check tile size and grid size.")
        # generate tiles
        print(x)
        print(y)
        # tiles: list[GeoDataFrame] = []
        # for x0 in x:
        #     for y0 in y:
        #         x1 = x0 + tile_size
        #         y1 = y0 + tile_size
        #         bbox = box(x0, y0, x1, y1)
        #         # Select features within the tile
        #         tile = self[self.geometry.within(bbox)].copy()
        #         if not tile.empty:
        #             tiles.append(tile)
        # self._tiles = tiles  # type: ignore
        return self

    def apply_to_tiles(self, func: Callable[[GeoDataFrame], GeoDataFrame]) -> TileMixin:
        if hasattr(self, "_tiles"):
            self._tiles = [func(tile) for tile in self._tiles]  # type: ignore
        return self

    def merge(self) -> GeoDataFrame:
        if hasattr(self, "_tiles") and self._tiles:
            merged = pd.concat(self._tiles, ignore_index=True)  # type: ignore
            result = gpd.GeoDataFrame(merged)
            # Clean up tiles
            delattr(self, "_tiles")
            return result
        # If no tiles present, return self casted as GeoDataFrame
        return gpd.GeoDataFrame(self)


if __name__ == "__main__":
    print("read")
    gdf = gpd.read_file("example_data/bnd_oa_00_2024_2Q.shp")
    tgdf = TileMixin(gdf)
    print("end")
