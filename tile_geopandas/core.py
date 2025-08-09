from __future__ import annotations

import geopandas as gpd
from typing import Optional, Callable
from geopandas import GeoDataFrame
from numpy import uint64
from typeguard import typechecked
from uuid import uuid4

class TileDataFrame():

    @typechecked
    def __init__(self, gdf: GeoDataFrame, id: Optional[str] = None):
        self.id: Optional[str] = id
        self.parallel: uint64 = 1
        self.tiles: Optional[GeoDataFrame] = None
        self.process: Callable[TileDataFrame] = lambda x: x
    
    @property
    def _constructor(self):
        return TileDataFrame

    @typechecked
    def set_id_column(self):
        if self.id in self.columns:
            is_unique = len(self) == len(self[self.id].unique())
            if not is_unique:
                raise ValueError("ID column is not unique")
            return self
        elif self.id is None:
            self.id = "id"
            while self.id not in self.columns:
                self.id = str(uuid4())
        elif self.id not in self.
        
        el
            
        if self.id not in self.columns:
            self["id"] = range(len(self))
    
    @typechecked
    def set_parallel(self, parallel: int) -> TileDataFrame:
        self.parallel = parallel
        return self
    
    @typechecked
    def add_process(self, func: Callable[GeoDataFrame]) -> GeoDataFrame:
        self.process = lambda func: process(func) 
        return self
    
    @typechecked
    def process(self) -> GeoDataFrame:
        return self.process(self)
    
    class TileMixin(gpd.GeoDataFrame):
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
    print("=" * 100)
    gdf = gpd.read_file("example_data/bnd_oa_00_2024_2Q.shp")
    tgdf = TileDataFrame(id="id", data=gdf)

    print("=" * 100)
    print(type(tgdf))
    print(tgdf.head())
    print("end")
