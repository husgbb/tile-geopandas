# tile-geopandas

Tile-GeoPandas: Parallel geospatial data processing using tile-splitting patterns.

## Installation

```bash
pip install tile-geopandas
```

## Usage

```python
from tile_geopandas.core import TileGeoDataFrame

# Load data
tdf = TileGeoDataFrame.from_file('path/to/data.shp')

# Split into tiles of 1 degree, apply a function, and merge
result = (
    tdf
    .tile(1.0)
    .apply_to_tiles(lambda gdf: gdf.buffer(0.01))
    .merge()
)
```

## License

MIT
