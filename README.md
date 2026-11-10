# ZIP DTM -> STL (Jupyter + uv)

Free/open-source Python workflow for Norwegian elevation ZIP exports (such as those from Høydedata.no). **This package reads tiled GeoTIFF DTM files under `*/data/dtm/*.tif`**. It does not treat SOSI files, metadata overviews, or digital surface models as terrain.

## Quick start (Fedora / Linux)

1. Install `uv` if you don't already have it: https://docs.astral.sh/uv/getting-started/installation/
2. Unzip this PROJECT archive into a folder.
3. Put your terrain ZIP (for example `eksport_1811757_20261105.zip`) into the **`input/`** folder. Do not unzip the terrain archive.
4. Open a terminal in the project folder and run:

   ```bash
   uv sync
   uv run jupyter lab
   ```

5. Open `01_terrain_to_stl.ipynb` and run cells top to bottom. Files are saved in `output/`.

You can use an existing Python environment with `pip install -e .` instead of uv.

## What this notebook reports

- Input DTM tiles and acquisition name (data from distinct acquisitions remain separate)
- Native raster cell size in metres (not the same as measurement accuracy)
- CRS and bounding box in projected metre coordinates
- Bounding rectangle width × length in metres, and min/max elevation in metres
- Real elevation range (maximum minus minimum, **not** absolute altitude)
- Coverage and elevation range of your selected crop, without need to mosaic the entire project
- Projected horizontal print size, vertical exaggeration, printed model height, approximate scale, STL face count, water-tightness

### What these dimensions mean

The reported dimensions for an acquisition are the **bounding box around all its tiles**, not necessarily an entirely filled rectangle. For each crop, the notebook uses your selected rectangular area. Ocean may have an elevation of zero; zero is NOT considered missing. Masked rasters/true no-data cells are treated separately.

## First recommended example for the attached export

```
ZIP: eksport_1811757_20261105.zip
Source: NDH Vågan-Hadsel sør 2pkt 2017
BBox (UTM 33N metres): (468800, 7563000, 469600, 7563600)
Raw cell size: 0.5 m
Crop dimensions: 800 m × 600 m
Output grid: 2 m (fast initial preview)
```

After inspecting the first model, try `OUTPUT_RESOLUTION_M = 0.5` to keep the full native grid resolution. A native-resolution STL for an 800 m × 600 m region can contain **millions of triangles**, so allow enough memory and disk space and increase the configured mesh face limit if necessary.

## Parameters

Edit in the notebook:

- `ZIP_PATH`: zipped elevation export, or use automatic detection of `input/*.zip`
- `SOURCE`: acquisition to use (important when tiles from different years overlap)
- `BBOX`: `(min_x, min_y, max_x, max_y)` in the dataset CRS, in metres
- `OUTPUT_RESOLUTION_M`: grid spacing for STL; must be at least the native input resolution
- `MODEL_WIDTH_MM`: physical printed width of the STL
- `VERTICAL_EXAGGERATION`: 1.0 = true proportions; 2.0 doubles the relief only
- `BASE_MM`: thickness of flat supporting base
- `FILL_HOLES`: by default **False**; if enabled, missing cells are filled with nearest available values, which may be false terrain
- `MAX_CELLS`, `MAX_FACES`: protection from accidentally creating enormous models

The program produces a **solid watertight rectangular STL** with a bottom and walls, ready for inspection in Blender, PrusaSlicer, or another mesh/slicer tool. STL uses numerical coordinate units but does not tag them; here all model coordinates are **millimetres**.

### Important limitations

- Geometry is sampled at **cell centres**, scaled to the chosen requested ground footprint. The representative outer edge is the declared rectangular crop extent, so nominal width and length equal the specified bbox, up to output-grid rounding.
- Height extrema are taken from pixel values. Elevation-model grid spacing is not the same as surveying accuracy.
- Overlaps between different collection years are **not** silently mixed: choose the `SOURCE` explicitly. Within one acquisition, the first available tile wins in overlaps.
- Large crops can be hundreds of millions of cells: start with 2–5 m output resolution or a smaller bbox.
- This is a DTM (bare-ground) workflow; it does not overlay Norgeskart textures, roads or vegetation.
- The notebook can optionally read local N50 File Geodatabase data from `n50_input/Basisdata_25833_N50Kartdata_FGDB.gdb`. It previews all N50 features and sea/lake polygons. Set `ADD_LAKES = True` to draw the original N50 vector boundaries directly over the terrain preview.
- STL does not support textures or material colours. With `ADD_LAKES = True`, the notebook additionally exports a watertight water-overlay STL. Its coastline comes from the N50 polygons (not DTM cells), while its height is sampled from the terrain. Import the two STLs together in a multi-material slicer and assign a separate filament to the water part.
- Do not use filled/interpolated regions as authoritative height data.
