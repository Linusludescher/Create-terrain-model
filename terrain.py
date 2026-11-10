"""Read zipped Kartverket DTM GeoTIFF tiles, validate, crop and export printable STL.

Coordinates are assumed projected in metres. No interpolation is used unless the
requested output grid is coarser than the original data. All heights are metres.
"""
from __future__ import annotations

import math
from pathlib import Path
from zipfile import ZipFile
from contextlib import ExitStack

import numpy as np
import rasterio
from rasterio.io import MemoryFile
from rasterio.merge import merge
from rasterio.enums import Resampling
import trimesh


def list_tiles(zip_path):
    """Read ONLY terrain DTM tiles; exclude overview/punktetthet/DSM GeoTIFFs."""
    with ZipFile(zip_path) as z:
        names = [n for n in z.namelist() if n.lower().endswith('.tif')
                 and '/data/dtm/' in n.lower() and not n.startswith('__MACOSX')]
        if not names:
            raise ValueError('Fant ingen GeoTIFF DTM-filer under data/dtm/ i ZIP-filen.')
        return sorted(names)


def _source(tile):
    return tile.split('/')[0]


def scan_archive(zip_path, read_heights=True):
    """Return per-tile metadata and global/group summary; stream one tile at a time."""
    rows = []
    with ZipFile(zip_path) as z:
        for name in list_tiles(zip_path):
            with MemoryFile(z.read(name)) as mem, mem.open() as ds:
                if not ds.crs or not ds.crs.is_projected:
                    raise ValueError(f'Projected CRS expected, got {ds.crs}: {name}')
                unit = ds.crs.linear_units
                if unit.lower() not in ('metre', 'meter', 'metres', 'meters', 'm'):
                    raise ValueError(f'Expected metre units, got {unit}: {name}')
                b = ds.bounds
                stat = {'source': _source(name), 'file': name, 'crs': str(ds.crs),
                        'resolution_x_m': abs(ds.res[0]), 'resolution_y_m': abs(ds.res[1]),
                        'columns': ds.width, 'rows': ds.height,
                        'min_x': b.left, 'min_y': b.bottom, 'max_x': b.right, 'max_y': b.top,
                        'valid_pixels': None, 'min_elevation_m': None, 'max_elevation_m': None}
                if read_heights:
                    a = ds.read(1, masked=True)
                    valid = a.compressed()
                    # Important: zero is valid height (sea level), not missing data.
                    stat['valid_pixels'] = int(valid.size)
                    if valid.size:
                        stat['min_elevation_m'] = float(valid.min())
                        stat['max_elevation_m'] = float(valid.max())
                rows.append(stat)
    return rows


def summarize(rows, source=None):
    matching = [r for r in rows if source is None or r['source'] == source]
    if not matching:
        raise ValueError(f'No tiles matching source {source!r}')
    crs = {r['crs'] for r in matching}
    if len(crs) != 1:
        raise ValueError(f'Multiple CRSs in selection: {crs}')
    xmin=min(r['min_x'] for r in matching); xmax=max(r['max_x'] for r in matching)
    ymin=min(r['min_y'] for r in matching); ymax=max(r['max_y'] for r in matching)
    mins=[r['min_elevation_m'] for r in matching if r['min_elevation_m'] is not None]
    maxs=[r['max_elevation_m'] for r in matching if r['max_elevation_m'] is not None]
    return {'tiles':len(matching), 'crs':next(iter(crs)),
            'raw_resolution_m': sorted(set((r['resolution_x_m'],r['resolution_y_m']) for r in matching)),
            'bbox_m':(xmin,ymin,xmax,ymax), 'bounding_width_m':xmax-xmin,
            'bounding_length_m':ymax-ymin, 'min_height_m':min(mins) if mins else None,
            'max_height_m':max(maxs) if maxs else None,
            'vertical_range_m':max(maxs)-min(mins) if mins and maxs else None,
            'note':'Bounding box includes gaps and overlaps; it is NOT necessarily continuous terrain coverage.'}


def crop_dem(zip_path, bbox, resolution_m=2.0, source=None, max_cells=8_000_000,
             fill_holes=False):
    """Read overlapping tiles directly from ZIP and combine into requested DEM.

    bbox=(xmin,ymin,xmax,ymax) in the GeoTIFF projected metre CRS.
    If source is None, only accept data from one collection in the crop.
    """
    xmin,ymin,xmax,ymax = map(float,bbox)
    if not (xmax>xmin and ymax>ymin and resolution_m>0):
        raise ValueError('Invalid bbox/resolution')
    ncols=math.ceil((xmax-xmin)/resolution_m)
    nrows=math.ceil((ymax-ymin)/resolution_m)
    if ncols*nrows>max_cells:
        raise ValueError(f'{ncols*nrows:,} cells exceeds max_cells={max_cells:,}; select a smaller bbox or larger resolution_m')
    with ZipFile(zip_path) as z:
        chosen=[]
        # Read individual tile headers only when selecting the crop.
        for name in list_tiles(zip_path):
            if source is not None and _source(name)!=source:
                continue
            with MemoryFile(z.read(name)) as mem, mem.open() as ds:
                b=ds.bounds
                if b.right>xmin and b.left<xmax and b.top>ymin and b.bottom<ymax:
                    chosen.append(name)
        if not chosen:
            raise ValueError('No DTM tiles intersect requested bbox/source.')
        sources={_source(n) for n in chosen}
        if len(sources)>1:
            raise ValueError(f'Multiple overlapping acquisitions {sources}; set source to select one explicitly.')
        with ExitStack() as stack:
            datasets=[]
            for name in chosen:
                mem=stack.enter_context(MemoryFile(z.read(name)))
                datasets.append(stack.enter_context(mem.open()))
            crs={str(ds.crs) for ds in datasets}
            if len(crs)>1:
                raise ValueError(f'Mixed CRSs: {crs}')
            raw_resolution=min(abs(ds.res[0]) for ds in datasets)
            if resolution_m < raw_resolution-1e-7:
                raise ValueError(f'Output resolution {resolution_m} m is finer than raw {raw_resolution} m. Oversampling is disallowed.')
            # first source used for overlapping pixels, no averaging across acquisitions.
            arr,transform=merge(datasets,bounds=bbox,res=(resolution_m,resolution_m),
                                masked=True,method='first',resampling=Resampling.bilinear)
            data=arr[0].astype('float32')
            if not np.ma.isMaskedArray(data):
                data=np.ma.masked_invalid(data)
            mask=np.ma.getmaskarray(data)|(~np.isfinite(np.asarray(data)))
            coverage=1.0-float(mask.mean())
            if mask.all():
                raise ValueError('Selection has no valid terrain samples')
            if mask.any():
                if not fill_holes:
                    raise ValueError(f'Crop includes {mask.mean():.1%} missing terrain cells. Choose a fully covered bbox or set fill_holes=True (estimated values).')
                from scipy.ndimage import distance_transform_edt
                indices=distance_transform_edt(mask,return_distances=False,return_indices=True)
                data=np.asarray(data)[tuple(indices)]
            else:
                data=np.asarray(data)
            return {'elevation_m':np.asarray(data,dtype='float32'), 'transform':transform,
                    'crs':str(datasets[0].crs), 'source':next(iter(sources)), 'coverage':coverage,
                    'raw_resolution_m':raw_resolution, 'output_resolution_m':resolution_m,
                    'tile_count':len(chosen), 'bbox_requested_m':tuple(bbox),
                    'bbox_actual_m':(transform.c,transform.f-data.shape[0]*resolution_m,
                                     transform.c+data.shape[1]*resolution_m,transform.f),
                    'filled_holes':bool(mask.any())}


def make_mesh(dem, model_width_mm=180.0, vertical_exaggeration=1.0,
              base_mm=3.0, max_faces=4_000_000):
    """Watertight rectangular solid. Height units and physical scale explicit."""
    h=np.asarray(dem['elevation_m'],dtype='float64')
    if h.ndim!=2 or min(h.shape)<2 or not np.isfinite(h).all():
        raise ValueError('Need a finite 2D grid of at least 2x2')
    nr,nc=h.shape
    faces_expected=4*(nr-1)*(nc-1)+4*(nr+nc-2)+2
    if faces_expected>max_faces:
        raise ValueError(f'{faces_expected:,} STL faces exceeds max_faces={max_faces:,}. Use a coarser resolution.')
    if not (model_width_mm>0 and vertical_exaggeration>0 and base_mm>0):
        raise ValueError('Model width, exaggeration and base must be positive')
    res=float(dem['output_resolution_m'])
    width_m=nc*res; length_m=nr*res
    # Pixel-centre geometry: represented domain width=(nc-1)*res, but nominal
    # bounding area is nc*res. Report both clearly to avoid scale ambiguity.
    horizontal_scale=model_width_mm/width_m
    model_length_mm=length_m*horizontal_scale
    x=np.linspace(0,model_width_mm,nc)
    y=np.linspace(model_length_mm,0,nr)  # raster row zero is north
    xx,yy=np.meshgrid(x,y)
    zz=base_mm+(h-h.min())*horizontal_scale*vertical_exaggeration
    upper=np.column_stack((xx.ravel(),yy.ravel(),zz.ravel()))
    bottom=np.column_stack((xx.ravel(),yy.ravel(),np.zeros(nr*nc)))
    vert=np.vstack([upper,bottom]); N=nr*nc
    rr=np.arange(nr-1)[:,None];cc=np.arange(nc-1)[None,:]
    a=(rr*nc+cc).ravel();b=a+1;c=a+nc;d=c+1
    top=np.vstack([np.stack([a,c,b],axis=1),np.stack([b,c,d],axis=1)])
    # Bottom grid has matching tessellation and opposite winding.
    bot=np.vstack([np.stack([a+N,b+N,c+N],axis=1),np.stack([b+N,d+N,c+N],axis=1)])
    boundary=np.concatenate([np.arange(nc),
                             np.arange(1,nr)*nc+(nc-1),
                             (nr-1)*nc+np.arange(nc-2,-1,-1),
                             np.arange(nr-2,0,-1)*nc])
    bnext=np.roll(boundary,-1)
    sides=np.vstack([np.stack([boundary,bnext,boundary+N],axis=1),
                     np.stack([bnext,bnext+N,boundary+N],axis=1)])
    faces=np.vstack([top,bot,sides])
    mesh=trimesh.Trimesh(vertices=vert,faces=faces,process=False)
    trimesh.repair.fix_normals(mesh)
    info={'terrain_width_m':width_m,'terrain_length_m':length_m,
          'raw_height_min_m':float(h.min()),'raw_height_max_m':float(h.max()),
          'real_height_range_m':float(np.ptp(h)), 'model_width_mm':model_width_mm,
          'model_length_mm':model_length_mm,'model_total_height_mm':float(zz.max()),
          'vertical_exaggeration':vertical_exaggeration,'horizontal_scale_ratio':1/(horizontal_scale/1000),
          'triangles':len(mesh.faces),'watertight':bool(mesh.is_watertight),
          'winding_consistent':bool(mesh.is_winding_consistent)}
    if not mesh.is_watertight:
        raise RuntimeError('Generated STL mesh is not watertight')
    return mesh,info


def export_stl(dem, output_path, **model_parameters):
    mesh,info=make_mesh(dem,**model_parameters)
    output_path=Path(output_path);output_path.parent.mkdir(parents=True,exist_ok=True)
    mesh.export(str(output_path), file_type='stl')
    info['path']=str(output_path)
    info['file_size_mb']=round(output_path.stat().st_size/1e6,2)
    return info
