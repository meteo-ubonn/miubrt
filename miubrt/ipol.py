# Copyright (c) 2023, miub developers.
# Distributed under the MIT License. See LICENSE.txt for more info.

import numpy as np
import wradlib as wrl
import xarray as xr
from osgeo import osr


def get_target_grid(da, res):
    if isinstance(res, tuple):
        res_x = res[0]
        res_y = res[1]
    else:
        res_x = res
        res_y = res

    proj = osr.SpatialReference()
    proj.ImportFromWkt(da["spatial_ref"].attrs["crs_wkt"])
    x0, y0 = wrl.georef.reproject(
        da.longitude.values, da.latitude.values, projection_target=proj
    )

    res_x0 = abs(da.x.max().values - da.x.min().values) / res_x
    num_pixels_x = int(res_x0 - res_x0 % 2) + 2
    res_y0 = abs(da.y.max().values - da.y.min().values) / res_y
    num_pixels_y = int(res_y0 - res_y0 % 2) + 2
    x_range = num_pixels_x // 2 * res_x
    y_range = num_pixels_y // 2 * res_y
    xgrid = (
        np.linspace(x0 - x_range, x0 + x_range, num_pixels_x, endpoint=False)
        + res_x / 2
    )
    ygrid = (
        np.linspace(y0 - y_range, y0 + y_range, num_pixels_y, endpoint=False)
        + res_y / 2
    )
    grid_xy_raw = np.meshgrid(xgrid, ygrid)
    grid_xy_grid = np.dstack((grid_xy_raw[0], grid_xy_raw[1]))
    return xgrid, ygrid, grid_xy_grid


def get_target_coordinates(grid):
    grid_xy = np.stack((grid[..., 0].ravel(), grid[..., 1].ravel()), axis=-1)
    return grid_xy


def get_source_coordinates(da):
    xy = np.stack((da.x.values.ravel(), da.y.values.ravel()), axis=-1)
    return xy


def create_dataarray(data, attrs, xgrid, ygrid):
    attrs.pop("_Undetect", None)
    # create x,y coordinate DataArrays
    x = xr.DataArray(
        xgrid,
        dims=["x"],
        attrs=dict(
            standard_name="projection_x_coordinate",
            long_name="x coordinate of projection",
            units="m",
        ),
    )
    y = xr.DataArray(
        ygrid,
        dims=["y"],
        attrs=dict(
            standard_name="projection_y_coordinate",
            long_name="y coordinate of projection",
            units="m",
        ),
    )
    # create data DataArray
    data = xr.DataArray(data, dims=["y", "x"], coords={"x": x, "y": y}, attrs=attrs)
    return data


def create_dataarray1(ds, data, moment):
    # get moment attrributes
    mom_attrs = ds[moment].attrs
    # remove _Undetect since it's not in the original data
    mom_attrs.pop("_Undetect", None)
    # create moment DataArray
    mom = xr.DataArray(data.astype(np.float32), dims=["y", "x"], attrs=mom_attrs)
    # fix encoding
    mom.encoding = dict(zlib=True, _FillValue=0.0, least_significant_digit=2)
    return mom


def create_dataarray_attrs(data, attrs):
    # get moment attrributes
    # mom_attrs = ds[moment].attrs
    # remove _Undetect since it's not in the original data
    attrs = attrs.copy()
    attrs.pop("_Undetect", None)
    # create moment DataArray
    mom = xr.DataArray(data.astype(np.float32), dims=["y", "x"], attrs=attrs)
    # fix encoding
    mom.encoding = dict(zlib=True, _FillValue=0.0, least_significant_digit=2)
    return mom


def polar_to_grid_array(da, proj, nb_pixels, method):
    """Polar to Cartesian Interpolation"""
    # georeference single sweep
    da = da.pipe(wrl.georef.georeference_dataset, proj=proj)
    # get source coordinates
    src = get_source_coordinates(da)
    # create target grid
    xgrid, ygrid, trg_grid = get_target_grid(da, nb_pixels)
    # get target coordinates
    trg = get_target_coordinates(trg_grid)
    # setup interpolator
    if method == "nearest":
        ip_ = wrl.ipol.Nearest(src, trg)
    elif method == "linear":
        ip_ = wrl.ipol.Linear(src, trg)
    elif method == "idw":
        ip_ = wrl.ipol.Idw(src, trg)
    else:
        raise ValueError(f"unknown interpolation method {method}")

    res = ip_(da.values.ravel()).reshape((len(ygrid), len(xgrid)))
    return create_dataarray(res, da.attrs, xgrid, ygrid)


def create_dataset(moments, xgrid, ygrid):
    """Create xarray dataset"""
    # create x,y DataArrays
    x = xr.DataArray(
        xgrid,
        dims=["x"],
        attrs=dict(
            standard_name="projection_x_coordinate",
            long_name="x coordinate of projection",
            units="m",
        ),
    )
    y = xr.DataArray(
        ygrid,
        dims=["y"],
        attrs=dict(
            standard_name="projection_y_coordinate",
            long_name="y coordinate of projection",
            units="m",
        ),
    )
    data_vars = dict()
    data_vars.update(moments)
    # data_vars.update({"time": ds.reset_coords().time})

    # create Dataset
    ds_out = xr.Dataset(
        data_vars=data_vars,
        coords={"x": x, "y": y},
    )
    # apply CF Conventions
    ds_out.attrs["Conventions"] = "CF-1.5"
    return ds_out


def to_cartesian(ds, proj=None, res=2000, method="idw"):
    """Georeference Dataset, Grid/Project moments to raster dataset with EPSG,
    output to netcdf

    Parameter
    ---------
    proj : GDAL OSR projection
        projection to use for cartesian interpolation, if None defaults to dataset spatial_ref (if available) or aeqd.
    res : int | float
        resolution in projection units
    method : str
        string describing interpolation method
    """
    if "time" in ds.dims:
        raise TypeError("Only works for one sweep, please select a single timestep.")

    has_spatial_ref = "spatial_ref" in ds.coords
    if proj is None:
        if has_spatial_ref:
            proj = osr.SpatialReference()
            proj.ImportFromWkt(ds["spatial_ref"].attrs["crs_wkt"])
        else:
            proj = wrl.georef.create_osr(
                "aeqd", lat_0=ds.latitude.values, lon_0=ds.longitude.values
            )

    if not has_spatial_ref:
        ds = ds.pipe(wrl.georef.georeference_dataset, proj=proj).rio.write_crs(
            proj.ExportToWkt()
        )

    # get source coordinates
    src = get_source_coordinates(ds)
    # create target grid
    xgrid, ygrid, trg_grid = get_target_grid(ds, res=res)
    # get target coordinates
    trg = get_target_coordinates(trg_grid)
    # setup interpolator
    if method == "nearest":
        ip_ = wrl.ipol.Nearest(src, trg)
    elif method == "linear":
        ip_ = wrl.ipol.Linear(src, trg)
    elif method == "idw":
        ip_ = wrl.ipol.Idw(src, trg)
    else:
        raise ValueError(f"unknown interpolation method {method}")

    # interpolate data
    moment_dict = {k: v for k, v in ds.data_vars.items()}
    for name, mom in ds.data_vars.items():
        data = ip_(mom.values.ravel()).reshape((len(ygrid), len(xgrid)))
        moment_dict[name] = create_dataarray_attrs(data, mom.attrs)
        ds_out = create_dataset(moment_dict, xgrid, ygrid)

    # mask values outside range
    maxrange = ds.range.max().values
    x, y = ds_out.x, ds_out.y
    ds_out = ds_out.where((x**2 + y**2) ** 0.5 < maxrange)

    # write current spatial_ref
    ds_out = ds_out.rio.write_crs(proj.ExportToWkt())

    # apply coordinates
    ds_out = ds_out.assign_coords(ds.time.coords)
    return ds_out
