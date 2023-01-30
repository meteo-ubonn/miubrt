#!/usr/bin/env python
# Copyright (c) 2023, miub developers.
# Distributed under the MIT License. See LICENSE.txt for more info.

import datetime as dt
import os

import numpy as np
import wradlib as wrl
import xarray as xr

from .util import create_xpol_filelist, create_cpol_filelist


def get_xpol_path(inpath=None, start_time=None, loc="boxpol"):
    """Create path of BoXPol/JuXPol radar data files

    Parameter
    ---------
    path : str
        Path to root folder of radar data,
        defaults to None (/automount/radar/scans or /automount/radar-archiv/scans)
    start_time : datetime.datetime
        datetime - object to select correct folder
    loc : str
        "boxpol" or "juxpol" (not case sensitive)

    Return
    ------
    radar_path : str
        Path to radar data

    """
    if start_time is None:
        start_time = dt.datetime.today()
    loc = "" if loc.lower()[0:2] == "bo" else "_juelich"
    if inpath is None:
        ins = "-archiv" if start_time < dt.datetime(2015, 1, 1) else ""
        inpath = f"/automount/radar{ins}/scans{loc}"
        if not os.path.exists(inpath):
            # try MIUBRT_DATA
            if os.environ.get("MIUBRT_DATA", None) is not None:
                inpath = os.environ["MIUBRT_DATA"]
            else:
                raise ValueError(
                    f"Data Path: '{inpath}' not found and 'MIUBRT_DATA' environment variable not set."
                )

    radar_path = os.path.join(inpath, "{0}/{0}-{1:02}/{0}-{1:02d}-{2:02d}")
    return radar_path.format(start_time.year, start_time.month, start_time.day)


def get_cpol_path(path=None, date=None, loc=None, sname=None, elev_id=None):
    """Create path of DWD radar data files.

    Parameter
    ---------
    path : str
        Path to root folder of radar data,
        defaults to None (/automount/realpep/upload/RealPEP-SPP/DWD-CBand)
    start_time : datetime.datetime
        datetime - object to select correct folder
    loc : str
        3 char DWD radar abbrevation (eg. "ESS", not case sensitive)

    Return
    ------
    radar_path : str
        Path to radar data
    """
    if path is None:
        path = f'/automount/realpep/upload/RealPEP-SPP/DWD-CBand'
    if date is None:
        raise ValueError("need date")
    if loc is None:
        raise ValueError("need loc")
    date = f'{date.year}/{date.year}-{date.month:02d}/{date.year}-{date.month:02d}-{date.day:02d}'
    radar_path = os.path.join(path, date)
    radar_path = os.path.join(radar_path, f"{loc.lower()}")
    radar_path += f"/{sname}/{elev_id:02d}/ras07-{sname}_sweeph5onem_"
    return radar_path


def to_netcdf(ds, fpath, engine="netcdf4"):
    """Export Dataset to NetCDF file"""
    delayed_obj = ds.to_netcdf(fpath, engine=engine, compute=False)
    from dask.diagnostics import ProgressBar

    with ProgressBar():
        delayed_obj.compute()


def save_netcdf_dataset(ds, path, post=None, engine="netcdf4"):
    """Save Dataset to NetCDF file"""
    elevation = np.median(ds.elevation)
    t0 = ds.time[0].values.astype("M8[s]").astype("O")
    t0 = t0.replace(second=0, microsecond=0, minute=((t0.minute // 5) * 5))
    t1 = ds.time[-1].values.astype("M8[s]").astype("O") + dt.timedelta(minutes=5)
    t1 = t1.replace(second=0, microsecond=0, minute=((t1.minute // 5) * 5))
    fname = [
        f"{ds.location.values}",
        f"{int(round(elevation, 1) * 10):03d}",
        f"{t0:%Y%m%d}",
        f"{t0:%H%M}",
        f"{t1:%H%M}",
    ]
    if post is not None:
        fname.append(post)
    fname = "_".join(fname) + ".nc"
    ofname = os.path.join(path, fname)
    ds.pipe(to_netcdf, ofname, engine=engine)
    return ofname


def create_netcdf_dataset(
    location,
    name,
    start_time,
    end_time,
    sweep,
    inpath=None,
    outpath="",
    chunks={},
    moments=None,
    engine="h5netcdf",
):
    """Create NetCDF file from radar data

    This function automagically creates an xarray Dataset conssisting of
    georeferenced radar data.

    Parameter
    ---------
    location : str
        "boxpol" or "juxpol" or 3-char DWD radar string
    name : str
        Sweep name (eg. "n_ppi_010deg" for XBand or "00" for DWD)
    start_time : datetime.datetime
        datetime - object to select start time
    end_time : datetime.datetime
        datetime - object to select end time
    sweep : int
        If multiple sweeps in one file, select sweep number,
        defaults to 0

    Return
    ------
    ds : xarray.Dataset
        Dataset of radar sweep data.
    """
    if "xpol" in location.lower():
        radar_path = get_xpol_path(inpath=inpath, start_time=start_time, loc=location)
        print("radar path:", radar_path)
        read_engine = wrl.io.backends.GamicBackendEntrypoint
        group = f"scan{sweep}"
        file_path = os.path.join(radar_path, name)
    else:
        file_path = get_cpol_path(date=start_time, loc=location, sname=name, elev_id=sweep)
        print("radar path:", radar_path)
        read_engine = wrl.io.backends.OdimBackendEntrypoint
        group = f"dataset{sweep+1}"
    if "xpol" in location.lower():
        file_obj = list(
            create_xpol_filelist(os.path.join(file_path, "*"), start_time, end_time))
        ds = wrl.io.open_radar_mfdataset(file_obj, engine=read_engine, group=group)
    else:
        file_obj = list(
            create_cpol_filelist(file_path, start_time, end_time, moments))
        cpol_names = np.array(file_obj).reshape((len(moments), -1)).T.tolist()
        ds = xr.open_mfdataset(cpol_names, concat_dim=["time", None], combine="nested",
                               engine="wradlib-odim", group="dataset1")
    ds = ds.assign_coords({"location": location})
    # ds = ds.chunk({"time": 24})
    ofname = save_netcdf_dataset(ds, outpath, engine=engine)
    del ds
    return os.path.abspath(ofname)
