#!/usr/bin/env python
# Copyright (c) 2023, miub developers.
# Distributed under the MIT License. See LICENSE.txt for more info.

import datetime as dt
import os

import numpy as np
import wradlib as wrl

from .util import create_filelist


def get_xpol_path(inpath=None, start_time=None, loc="boxpol"):
    """Create path of BoXPol/JuXPol radar data files

    Parameter
    ---------
    path : str
        Path to root folder of radar data, Defaults to None (/auotomount)

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
        f"{ds.location.values} ",
        f"{int(round(elevation, 1) * 10)}",
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
    engine="h5netcdf",
):
    """Create NetCDF file from radar data"""
    radar_path = get_xpol_path(inpath=inpath, start_time=start_time, loc=location)
    file_path = os.path.join(radar_path, name)
    file_obj = list(create_filelist(os.path.join(file_path, "*"), start_time, end_time))
    vol = wrl.io.open_odim(file_obj, loader="h5py", flavour="GAMIC", chunks=chunks)
    ds = vol[sweep].data
    ds = ds.assign_coords({"location": location})
    # ds = ds.chunk({"time": 24})
    ofname = save_netcdf_dataset(ds, outpath, engine=engine)
    del vol
    return os.path.abspath(ofname)
