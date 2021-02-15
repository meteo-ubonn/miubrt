#!/usr/bin/env python
# Copyright (c) 2019-2020, miub developers.
# Distributed under the MIT License. See LICENSE.txt for more info.

import datetime as dt
import numpy as np
import xarray as xr



def calib_mask(ds,i, par):
    """Calibrate Polarimetric Moments"""
    height = ds['z']
    if (i==0):
        zh    = ds['DBZH'].where((ds.RHOHV > par.rhohv_th) & (ds.PHIDP < par.phidp_th)
                              & (ds.DBZH < par.dbz_th) & (ds.DBZH > par.dbz_thm)
                              & (height < par.z_melt) & (height > par.z_noise))
        zdr    = ds['ZDR'].where((ds.RHOHV > par.rhohv_th) & (ds.PHIDP < par.phidp_th)
                              & (ds.DBZH < par.dbz_th) & (ds.DBZH > par.dbz_thm)
                              & (height < par.z_melt) & (height > par.z_noise))
        kdp    = ds['KDP'].where((ds.RHOHV > par.rhohv_th) & (ds.PHIDP < par.phidp_th)
                              & (ds.DBZH < par.dbz_th) & (ds.DBZH > par.dbz_thm)
                              & (height < par.z_melt) & (height > par.z_noise))
        rhohv    = ds['RHOHV'].where((ds.RHOHV > par.rhohv_th) & (ds.PHIDP < par.phidp_th)
                              & (ds.DBZH < par.dbz_th) & (ds.DBZH > par.dbz_thm)
                              & (height < par.z_melt) & (height > par.z_noise))
        phidp    = ds['PHIDP'].where((ds.RHOHV > par.rhohv_th) & (ds.PHIDP < par.phidp_th)
                              & (ds.DBZH < par.dbz_th) & (ds.DBZH > par.dbz_thm)
                              & (height < par.z_melt) & (height > par.z_noise))                         
    elif (i==1):
        zh    = ds['DBZH']
        zdr   = ds['ZDR']
        kdp   = ds['KDP']
        rhohv   = ds['RHOHV']
        phidp   = ds['PHIDP']
    
    return zh, zdr, kdp, rhohv, phidp, height


def phase_offset(ds,rhohv_min):
    """Phase offset via phase histogram"""
    phioff = ds.PHIDP.copy().where((ds.RHOHV>=rhohv_min) & (ds.DBZH>=0))
    #create binary array
    phib   = xr.where(np.isnan(phioff), 0, 1) 
    #calculuate rolling sum
    phib_sum = phib.rolling(range=10, center=True).sum(skipna=True)
    #find first maximum
    start_range = phib_sum.idxmax(dim="range")
    #calculate offset per ray
    # range in m starting from first valid bin
    rng = 10000
    # add range 
    stop_range = start_range + rng

    # get phase values in specified range
    off = phioff.where((phioff.range >= start_range) & (phioff.range <= stop_range))
    
    # calculate nan median over range
    off = off.median(dim='range', skipna=True)
    
    #hist = off.hvplot.line(groupby="time", x="azimuth",
    #                      #xlim=(100, 160), 
    #                      ylim=(120, 150),
    #                     )
    #hist
    #calculate offset per sweep
    from dask.diagnostics import ProgressBar
    with ProgressBar():
        # calculate median per sweep
        phi_offset1 = off.load().median(dim='azimuth', skipna=True)

    offset = phi_offset1.copy()
    # Fill start of radial with offset median
    phi_fix = ds.PHIDP.where(ds.PHIDP.range > start_range).fillna(offset)
    
    phi_fix = phi_fix - offset
    
    #dealias
    phi_fix = xr.where((phi_fix <= (offset-180)), phi_fix + 360., phi_fix)
    phi_fix = phi_fix.where(ds.VRADH > ds.VRADH.min())
    return(phi_fix)

def phidp_filter(phidpR,rhohv,len1,len2,na,nrange,az):
    nfin = np.empty([na, nrange])
    phib   = xr.where(np.isnan(phidpR), 0, 1) 
    lphi = 0
    lconf = len2
    lra = 0
    for ri in range(nrange):
        up = np.minimum(ri+len1,nrange-2)
        down = np.maximum(ri-len1,0)
        phit = phib[az,down:up]
        indphit = np.where(phit ==1)
        phis = np.sum(phit)
        if (phis>0):
            nfin[az,ri] = np.sum(rhohv.values[az,down+indphit],axis=1)
        if (phib[az,ri]==1):
            if ((phidpR[az,ri]<lphi-1) | (phidpR[az,ri]>lphi+5)):
                if (nfin[az,ri]>lconf):
                    lphi = phidpR[az,ri]
                    lconf = nfin[az,ri]
                    phidpR[az,ri] = np.nan
                else:
                    phidpR[az,ri] = np.nan
            else:
                lphi = phidpR[az,ri]
                lconf = nfin[az,ri]
                lra   = ri     
    return(phidpR)