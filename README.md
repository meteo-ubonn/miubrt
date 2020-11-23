# README #

This is a collection of tools for use with the X-band radar data from BoXPol/JuXPol radar and DWD C-Band radar network.
All processing are done in their native polar coordinates.

### Documentation ###

# BoXPol data #
Data Format    : mvol ~ hdf5 format, ODIM_H5 (GAMIC users)
Data Structure : PPI  ~ 288 files for 1 elevation angle stored in named elevation folders
               : RHI  ~ 1 azimuth angle (LACROS)

# JuXPol Data #
Data Format    : mvol ~ hdf5 format, ODIM_H5 (GAMIC users)
Data Structure : PPI  ~ 288 files for 10 elevation angle stored in one common folder
               : RHI  ~ 1 azimuth angle


# Radar data processing flowchart (X-band) #

1. Process netcdf files in native polar-cordinates for all PPI sweeps and RHI
2. Apply calibration to Zh and Zdr
3. Remove ground clutter and non-meteorological data
4. Attenuation and Differential attenuation correction
5. Estimate KDP
